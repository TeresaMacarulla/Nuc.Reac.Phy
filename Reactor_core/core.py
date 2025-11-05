import numpy as np
import glob
import matplotlib.pyplot as plt
import openmc
import math
from pathlib import Path
from pins import build_pin_universes
import openmc.deplete as deplete

def build_assembly_universe(assembly_x=17, assembly_y=17, assembly_z=150, pinx=1.4, piny=1.4, pinz=2, with_ctrl_3x3=True):
    
    # Use pin-cells
    f, c, w, materials, light_water = build_pin_universes()

    # 3x3 central map for control rows
    cx, cy = assembly_x//2, assembly_y//2
    ctrl_i = {cx-1, cx, cx+1}
    ctrl_j = {cy-1, cy, cy+1}

    universii = np.empty((assembly_z, assembly_x, assembly_y), dtype=object)
    for z in range(assembly_z):
        for i in range(assembly_x):
            for j in range(assembly_y):
                universii[z,i,j] = (c if with_ctrl_3x3 and (i in ctrl_i and j in ctrl_j) else f)

    lat = openmc.RectLattice(name="fuel_assembly")
    lat.lower_left = (-assembly_x*pinx/2, -assembly_y*piny/2, -assembly_z*pinz/2)
    lat.pitch = (pinx, piny, pinz)
    lat.universes = universii.tolist()

    # Outside pin grid, water 
    water_cell = openmc.Cell(fill=light_water)
    outer_uni = openmc.Universe(cells=[water_cell])
    lat.outer = outer_uni

    # Save the assembly into a Universe
    lat_region = +openmc.ZPlane(-assembly_z*pinz/2) & -openmc.ZPlane(+assembly_z*pinz/2) & \
                 openmc.model.rectangular_prism(width=assembly_x*pinx, height=assembly_y*piny, boundary_type='transmission')
    lat_cell = openmc.Cell(fill=lat, region=lat_region)
    assm_uni = openmc.Universe(cells=[lat_cell])

    return assm_uni, materials, light_water


def build_core():
    # 1) Assembly woth no control rows to start with criticality
    assembly_uni, materials, light_water = build_assembly_universe(with_ctrl_3x3=False)

    # 2) (Optional) 3x3 rows for insertion
    ctrl_assembly_uni, *_ = build_assembly_universe(with_ctrl_3x3=True)

    # 3) Assembly lattice (e.g., 9x9)
    nx = ny = 9
    pitch_assm = 23.8  # 17 * 1.4 cm (from assembly.py)
    latt = openmc.RectLattice(name="core_lattice")
    latt.lower_left = (-nx*pitch_assm/2, -ny*pitch_assm/2, -150*2/2)  # high
    latt.pitch = (pitch_assm, pitch_assm, 2)  # 2 cm for each z surface
    core_cells = np.full((150, nx, ny), assembly_uni, dtype=object) # fill with fuel assemblies

    # Place some assemblies with control rows 
    #mid = nx//2
    #for (i,j) in [(mid,mid), (mid-2,mid), (mid+2,mid), (mid,mid-2), (mid,mid+2)]:
        #core_cells[:, i, j] = ctrl_assembly_uni
    core_cells[:, 1:8, 1:8] = ctrl_assembly_uni

    latt.universes = core_cells.tolist()

    # Outside assemblies grid, water
    outer_uni = openmc.Universe(cells=[openmc.Cell(fill=light_water)])
    latt.outer = outer_uni

    # 4) Circular section of the core + reflectors and vessel
    hz = 150*2/2
    zmin = openmc.ZPlane(-hz, boundary_type='transmission')
    zmax = openmc.ZPlane(+hz, boundary_type='transmission')

    core_radius = 105.0   # ~ 9*23.8/2 cm ≈ 107 cm
    r_core = openmc.ZCylinder(r=core_radius, boundary_type='transmission')

    # reflectors/vessel
    r_reflector = openmc.ZCylinder(r=core_radius+25.0, boundary_type='transmission')  # ~25 cm water
    r_vessel    = openmc.ZCylinder(r=core_radius+25.0+10.0, boundary_type='vacuum')   # ~10 cm steel + vacuum outside

    # materials
    steel = openmc.Material(name="Steel"); steel.add_element("Fe", 1.0); steel.set_density("g/cm3", 7.9)
    mats = openmc.Materials(materials + [steel]); mats.export_to_xml()

    # cells by region
    core_cell    = openmc.Cell(fill=latt,        region = (-r_core) & +zmin & -zmax)
    refl_cell    = openmc.Cell(fill=light_water, region = (+r_core & -r_reflector) & +zmin & -zmax)
    vessel_cell  = openmc.Cell(fill=steel,       region = (+r_reflector & -r_vessel) & +zmin & -zmax)
    outside_cell = openmc.Cell(region = (+r_vessel) | -zmin | +zmax) 

    root = openmc.Universe(cells=[core_cell, refl_cell, vessel_cell, outside_cell])

    geom = openmc.Geometry(root)

    # Make sure that materials and geometry are updated
    all_mats = list(geom.get_all_materials().values())
    openmc.Materials(all_mats).export_to_xml()

    geom.export_to_xml()

    # Simulation parameters
    settings = openmc.Settings()
    settings.batches = 80 #150
    settings.inactive = 20 #30
    settings.particles = 10000 #30000
    bounds = [-core_radius, -core_radius, -hz, core_radius, core_radius, hz]
    src = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    settings.source = openmc.IndependentSource(space=src)
    settings.export_to_xml()

    return geom

def plot_core_xy(filename="core_xy.png", width_cm=260.0, pixels=1400):

    plot = openmc.Plot()
    plot.filename = filename.rsplit(".", 1)[0]  
    plot.basis = 'xy'
    plot.origin = (0.0, 0.0, 0.0)     
    plot.width = (width_cm, width_cm) 
    plot.pixels = (pixels, pixels)
    plot.color_by = 'material'        

    plots = openmc.Plots([plot])
    plots.export_to_xml()

    openmc.plot_geometry()
    print(f"Plot generated: {filename}")

def core_flux_tallies(core_radius, hz, nx=200, ny=200, nz=1):

    # Mesh covering the core + reflector + vessel
    mesh = openmc.RegularMesh()
    mesh.dimension   = [nx, ny, nz]
    mesh.lower_left  = [-core_radius, -core_radius, -hz]
    mesh.upper_right = [ core_radius,  core_radius,  hz]
    mf = openmc.MeshFilter(mesh)

    tallies = openmc.Tallies()

    # Total flux
    t_flux = openmc.Tally(name="core_flux_mesh")
    t_flux.filters = [mf]
    t_flux.scores  = ["flux"]
    tallies.append(t_flux)
    tallies.export_to_xml()

def plot_flux(tally_name="core_flux_mesh", out_pdf="flux_XY.pdf"):
    
    # Localiza el último statepoint
    sp_files = sorted(glob.glob("statepoint.*.h5"))
    if not sp_files:
        raise RuntimeError("No statepoint.*.h5 found. Run the simulation first.")
    sp_path = sp_files[-1]

    with openmc.StatePoint(sp_path) as sp:
        # Recupera el tally por NOMBRE (no por índice)
        tally = sp.get_tally(name=tally_name)

        # La primera filter debe ser el MeshFilter
        mesh = tally.filters[0].mesh
        nx, ny, nz = mesh.dimension

        # mean es un vector; dale forma según la malla
        data = tally.mean.ravel(order='C')
        data3d = data.reshape((nx, ny, nz), order='C')  # (x,y,z)

        # Integra en z (si nz=1, es lo mismo)
        img = data3d.sum(axis=2).T  # transpuesta para que Y vaya vertical

    # Plot absoluto
    plt.figure(figsize=(6,5))
    im = plt.imshow(img, origin='lower', interpolation='nearest',
                    extent=[mesh.lower_left[0], mesh.upper_right[0],
                            mesh.lower_left[1], mesh.upper_right[1]])
    plt.colorbar(label="Flux")
    plt.xlabel("x [cm]"); plt.ylabel("y [cm]")
    plt.tight_layout()
    plt.savefig(out_pdf, dpi=200)
    plt.close()

def calc_fuel_volume_cm3(nx=9, ny=9, ax=17, ay=17, n_ctrl_assemblies=49, ctrl_ax=3, ctrl_ay=3, r_fuel_cm=0.40, nz_layers=150, dz_cm=2.0, nonfuel_pins_per_assembly=0):

    n_assemblies = nx * ny
    pins_per_assembly = ax * ay
    total_pins_all = n_assemblies * pins_per_assembly
    total_ctrl_pins = n_ctrl_assemblies * (ctrl_ax * ctrl_ay)
    total_other_nonfuel = n_assemblies * nonfuel_pins_per_assembly

    n_fuel_pins = total_pins_all - total_ctrl_pins - total_other_nonfuel

    H_core_cm = nz_layers * dz_cm
    pin_area = math.pi * (r_fuel_cm ** 2)  # cm^2
    V_fuel_cm3 = n_fuel_pins * pin_area * H_core_cm
    return V_fuel_cm3

def set_uo2_volume(geom_or_materials, volume_cm3):

    # Obtener lista de materiales
    if isinstance(geom_or_materials, openmc.Geometry):
        mats = geom_or_materials.get_all_materials().values()
    elif isinstance(geom_or_materials, openmc.Materials):
        mats = geom_or_materials
    elif isinstance(geom_or_materials, (list, tuple)):
        mats = geom_or_materials
    else:
        raise TypeError("Provide openmc.Geometry, openmc.Materials or materials list.")

    # Search "Uranium dioxide"
    target1 = None
    target2 = None
    for m in mats:
        if (m.name or "") == "Uranium dioxide":
            target1 = m
            break
    for m in mats:
        if (m.name or "") == "Uranium dioxide":
            target2 = m

    target1.volume = float(volume_cm3)
    target2.volume = float(volume_cm3)
    return target1, target2


if __name__ == "__main__":

    print(" ")
    print("-------------------------------------------------------------------------------")
    print(" ")
    print("1: Plot an xy plane of the core.")
    print("2: Run the simulation.")
    print("3: Run the simulation fixing the fuel temperature.")
    print("4: Plot neutron flux distribution.")
    print("5: Run the reactor for 6 months at a power level of 1200 MW thermal.")
    number = float(input("Choose a number: "))

    if number == 1:
        # Build and export the core model
        geom = build_core()
        plot_core_xy(filename="core_xy.png", width_cm=280.0, pixels=1400)
        print("See the core in: <<core_xy.png>> inside this folder")

    elif number == 2:

        geom = build_core()
        openmc.run()

    elif number == 3:
        Tf = float(input("Choose a fuel temperature [K]: "))
        geom = build_core()

        for cell in geom.get_all_cells().values():
            if isinstance(cell.fill, openmc.Material) and cell.fill.name == "Uranium dioxide":
                cell.temperature = Tf
            geom.export_to_xml()

        settings = openmc.Settings()
        settings.temperature = {"default": 600, "method": "interpolation"}
        settings.batches, settings.inactive, settings.particles = 150, 30, 30000
        settings.source = openmc.IndependentSource(
            space=openmc.stats.Box([-105.0,-105.0,-150*2/2],[105.0,105.0,150*2/2], only_fissionable=True))
        settings.export_to_xml()

        openmc.run() 

    elif number == 4:
        geom = build_core()

        # Create tallies 
        core_flux_tallies(core_radius=140.0, hz=150*2/2, nx=220, ny=220, nz=1)
        openmc.run()
        plot_flux()
        print("See the core in: <<flux_XY_relative.pdf>> inside this folder")

    elif number == 5:
        
        openmc.reset_auto_ids()
        geom = build_core()

        Vfuel = calc_fuel_volume_cm3()

        # Assign the volume to the UO2 material.
        n_marked = set_uo2_volume(geom, Vfuel)

        # Re-exports materials so that they remain persistent in XML.
        openmc.Materials(list(geom.get_all_materials().values())).export_to_xml()

        # Designating the model to run depletion on
        reactor_model = openmc.Model()
        reactor_model.materials = openmc.Materials.from_xml('materials.xml')
        reactor_model.geometry  = openmc.Geometry.from_xml('geometry.xml')
        reactor_model.settings  = openmc.Settings.from_xml('settings.xml')
        reactor_model.tallies = openmc.Tallies.from_xml('tallies.xml')

        #Designating the depletion chain being used for our depletion calculation
        chain = openmc.deplete.Chain.from_xml("./chain_endfb80_pwr.xml")
        operator = openmc.deplete.CoupledOperator(reactor_model, "./chain_endfb80_pwr.xml")
        
        # Designating power settings in watts
        power = 1200e6 
        # Number of timesteps and their length in this case 6 steps of one month length in seconds.
        time_steps = [30*24*60*60]*6

        # Running the depletion calculation
        integrator = openmc.deplete.PredictorIntegrator(operator, time_steps, power)
        integrator.integrate()
        
        """
        # Leer resultados
        res = deplete.Results("depletion_results.h5")

        # API puede ser get_eigenvalue() o get_keff() según versión
        try:
            times, keff, keff_sd = res.get_eigenvalue()          # tiempos [s], k mean, k std
        except AttributeError:
            times, (keff, keff_sd) = res.get_keff(time_units='s')  # alternativa

        # Pasar a días para el plot
        t_days = np.array(times, dtype=float)/86400.0

        # Graficar
        plt.figure()
        plt.errorbar(t_days, keff, yerr=keff_sd, marker='o', linestyle='-')
        plt.xlabel("Time [days]")
        plt.ylabel(r"$k_{\mathrm{eff}}$")
        plt.title(r"$k_{\mathrm{eff}}$ during first 6 months at 1200 MW$_\mathrm{th}$")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("keff_vs_time_burnup.pdf", dpi=160)
        plt.show()
        """

    else:
        print("Not an option.")





