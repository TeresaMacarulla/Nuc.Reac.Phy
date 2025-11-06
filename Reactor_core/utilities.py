import numpy as np
import glob
import matplotlib.pyplot as plt
import openmc
import math
from pathlib import Path
from pins import build_pin_universes

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
    # 1) Assembly with no control rows to start with criticality
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
    core_cells[:, 1:8, 1:8] = assembly_uni #ctrl_assembly_uni

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
    settings.batches = 150
    settings.inactive = 30
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
