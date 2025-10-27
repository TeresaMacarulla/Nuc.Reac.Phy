import numpy as np
import openmc
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

    # example: place 4 assemblies with control rows 
    mid = nx//2
    #for (i,j) in [(mid,mid), (mid-2,mid), (mid+2,mid), (mid,mid-2), (mid,mid+2)]:
    #for (i,j) in [(mid,mid), (mid-2,mid), (mid+2,mid), (mid,mid-2), (mid,mid+2), (mid+2,mid+2), (mid-2,mid+2), (mid+2,mid-2), (mid-2,mid-2), (mid,mid-4), (mid,mid+4), (mid-4,mid), (mid+4,mid)]:
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
    #r_reflector = openmc.ZCylinder(r=core_radius+10.0, boundary_type='transmission')  # ~10 cm water
    r_vessel    = openmc.ZCylinder(r=core_radius+25.0+10.0, boundary_type='vacuum')   # ~10 cm steel + vacuum outside
    #r_vessel    = openmc.ZCylinder(r=core_radius+10.0+25.0, boundary_type='vacuum')   # ~20 cm steel + vacuum outside

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
    settings.particles = 30000
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

if __name__ == "__main__":
    # Build and export the core model
    geom = build_core()

    plot_core_xy(filename="core_xy.png", width_cm=280.0, pixels=1400)

    # (Optional) run the simulation after plotting
    openmc.run()

    #k-effective (Collision)     = 1.25289 +/- 0.00055
    #k-effective (Track-length)  = 1.25364 +/- 0.00067
    #k-effective (Absorption)    = 1.25337 +/- 0.00055
    #Combined k-effective        = 1.25321 +/- 0.00043
    #Leakage Fraction            = 0.00420 +/- 0.00005

