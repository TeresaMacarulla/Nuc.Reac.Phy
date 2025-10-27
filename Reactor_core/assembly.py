# assembly.py
import numpy as np
import openmc

from pins import build_pin_universes  

# --- Constructing geometry ---
f, c, w, materials, light_water = build_pin_universes()

# Export materials to XML for this run
openmc.Materials(materials).export_to_xml()

# Dimensions of the pincells / assembly
pinx = 1.4
piny = 1.4
pinz = 2              # cm per layer in z (pitch in z)
assembly_x = 17
assembly_y = 17
assembly_z = 150      # number of layers (for plotting, professor's value)

cx = assembly_x // 2  
cy = assembly_y // 2
ctrl_i = [cx-1, cx, cx+1]
ctrl_j = [cy-1, cy, cy+1]

# Making array for universes (Nz, Nx, Ny) as object array
universii = np.empty((assembly_z, assembly_x, assembly_y), dtype=object)

# Fill: control pins at center, otherwise fuel pins
for z in range(assembly_z):
    for i in range(assembly_x):
        for j in range(assembly_y):
            if (i in ctrl_i) and (j in ctrl_j):
                universii[z, i, j] = c   # Place controlpins
            else:
                universii[z, i, j] = f   # Place fuel pins

# Making the lattice 
f_assembly = openmc.RectLattice(name="fuel_assembly")
f_assembly.lower_left = (-assembly_x*abs(pinx)/2, -assembly_y*abs(piny)/2, -assembly_z*pinz/2)
f_assembly.pitch = (pinx, piny, pinz)
f_assembly.universes = universii.tolist()  

# Making outer universe 
water_cell = openmc.Cell(fill=light_water)
outer_uni = openmc.Universe(cells=[water_cell])
f_assembly.outer = outer_uni

# Putting the lattice into a simulation 
outedge = "reflective"  # Infinite reactor or singular assembly
outcap = openmc.ZPlane(assembly_z, boundary_type=outedge)
outbot = openmc.ZPlane(-assembly_z, boundary_type=outedge)
outcell_xy = openmc.model.rectangular_prism(width=23.8, height=23.8, boundary_type=outedge)
out_region = outcell_xy & +outbot & -outcap

f_cell = openmc.Cell(fill=f_assembly, region=out_region)  # making cell for assembly
assembly_uni = openmc.Universe(cells=[f_cell])            # Putting the assembly into a universe

geometry = openmc.Geometry()
geometry.root_universe = assembly_uni
geometry.export_to_xml()

if __name__ == "__main__":

    b4c_mat = next(m for m in materials if m.name == "B4C")

    plot = openmc.Plot()
    plot.filename = "assembly_xy"
    plot.basis = "xy"
    plot.width = (assembly_x*pinx, assembly_y*piny)
    plot.pixels = (900, 900)

    plot.color_by = "material"   

    DARK_BLUE = (11, 61, 145)
    WHITE     = (255, 255, 255)
    plot.colors = {
        light_water:DARK_BLUE,   
        b4c_mat:WHITE,   
        # next(m for m in materials if m.name=="Uranium dioxide"): "#7f8c00",
        # next(m for m in materials if m.name=="Zircaloy"):        "#808080",
        # next(m for m in materials if m.name=="Helium"):          "#8b4513",
    }

    plot.background = "white"

    plots = openmc.Plots([plot])
    plots.export_to_xml()
    openmc.plot_geometry()
    print("assembly.py: plot generated")
