import openmc
import numpy as np
import matplotlib.pyplot as plt

# EXERCISE A---------------------------------------------------------------------------

#  UO2
urox = openmc.Material(name = "Uranium dioxide")
urox.add_nuclide("U235", 0.07)
urox.add_nuclide("U238", 0.93)
urox.add_element("O", 2.0)
urox.set_density("g/cm3",10.97)
# Enhanced scattering cross sections thermal regions
urox.add_s_alpha_beta("c_O_in_UO2")

# Heavy water (D2O)
heavy_water = openmc.Material(name="Heavy water")
heavy_water.add_nuclide("H2", 2.0)      
heavy_water.add_element("O", 1.0)
heavy_water.set_density("g/cm3", 1.105) 
heavy_water.add_s_alpha_beta("c_D_in_D2O")

# Graphite 
graphite = openmc.Material(name="Graphite")
graphite.add_element("C", 1.0)
graphite.set_density("g/cm3", 1.70)    
graphite.add_s_alpha_beta("c_Graphite")      


# Helium 
helium = openmc.Material(name="Helium")
helium.add_nuclide("He4", 1.0)
helium.set_density("g/cm3", 0.0001786)  

# Zircaloy (cladding)
zircaloy = openmc.Material(name="Zircaloy")
zircaloy.add_element("Zr", 1.0)
zircaloy.set_density("g/cm3", 6.52)  

materials = openmc.Materials([urox, helium, zircaloy, graphite])
materials.export_to_xml()

# EXERCISE B-----------------------------------------------------------------------------------------------

# Box of size 1.4 cm x 1.4 cm x 2 cm with “vacuum” or "reflective" boundary condition
box = openmc.model.RectangularPrism(width=2, height=2, boundary_type='reflective')
zmin = openmc.ZPlane(-1.0, boundary_type='reflective')
zmax = openmc.ZPlane(+1.0, boundary_type='reflective')

# Concentric cylinders: fuel, cavity, coating
c_fuel = openmc.ZCylinder(r=0.40, boundary_type='transmission') # 4.0 mm
c_gap = openmc.ZCylinder(r=0.42, boundary_type='transmission') # 4.2 mm
c_clad = openmc.ZCylinder(r=0.43, boundary_type='transmission') # 4.5 mm

box_outside_cylinder = (-box)&(+zmin)&(-zmax)#&(~c_fuel)
volume_c_fuel=(-c_fuel)&(+zmin)&(-zmax)
volume_c_gap=(-c_gap)&(+zmin)&(-zmax)
volume_c_clad=(-c_clad)&(+zmin)&(-zmax)

fuel_region = (-c_fuel)&box_outside_cylinder # inside the fuel
gap_region = (+c_fuel)&(-c_gap)&box_outside_cylinder # between fuel and gap
clad_region = (+c_gap)&(-c_clad)&box_outside_cylinder # between gap and clad
water_region = (+c_clad)&box_outside_cylinder # outside clad, inside the box

fuel_cell = openmc.Cell(name='fuel',fill=urox,region=fuel_region)
gap_cell = openmc.Cell(name='gap',fill=helium,region=gap_region)
clad_cell = openmc.Cell(name='clad',fill=zircaloy,region=clad_region)
mod_cell = openmc.Cell(name='water',fill=graphite,region=water_region)

# Universe and Geometry
root = openmc.Universe(cells=[fuel_cell, gap_cell, clad_cell, mod_cell])
geom = openmc.Geometry(root)
geom.export_to_xml()

# Plots
colors = {urox: 'green', helium: 'brown', zircaloy: 'gray', graphite: 'white'}

# Plot XY
plot_xy           = openmc.Plot()
plot_xy.filename  = 'pin_xy'
plot_xy.basis     = 'xy'
plot_xy.width     = (1.4, 1.4)
plot_xy.pixels    = (800, 800)
plot_xy.color_by  = 'material'
plot_xy.colors    = colors

# Plot XZ
plot_xz           = openmc.Plot()
plot_xz.filename  = 'pin_xz'
plot_xz.basis     = 'xz'
plot_xz.width     = (1.4, 2.0)    # width x, height z
plot_xz.pixels    = (800, 800)
plot_xz.color_by  = 'material'
plot_xz.colors    = colors

plots = openmc.Plots([plot_xy, plot_xz])
plots.export_to_xml()

# Plot generation
openmc.plot_geometry()

# EXERCISE C-----------------------------------------------------------------------------------------------

settings = openmc.Settings()
settings.batches = 100 # How many iterations of neutron generations
settings.inactive = 10 # How many generations used to determine sample spots
settings.particles = 10000 # Number of neutrons in each generation

materials = openmc.Materials([graphite, urox, helium, zircaloy])
materials.export_to_xml()

# The following section tells OpenMC where to sample particle sites.
bounds = [-0.7,-0.7,-1,0.7,0.7,1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],only_fissionable=True)
settings.source = openmc.IndependentSource(space=uniform_dist)
settings.export_to_xml()

if __name__ == "__main__":
 openmc.run(threads = 4) # Threads set for heplab primarily

# EXERCISE E-----------------------------------------------------------------------------------------------

mesh = openmc.RegularMesh() #Making a mesh instance
mesh.dimension = [100,100] #Setting the bins of the mesh to 100 x 100 bins
mesh.lower_left= [-0.7,-0.7] #Setting lower left coordinates
mesh.upper_right= [0.7,0.7] #Setting upper right coordinates

mesh_filter = openmc.MeshFilter(mesh) #Making a meshfilter
tallies_file = openmc.Tallies() #For our tallies to be included we have to put them into an input file
 
tally = openmc.Tally(name = "Tally") #Making a tally instance
tally.filters = [mesh_filter] #Designatingen the mesh for the tally
tally.scores = ["flux", 'prompt-nu-fission'] #Designating the objects to be tallied
tallies_file.append(tally) #Appending the tally to our tallies file

tallies_file.export_to_xml()
openmc.run()

# EXERCISE D-----------------------------------------------------------------------------------------------

# Define an energy filter (example: from 1e-5 eV to 20 MeV, 100 bins log-spaced)
energy_bins = np.logspace(-5, 7, 200)   # Adjust the number of bins 
energy_filter = openmc.EnergyFilter(energy_bins)

tally2 = openmc.Tally(name = "Energy") #Making a new tally only looking at flux
tally2.filters = [energy_filter] #Designating the energy filter
tally2.scores = ["flux"]
tallies_file.append(tally2)

tallies_file.export_to_xml()
openmc.run()

# K values with Heavy water
#k-effective (Collision)     = 0.86895 +/- 0.00137
#k-effective (Track-length)  = 0.86909 +/- 0.00150
#k-effective (Absorption)    = 0.87020 +/- 0.00110
#Combined k-effective        = 0.86985 +/- 0.00106
#Leakage Fraction            = 0.00000 +/- 0.00000

# K values with Heavy water and gemotery changes fuel 0.38, gap 0.44, clad 0.47
#k-effective (Collision)     = 0.88454 +/- 0.00147
#k-effective (Track-length)  = 0.88347 +/- 0.00157
#k-effective (Absorption)    = 0.88561 +/- 0.00120
#Combined k-effective        = 0.88536 +/- 0.00118
#Leakage Fraction            = 0.00000 +/- 0.00000

# K values with Graphite and original geometry
#k-effective (Collision) = 0.52593 +/- 0.00070 
#k-effective (Track-length) = 0.52587 +/- 0.00075 
#k-effective (Absorption) = 0.52501 +/- 0.00089 
#Combined k-effective = 0.52558 +/- 0.00063 
#Leakage Fraction = 0.00000 +/- 0.00000

# K values with Graphite, 4% enrichemnt, fuel 0.40, gap 0.42, clad 0.44
#k-effective (Collision)     = 0.71252 +/- 0.00086
#k-effective (Track-length)  = 0.71271 +/- 0.00095
#k-effective (Absorption)    = 0.71242 +/- 0.00091
#Combined k-effective        = 0.71252 +/- 0.00069
#Leakage Fraction            = 0.00000 +/- 0.00000

# K values with Graphite, 6% enrichemnt, fuel 0.40, gap 0.42, clad 0.43, external box 1.6x1.6
#k-effective (Collision)     = 0.86151 +/- 0.00100
#k-effective (Track-length)  = 0.86098 +/- 0.00102
#k-effective (Absorption)    = 0.86430 +/- 0.00101
#Combined k-effective        = 0.86274 +/- 0.00088
#Leakage Fraction            = 0.00000 +/- 0.00000

# K values with Graphite, 6% enrichemnt, fuel 0.40, gap 0.42, clad 0.43, external box 2x2
#k-effective (Collision)     = 0.93768 +/- 0.00114
#k-effective (Track-length)  = 0.93835 +/- 0.00142
#k-effective (Absorption)    = 0.93789 +/- 0.00093
#Combined k-effective        = 0.93773 +/- 0.00081
#Leakage Fraction            = 0.00000 +/- 0.00000

# K values with Graphite, 7% enrichemnt, fuel 0.40, gap 0.42, clad 0.43, external box 2x2
#k-effective (Collision)     = 0.97525 +/- 0.00127
#k-effective (Track-length)  = 0.97510 +/- 0.00147
#k-effective (Absorption)    = 0.97484 +/- 0.00095
#Combined k-effective        = 0.97496 +/- 0.00090
#Leakage Fraction            = 0.00000 +/- 0.00000