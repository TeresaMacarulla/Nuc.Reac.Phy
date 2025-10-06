import openmc
import numpy as np
import matplotlib.pyplot as plt

# EXERCISE A---------------------------------------------------------------------------

#  UO2
urox = openmc.Material(name = "Uranium dioxide")
urox.add_nuclide("U235", 0.10)
urox.add_nuclide("U238", 0.90)
urox.add_element("O", 2.0)
urox.set_density("g/cm3",10.97)
# Enhanced scattering cross sections thermal regions
urox.add_s_alpha_beta("c_O_in_UO2")

#PuO2
puox = openmc.Material(name="Plutonium dioxide")
puox.add_nuclide("Pu238", 0.02)  # ~ 2%
puox.add_nuclide("Pu239", 0.50)  # ~56%
puox.add_nuclide("Pu240", 0.40)  # ~24%
puox.add_nuclide("Pu241", 0.02)  # ~12%
puox.add_nuclide("Pu242", 0.06)  # ~ 6%
puox.add_element("O", 2.0)
puox.set_density("g/cm3", 11.45)

# Heavy water (D2O)
heavy_water = openmc.Material(name="Heavy water")
heavy_water.add_nuclide("H2", 2.0)      
heavy_water.add_element("O", 1.0)
heavy_water.set_density("g/cm3", 1.105) 
heavy_water.add_s_alpha_beta("c_D_in_D2O")

# Light water (H2O)
light_water = openmc.Material(name="Light water")
light_water.add_nuclide("H1", 2.0)
light_water.add_element("O", 1.0)
light_water.set_density("g/cm3", 1.0)
# Enhanced scattering cross sections thermal regions
light_water.add_s_alpha_beta("c_H_in_H2O")

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

materials = openmc.Materials([puox, helium, zircaloy, light_water])
materials.export_to_xml()

# EXERCISE B-----------------------------------------------------------------------------------------------

# Box of size 1.4 cm x 1.4 cm x 2 cm with “vacuum” or "reflective" boundary condition
box = openmc.model.RectangularPrism(width=1.4, height=1.4, boundary_type='reflective')
zmin = openmc.ZPlane(-1.0, boundary_type='reflective')
zmax = openmc.ZPlane(+1.0, boundary_type='reflective')

# Concentric cylinders: fuel, cavity, coating
c_fuel = openmc.ZCylinder(r=0.40, boundary_type='transmission') # 4.0 mm
c_gap = openmc.ZCylinder(r=0.42, boundary_type='transmission') # 4.2 mm
c_clad = openmc.ZCylinder(r=0.45, boundary_type='transmission') # 4.5 mm

box_outside_cylinder = (-box)&(+zmin)&(-zmax)#&(~c_fuel)
volume_c_fuel=(-c_fuel)&(+zmin)&(-zmax)
volume_c_gap=(-c_gap)&(+zmin)&(-zmax)
volume_c_clad=(-c_clad)&(+zmin)&(-zmax)

fuel_region = (-c_fuel)&box_outside_cylinder # inside the fuel
gap_region = (+c_fuel)&(-c_gap)&box_outside_cylinder # between fuel and gap
clad_region = (+c_gap)&(-c_clad)&box_outside_cylinder # between gap and clad
water_region = (+c_clad)&box_outside_cylinder # outside clad, inside the box

fuel_cell = openmc.Cell(name='fuel',fill=puox,region=fuel_region)
gap_cell = openmc.Cell(name='gap',fill=helium,region=gap_region)
clad_cell = openmc.Cell(name='clad',fill=zircaloy,region=clad_region)
mod_cell = openmc.Cell(name='water',fill=light_water,region=water_region)

# Universe and Geometry
root = openmc.Universe(cells=[fuel_cell, gap_cell, clad_cell, mod_cell])
geom = openmc.Geometry(root)
geom.export_to_xml()

# Plots
#colors = {urox: 'green', helium: 'brown', zircaloy: 'gray', light_water: 'white'}

# Plot XY
#plot_xy           = openmc.Plot()
#plot_xy.filename  = 'pin_xy'
#plot_xy.basis     = 'xy'
#plot_xy.width     = (1.4, 1.4)
#plot_xy.pixels    = (800, 800)
#plot_xy.color_by  = 'material'
#plot_xy.colors    = colors

# Plot XZ
#plot_xz           = openmc.Plot()
#plot_xz.filename  = 'pin_xz'
#plot_xz.basis     = 'xz'
#plot_xz.width     = (1.4, 2.0)    # width x, height z
#plot_xz.pixels    = (800, 800)
#plot_xz.color_by  = 'material'
#plot_xz.colors    = colors

#plots = openmc.Plots([plot_xy, plot_xz])
#plots.export_to_xml()

# Plot generation
#openmc.plot_geometry()

# EXERCISE C-----------------------------------------------------------------------------------------------

settings = openmc.Settings()
settings.batches = 100 # How many iterations of neutron generations
settings.inactive = 10 # How many generations used to determine sample spots
settings.particles = 10000 # Number of neutrons in each generation

materials = openmc.Materials([light_water, puox, helium, zircaloy])
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

# K values UO2 1%
# k-effective (Collision)     = 1.02685 +/- 0.00111
# k-effective (Track-length)  = 1.02677 +/- 0.00121
# k-effective (Absorption)    = 1.02586 +/- 0.00121
# Combined k-effective        = 1.02636 +/- 0.00084
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values UO2 2%
# k-effective (Collision)     = 1.29182 +/- 0.00128
# k-effective (Track-length)  = 1.29260 +/- 0.00151
# k-effective (Absorption)    = 1.29092 +/- 0.00107
# Combined k-effective        = 1.29147 +/- 0.00086
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values UO2 3%
# k-effective (Collision)     = 1.41906 +/- 0.00132
# k-effective (Track-length)  = 1.41967 +/- 0.00178
# k-effective (Absorption)    = 1.41642 +/- 0.00099
# Combined k-effective        = 1.41709 +/- 0.00098
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values UO2 4%
# k-effective (Collision)     = 1.48799 +/- 0.00142
# k-effective (Track-length)  = 1.48651 +/- 0.00187
# k-effective (Absorption)    = 1.48883 +/- 0.00099
# Combined k-effective        = 1.48860 +/- 0.00094
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values UO2 5%
# k-effective (Collision)     = 1.53451 +/- 0.00126
# k-effective (Track-length)  = 1.53325 +/- 0.00190
# k-effective (Absorption)    = 1.53510 +/- 0.00095
# Combined k-effective        = 1.53508 +/- 0.00085
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values UO2 6%
# k-effective (Collision)     = 1.56792 +/- 0.00138
# k-effective (Track-length)  = 1.56936 +/- 0.00197
# k-effective (Absorption)    = 1.56693 +/- 0.00106
# Combined k-effective        = 1.56730 +/- 0.00101
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values UO2 7%
# k-effective (Collision)     = 1.59278 +/- 0.00146
# k-effective (Track-length)  = 1.59645 +/- 0.00182
# k-effective (Absorption)    = 1.59233 +/- 0.00090
# Combined k-effective        = 1.59281 +/- 0.00094
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values UO2 8%
# k-effective (Collision)     = 1.60879 +/- 0.00141
# k-effective (Track-length)  = 1.60783 +/- 0.00211
# k-effective (Absorption)    = 1.61083 +/- 0.00085
# Combined k-effective        = 1.61049 +/- 0.00083
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values UO2 9%
# k-effective (Collision)     = 1.62587 +/- 0.00133
# k-effective (Track-length)  = 1.62672 +/- 0.00190
# k-effective (Absorption)    = 1.62552 +/- 0.00087
# Combined k-effective        = 1.62556 +/- 0.00084
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values UO2 10%
# k-effective (Collision)     = 1.63617 +/- 0.00155
# k-effective (Track-length)  = 1.63531 +/- 0.00215
# k-effective (Absorption)    = 1.63698 +/- 0.00095
# Combined k-effective        = 1.63679 +/- 0.00093
# Leakage Fraction            = 0.00000 +/- 0.00000



# K values PuO2 5% ("Pu238", 0.02) ("Pu239", 0.05) ("Pu240", 0.85) ("Pu241", 0.02) ("Pu242", 0.06)
# k-effective (Collision)     = 0.84400 +/- 0.00075
# k-effective (Track-length)  = 0.84429 +/- 0.00081
# k-effective (Absorption)    = 0.84259 +/- 0.00110
# Combined k-effective        = 0.84371 +/- 0.00058
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values PuO2 10%  ("Pu238", 0.02) ("Pu239", 0.10) ("Pu240", 0.80) ("Pu241", 0.02) ("Pu242", 0.06) 
# k-effective (Collision)     = 1.00066 +/- 0.00083
# k-effective (Track-length)  = 1.00229 +/- 0.00094
# k-effective (Absorption)    = 1.00274 +/- 0.00122
# Combined k-effective        = 1.00178 +/- 0.00068
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values PuO2 15% ("Pu238", 0.02) ("Pu239", 0.15) ("Pu240", 0.75) ("Pu241", 0.02) ("Pu242", 0.06)
# k-effective (Collision)     = 1.12026 +/- 0.00077
# k-effective (Track-length)  = 1.12024 +/- 0.00116
# k-effective (Absorption)    = 1.11888 +/- 0.00148
# Combined k-effective        = 1.12001 +/- 0.00073
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values PuO2 20% ("Pu238", 0.02) ("Pu239", 0.20) ("Pu240", 0.70) ("Pu241", 0.02) ("Pu242", 0.06)
# k-effective (Collision)     = 1.21523 +/- 0.00093
# k-effective (Track-length)  = 1.21500 +/- 0.00113
# k-effective (Absorption)    = 1.21445 +/- 0.00124
# Combined k-effective        = 1.21490 +/- 0.00074
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values PuO2 25% ("Pu238", 0.02) ("Pu239", 0.25) ("Pu240", 0.65) ("Pu241", 0.02) ("Pu242", 0.06)
# k-effective (Collision)     = 1.29242 +/- 0.00098
# k-effective (Track-length)  = 1.29247 +/- 0.00143
# k-effective (Absorption)    = 1.29498 +/- 0.00121
# Combined k-effective        = 1.29340 +/- 0.00082
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values PuO2 30% ("Pu238", 0.02) ("Pu239", 0.30) ("Pu240", 0.60) ("Pu241", 0.02) ("Pu242", 0.06)
# k-effective (Collision)     = 1.36191 +/- 0.00100
# k-effective (Track-length)  = 1.36143 +/- 0.00140
# k-effective (Absorption)    = 1.36296 +/- 0.00121
# Combined k-effective        = 1.36223 +/- 0.00079
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values PuO2 35% ("Pu238", 0.02) ("Pu239", 0.35) ("Pu240", 0.55) ("Pu241", 0.02) ("Pu242", 0.06)
# k-effective (Collision)     = 1.42192 +/- 0.00100
# k-effective (Track-length)  = 1.42137 +/- 0.00156
# k-effective (Absorption)    = 1.42358 +/- 0.00119
# Combined k-effective        = 1.42255 +/- 0.00082
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values PuO2 40% ("Pu238", 0.02) ("Pu239", 0.40) ("Pu240", 0.50) ("Pu241", 0.02) ("Pu242", 0.06)
# k-effective (Collision)     = 1.47407 +/- 0.00123
# k-effective (Track-length)  = 1.47536 +/- 0.00168
# k-effective (Absorption)    = 1.47610 +/- 0.00124
# Combined k-effective        = 1.47519 +/- 0.00089
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values PuO2 45% ("Pu238", 0.02) ("Pu239", 0.45) ("Pu240", 0.45) ("Pu241", 0.02) ("Pu242", 0.06)
# k-effective (Collision)     = 1.52553 +/- 0.00120
# k-effective (Track-length)  = 1.52360 +/- 0.00165
# k-effective (Absorption)    = 1.52653 +/- 0.00102
# Combined k-effective        = 1.52598 +/- 0.00089
# Leakage Fraction            = 0.00000 +/- 0.00000

# K values PuO2 50% ("Pu238", 0.02) ("Pu239", 0.50) ("Pu240", 0.40) ("Pu241", 0.02) ("Pu242", 0.06)
# k-effective (Collision)     = 1.57334 +/- 0.00103
# k-effective (Track-length)  = 1.57367 +/- 0.00163
# k-effective (Absorption)    = 1.57386 +/- 0.00098
# Combined k-effective        = 1.57364 +/- 0.00071
# Leakage Fraction            = 0.00000 +/- 0.00000


