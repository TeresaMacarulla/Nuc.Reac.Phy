# conda install -c conda-forge openmc-plotter
# python -m pip install --upgrade "git+https://github.com/openmc-dev/plotter.git@develop"
# openmc-plotter
import openmc
import numpy as np
import matplotlib.pyplot as plt
from utilities import build_assembly_universe

# ============================================================
# Existing objects
# ============================================================

# assembly_universe
# water
# stainless_steel
# vessel_steel
assembly_uni, materials, light_water = build_assembly_universe(with_ctrl_3x3=True)
steel = next(mat for mat in materials if mat.name == "Steel")
steel2 = next(mat for mat in materials if mat.name == "Steel2")

# ============================================================
# Water universe
# ============================================================

water_cell = openmc.Cell(fill=light_water)
water_universe = openmc.Universe(name="Water universe", cells=[water_cell])

# ============================================================
# Core lattice
# ============================================================

pin_pitch = 1.4
pins_per_assembly = 17
assembly_pitch = pin_pitch * pins_per_assembly   # 23.8 cm
N = 15
core_array = np.full((N, N), water_universe, dtype=openmc.Universe)
center = (N - 1) / 2

core_array[5:10, 1] = assembly_uni
core_array[4:11, 2] = assembly_uni
core_array[3:12, 3] = assembly_uni
core_array[2:13, 4] = assembly_uni
core_array[1:14, 5:10] = assembly_uni
core_array[2:13, 10] = assembly_uni
core_array[3:12, 11] = assembly_uni
core_array[4:11, 12] = assembly_uni
core_array[5:10, 13] = assembly_uni

core_lattice = openmc.RectLattice(name="Core lattice")
core_lattice.pitch = (assembly_pitch, assembly_pitch)
core_lattice.lower_left = (-N * assembly_pitch / 2, -N * assembly_pitch / 2)
core_lattice.universes = core_array
core_lattice.outer = water_universe

# ============================================================
# Radial reactor structures
# ============================================================

barrel_inner = openmc.ZCylinder(r=180.0)
barrel_outer = openmc.ZCylinder(r=190.0)
vessel_inner = openmc.ZCylinder(r=190.0)
vessel_outer = openmc.ZCylinder(r=220.0, boundary_type="vacuum")

# ============================================================
# Axial boundaries
# ============================================================

core_height = 300.0
z_bottom = openmc.ZPlane(z0=-core_height / 2, boundary_type="vacuum")
z_top = openmc.ZPlane(z0=core_height / 2, boundary_type="vacuum")
axial = +z_bottom & -z_top

# ============================================================
# Cells
# ============================================================

core_cell = openmc.Cell(name="Core and moderator", fill=core_lattice, region=-barrel_inner & axial)
barrel_cell = openmc.Cell(name="Core barrel", fill=steel, region=+barrel_inner & -barrel_outer & axial)
vessel_cell = openmc.Cell(name="Pressure vessel", fill=steel2, region=+vessel_inner & -vessel_outer & axial)

# ============================================================
# Root universe / geometry
# ============================================================

root = openmc.Universe(name="Reactor", cells=[core_cell, barrel_cell, vessel_cell])
geometry = openmc.Geometry(root)
geometry.export_to_xml()

openmc.run()