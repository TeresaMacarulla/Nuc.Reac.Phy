# pins.py
import openmc

def build_pin_universes():
    # ----------------- Materials -----------------
    # UO2 
    urox = openmc.Material(name="Uranium dioxide")
    urox.add_nuclide("U235", 0.012)
    urox.add_nuclide("U238", 0.988)
    urox.add_element("O", 2.0)
    urox.set_density("g/cm3", 10.97)
    urox.add_s_alpha_beta("c_O_in_UO2")

    # Light water (H2O)
    light_water = openmc.Material(name="Light water")
    light_water.add_nuclide("H1", 2.0)
    light_water.add_element("O", 1.0)
    light_water.set_density("g/cm3", 1.0)
    light_water.add_s_alpha_beta("c_H_in_H2O")

    # Helium gap
    helium = openmc.Material(name="Helium")
    helium.add_nuclide("He4", 1.0)
    helium.set_density("g/cm3", 0.0001786)

    # Zircaloy cladding
    zircaloy = openmc.Material(name="Zircaloy")
    zircaloy.add_element("Zr", 1.0)
    zircaloy.set_density("g/cm3", 6.52)

    # --- B4C (control rod) as a single material (natural boron) ---
    b4c = openmc.Material(name="B4C")
    # natural boron: ~19.9% B-10, 80.1% B-11
    b4c.add_nuclide("B10", 0.199 * 4.0)   # 4 boron atoms per formula unit
    b4c.add_nuclide("B11", 0.801 * 4.0)
    b4c.add_element("C", 1.0)             # 1 carbon atom per formula unit
    b4c.set_density("g/cm3", 2.52)


    # Keep a list to export later from assembly.py
    materials = [urox, light_water, helium, zircaloy, b4c]

    # ----------------- Pin geometry (Part 1) -----------------
    # Box 1.4x1.4 cm in XY, height 2 cm is handled by assembly; here just local z planes
    box = openmc.model.RectangularPrism(width=1.4, height=1.4, boundary_type="transmission")
    zmin = openmc.ZPlane(-1.0, boundary_type="transmission")
    zmax = openmc.ZPlane(+1.0, boundary_type="transmission")

    # Concentric cylinders (r in cm)
    c_fuel = openmc.ZCylinder(r=0.40, boundary_type="transmission")
    c_gap  = openmc.ZCylinder(r=0.42, boundary_type="transmission")
    c_clad = openmc.ZCylinder(r=0.45, boundary_type="transmission")

    # Regions
    box_outside_cyl = (-box) & (+zmin) & (-zmax)
    fuel_region = (-c_fuel) & box_outside_cyl
    gap_region  = (+c_fuel) & (-c_gap) & box_outside_cyl
    clad_region = (+c_gap)  & (-c_clad) & box_outside_cyl
    water_region= (+c_clad) & box_outside_cyl

    # -------- Universe: fuel pin (UO2) --------
    fuel_cell = openmc.Cell(name="fuel", fill=urox, region=fuel_region)
    gap_cell  = openmc.Cell(name="gap",  fill=helium, region=gap_region)
    clad_cell = openmc.Cell(name="clad", fill=zircaloy, region=clad_region)
    mod_cell  = openmc.Cell(name="mod",  fill=light_water, region=water_region)
    pin_universe = openmc.Universe(cells=[fuel_cell, gap_cell, clad_cell, mod_cell])

    # -------- Universe: control-rod pin (B4C in the fuel region) --------
    fuel_cell_c = openmc.Cell(name="fuel_ctrlpin", fill=b4c,       region=fuel_region)
    gap_cell_c  = openmc.Cell(name="gap_ctrlpin",  fill=helium,    region=gap_region)
    clad_cell_c = openmc.Cell(name="clad_ctrlpin", fill=zircaloy,  region=clad_region)
    mod_cell_c  = openmc.Cell(name="mod_ctrlpin",  fill=light_water, region=water_region)
    controlpin_universe = openmc.Universe(cells=[fuel_cell_c, gap_cell_c, clad_cell_c, mod_cell_c])

    # -------- Universe: "water pin" (fuel region replaced by water) --------
    water_fuel_cell = openmc.Cell(name="fuel_as_water", fill=light_water, region=fuel_region)
    water_universe = openmc.Universe(cells=[water_fuel_cell, gap_cell, clad_cell, mod_cell])

    return pin_universe, controlpin_universe, water_universe, materials, light_water

if __name__ == "__main__":
    f, c, w, mats, lw = build_pin_universes()
    openmc.Materials(mats).export_to_xml()
    print("pins.py: pin universes ready.")

