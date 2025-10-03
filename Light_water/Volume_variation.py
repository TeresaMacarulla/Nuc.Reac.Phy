#!/usr/bin/env python3
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import openmc

# ----------------- User controls -----------------
# Pitch values (cm) for the rectangular prism (width=height=pitch)
PITCHES = [1.2, 1.3, 1.4, 1.6, 1.8, 2.0]

# Fuel / gap / clad radii (cm) — keep as in your base case
R_FUEL = 0.40
R_GAP  = 0.42
R_CLAD = 0.45

# Axial planes (cm) — same as your model: total height = 2.0 cm
ZMIN, ZMAX = -1.0, 1.0
HEIGHT = ZMAX - ZMIN

# Simulation settings (match your script)
BATCHES   = 100
INACTIVE  = 10
PARTICLES = 10000
THREADS   = 4

# Thermal cutoff (eV) for the "thermal utilization" tallies
E_TH = 1.0
# -------------------------------------------------


def get_material(mats, name):
    """Return a material by its .name from a Materials object."""
    for m in mats:
        if m.name == name:
            return m
    raise KeyError(f"Material named '{name}' not found in materials.xml")


def build_model_for_pitch(pitch, mats):
    """
    Build materials/geometry/settings/tallies for a given pitch.
    Returns (geometry, settings, tallies) and the fuel_cell object (for filters).
    """
    # --- Materials: reuse from materials.xml (same as your base script) ---
    urox       = get_material(mats, "Uranium dioxide")
    helium     = get_material(mats, "Helium")
    zircaloy   = get_material(mats, "Zircaloy")
    light_water= get_material(mats, "Light water")

    # --- Surfaces / Regions (same BCs: reflective on box and Z-planes) ---
    box  = openmc.model.RectangularPrism(width=pitch, height=pitch, boundary_type='reflective')
    zmin = openmc.ZPlane(ZMIN, boundary_type='reflective')
    zmax = openmc.ZPlane(ZMAX, boundary_type='reflective')

    c_fuel = openmc.ZCylinder(r=R_FUEL, boundary_type='transmission')
    c_gap  = openmc.ZCylinder(r=R_GAP,  boundary_type='transmission')
    c_clad = openmc.ZCylinder(r=R_CLAD, boundary_type='transmission')

    box_outside_cyl = (-box) & (+zmin) & (-zmax)
    fuel_region = (-c_fuel) & box_outside_cyl
    gap_region  = (+c_fuel) & (-c_gap)  & box_outside_cyl
    clad_region = (+c_gap)  & (-c_clad) & box_outside_cyl
    mod_region  = (+c_clad) & box_outside_cyl

    fuel_cell = openmc.Cell(name='fuel', fill=urox,     region=fuel_region)
    gap_cell  = openmc.Cell(name='gap',  fill=helium,   region=gap_region)
    clad_cell = openmc.Cell(name='clad', fill=zircaloy, region=clad_region)
    mod_cell  = openmc.Cell(name='water',fill=light_water, region=mod_region)

    univ = openmc.Universe(cells=[fuel_cell, gap_cell, clad_cell, mod_cell])
    geometry = openmc.Geometry(univ)

    # --- Settings (match your base script, but source box matches new pitch) ---
    settings = openmc.Settings()
    settings.batches   = BATCHES
    settings.inactive  = INACTIVE
    settings.particles = PARTICLES

    bounds = [-pitch/2, -pitch/2, ZMIN, pitch/2, pitch/2, ZMAX]
    src = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    settings.source = openmc.IndependentSource(space=src)

    # --- Tallies for thermal utilization f = <Σaφ>_T^F / <Σaφ>_T ---
    thermal_E = openmc.EnergyFilter([0.0, E_TH])
    fuel_only = openmc.CellFilter([fuel_cell])

    t_num = openmc.Tally(name="therm_abs_in_fuel")
    t_num.filters = [fuel_only, thermal_E]
    t_num.scores  = ["absorption"]

    t_den = openmc.Tally(name="therm_abs_total")
    t_den.filters = [thermal_E]
    t_den.scores  = ["absorption"]

    tallies = openmc.Tallies([t_num, t_den])

    return geometry, settings, tallies, fuel_cell


def run_one_pitch(pitch, mats):
    """Build, export, run, compute f, and return (volume_cm3, f)."""
    geometry, settings, tallies, _fuel_cell = build_model_for_pitch(pitch, mats)

    # Export current model to XML (overwrite each iteration)
    mats.export_to_xml()
    geometry.export_to_xml()
    tallies.export_to_xml()
    settings.export_to_xml()

    print(f"\n--- Running OpenMC for pitch = {pitch:.3f} cm ---")
    openmc.run(threads=THREADS)

    # Read the latest statepoint and compute f
    sp_path = sorted(glob.glob("statepoint.*.h5"))[-1]
    sp = openmc.StatePoint(sp_path)

    num = sp.get_tally(name="therm_abs_in_fuel").mean.sum()  # <Σaφ>_T^F
    den = sp.get_tally(name="therm_abs_total").mean.sum()    # <Σaφ>_T
    f   = float(num/den) if den > 0 else float("nan")

    # Pin-cell volume = pitch * pitch * HEIGHT
    volume = pitch * pitch * HEIGHT

    print(f"  <Σaφ>_T^F = {num:.6e},  <Σaφ>_T = {den:.6e}  ->  f = {f:.6f}")
    print(f"  Pin-cell volume = {volume:.6f} cm^3")
    return volume, f


def main():
    # Load materials from the XML produced by your base script
    if not os.path.exists("materials.xml"):
        raise RuntimeError("materials.xml not found. Run your base script once to generate XMLs.")

    mats = openmc.Materials.from_xml("materials.xml")

    vols, fs = [], []
    for pitch in PITCHES:
        V, f = run_one_pitch(pitch, mats)
        vols.append(V); fs.append(f)

    # Plot f vs pin-cell volume and export to PDF
    plt.figure()
    plt.plot(vols, fs, marker='o', linewidth=1.5)
    plt.xlabel(r"Pin-cell volume $V$ [cm$^3$]")
    plt.ylabel(r"Thermal utilization $f$ (E $\leq$ 1 eV)")
    plt.title("Thermal utilization vs pin-cell volume (reflective BC)")
    plt.grid(True, which="both", alpha=0.35)
    plt.tight_layout()
    plt.savefig("f_vs_pincell_volume.pdf")
    plt.show()

    # Also dump a small CSV for convenience
    with open("f_vs_volume.csv", "w") as fh:
        fh.write("volume_cm3,f\n")
        for V, f in zip(vols, fs):
            fh.write(f"{V:.6f},{f:.8f}\n")
    print("\nSaved: f_vs_pincell_volume.pdf and f_vs_volume.csv")

if __name__ == "__main__":
    main()

