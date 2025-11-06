import numpy as np
import glob
import matplotlib.pyplot as plt
import openmc
import math
from pathlib import Path
from utilities import build_core
from utilities import plot_core_xy
from utilities import core_flux_tallies
from utilities import plot_flux
from utilities import calc_fuel_volume_cm3
from utilities import set_uo2_volume
import openmc.deplete as deplete

if __name__ == "__main__":

    print(" ")
    print("-------------------------------------------------------------------------------")
    print(" ")
    print("1: Plot an xy plane of the core.")
    print("2: Run the simulation.")
    print("3: Run the simulation fixing the fuel temperature.")
    print("4: Plot neutron flux distribution.")
    print("5: Run the reactor for 6 months at a power level of 1200 MW thermal.")
    print("6: See some results from 6 month operation")

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
        print("See the core flux in: <<flux_XY.pdf>> inside this folder")

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

    elif number == 6:

        """
        # Read the results
        res = deplete.Results("depletion_results.h5")

        # get_eigenvalue() or get_keff() 
        try:
            times, keff, keff_sd = res.get_eigenvalue()          # time [s], k mean, k std
        except AttributeError:
            times, (keff, keff_sd) = res.get_keff(time_units='s')  # alternative

        # From seconds to days
        t_days = np.array(times, dtype=float)/86400.0

        # Plot
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





