import openmc
import numpy as np
import matplotlib.pyplot as plt

urox = openmc.Material(name = "Uranium dioxide")
#This makes an instance of material that we may modify at will
urox.add_nuclide("U235", 0.02)
#This above adds uranium 235 and sets the relative amount of atoms with the float.
urox.add_nuclide("U238", 0.98)
urox.add_element("O", 2.0)
urox.set_density("g/cm3",10.97)
#Theoretical density of UO2 at room temperature, change when needed
urox.add_s_alpha_beta("c_O_in_UO2")
#Enhanced scattering cross sections thermal regions