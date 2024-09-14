import numpy as np
from scipy.optimize import curve_fit

# Read data from file
data = np.loadtxt('EvsV')
volumes = data[:, 0]  # Total volume of the unit cell in Å³
energies = data[:, 1]  # Total energy of the unit cell in eV

# Birch-Murnaghan EOS
def birch_murnaghan(V, E0, V0, B0, B0_prime):
    eta = (V0/V)**(2/3)
    return E0 + (9/16) * V0 * B0 * ((eta-1)**3 * B0_prime + (eta-1)**2 * (6 - 4*eta))


# Initial guesses for E0, V0, B0, B0_prime
initial_guesses = [min(energies), volumes[np.argmin(energies)], 1, 4]


# Curve fitting to obtain E0, V0, B0, B0'
params, covariance = curve_fit(birch_murnaghan, volumes, energies, p0=initial_guesses)

# Extracting the fitted parameters
E0, V0, B0, B0_prime = params

# B0 is in units of energy per unit volume. Convert to GPa if necessary
# For example, if your energies are in eV and volumes in Å^3, convert B0 to GPa:
B0_GPa = B0 * 160.21766208

# Print the results
print(f"Equilibrium volume (V0): {V0} Å^3")
print(f"Minimum energy (E0): {E0} eV")
print(f"Bulk modulus (B0): {B0_GPa} GPa")
print(f"Pressure derivative of bulk modulus (B0'): {B0_prime}")
