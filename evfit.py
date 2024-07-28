# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

# Define the fourth-order Birch-Murnaghan equation of state for energy
def birch_murnaghan_eos(V, E0, V0, B0, B0_prime, B0_double_prime):
    term1 = (V0 / V)**(2/3) - 1
    term2 = term1**2
    term3 = term1**3
    return (E0 + (9 * V0 * B0 / 16) * 
            (term3 * B0_prime + term2 * (6 - 4 * (V0 / V)**(2/3)) + 
            0.5 * term3 * (B0_double_prime - B0_prime + (4/3) * B0_prime * (B0_prime - 4))))

# Objective function for optimization
def objective_function(params, volume, energy):
    E0, V0, B0, B0_prime, B0_double_prime = params
    
    # Ensure V0 is positive
    if V0 <= 0:
        return np.inf
    
    # Calculate energy using Birch-Murnaghan EOS
    energy_calc = birch_murnaghan_eos(volume, E0, V0, B0, B0_prime, B0_double_prime)
    
    # Calculate the sum of squared residuals
    residuals = energy - energy_calc
    return np.sum(residuals**2)

# Load data from a text file
def load_data_output(filename):
    try:
        data = np.loadtxt(filename, skiprows=1)  # Skip the header row
        volume = data[:, 0]  # First column is volume
        energy = data[:, 1]   # Second column is energy
        return volume, energy  # Return volume and energy
    except Exception as e:
        print(f"Error loading data: {e}")
        return None, None


# Load data from a text file
def load_data(filename):
    try:
        data = np.loadtxt(filename)  # Load data without skipping rows
        volume = data[:, 0]  # First column is volume
        energy = data[:, 1]   # Second column is energy
        return volume, energy  # Return volume and energy
    except Exception as e:
        print(f"Error loading data: {e}")
        return None, None

# Read the number of atoms from the POSCAR file
def read_number_of_atoms(poscar_file):
    try:
        with open(poscar_file, 'r') as f:
            lines = f.readlines()
            # The number of atoms is in line 7 (index 6)
            atom_counts = lines[6].strip().split()
            total_atoms = sum(int(count) for count in atom_counts)
            return total_atoms
    except Exception as e:
        print(f"Error reading POSCAR file: {e}")
        return None

# Fit the data
def fit_data(volume, energy):
    # Filter out non-positive volume values
    mask = volume > 0
    volume_filtered = volume[mask]
    energy_filtered = energy[mask]
    
    if len(volume_filtered) == 0 or len(energy_filtered) == 0:
        raise ValueError("No valid volume and energy data available for fitting.")
    
    # Initial guess for parameters: E0, V0, B0, B0_prime, B0_double_prime
    initial_guess = [min(energy_filtered), np.mean(volume_filtered), 1.0, 4.0, 0.0]
    
    # Perform optimization using scipy.optimize.minimize
    res = opt.minimize(objective_function, initial_guess, args=(volume_filtered, energy_filtered),
                       bounds=[(None, None), (1e-6, None), (1e-6, None), (1e-6, None), (None, None)])
    
    return res.x

# Write original data to a separate output file
def write_original_data_output(filename, volume, energy, total_atoms):
    with open(filename, 'w') as f:
        f.write("Volume\tEnergy\n")
        for v, e in zip(volume * total_atoms, energy):
            f.write(f"{v:.6f}\t{e:.6f}\n")

# Write fitted parameters to a separate output file
def write_fit_parameters_output(filename, params, total_atoms):
    with open(filename, 'w') as f:
        f.write("Parameter\tValue\tUnit\n")
        f.write(f"E0\t{params[0]:.6f}\teV\n")
        f.write(f"V0\t{params[1] * total_atoms:.6f}\t\u00C5^3\n")  # Multiply V0 by total_atoms
        f.write(f"B0\t{params[2]:.6f}\tGPa\n")
        f.write(f"B0_prime\t{params[3]:.6f}\t-\n")
        f.write(f"B0_double_prime\t{params[4]:.6f}\tGPa^-1\n")


# Write fitted curve data to a separate output file
def write_fit_curve_output(filename, volume_fit, energy_fit, total_atoms):
    with open(filename, 'w') as f:
        f.write("Volume_Fit\tEnergy_Fit\n")
        for v, e in zip(volume_fit * total_atoms, energy_fit):
            f.write(f"{v:.6f}\t{e:.6f}\n")
            
# Function to plot the results from output files
def plot_results(original_data_file, fit_curve_file, total_atoms, params):
    # Load original data
    volume_orig, energy_orig = load_data_output(original_data_file)
    
    # Load fitted curve data
    volume_fit, energy_fit = load_data_output(fit_curve_file)
    
    plt.figure(figsize=(10, 6))
    plt.scatter(volume_orig * total_atoms, energy_orig, label='Original Data', color='blue', marker='o')
    plt.plot(volume_fit * total_atoms, energy_fit, label='Fitted Curve', color='red')
    plt.xlabel('Volume (Å³)')
    plt.ylabel('Energy (eV)')
    plt.title('Fitted Birch-Murnaghan Equation of State')
    plt.legend()
    plt.grid()

    # Prepare text for overlay
    textstr = '\n'.join((
        f"E0: {params[0]:.6f} eV",
        f"V0: {params[1] * total_atoms:.6f} Å³",
        f"B0: {params[2]:.6f} GPa",
        f"B0': {params[3]:.6f}",
        f"B0'': {params[4]:.6f} GPa^-1"
    ))

    # Positioning the text box in the center of the plot
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    plt.text(0.5, 0.5, textstr, transform=plt.gca().transAxes, fontsize=12,
             verticalalignment='center', horizontalalignment='center', bbox=props)

    plt.savefig('fit_plot.png', dpi=500)  # Save the plot as a PNG file
    #plt.show()  # Display the plot

# Main function
def main():
    input_file = 'EvsV'  # Your data file
    poscar_file = 'POSCAR'  # Your POSCAR file
    original_data_output_file = 'simulation_data.txt'
    fit_parameters_output_file = 'fit_parameters.txt'
    fit_curve_output_file = 'fit_curve_data.txt'
    
    # Load data
    volume, energy = load_data(input_file)
    
    if volume is None or energy is None:
        print("Failed to load data. Exiting.")
        return
    
    # Read the number of atoms from the POSCAR file
    total_atoms = read_number_of_atoms(poscar_file)
    if total_atoms is None:
        print("Failed to read number of atoms. Exiting.")
        return
    
    try:
        # Fit the data
        params = fit_data(volume, energy)
        
        # Generate fitted curve data points with a wider range
        volume_fit = np.linspace(0.6 * min(volume), 1.4 * max(volume), 5000)  # Wider range for extrapolation
        energy_fit = birch_murnaghan_eos(volume_fit, *params)
        
       # Write original data to output file
        write_original_data_output(original_data_output_file, volume, energy, total_atoms)
        
        # Write fitted parameters to output file
        write_fit_parameters_output(fit_parameters_output_file, params, total_atoms)
        
        # Write fitted curve data to output file
        write_fit_curve_output(fit_curve_output_file, volume_fit, energy_fit, total_atoms)
        
        # Plot the results
        plot_results(original_data_output_file, fit_curve_output_file, total_atoms, params)
        
        print("Original data written to", original_data_output_file)
        print("Fitted parameters written to", fit_parameters_output_file)
        print("Fitted curve data written to", fit_curve_output_file)
    except Exception as e:
        print(f"Error during fitting or writing output: {e}")

if __name__ == "__main__":
    main()

