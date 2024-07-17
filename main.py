import tkinter as tk
from tkinter import messagebox
import pint
import numpy as np
from time import sleep
from scipy.optimize import fsolve
# Create a unit registry
ureg = pint.UnitRegistry()

e = 1.6 * 10 **(-19) * ureg('coulomb')

rho = 1.293 * ureg('kg/m^3') # density of air
lambd = 6.73 * 10**-8 * ureg('meter') #mean free path
k = 1.38 * 10**-23 *ureg('J/K') # boltzmann constant
rho = 1.293 * ureg('kg/m^3') # density of air
lambd = 6.73 * 10**-8 * ureg('meter') #mean free path
k = 1.38 * 10**-23 *ureg('J/K') # boltzmann constant

def dynamic_viscosity_air(T):
    return (2.791 * 10**-7 * T**0.7355) * ureg('N*s/m^2')
mu = dynamic_viscosity_air(25+273.15)
def cunningham_coefficient(dp):
    mfp = 68 * ureg('nanometer')
    return 1 + mfp/dp * (2.514 + 0.800 * np.exp(-0.55* (dp/mfp)))

def electrical_mobility(dp, C):
    return (e*C)/(3*np.pi*mu*dp)



# Function to handle the submit button click
def submit():

    try:

        valid_units = {
            'DMA voltage': ['kilogram * meter ** 2 / ampere / second ** 3'],
            'DMA L': ['meter'],
            'Sheath flowrate': ['meter ** 3 / second'],
            'DMA r1': ['meter'],
            'DMA r2': ['meter'],
        }
        # Read and parse input values with units
        inputs = {
            'DMA voltage': ureg(entry1.get()),
            'Sheath flowrate': ureg(entry2.get()),
            'DMA r1': ureg(entry3.get()),
            'DMA r2': ureg(entry4.get()),
            'DMA L': ureg(entry5.get())
        }
        for key, quantity in inputs.items():
            try:
                unit = str(quantity.to_base_units().units)
                if unit not in valid_units[key]:
                    raise pint.errors.UndefinedUnitError(f"Invalid unit '{unit}' for {key}")
            except (ValueError, AttributeError, pint.errors.UndefinedUnitError) as e:
                messagebox.showerror("Invalid input",
                                      f"Please enter valid floating point numbers with units. Error: {e}")
                break
        
        K = (inputs['Sheath flowrate']*np.log(inputs['DMA r2']/inputs['DMA r1']))/(2*np.pi*inputs['DMA L'])


        e = 1.6 * 10 **(-19) * ureg('coulomb')
        A = ((e/(3*np.pi*mu)) / K)

        dma_voltage = inputs['DMA voltage']
        i = 0
        old = 300 * ureg('nanometer')
        new = 1 *ureg('nanometer')

        while True:
            new=A*dma_voltage*cunningham_coefficient(old)
            iteration_var.set(f"Iteration: {i}, D_p: {new.to('nanometer'):.2f}")
            root.update_idletasks()  # Update the GUI
            if np.abs(old.to('nanometer').magnitude-new.to('nanometer').magnitude) <= 0.01:
                break
            if i >= 10000:
                messagebox.showerror('Convergence error', 'Selected parameters did not converge')
                break
            i += 1
            old = new


        
        # Show the messagebox with formatted output
        messagebox.showinfo("Calculation",f"D_p = {new.to('nanometer')},\nIterations: {i}")

    except (ValueError, pint.errors.UndefinedUnitError) as e:
        messagebox.showerror("Invalid Input", f"Please enter valid floating point numbers with units. Error: {e}")

# Create the main window
root = tk.Tk()
root.title("Cunningham coefficient calculator")

# Create labels and entry widgets
tk.Label(root, text="Enter DMA voltage to calculate particle diameter:").grid(row=0, columnspan=2)

tk.Label(root, text="DMA voltage (e.g., 1000V):").grid(row=1, column=0)
entry1 = tk.Entry(root)
entry1.grid(row=1, column=1)

tk.Label(root, text="Sheath flowrate (e.g., 10 L/min):").grid(row=2, column=0)
entry2= tk.Entry(root)
entry2.grid(row=2, column=1)

tk.Label(root, text="DMA r1 (e.g., 3cm):").grid(row=3, column=0)
entry3 = tk.Entry(root)
entry3.grid(row=3, column=1)

tk.Label(root, text="DMA r2 (e.g., 1cm):").grid(row=4, column=0)
entry4 = tk.Entry(root)
entry4.grid(row=4, column=1)

tk.Label(root, text="DMA effective length (e.g., 10cm):").grid(row=5, column=0)
entry5 = tk.Entry(root)
entry5.grid(row=5, column=1)

iteration_var = tk.StringVar()
iteration_label = tk.Label(root, textvariable=iteration_var)
iteration_label.grid(row=6, columnspan=2)

# Create submit and cancel buttons
tk.Button(root, text="Calculate", command=submit).grid(row=7, column=0)
tk.Button(root, text="Exit", command=root.quit).grid(row=7, column=1)

# Run the main event loop
root.mainloop()
