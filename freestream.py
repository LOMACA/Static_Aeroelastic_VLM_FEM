###############################################################

# Freestream conditions module

###############################################################

import tkinter as tk
from tkinter import ttk
import numpy as np

class Freestream(ttk.Frame):

    """
    Class to handle the user input for freestream conditions including velocity, density, and angle of attack.
    Provides methods for retrieving and submitting user input values.
    """

    def __init__(self, parent, *args, **kwargs):

        """
        Initializes the Freestream class UI elements for velocity, density, and angle of attack.
        """

        super().__init__(parent, *args, **kwargs)

        self.submitted_data = None

        ttk.Label(self, text="Freestream Velocity (m/s):", font=("Arial", 12)).grid(row=0, column=0, padx=5, pady=5,
                                                                                    sticky="w")
        self.velocity_entry = ttk.Entry(self)
        self.velocity_entry.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

        ttk.Label(self, text="Freestream Density (kg/m^3):", font=("Arial", 12)).grid(row=1, column=0, padx=5, pady=5,
                                                                                      sticky="w")
        self.density_entry = ttk.Entry(self)
        self.density_entry.grid(row=1, column=1, padx=5, pady=5, sticky="ew")

        ttk.Label(self, text="Angle of Attack (degrees):", font=("Arial", 12)).grid(row=2, column=0, padx=5, pady=5,
                                                                                    sticky="w")
        self.aoa_entry = ttk.Entry(self)
        self.aoa_entry.grid(row=2, column=1, padx=5, pady=5, sticky="ew")

        submit_button = ttk.Button(self, text="Submit", command=self.submit_freestream_inputs)
        submit_button.grid(row=3, column=0, columnspan=2, padx=5, pady=10, sticky="ew")

        self.status_label = ttk.Label(self, text="", font=("Arial", 10, 'italic'))
        self.status_label.grid(row=4, column=0, columnspan=2, padx=5, pady=5)

    def get_freestream_inputs(self):

        """
        Retrieves and converts the user inputs for velocity, density, and angle of attack (AoA).
        Converts AoA to radians.

        Returns:
            dict: Freestream input values in SI units, or None if invalid input is detected.
        """

        try:
            velocity = float(self.velocity_entry.get())
            density = float(self.density_entry.get())
            alpha_deg = float(self.aoa_entry.get())
            alpha_rad = np.radians(alpha_deg)

            return {
                'velocity': velocity,
                'density': density,
                'alpha': alpha_rad
            }
        except ValueError as e:
            print(f"Invalid input for freestream conditions: {e}")
            return None

    def submit_freestream_inputs(self):

        """
        Retrieves the freestream inputs and updates the status label to indicate success or failure.
        Sets `self.submitted_data` with the input dictionary if inputs are valid.
        """

        freestream_data = self.get_freestream_inputs()
        if freestream_data is not None:
            self.submitted_data = freestream_data
            self.status_label.config(text="Freestream inputs submitted successfully!", foreground="green")
        else:
            self.status_label.config(text="Invalid inputs. Please check your values.", foreground="red")

    def set_freestream_inputs(self, input_values):

        """
        Sets the input fields for velocity, density, and angle of attack from a provided dictionary.

        """

        if 'velocity' in input_values:
            self.velocity_entry.delete(0, tk.END)
            self.velocity_entry.insert(0, input_values['velocity'])

        if 'density' in input_values:
            self.density_entry.delete(0, tk.END)
            self.density_entry.insert(0, input_values['density'])

        if 'aoa' in input_values:
            self.aoa_entry.delete(0, tk.END)
            self.aoa_entry.insert(0, input_values['aoa'])

    def get_submitted_data(self):

        """
        Returns the submitted freestream data after calling `submit_freestream_inputs`.

        Returns:
            dict: Contains velocity, density, and AoA in radians, or None if not yet submitted.
        """

        return self.submitted_data

    def get_input_values(self):

        """
        Retrieves raw user input values for velocity, density, and angle of attack (in degrees).
        Returns:
            dict: Freestream input values without conversions, or None if invalid.
        """

        try:
            velocity = float(self.velocity_entry.get())
            density = float(self.density_entry.get())
            alpha_deg = float(self.aoa_entry.get())
            return {
                'velocity': velocity,
                'density': density,
                'aoa': alpha_deg
            }
        except ValueError as e:
            print(f"Error getting freestream input values: {e}")
            return None

    def set_input_values(self, input_values):

        """
        Sets the input fields based on provided values for velocity, density, and AoA in degrees.
        """

        if 'velocity' in input_values:
            self.velocity_entry.delete(0, tk.END)
            self.velocity_entry.insert(0, input_values['velocity'])

        if 'density' in input_values:
            self.density_entry.delete(0, tk.END)
            self.density_entry.insert(0, input_values['density'])

        if 'aoa' in input_values:
            self.aoa_entry.delete(0, tk.END)
            self.aoa_entry.insert(0, input_values['aoa'])
