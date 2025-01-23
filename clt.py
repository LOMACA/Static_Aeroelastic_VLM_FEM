###########################################################################################################

# Classical Laminate Theory - According to Schürmann in Konstruieren mit Faser-Kunststoffverbunden

###########################################################################################################

import tkinter as tk
from tkinter import ttk
import numpy as np

class CLT(ttk.Frame):

    def __init__(self, parent, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        self.create_input_widgets()
        self.create_output_field()
        self.create_buttons()

    def create_input_widgets(self):

        """
        Create the input widgets for laminate properties.
        """

        input_frame = ttk.LabelFrame(self, text="Input Parameters")
        input_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        ttk.Label(input_frame, text="Layer Thickness (mm):").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.layer_thickness_entry = ttk.Entry(input_frame)
        self.layer_thickness_entry.grid(row=0, column=1, padx=5, pady=5, sticky="w")

        ttk.Label(input_frame, text="Fiber Volume Fraction (%):").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.fiber_volume_fraction_entry = ttk.Entry(input_frame)
        self.fiber_volume_fraction_entry.grid(row=1, column=1, padx=5, pady=5, sticky="w")

        ttk.Label(input_frame, text="Total Thickness (mm):").grid(row=2, column=0, padx=5, pady=5, sticky="w")
        self.root_width_entry = ttk.Entry(input_frame)
        self.root_width_entry.grid(row=2, column=1, padx=5, pady=5, sticky="w")

        ttk.Label(input_frame, text="Layer Orientations (comma-separated):").grid(row=3, column=0, padx=5, pady=5, sticky="w")
        self.layer_orientations_entry = ttk.Entry(input_frame)
        self.layer_orientations_entry.grid(row=3, column=1, padx=5, pady=5, sticky="w")

        ttk.Label(input_frame, text="Fiber Type:").grid(row=4, column=0, padx=5, pady=5, sticky="w")
        self.fiber_type_var = tk.StringVar(value="HT")
        fiber_type_dropdown = ttk.Combobox(input_frame, textvariable=self.fiber_type_var, values=["HT", "HM"])
        fiber_type_dropdown.grid(row=4, column=1, padx=5, pady=5, sticky="w")

    def create_output_field(self):

        """
        Create the output field for displaying results.
        """

        self.output_frame = ttk.LabelFrame(self, text="Results")
        self.output_frame.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")

        self.output_field = tk.Text(self.output_frame, wrap=tk.WORD, width=50, height=15)
        self.output_field.grid(row=0, column=0, padx=5, pady=5)

    def create_buttons(self):

        """
        Create buttons for calculation and clearing inputs.
        """

        button_frame = ttk.Frame(self)
        button_frame.grid(row=2, column=0, padx=10, pady=10, sticky="ew")

        calculate_button = ttk.Button(button_frame, text="Calculate Laminate Properties", command=self.calculate_and_display_laminate_properties)
        calculate_button.grid(row=0, column=0, padx=5, pady=5, sticky="ew")

        clear_button = ttk.Button(button_frame, text="Clear Inputs", command=self.clear_inputs)
        clear_button.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

    def assign_fiber_properties(self):

        """
        Assign fiber material properties based on selected fiber type (HT or HM).

        Returns: (E_parallel, E_perpendicular, G_parallel_perpendicular,
                    nu_parallel_perpendicular, rho_f) as tuple
        """

        fiber_type = self.fiber_type_var.get()

        if fiber_type == "HT":
            E_parallel = 230000 * 1e6
            E_perpendicular = 28000 * 1e6
            G_parallel_perpendicular = 50000 * 1e6
            nu_parallel_perpendicular = 0.23
            rho_f = 1.74 * 1000
        elif fiber_type == "HM":
            E_parallel = 392000 * 1e6
            E_perpendicular = 15200 * 1e6
            G_parallel_perpendicular = 28600 * 1e6
            nu_parallel_perpendicular = 0.2
            rho_f = 1.81 * 1000

        else:
            raise ValueError("Invalid fiber type selected.")

        return E_parallel, E_perpendicular, G_parallel_perpendicular, nu_parallel_perpendicular, rho_f

        ### CALCULATE SINGLE LAYER PROPERTIES ###

    def calculate_youngs_modulus_parallel(self, E_f_parallel, fiber_volume_fraction, E_m):

        """
        Calculate the longitudinal Young's modulus of the composite material.

        Arguments:
            E_f_parallel: Longitudinal Young's modulus of fiber.
            fiber_volume_fraction: Fiber volume fraction.
            E_m: Young's modulus of the matrix.

        Returns: Effective Young's modulus in the parallel direction.
        """

        return E_f_parallel * fiber_volume_fraction + E_m * (1 - fiber_volume_fraction)

    def calculate_youngs_modulus_orthogonal(self, E_m, v_m, fiber_volume_fraction, E_f_perp):

        """
        Calculate Young's modulus orthogonal to fiber direction.

        Arguments:
            E_m: Matrix Young's modulus.
            v_m: Matrix Poisson's ratio.
            fiber_volume_fraction: Fiber volume fraction.
            E_f_perp: Fiber Young's modulus orthogonal to fiber direction.

        Returns: Effective Young's modulus in the orthogonal direction.
        """

        return (E_m / (1 - v_m ** 2)) * (1 + 0.85 * fiber_volume_fraction ** 2) / \
            ((1 - fiber_volume_fraction) ** 1.25 + (E_m / ((1 - v_m ** 2) * E_f_perp)) * fiber_volume_fraction)

    def calculate_shear_modulus(self, fiber_volume_fraction, G_f):

        """
        Calculate shear modulus for the composite material.

        Arguments:
            fiber_volume_fraction: Fiber volume fraction.
            G_f: Shear modulus of fiber.

        Returns: Effective shear modulus.
        """

        G_m = 1019 * 1e6
        return G_m * (1 + 0.4 * fiber_volume_fraction ** 0.5) / \
            ((1 - fiber_volume_fraction) ** 1.45 + (G_m / G_f) * fiber_volume_fraction)

    def calculate_poissons_ratio(self, fiber_volume_fraction, v_f, v_m):

        """
        Calculate Poisson's ratio for the composite material.

        Arguments:
            fiber_volume_fraction: Fiber volume fraction.
            v_f: Poisson's ratio of fiber.
            v_m: Poisson's ratio of matrix.

        Returns: Effective Poisson's ratio.
        """

        return fiber_volume_fraction * v_f + (1 - fiber_volume_fraction) * v_m


    def calculate_effective_moduli(self, A_matrix, t):

        """
        Calculate effective laminate parameters from stress-based stiffness matrix of the laminate.

        Arguments:
            A_matrix: A matrix from the laminate theory.
            t: Total thickness of the laminate.

        Returns: (E_x, E_y, G_xy, nu_xy) Effective moduli and Poisson's ratio.
        """

        A_inv = np.linalg.inv(A_matrix)

        A_inv_11 = A_inv[0, 0]
        A_inv_22 = A_inv[1, 1]
        A_inv_66 = A_inv[2, 2]
        A_inv_12 = A_inv[0, 1]

        E_x = 1 / (A_inv_11 * t)
        E_y = 1 / (A_inv_22 * t)
        G_xy = 1 / (A_inv_66 * t)

        nu_xy = -A_inv_12 / A_inv_22
        nu_yx = -A_inv_12 / A_inv_11

        return E_x, E_y, G_xy, nu_xy

    def create_symmetric_matrices(self, num_layers_half, num_layers_total, layer_orientations, E_parallel,
                                  E_perpendicular, G_parallel_perpendicular, nu_parallel_perpendicular):

        """
        Create symmetric stiffness matrices for each layer orientation to represent laminate properties.

        Arguments:
            num_layers_half: Half of the total layers (used for symmetry).
            num_layers_total: Total number of layers.
            layer_orientations: List of layer orientations (angles).
            E_parallel: Effective Young's modulus parallel to fiber.
            E_perpendicular: Effective Young's modulus perpendicular to fiber.
            G_parallel_perpendicular: Effective shear modulus.
            nu_parallel_perpendicular: Effective Poisson's ratio.

        Returns: List of transformed Q matrices for each layer.
        """

        Q_matrices = []

        for i in range(num_layers_half):
            orientation = layer_orientations[i % len(layer_orientations)]
            Q_bar = self.calculate_transformed_Q_matrix(E_parallel, E_perpendicular, G_parallel_perpendicular,
                                                        nu_parallel_perpendicular, orientation)
            Q_matrices.append(Q_bar)

        if num_layers_half * 2 < num_layers_total:
            orientation = layer_orientations[0]
            Q_bar = self.calculate_transformed_Q_matrix(E_parallel, E_perpendicular, G_parallel_perpendicular,
                                                        nu_parallel_perpendicular, orientation)
            Q_matrices.append(Q_bar)

        for i in range(num_layers_half - 1, -1, -1):
            Q_matrices.append(Q_matrices[i])

        return Q_matrices

    def calculate_transformed_Q_matrix(self, E_parallel, E_perpendicular, G_parallel_perpendicular,
                                       nu_parallel_perpendicular,
                                       alpha):

        """
        Transform stiffness matrix based on fiber orientation angle.

        Arguments:
            E_parallel: Young's modulus parallel to fiber.
            E_perpendicular: Young's modulus perpendicular to fiber.
            G_parallel_perpendicular: Shear modulus.
            nu_parallel_perpendicular: Poisson's ratio.
            alpha: Fiber orientation angle in degrees.

        Returns: Transformed stiffness matrix Q.
        """

        alpha_rad = np.radians(alpha)

        Q11 = E_parallel / (1 - nu_parallel_perpendicular ** 2)
        Q22 = E_perpendicular / (1 - nu_parallel_perpendicular ** 2)
        Q12 = nu_parallel_perpendicular * E_perpendicular / (1 - nu_parallel_perpendicular ** 2)
        Q66 = G_parallel_perpendicular

        cos_alpha = np.cos(alpha_rad)
        sin_alpha = np.sin(alpha_rad)
        cos2_alpha = cos_alpha ** 2
        sin2_alpha = sin_alpha ** 2
        cos4_alpha = cos2_alpha ** 2
        sin4_alpha = sin2_alpha ** 2
        sin2_alpha_2 = np.sin(2 * alpha_rad) ** 2

        Q_bar11 = Q11 * cos4_alpha + Q22 * sin4_alpha + 0.5 * (Q12 + 2 * Q66) * sin2_alpha_2
        Q_bar22 = Q11 * sin4_alpha + Q22 * cos4_alpha + 0.5 * (Q12 + 2 * Q66) * sin2_alpha_2
        Q_bar12 = Q12 + 0.25 * (Q11 + Q22 - 2 * Q12 - 4 * Q66) * sin2_alpha_2
        Q_bar66 = Q66 + 0.25 * (Q11 + Q22 - 2 * Q12 - 4 * Q66) * sin2_alpha_2
        Q_bar16 = -0.5 * ((Q11 + Q22 - 2 * Q12 - 4 * Q66) * sin2_alpha - (Q11 - Q12 - 2 * Q66)) * np.sin(2 * alpha_rad)
        Q_bar26 = -0.5 * ((Q22 - Q12 - 2 * Q66) - (Q11 + Q22 - 2 * Q12 - 4 * Q66) * sin2_alpha) * np.sin(2 * alpha_rad)

        Q_bar = np.array([
            [Q_bar11, Q_bar12, Q_bar16],
            [Q_bar12, Q_bar22, Q_bar26],
            [Q_bar16, Q_bar26, Q_bar66]
        ])

        return Q_bar

    def calculate_transformed_S_matrix(self, E_parallel, E_perpendicular, G_parallel_perpendicular,
                                       nu_parallel_perpendicular,
                                       alpha):

        """
        Transform compliance matrix based on fiber orientation angle.

        Arguments:
            E_parallel: Young's modulus parallel to fiber.
            E_perpendicular: Young's modulus perpendicular to fiber.
            G_parallel_perpendicular: Shear modulus.
            nu_parallel_perpendicular: Poisson's ratio.
            alpha: Fiber orientation angle in degrees.

        Returns: Transformed S matrix.
        """

        alpha_rad = np.radians(alpha)

        S11 = 1 / E_parallel
        S22 = 1 / E_perpendicular
        S66 = 1 / G_parallel_perpendicular
        S12 = nu_parallel_perpendicular / E_parallel

        cos_alpha = np.cos(alpha_rad)
        sin_alpha = np.sin(alpha_rad)
        cos2_alpha = cos_alpha ** 2
        sin2_alpha = sin_alpha ** 2
        sin3_alpha = sin2_alpha * sin_alpha
        cos3_alpha = cos2_alpha * cos_alpha
        cos4_alpha = cos2_alpha ** 2
        sin4_alpha = sin2_alpha ** 2
        sin2_alpha_2 = np.sin(2 * alpha_rad) ** 2
        cos2_alpha_2 = np.cos(2 * alpha_rad) ** 2

        S_bar11 = cos4_alpha / E_parallel + sin4_alpha / E_perpendicular + 0.25 * (S66 - 2 * S12) * sin2_alpha_2
        S_bar22 = sin4_alpha / E_parallel + cos4_alpha / E_perpendicular + 0.25 * (S66 - 2 * S12) * sin2_alpha_2
        S_bar66 = cos2_alpha_2 / G_parallel_perpendicular + (S11 + S22 + 2 * S12) * sin2_alpha_2
        S_bar12 = 0.25 * (S11 + S22 - S66) * sin2_alpha_2 - S12 * (sin4_alpha + cos4_alpha)
        S_bar16 = -(2 / E_perpendicular + 2 * S12 - S66) * sin3_alpha * cos_alpha + (
                2 / E_parallel + 2 * S12 - S66) * cos3_alpha * sin_alpha
        S_bar26 = -(2 / E_perpendicular + 2 * S12 - S66) * cos3_alpha * sin_alpha + (
                2 / E_parallel + 2 * S12 - S66) * sin3_alpha * cos_alpha

        S_bar = np.array([
            [S_bar11, S_bar12, S_bar16],
            [S_bar12, S_bar22, S_bar26],
            [S_bar16, S_bar26, S_bar66]
        ])

        return S_bar

    def calculate_A_matrix(self, Q_matrices, layer_thickness):

        """
        Calculate stress-transformed stiffness matrix for the laminate by integrating Q matrices across thickness.

        Parameters:
            Q_matrices: List of transformed Q matrices for each layer.
            layer_thickness: Thickness of each layer.

        Returns: stress-transformed stiffness matrix A
        """

        A_matrix = np.zeros((3, 3))

        for Q_bar in Q_matrices:
            A_matrix[0, 0] += (Q_bar[0, 0] * layer_thickness)
            A_matrix[0, 1] += (Q_bar[0, 1] * layer_thickness)
            A_matrix[0, 2] += (Q_bar[0, 2] * layer_thickness)
            A_matrix[1, 0] += (Q_bar[1, 0] * layer_thickness)
            A_matrix[1, 1] += (Q_bar[1, 1] * layer_thickness)
            A_matrix[1, 2] += (Q_bar[1, 2] * layer_thickness)
            A_matrix[2, 0] += (Q_bar[2, 0] * layer_thickness)
            A_matrix[2, 1] += (Q_bar[2, 1] * layer_thickness)
            A_matrix[2, 2] += (Q_bar[2, 2] * layer_thickness)

        return A_matrix

    def calculate_and_display_laminate_properties(self):

        """
        Outputs the laminate properties to the GUI text field.
        """

        try:

            layer_thickness = float(self.layer_thickness_entry.get()) / 1000
            fiber_volume_fraction = float(self.fiber_volume_fraction_entry.get()) / 100
            total_thickness = float(self.root_width_entry.get()) / 1000
            layer_orientations = list(map(float, self.layer_orientations_entry.get().split(',')))

            E_f_parallel, E_f_perp, G_f, nu_f, rho_f = self.assign_fiber_properties()
            E_m = 3150 * 1e6
            v_m = 0.37

            E_parallel = self.calculate_youngs_modulus_parallel(E_f_parallel, fiber_volume_fraction, E_m)
            E_perpendicular = self.calculate_youngs_modulus_orthogonal(E_m, v_m, fiber_volume_fraction, E_f_perp)
            G_parallel_perpendicular = self.calculate_shear_modulus(fiber_volume_fraction, G_f)
            nu_parallel_perpendicular = self.calculate_poissons_ratio(fiber_volume_fraction, nu_f, v_m)

            num_layers = int(total_thickness / layer_thickness)
            num_layers_half = num_layers // 2

            Q_matrices = self.create_symmetric_matrices(num_layers_half, num_layers, layer_orientations, E_parallel,
                                                        E_perpendicular, G_parallel_perpendicular,
                                                        nu_parallel_perpendicular)

            A_matrix = self.calculate_A_matrix(Q_matrices, layer_thickness)

            self.E_x, self.E_y, self.G_xy, self.nu_xy = self.calculate_effective_moduli(A_matrix, total_thickness)

            self.output_field.delete(1.0, tk.END)
            self.output_field.insert(tk.END, "Laminate Properties:\n")
            self.output_field.insert(tk.END, f"Stiffness Matrix:\n{A_matrix}\n\n")
            self.output_field.insert(tk.END, f"Effective Moduli:\n")
            self.output_field.insert(tk.END, f"E_x: {self.E_x / 1e6:.2f} MPa\n")
            self.output_field.insert(tk.END, f"E_y: {self.E_y / 1e6:.2f} MPa\n")
            self.output_field.insert(tk.END, f"G_xy: {self.G_xy / 1e6:.2f} MPa\n")
            self.output_field.insert(tk.END, f"nu_xy: {self.nu_xy:.4f}\n")

        except Exception as e:
            self.output_field.insert(tk.END, f"Error: {str(e)}\n")

    def get_laminate_properties(self):

        """
        Returns the computed laminate properties including stiffness and effective moduli.
        """

        return {
            'E_x': self.E_x,
            'E_y': self.E_y,
            'G_xy': self.G_xy,
            'nu_xy': self.nu_xy
        }

    def clear_inputs(self):

        """
        Clear all input fields.
        """

        self.layer_thickness_entry.delete(0, tk.END)
        self.fiber_volume_fraction_entry.delete(0, tk.END)
        self.root_width_entry.delete(0, tk.END)
        self.layer_orientations_entry.delete(0, tk.END)
        self.output_field.delete("1.0", tk.END)

    def get_input_values(self):

        """
        Retrieve input values from GUI fields for CLT parameters.

        Returns: Dictionary containing the CLT input values.
        """

        try:
            layer_thickness = float(self.layer_thickness_entry.get())
            fiber_volume_fraction = float(self.fiber_volume_fraction_entry.get())
            total_thickness = float(self.root_width_entry.get())
            layer_orientations = list(
                map(float, self.layer_orientations_entry.get().split(',')))
            fiber_type = self.fiber_type_var.get()

            return {
                'layer_thickness': layer_thickness,
                'fiber_volume_fraction': fiber_volume_fraction,
                'total_thickness': total_thickness,
                'layer_orientations': layer_orientations,
                'fiber_type': fiber_type,
            }
        except ValueError as e:
            print(f"Error getting CLT input values: {e}")
            return None

    def set_input_values(self, input_values):

        """
        Set input values in GUI fields for CLT parameters.

        Parameters:
            input_values: Dictionary containing CLT input values.
        """

        if 'layer_thickness' in input_values:
            self.layer_thickness_entry.delete(0, tk.END)
            self.layer_thickness_entry.insert(0, input_values['layer_thickness'])

        if 'fiber_volume_fraction' in input_values:
            self.fiber_volume_fraction_entry.delete(0, tk.END)
            self.fiber_volume_fraction_entry.insert(0, input_values['fiber_volume_fraction'])

        if 'total_thickness' in input_values:
            self.root_width_entry.delete(0, tk.END)
            self.root_width_entry.insert(0, input_values['total_thickness'])

        if 'layer_orientations' in input_values:
            orientations = ', '.join(map(str, input_values['layer_orientations']))
            self.layer_orientations_entry.delete(0, tk.END)
            self.layer_orientations_entry.insert(0, orientations)

        if 'fiber_type' in input_values:
            self.fiber_type_var.set(input_values['fiber_type'])


