#################################################################################

# FEM module with fixed cross-sectional properties input by the user

#################################################################################

import numpy as np
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog
#from matplotlib.figure import Figure
#import pandas as pd

class FEMAnalysisAlt(ttk.Frame):

    """
    Finite Element Analysis class for a 3D Euler-Bernoulli Beam integrated into the GUI.
    This class initializes the beam properties, computes the element stiffness matrix,
    assembles the global stiffness matrix, applies aerodynamic loads, and solves for
    displacements.
    """

    def __init__(self, parent, freestream, geometry, vlm, clt,  *args, **kwargs):

        super().__init__(parent, *args, **kwargs)
        self.freestream = freestream
        self.geometry = geometry
        self.vlm = vlm
        self.clt = clt

        self.entry_frame = ttk.Frame(self)
        self.entry_frame.grid(row=0, column=0, columnspan=2, pady=5, sticky="n")

        ttk.Label(self.entry_frame, text="Beam Position (% of Chord Length):").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self.beam_pos_entry = ttk.Entry(self.entry_frame)
        self.beam_pos_entry.grid(row=0, column=1, padx=5, pady=5, sticky="w")

        ttk.Label(self.entry_frame, text="Axial Rigidity (EA):").grid(row=1, column=0, padx=5, pady=5, sticky="w")
        self.EA_entry = ttk.Entry(self.entry_frame)
        self.EA_entry.grid(row=1, column=1, padx=5, pady=5, sticky="w")

        ttk.Label(self.entry_frame, text="Torsional Rigidity (GJ):").grid(row=2, column=0, padx=5, pady=5, sticky="w")
        self.GJ_entry = ttk.Entry(self.entry_frame)
        self.GJ_entry.grid(row=2, column=1, padx=5, pady=5, sticky="w")

        ttk.Label(self.entry_frame, text="Flexural Rigidity (EIy):").grid(row=3, column=0, padx=5, pady=5, sticky="w")
        self.EIy_entry = ttk.Entry(self.entry_frame)
        self.EIy_entry.grid(row=3, column=1, padx=5, pady=5, sticky="w")

        ttk.Label(self.entry_frame, text="Flexural Rigidity (EIz):").grid(row=4, column=0, padx=5, pady=5, sticky="w")
        self.EIz_entry = ttk.Entry(self.entry_frame)
        self.EIz_entry.grid(row=4, column=1, padx=5, pady=5, sticky="w")

        run_button = ttk.Button(self.entry_frame, text="Run FE Analysis", command=self.solve_displacements)
        run_button.grid(row=6, column=1, padx=5, pady=5, sticky="ew")

        ttk.Button(self.entry_frame, text="View Structural Mesh", command=self.update_plot_with_geometry).grid(
            row=6, column=0, padx=5, pady=5, sticky="ew")

        self.material_source = tk.StringVar(value="manual")

        self.material_frame = ttk.LabelFrame(self, text="Material Properties")
        self.material_frame.grid(row=7, column=0, padx=10, pady=10, sticky="nsew")

        ttk.Radiobutton(
            self.material_frame, text="Manual Input", variable=self.material_source, value="manual",
            command=self.toggle_material_input
        ).grid(row=0, column=0, padx=5, pady=5, sticky="w")
        ttk.Radiobutton(
            self.material_frame, text="Import from Classical Laminate Theory", variable=self.material_source, value="clt",
            command=self.toggle_material_input
        ).grid(row=0, column=1, padx=5, pady=5, sticky="w")

        self.manual_material_frame = ttk.LabelFrame(self, text="Manual Material Properties")
        self.manual_material_frame.grid(row=8, column=0, padx=10, pady=10, sticky="nsew")

        ttk.Label(self.manual_material_frame, text="Young's Modulus (E):").grid(row=0, column=0, padx=5, pady=5,
                                                                                     sticky="w")
        self.youngs_modulus_entry = ttk.Entry(self.manual_material_frame)
        self.youngs_modulus_entry.grid(row=0, column=1, padx=5, pady=5, sticky="w")

        ttk.Label(self.manual_material_frame, text="Poisson's Ratio:").grid(row=1, column=0, padx=5, pady=5,
                                                                                 sticky="w")
        self.poisson_ratio_entry = ttk.Entry(self.manual_material_frame)
        self.poisson_ratio_entry.grid(row=1, column=1, padx=5, pady=5, sticky="w")

        ttk.Label(self.manual_material_frame, text="Shear Modulus (G):").grid(row=2, column=0, padx=5, pady=5,
                                                                                   sticky="w")
        self.shear_modulus_entry = ttk.Entry(self.manual_material_frame)
        self.shear_modulus_entry.grid(row=2, column=1, padx=5, pady=5, sticky="w")

        strain_button = ttk.Button(self.entry_frame, text="Calculate Strains", command=self.calculate_strains)
        strain_button.grid(row=8, column=0, padx=5, pady=5, sticky="ew")

        stress_button = ttk.Button(self.entry_frame, text="Calculate Stresses", command=self.calculate_stresses)
        stress_button.grid(row=8, column=1, padx=5, pady=5, sticky="ew")

        self.button_frame = ttk.Frame(self)
        self.button_frame.grid(row=0, column=6, columnspan=2, pady=5, sticky="nsew")

        ttk.Button(self.button_frame, text="View Displaced Geometry", command=self.update_plot_with_displacement).grid(row=0, column=0, pady=5)
        ttk.Button(self.button_frame, text="Plot Z-Displacement", command=self.update_plot_with_z_displacement).grid(row=1, column=0, pady=5)
        ttk.Button(self.button_frame, text="Plot X-Displacement", command=self.update_plot_with_x_displacement).grid(
            row=2, column=0, pady=5)
        ttk.Button(self.button_frame, text="Plot Twist", command=self.update_plot_with_twist).grid(
            row=3, column=0, pady=5)
        ttk.Button(self.button_frame, text="Plot Strains", command=self.update_plot_with_strains).grid(row=4,
                                                                                                           column=0,
                                                                                                           pady=5)
        ttk.Button(self.button_frame, text="Plot Stresses", command=self.update_plot_with_stresses).grid(row=5,
                                                                                                             column=0,
                                                                                                             pady=5)

        self.export_option = tk.BooleanVar()
        self.export_checkbox = ttk.Checkbutton(self.button_frame, text="Export Plots", variable=self.export_option, command=self.toggle_export_folder)
        self.export_checkbox.grid(row=6, column=0, pady=5)
        self.select_folder_button = ttk.Button(self.button_frame, text="Select Export Folder",
                                               command=self.select_export_folder)
        self.select_folder_button.grid(row=7, column=0, pady=5)
        self.select_folder_button.config(state=tk.DISABLED)

        self.plot_frame = ttk.Frame(self)
        self.plot_frame.grid(row=0, column=4, columnspan=2, padx=5, pady=5, sticky="nsew")
        self.figure = plt.figure(figsize=(8, 7))
        self.canvas = FigureCanvasTkAgg(self.figure, self.plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH)

        self.output_field = tk.Text(self.entry_frame, width=80, height=20)
        self.output_field.grid(row=9, column=0, columnspan=2, padx=10, pady=10)

        self.export_folder = None

    def toggle_material_input(self):

        """
        Toggles between manual input and importing material properties from CLT.
        """

        if self.material_source.get() == "manual":
            self.manual_material_frame.grid(row=10, column=0, padx=10, pady=10, sticky="nsew")
        else:
            self.manual_material_frame.grid_forget()

    def select_export_folder(self):

        """
        Opens a dialog to select the folder where exported plots will be saved.
        """

        folder = filedialog.askdirectory(title="Select Folder for Exported Plots")
        if folder:
            self.export_folder = folder
            self.output_field.insert(tk.END, f"Export folder selected: {folder}\n")

    def toggle_export_folder(self):

        """
        Enables or disables the folder selection button based on the checkbox state.
        """

        if self.export_option.get():
            self.select_folder_button.config(state=tk.NORMAL)
        else:
            self.select_folder_button.config(state=tk.DISABLED)

    def update_plot(self):

        self.figure.clf()
        self.current_plot_function()
        self.canvas.draw()

    def update_plot_with_geometry(self):

        """
        Updates the plot with the structural mesh.
        """

        self.current_plot_function = self.build_structural_geometry
        self.update_plot()

    def update_plot_with_displacement(self):

        """
        Updates the plot with the total displacement.
        """

        self.current_plot_function = self.plot_displaced_structure
        self.update_plot()

    def update_plot_with_z_displacement(self):

        """
        Updates the plot with the z-displacement distribution.
        """

        self.current_plot_function = self.plot_z_displacement
        self.update_plot()

    def update_plot_with_x_displacement(self):

        """
        Updates the plot with the x-displacement distribution.
        """

        self.current_plot_function = self.plot_x_displacement
        self.update_plot()

    def update_plot_with_twist(self):

        """
        Updates the plot with the twist distribution.
        """

        self.current_plot_function = self.plot_twist
        self.update_plot()

    def update_plot_with_strains(self):

        """
        Updates the plot with the strain distribution.
        """

        self.current_plot_function = self.plot_strains
        self.update_plot()

    def update_plot_with_stresses(self):

        """
        Updates the plot with the stress distribution.
        """

        self.current_plot_function = self.plot_stresses
        self.update_plot()

    def initialize_beam_properties(self):

        """
        Initializes the beam properties based on the spanwise panel configuration from the Geometry class.
        The function retrieves the number of spanwise panels and computes properties like element length.
        """

        try:
            self.EA = float(self.EA_entry.get())
            self.GJ = float(self.GJ_entry.get())
            self.EIy = float(self.EIy_entry.get())
            self.EIz = float(self.EIz_entry.get())
            num_spanwise_panels = int(self.geometry.n_panels_span_entry.get())
            self.length = float(self.geometry.span_entry.get()) / 2
        except ValueError:
            self.output_field.insert(tk.END, "Error: Please provide valid inputs in the Geometry tab.\n")
            return False

        self.num_elements = num_spanwise_panels - 1
        self.num_nodes = self.num_elements + 1

        return True

    def beam_element_stiffness_matrix(self, l):

        """
        Computes the element stiffness matrix for a beam element based on its geometry and material properties.

        Arguments:
            l: Width of the panel (element length).
        """

        self.initialize_beam_properties()

        L = self.length

        EA = self.EA
        EIy = self.EIy
        EIz = self.EIz
        GJ = self.GJ

        k = np.zeros((12, 12))

        k_1 = EA / L
        k_2_z = EIz / L ** 3
        k_2_y = EIy / L ** 3
        k_3 = GJ / L

        k[0, 0] = k[6, 6] = k_1
        k[0, 6] = k[6, 0] = -k_1

        k[3, 3] = k[9, 9] = k_3
        k[3, 9] = k[9, 3] = -k_3

        k[1, 1] = k[7, 7] = 12 * k_2_z
        k[1, 7] = k[7, 1] = -12 * k_2_z
        k[1, 5] = k[5, 1] = 6 * k_2_z * l
        k[1, 11] = k[11, 1] = 6 * k_2_z * l
        k[5, 7] = k[7, 5] = -6 * k_2_z * l
        k[5, 11] = k[11, 5] = 2 * k_2_z * l ** 2
        k[5, 5] = k[11, 11] = 4 * k_2_z * l ** 2
        k[7, 11] = k[11, 7] = -6 * k_2_z * l
        k[8, 10] = k[10, 8] = 6 * k_2_y * l

        k[2, 2] = k[8, 8] = 12 * k_2_y
        k[2, 8] = k[8, 2] = -12 * k_2_y
        k[2, 4] = k[4, 2] = -6 * k_2_y * l
        k[2, 10] = k[10, 2] = -6 * k_2_y * l
        k[4, 8] = k[8, 4] = 6 * k_2_y * l
        k[4, 10] = k[10, 4] = 2 * k_2_y * l ** 2
        k[4, 4] = k[10, 10] = 4 * k_2_y * l ** 2

        return k

    def get_local_chords(self, y_stations):

        """
        Calculate local chord lengths at each spanwise station for the wing.

        Arguments:
            y_stations: List of spanwise positions.

        Returns: Array of chord lengths at each y_station.
        """

        root_chord = float(self.geometry.root_chord_entry.get())
        tip_chord = float(self.geometry.tip_chord_entry.get())
        span = float(self.geometry.span_entry.get()) / 2

        local_chords = root_chord - (root_chord - tip_chord) * (np.array(y_stations) / span)

        return local_chords

    def build_structural_geometry(self):

        """
        Builds a structural mesh overlaying the aerodynamic mesh.
        Each aerodynamic panel defines structural nodes on its edges,
        and CP is at the center of the panel at 25% of the chord.
        """

        wing_geometry_data = self.geometry.build_aerodynamic_mesh()
        panel_coordinates = wing_geometry_data['panel_coordinates']

        dihedral_angle_deg = float(self.geometry.dihedral_entry.get())
        dihedral_angle_rad = np.radians(dihedral_angle_deg)
        sweep_angle_deg = float(self.geometry.sweep_entry.get())
        sweep_angle_rad = np.radians(sweep_angle_deg)
        num_chordwise_panels = int(self.geometry.n_panels_chord_entry.get())

        elastic_axis_offset = (float(self.beam_pos_entry.get()) / 100) - 0.25

        root_chord = float(self.geometry.root_chord_entry.get())
        tip_chord = float(self.geometry.tip_chord_entry.get())
        span = float(self.geometry.span_entry.get()) / 2

        structural_nodes = []
        elastic_axis_lines = []
        cp_positions = []
        vectors_to_cp = []

        for i in range(0, len(panel_coordinates) - num_chordwise_panels, num_chordwise_panels):

            y_current = panel_coordinates[i][0][1]
            y_next = panel_coordinates[i + num_chordwise_panels][0][1]

            z_current = y_current * np.tan(dihedral_angle_rad)
            z_next = y_next * np.tan(dihedral_angle_rad)

            chord_current = root_chord - (root_chord - tip_chord) * (y_current / span)
            chord_next = root_chord - (root_chord - tip_chord) * (y_next / span)
            x_leading_current = y_current * np.tan(sweep_angle_rad)
            x_leading_next = y_next * np.tan(sweep_angle_rad)

            left_node = [
                x_leading_current + (0.25 + elastic_axis_offset) * chord_current,
                y_current,
                z_current,
            ]
            right_node = [
                x_leading_next + (0.25 + elastic_axis_offset) * chord_next,
                y_next,
                z_next,
            ]
            structural_nodes.append(left_node)
            structural_nodes.append(right_node)

            elastic_axis_lines.append((left_node, right_node))

            cp_y = (y_current + y_next) / 2
            cp_z = (z_current + z_next) / 2
            cp_chord = (chord_current + chord_next) / 2
            cp_x = ((x_leading_current + x_leading_next) / 2) + 0.25 * cp_chord
            cp_positions.append([cp_x, cp_y, cp_z])

        last_panel = panel_coordinates[-num_chordwise_panels]

        y_current = last_panel[1][1]
        y_next = last_panel[2][1]

        z_current = y_current * np.tan(dihedral_angle_rad)
        z_next = y_next * np.tan(dihedral_angle_rad)

        chord_current = root_chord - (root_chord - tip_chord) * (y_current / span)
        chord_next = root_chord - (root_chord - tip_chord) * (y_next / span)
        x_leading_current = y_current * np.tan(sweep_angle_rad)
        x_leading_next = y_next * np.tan(sweep_angle_rad)

        left_node_last = [
            x_leading_current + (0.25 + elastic_axis_offset) * chord_current,
            y_current,
            z_current,
        ]
        right_node_last = [
            x_leading_next + (0.25 + elastic_axis_offset) * chord_next,
            y_next,
            z_next,
        ]
        structural_nodes.append(left_node_last)
        structural_nodes.append(right_node_last)

        elastic_axis_lines.append((left_node_last, right_node_last))

        cp_y_last = (y_current + y_next) / 2
        cp_z_last = (z_current + z_next) / 2
        cp_chord_last = (chord_current + chord_next) / 2
        cp_x_last = ((x_leading_current + x_leading_next) / 2) + 0.25 * cp_chord_last
        cp_positions.append([cp_x_last, cp_y_last, cp_z_last])

        structural_nodes = np.array(structural_nodes)
        structural_nodes = np.unique(structural_nodes, axis=0)
        cp_positions = np.array(cp_positions)

        for panel in panel_coordinates:
            for i, corner in enumerate(panel):
                x, y, z = corner
                z += y * np.tan(dihedral_angle_rad)
                panel[i] = (x, y, z)


        for i in range(len(cp_positions)):
            cp = cp_positions[i]

            left_node = structural_nodes[i]
            right_node = structural_nodes[i + 1]

            vector_left_to_cp = [
                cp[0] - left_node[0],
                cp[1] - left_node[1],
                cp[2] - left_node[2],
            ]
            vector_right_to_cp = [
                cp[0] - right_node[0],
                cp[1] - right_node[1],
                cp[2] - right_node[2],
            ]

            vectors_to_cp.append((vector_left_to_cp, vector_right_to_cp))

        vectors_to_cp = np.array(vectors_to_cp)

        self.figure.clf()
        ax = self.figure.add_subplot(111, projection='3d')

        for i, panel in enumerate(panel_coordinates):
            x_coords, y_coords, z_coords = zip(*panel)
            ax.plot(
                [x_coords[0], x_coords[1], x_coords[2], x_coords[3], x_coords[0]],
                [y_coords[0], y_coords[1], y_coords[2], y_coords[3], y_coords[0]],
                color='blue', label="Aerodynamic Panels" if i == 0 else ""
            )

        for i in range(len(structural_nodes) - 1):
            ax.plot(
                [structural_nodes[i][0], structural_nodes[i + 1][0]],
                [structural_nodes[i][1], structural_nodes[i + 1][1]],
                [structural_nodes[i][2], structural_nodes[i + 1][2]],
                color='green', linestyle="--", label="Beam" if i == 0 else ""
            )

        for i, vectors in enumerate(vectors_to_cp):
            cp = cp_positions[i]
            left_node = structural_nodes[i]
            right_node = structural_nodes[i + 1]

            vector_left_endpoint = [
                left_node[0] + vectors[0][0],
                left_node[1] + vectors[0][1],
                left_node[2] + vectors[0][2],
            ]
            vector_right_endpoint = [
                right_node[0] + vectors[1][0],
                right_node[1] + vectors[1][1],
                right_node[2] + vectors[1][2],
            ]

            ax.plot(
                [left_node[0], vector_left_endpoint[0]],
                [left_node[1], vector_left_endpoint[1]],
                [left_node[2], vector_left_endpoint[2]],
                color='orange',
                label='Vector to CP (Left)' if i == 0 else ""
            )

            ax.plot(
                [right_node[0], vector_right_endpoint[0]],
                [right_node[1], vector_right_endpoint[1]],
                [right_node[2], vector_right_endpoint[2]],
                color='purple',
                label='Vector to CP (Right)' if i == 0 else "")

        ax.scatter(structural_nodes[:, 0], structural_nodes[:, 1], structural_nodes[:, 2], color='black', label="Nodes")
        ax.scatter(cp_positions[:, 0], cp_positions[:, 1], cp_positions[:, 2], color='magenta', label="CP Positions")

        ax.set_title("Structural Mesh")
        ax.set_xlabel("X (m)")
        ax.set_ylabel("Y (m)")
        ax.set_zlabel("Z (m)")
        ax.legend()

        self.canvas.draw()

        self.structural_nodes = structural_nodes
        self.cp_positions = cp_positions
        self.vectors_to_cp = vectors_to_cp

    def transform_stiffness_matrix(self, stiffness_matrix, dihedral_angle, sweep_angle):

        """
        Transform the stiffness matrix to the global coordinate system.

        Arguments:
            Local stiffness matrix.
            Dihedral angle.
            Sweep angle.
        """

        phi = np.radians(dihedral_angle)
        Lambda = np.radians(sweep_angle)

        T_sweep = np.array([
            [np.cos(Lambda), 0, np.sin(Lambda)],
            [0, 1, 0],
            [-np.sin(Lambda), 0, np.cos(Lambda)]
        ])

        T_dihedral = np.array([
            [1, 0, 0],
            [0, np.cos(phi), np.sin(phi)],
            [0, -np.sin(phi), np.cos(phi)]
        ])

        T = T_sweep @ T_dihedral

        T_expanded = np.zeros((12, 12))

        for i in range(4):

            T_expanded[i * 3:(i + 1) * 3, i * 3:(i + 1) * 3] = T

        transformed_matrix = T_expanded.T @ stiffness_matrix @ T_expanded

        return transformed_matrix

    def assemble_global_stiffness(self):

        """
        Assembles the global stiffness matrix for the FEM beam model.
        """

        self.build_structural_geometry()
        self.initialize_beam_properties()

        dihedral_angle = np.radians(float(self.geometry.dihedral_entry.get()))
        sweep_angle = np.radians(float(self.geometry.sweep_entry.get()))

        structural_nodes = self.structural_nodes
        n_nodes = len(structural_nodes)
        n_elements = n_nodes - 1

        total_dofs = 6 * n_nodes

        global_stiffness_matrix = np.zeros((total_dofs, total_dofs))

        for i in range(n_elements):

            node1, node2 = structural_nodes[i], structural_nodes[i + 1]
            element_length = np.linalg.norm(np.array(node2) - np.array(node1))

            element_stiffness_matrix = self.beam_element_stiffness_matrix(element_length)

            element_stiffness_matrix = self.transform_stiffness_matrix(
                element_stiffness_matrix, dihedral_angle, sweep_angle
            )

            element_dofs = [
                6 * i, 6 * i + 1, 6 * i + 2, 6 * i + 3, 6 * i + 4, 6 * i + 5,
                6 * (i + 1), 6 * (i + 1) + 1, 6 * (i + 1) + 2, 6 * (i + 1) + 3, 6 * (i + 1) + 4, 6 * (i + 1) + 5
            ]

            for row_local, row_global in enumerate(element_dofs):
                for col_local, col_global in enumerate(element_dofs):
                    global_stiffness_matrix[row_global, col_global] += element_stiffness_matrix[row_local, col_local]

        return global_stiffness_matrix

    def assemble_global_load_vector(self):

        """
        Assembles the global load vector for the FEM beam model using aerodynamic forces
        and optionally considers the wing's weight if gravity is enabled.
        """

        self.initialize_beam_properties()
        self.build_structural_geometry()

        phi = np.radians(float(self.geometry.dihedral_entry.get()))
        Lambda = np.radians(float(self.geometry.sweep_entry.get()))

        forces, _ = self.vlm.get_forces_bertin_vlm()
        lift_distribution = forces['spanwise_lift_distribution']
        drag_distribution = forces['spanwise_drag_distribution']

        structural_nodes = self.structural_nodes
        cp_positions = self.cp_positions
        r_cp_vectors = self.vectors_to_cp
        n_nodes = len(structural_nodes)
        n_elements = n_nodes - 1

        total_dofs = 6 * n_nodes
        global_load_vector = np.zeros(total_dofs)

        T_sweep = np.array([
            [np.cos(Lambda), 0, np.sin(Lambda)],
            [0, 1, 0],
            [-np.sin(Lambda), 0, np.cos(Lambda)]
        ])

        T_dihedral = np.array([
            [1, 0, 0],
            [0, np.cos(phi), np.sin(phi)],
            [0, -np.sin(phi), np.cos(phi)]
        ])

        T = T_sweep @ T_dihedral


        for i in range(n_elements):

            left_node = structural_nodes[i]
            right_node = structural_nodes[i + 1]
            cp = cp_positions[i]
            r_cp_left, r_cp_right = r_cp_vectors[i]

            lift = lift_distribution[i]
            drag = drag_distribution[i]
            force_local = np.array([drag, 0, lift])

            force_global = T @ force_local

            total_force = force_global

            force_left = 0.5 * total_force
            force_right = 0.5 * total_force

            moment_left = 0.5 * np.cross(r_cp_left, total_force)
            moment_right = 0.5 * np.cross(r_cp_right, total_force)

            left_dofs = [6 * i + j for j in range(6)]
            right_dofs = [6 * (i + 1) + j for j in range(6)]

            global_load_vector[left_dofs[0:3]] += force_left
            global_load_vector[left_dofs[3:6]] += moment_left

            global_load_vector[right_dofs[0:3]] += force_right
            global_load_vector[right_dofs[3:6]] += moment_right

        return global_load_vector

    def solve_displacements(self):

        """
        Solves the FEM system Ku = f for displacements and plots the deformed wing structure
        and the z-displacement over the span.
        """

        self.build_structural_geometry()
        self.initialize_beam_properties()

        global_stiffness_matrix = self.assemble_global_stiffness()
        global_load_vector = self.assemble_global_load_vector()

        fixed_dofs = [0, 1, 2, 3, 4, 5]

        for dof in fixed_dofs:
            global_stiffness_matrix[dof, :] = 0
            global_stiffness_matrix[:, dof] = 0
            global_stiffness_matrix[dof, dof] = 1
            global_load_vector[dof] = 0

        self.displacements = np.linalg.solve(global_stiffness_matrix, global_load_vector)

        displacements_str = "\n".join([f"DOF {i}: {value:.6f}" for i, value in enumerate(self.displacements)])

        self.output_field.delete("1.0", tk.END)

        self.output_field.insert(tk.END, "Displacements:\n")
        self.output_field.insert(tk.END, displacements_str)

        return self.displacements

    def calculate_strains(self):

        """
        Calculates the strain components for each beam element in the structural mesh.

        Returns: List of strain components for each element [[epsilon_xx, epsilon_yy, epsilon_zz, gamma_xy, gamma_xz, gamma_yz], ...]
        """

        self.strains = []

        if not hasattr(self, 'displacements'):
            raise ValueError("Displacements not calculated. Solve for displacements first.")

        n_nodes = len(self.structural_nodes)

        translational_displacements = self.displacements[:3*n_nodes].reshape(-1, 3)

        for i in range(len(self.structural_nodes) - 1):

            node1 = self.structural_nodes[i]
            node2 = self.structural_nodes[i + 1]

            disp1 = translational_displacements[i]
            disp2 = translational_displacements[i + 1]

            dx = node2[0] - node1[0]
            dy = node2[1] - node1[1]
            dz = node2[2] - node1[2]

            L = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
            if L == 0:
                raise ValueError(f"Element length is zero for element {i + 1}.")

            du = disp2[0] - disp1[0]
            dv = disp2[1] - disp1[1]
            dw = disp2[2] - disp1[2]
            epsilon_xx = du / dx if dx != 0 else 0
            epsilon_yy = dv / dy if dy != 0 else 0
            epsilon_zz = dw / dz if dz != 0 else 0
            gamma_xy = (du / dy if dy != 0 else 0) + (dv / dx if dx != 0 else 0)
            gamma_xz = (du / dz if dz != 0 else 0) + (dw / dx if dx != 0 else 0)
            gamma_yz = (dv / dz if dz != 0 else 0) + (dw / dy if dy != 0 else 0)

            self.strains.append([epsilon_xx, epsilon_yy, epsilon_zz, gamma_xy, gamma_xz, gamma_yz])

        strain_output = "Strain Components for Each Beam Element:\n"
        strain_output += "Element\t epsilon_xx\t epsilon_yy\t epsilon_zz\t gamma_xy\t gamma_xz\t gamma_yz\n"
        for idx, strain in enumerate(self.strains):
            strain_output += f"{idx + 1}\t {strain[0]:.6e}\t {strain[1]:.6e}\t {strain[2]:.6e}\t {strain[3]:.6e}\t {strain[4]:.6e}\t {strain[5]:.6e}\n"

        self.output_field.delete(1.0, tk.END)
        self.output_field.insert(tk.END, strain_output)

        return self.strains

    def fetch_material_properties(self):

        """
        Fetches material properties based on the selected source.
        """

        if self.material_source.get() == "manual":
            try:
                E = float(self.youngs_modulus_entry.get())
                nu = float(self.poisson_ratio_entry.get())
                G = float(self.shear_modulus_entry.get())
                return E, nu, G
            except ValueError:
                raise ValueError("Please enter valid material properties for manual input.")
        elif self.material_source.get() == "clt":
            laminate_properties = self.clt.get_laminate_properties()
            return laminate_properties['E_x'], laminate_properties['nu_xy'], laminate_properties['G_xy']

    def calculate_stresses(self):

        """
        Calculates the stress components for each beam element using the computed strains.
        """

        if not hasattr(self, 'strains'):
            raise ValueError("Strains not calculated. Solve for strains first.")

        try:

            E, nu, G = self.fetch_material_properties()

            D = np.array([
                [E / (1 - nu ** 2), nu * E / (1 - nu ** 2), 0],
                [nu * E / (1 - nu ** 2), E / (1 - nu ** 2), 0],
                [0, 0, G]
            ])

            self.stresses = []

            for strain in self.strains:
                epsilon = np.array([strain[0], strain[1], strain[3]])
                sigma = np.dot(D, epsilon)
                sigma_xx, sigma_yy, tau_xy = sigma
                self.stresses.append([sigma_xx, sigma_yy, tau_xy])

            stress_output = "Stress Components for Each Beam Element:\n"
            stress_output += "Element\t sigma_xx\t sigma_yy\t tau_xy\n"
            for idx, stress in enumerate(self.stresses):
                stress_output += f"{idx + 1}\t {stress[0]:.6e}\t {stress[1]:.6e}\t {stress[2]:.6e}\n"

            self.output_field.delete(1.0, tk.END)
            self.output_field.insert(tk.END, stress_output)

        except ValueError as e:
            self.output_field.delete(1.0, tk.END)
            self.output_field.insert(tk.END, f"Error: {str(e)}\n")

        return self.stresses

    def plot_strains(self):

        """
        Plots all strain components over the span of the wing.
        """

        if not hasattr(self, 'strains'):
            self.output_field.insert(tk.END, "Strains not calculated. Please calculate strains first.\n")
            return

        self.figure.clf()
        ax = self.figure.add_subplot(111)
        ax.set_title("Spanwise Strain Distribution")

        span_positions = [node[1] for node in self.structural_nodes[:-1]]

        epsilon_xx = [strain[0] for strain in self.strains]
        epsilon_yy = [strain[1] for strain in self.strains]
        gamma_xy = [strain[3] for strain in self.strains]

        ax.plot(span_positions, epsilon_xx, linestyle="-", label="epsilon_xx", color="blue")
        ax.plot(span_positions, epsilon_yy, linestyle="--", label="epsilon_yy", color="green")
        ax.plot(span_positions, gamma_xy, linestyle="-.", label="gamma_xy", color="orange")

        ax.set_xlabel("y (m)")
        ax.set_ylabel("Strain")
        ax.grid(True)
        ax.legend()

        self.canvas.draw()

        if self.export_option.get():
            if self.export_folder:
                filepath = f"{self.export_folder}/strain_distribution_all.png"
                self.figure.savefig(filepath)
                self.output_field.insert(tk.END, f"All strain components plot exported to: {filepath}\n")
            else:
                self.output_field.insert(tk.END, "No folder selected for exporting plots. Please select a folder.\n")

    def plot_stresses(self):

        """
        Plots all stress components over the span of the wing.
        """

        if not hasattr(self, 'stresses'):
            self.output_field.insert(tk.END, "Stresses not calculated. Please calculate stresses first.\n")
            return

        self.figure.clf()
        ax = self.figure.add_subplot(111)
        ax.set_title("Spanwise Stress Distribution")

        span_positions = [node[1] for node in self.structural_nodes[:-1]]

        sigma_xx = [stress[0] for stress in self.stresses]
        sigma_yy = [stress[1] for stress in self.stresses]
        tau_xy = [stress[2] for stress in self.stresses]

        ax.plot(span_positions, sigma_xx, linestyle="-", label="sigma_xx", color="blue")
        ax.plot(span_positions, sigma_yy, linestyle="--", label="sigma_yy", color="green")
        ax.plot(span_positions, tau_xy, linestyle="-.", label="tau_xy", color="orange")

        ax.set_xlabel("y (m)")
        ax.set_ylabel("Stress (MPa)")
        ax.grid(True)
        ax.legend()

        self.canvas.draw()

        if self.export_option.get():
            if self.export_folder:
                filepath = f"{self.export_folder}/stress_distribution_all.png"
                self.figure.savefig(filepath)
                self.output_field.insert(tk.END, f"All stress components plot exported to: {filepath}\n")
            else:
                self.output_field.insert(tk.END, "No folder selected for exporting plots. Please select a folder.\n")

    def plot_displaced_structure(self):

        """
        Plots the displaced structure on the existing canvas.
        """

        self.figure.clf()
        ax = self.figure.add_subplot(111, projection="3d")
        ax.set_title("Displaced Wing")

        x, y, z = zip(*self.structural_nodes)
        ax.plot(x, y, z, label="Original", color="blue", linestyle="-")

        n_nodes = len(self.structural_nodes)

        deformed_nodes = self.structural_nodes + self.displacements.reshape(n_nodes, 6)[:, :3]
        x_deformed, y_deformed, z_deformed = zip(*deformed_nodes)

        ax.plot(x_deformed, y_deformed, z_deformed, label="Deformed Structure", color="red")

        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.grid(True)
        ax.legend()

        self.canvas.draw()

    def plot_x_displacement(self):

        """
        Plots the Z-displacement along the span of the wing.
        """

        self.figure.clf()
        ax = self.figure.add_subplot(111)
        ax.set_title("Spanwise X-Displacement")

        self.u = self.displacements[0::6]
        span_positions = [node[1] for node in self.structural_nodes]

        ax.plot(span_positions, self.u, linestyle="-", label="VLM-FEM", color="blue")
        ax.set_xlabel("y (m)")
        ax.set_ylabel("x (m)")
        ax.grid(True)

        ax.legend()

        self.canvas.draw()

        if self.export_option.get():
            if self.export_folder:
                filepath = f"{self.export_folder}/x_displacement.png"
                self.figure.savefig(filepath)
                self.output_field.insert(tk.END, f"X-displacement plot exported to: {filepath}\n")
            else:
                self.output_field.insert(tk.END, "No folder selected for exporting plots. Please select a folder.\n")

    def plot_z_displacement(self):

        """
        Plots the Z-displacement along the span of the wing.
        """

        self.figure.clf()
        ax = self.figure.add_subplot(111)
        ax.set_title("Spanwise Z-Displacement")

        self.w = self.displacements[2::6]
        span_positions = [node[1] for node in self.structural_nodes]

        ax.plot(span_positions, self.w, linestyle="-", label="VLM-FEM", color="blue")
        ax.set_xlabel("y (m)")
        ax.set_ylabel("z (m)")
        ax.legend()
        ax.grid(True)

        self.canvas.draw()

        if self.export_option.get():
            if self.export_folder:
                filepath = f"{self.export_folder}/z_displacement.png"
                self.figure.savefig(filepath)
                self.output_field.insert(tk.END, f"Z-displacement plot exported to: {filepath}\n")
            else:
                self.output_field.insert(tk.END, "No folder selected for exporting plots. Please select a folder.\n")

    def transfer_displacements_to_cp(self):

        """
        Transfers the computed displacements from the structural nodes to the aerodynamic nodes.
        """

        displacements = self.displacements

        translational_displacements = displacements[0::6]
        translational_displacements = np.vstack((
            translational_displacements,
            displacements[1::6],
            displacements[2::6]
        )).T

        rotational_displacements = displacements[3::6]
        rotational_displacements = np.vstack((
            rotational_displacements,
            displacements[4::6],
            displacements[5::6]
        )).T

        aerodynamic_displacements = []

        for i in range(len(self.cp_positions)):
            r_left, r_right = self.vectors_to_cp[i]

            d_left = translational_displacements[i]
            d_right = translational_displacements[i + 1]
            theta_left = rotational_displacements[i]
            theta_right = rotational_displacements[i + 1]

            u_a = 0.5 * (
                    d_left + np.cross(theta_left, r_left) +
                    d_right + np.cross(theta_right, r_right)
            )

            aerodynamic_displacements.append(u_a)

        return np.array(aerodynamic_displacements)

    def plot_twist(self):

        """
        Computes and plots the twist distribution along the span of the wing
        using rotational displacements around the y-axis.
        """

        span_positions = [node[1] for node in self.structural_nodes]
        theta_y = self.displacements[4::6]

        aero_disp = self.transfer_displacements_to_cp()

        self.twist = theta_y * (180/np.pi)

        self.figure.clf()
        ax = self.figure.add_subplot(111)
        ax.set_title("Spanwise Twist Angle")
        ax.plot(span_positions, self.twist, linestyle="-", label="VLM-FEM", color="blue")
        ax.set_xlabel("Span Position y (m)")
        ax.set_ylabel("Twist Angle (deg)")
        ax.grid(True)
        ax.legend()
        self.canvas.draw()

        if self.export_option.get():
            if self.export_folder:
                filepath = f"{self.export_folder}/twist.png"
                self.figure.savefig(filepath)
                self.output_field.insert(tk.END, f"Twist plot exported to: {filepath}\n")
            else:
                self.output_field.insert(tk.END, "No folder selected for exporting plots. Please select a folder.\n")

    def get_input_values(self):

        """
        Retrieve input values from GUI fields for FEM parameters.

        Returns: Dictionary containing the FEM input values.
        """

        try:
            beam_pos = float(self.beam_pos_entry.get())
            EA = float(self.EA_entry.get())
            GJ = float(self.GJ_entry.get())
            EIy = float(self.EIy_entry.get())
            EIz = float(self.EIz_entry.get())
            E = float(self.youngs_modulus_entry.get())
            nu = float(self.poisson_ratio_entry.get())
            G = float(self.shear_modulus_entry.get())

            return {
                'shear_center': beam_pos,
                'EA': EA,
                'GJ': GJ,
                'EIy': EIy,
                'EIz': EIz,
                'E': E,
                'nu': nu,
                'G': G
            }
        except ValueError as e:
            print(f"Error getting FEM input values: {e}")
            return None

    def set_input_values(self, input_values):

        """
        Set input values in GUI fields for FEM parameters.

        Parameters:
            input_values: Dictionary containing FEM input values.
        """

        if 'shear_center' in input_values:
            self.beam_pos_entry.delete(0, tk.END)
            self.beam_pos_entry.insert(0, input_values['shear_center'])

        if 'EA' in input_values:
            self.EA_entry.delete(0, tk.END)
            self.EA_entry.insert(0, input_values['EA'])

        if 'GJ' in input_values:
            self.GJ_entry.delete(0, tk.END)
            self.GJ_entry.insert(0, input_values['GJ'])

        if 'EIy' in input_values:
            self.EIy_entry.delete(0, tk.END)
            self.EIy_entry.insert(0, input_values['EIy'])

        if 'EIz' in input_values:
            self.EIz_entry.delete(0, tk.END)
            self.EIz_entry.insert(0, input_values['EIz'])

        if 'E' in input_values:
            self.youngs_modulus_entry.delete(0, tk.END)
            self.youngs_modulus_entry.insert(0, input_values['E'])

        if 'nu' in input_values:
            self.poisson_ratio_entry.delete(0, tk.END)
            self.poisson_ratio_entry.insert(0, input_values['nu'])

        if 'G' in input_values:
            self.shear_modulus_entry.delete(0, tk.END)
            self.shear_modulus_entry.insert(0, input_values['G'])

    def get_results(self):

        """
        Returns the structural results computed by the FEM class, including displacements, twist, strains, and stresses.
        """

        spanwise_positions = [node[1] for node in self.structural_nodes]
        displacements = self.displacements
        twist_deg = self.twist

        n_nodes = len(spanwise_positions)
        n_elements = len(self.strains)
        node_strains = np.zeros((n_nodes, 6))
        node_stresses = np.zeros((n_nodes, 3))

        for i in range(n_elements):

            node_strains[i] += np.array(self.strains[i])
            node_strains[i + 1] += np.array(self.strains[i])

            node_stresses[i] += np.array(self.stresses[i])
            node_stresses[i + 1] += np.array(self.stresses[i])


        node_strains[1:-1] /= 2
        node_stresses[1:-1] /= 2

        return {
            "spanwise_positions": spanwise_positions,
            "displacements": displacements,
            "twist": twist_deg,
            "strains": node_strains.tolist(),
            "stresses": node_stresses.tolist(),
        }



