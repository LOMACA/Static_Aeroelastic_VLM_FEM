##################################################################

# VLM Implementation following Bertin's description

#################################################################

import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog

class VLM(ttk.Frame):

    """
    This class implements a VLM based on Bertin's description in "Aerodynamics for Engineers"
    """

    def __init__(self, parent, geometry, freestream, *args, **kwargs):

        """
        Initializes the VLM class, linking the geometry and freestream classes and setting up the GUI.
        """

        super().__init__(parent, *args, **kwargs)

        self.geometry = geometry
        self.freestream = freestream

        self.panel_coordinates = None
        self.vortex_points = None
        self.control_points = None
        self.velocity = None
        self.density = None
        self.alpha = None
        self.normal_vectors = None
        self.AIC_matrix = None
        self.Gamma = None
        self.rhs = None
        self.AIC_matrix_alt = None

        button_frame = tk.Frame(self)
        button_frame.grid(row=0, column=0, padx=5, pady=5, sticky="nsw")

        self.import_geometry_button = ttk.Button(button_frame, text="Import Geometry Data",
                                                 command=self.import_geometry_data)
        self.import_geometry_button.grid(row=0, column=0, padx=5, pady=5, sticky="w")

        self.import_freestream_button = ttk.Button(button_frame, text="Import Freestream Data",
                                                   command=self.import_freestream_data)
        self.import_freestream_button.grid(row=1, column=0, padx=5, pady=5, sticky="w")

        self.compute_freestream_button = ttk.Button(button_frame, text="Compute Velocity Vector",
                                                    command=self.compute_freestream_velocity)
        self.compute_freestream_button.grid(row=2, column=0, padx=5, pady=5, sticky="w")

        self.status_label = ttk.Label(button_frame, text="Status: Awaiting data import", font=("Arial", 10, 'italic'))
        self.status_label.grid(row=3, column=0, padx=5, pady=5, sticky="w")

        self.output_label = ttk.Label(button_frame, text="Status: Awaiting data import", font=("Arial", 10, 'italic'))
        self.output_label.grid(row=4, column=0, padx=5, pady=5, sticky="w")

        #self.compute_normals_button = ttk.Button(button_frame, text="Compute Normal Vectors",
                                                 #command=self.visualize_normal_vectors)
        #self.compute_normals_button.grid(row=5, column=0, padx=5, pady=5, sticky="w")

        self.compute_circulation_button = ttk.Button(button_frame, text="Compute Circulation",
                                                     command=self.compute_circulation)
        self.compute_circulation_button.grid(row=5, column=0, padx=5, pady=5, sticky="w")

        self.plot_circulation_button = ttk.Button(button_frame, text="Plot Circulation",
                                                  command=self.plot_circulation_over_span)
        self.plot_circulation_button.grid(row=6, column=0, padx=5, pady=5, sticky="w")

        self.compute_button = ttk.Button(button_frame, text="Compute Aerodynamic Forces",
                                         command=self.calculate_and_display_aerodynamics)
        self.compute_button.grid(row=7, column=0, padx=5, pady=5, sticky="w")

        self.plot_CL_button = ttk.Button(button_frame, text="Plot Spanwise Lift Coefficient",
                                         command=self.plot_panel_CL_over_span)
        self.plot_CL_button.grid(row=8, column=0, padx=5, pady=5, sticky="w")

        self.collect_force_button = ttk.Button(button_frame, text="Collect Aerodynamic Forces",
                                               command=self.get_forces_bertin_vlm)
        self.collect_force_button.grid(row=9, column=0, padx=5, pady=5, sticky="w")

        self.export_folder = None

        self.export_option = tk.BooleanVar()
        self.export_checkbox = ttk.Checkbutton(button_frame, text="Export Plots", variable=self.export_option,
                                               command=self.toggle_export_folder)
        self.export_checkbox.grid(row=10, column=0, padx=5, pady=5, sticky="w")

        self.select_folder_button = ttk.Button(button_frame, text="Select Export Folder",
                                               command=self.select_export_folder)
        self.select_folder_button.grid(row=11, column=0, padx=5, pady=5, sticky="w")
        self.select_folder_button.config(state=tk.DISABLED)

        text_frame = tk.Frame(self)
        text_frame.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

        self.scrollbar = tk.Scrollbar(text_frame)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.text_field = tk.Text(text_frame, height=55, width=80, yscrollcommand=self.scrollbar.set)
        self.text_field.pack(side=tk.LEFT, fill=tk.BOTH)

        self.scrollbar.config(command=self.text_field.yview)

        self.plot_frame = tk.Frame(self)
        self.plot_frame.grid(row=0, column=2, padx=5, pady=5, sticky="nsew")

        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=2)
        self.grid_rowconfigure(0, weight=1)

        self.figure = plt.figure(figsize=(7, 6))
        self.canvas = FigureCanvasTkAgg(self.figure, self.plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH)


    def import_geometry_data(self):

        """
        Imports the geometry data from the Geometry class. This includes panel coordinates, vortex points,
        control points, and trailing vortex points.
        """

        self.panel_coordinates = None
        self.vortex_points = None
        self.control_points = None

        wing_geometry_data = self.geometry.build_aerodynamic_mesh()

        if wing_geometry_data:

            self.panel_coordinates = wing_geometry_data.get('panel_coordinates', None)
            self.vortex_points = wing_geometry_data.get('vortex_points', None)
            self.control_points = wing_geometry_data.get('control_points', None)

            # if self.panel_coordinates is not None:
            #     mirrored_panel_coordinates = []
            #     for panel in self.panel_coordinates:
            #         mirrored_panel = [(x, -y, z) for (x, y, z) in panel]
            #         mirrored_panel_coordinates.append(mirrored_panel)
            #     self.panel_coordinates += mirrored_panel_coordinates

            # if self.vortex_points is not None:
            #     mirrored_vortex_points = [(x, -y, z) for (x, y, z) in self.vortex_points]
            #     self.vortex_points += mirrored_vortex_points

            # if self.control_points is not None:
            #     mirrored_control_points = [(x, -y, z) for (x, y, z) in self.control_points]
            #     self.control_points += mirrored_control_points

            if self.vortex_points and self.control_points:
                self.status_label.config(text="Geometry data imported successfully.", foreground="green")
            else:
                self.status_label.config(text="Error: Geometry data missing.", foreground="red")
        else:
            self.status_label.config(text="Error: Geometry data not available.", foreground="red")

    def import_freestream_data(self):

        """
        Imports freestream velocity, density, and angle of attack from the freestream class.
        """

        freestream_data = self.freestream.get_submitted_data()

        if freestream_data:

            self.velocity = freestream_data.get('velocity', None)
            self.density = freestream_data.get('density', None)
            self.alpha = freestream_data.get('alpha', None)

            if self.velocity is not None and self.density is not None and self.alpha is not None:
                self.status_label.config(text="Freestream data imported successfully!", foreground="green")
            else:
                self.status_label.config(text="Error: Freestream data missing.", foreground="red")
        else:
            self.status_label.config(text="Error: Freestream data not submitted.", foreground="red")

    def compute_freestream_velocity(self):

        """
        Computes the freestream velocity vector based on the given velocity magnitude and angle of attack.

        Returns: The freestream velocity vector.
        """

        if self.velocity is None or self.alpha is None:
            self.output_label.config(text="Please import freestream data first.")
            return

        V_inf = self.velocity
        alpha = self.alpha

        V_x = V_inf * np.cos(alpha)
        V_z = V_inf * np.sin(alpha)

        self.freestream_velocity = np.array([V_x, 0, V_z])
        #self.output_label.config(text=f"Freestream Velocity Vector: {self.freestream_velocity}")
        return self.freestream_velocity

    def compute_normal_vectors(self):

        """
        Computes normal vectors for each panel in the geometry, to be used in aerodynamic influence coefficient (AIC)
        calculations.

        Returns: List of normal vectors for each panel.
        """

        self.normal_vectors = []

        for panel in self.panel_coordinates:

            v1 = np.array(panel[1]) - np.array(panel[0])
            v2 = np.array(panel[3]) - np.array(panel[0])

            normal_vector = np.cross(v1, v2)

            normal_vector_norm = np.linalg.norm(normal_vector)

            if normal_vector_norm > 1e-8:
                normal_vector /= normal_vector_norm
            else:
                normal_vector = np.zeros_like(normal_vector)

            self.normal_vectors.append(normal_vector)

        return self.normal_vectors

    def visualize_normal_vectors(self):

        """
        Visualizes the normal vectors on a 3D plot overlaid on the wing panels.
        """

        self.normal_vectors = self.compute_normal_vectors()

        if self.normal_vectors is None:
            print("No panel or normal vector data available for visualization.")
            return

        plt.close('all')

        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')

        for i, panel in enumerate(self.panel_coordinates):
            x_coords, y_coords, z_coords = zip(*panel)

            x_coords += (x_coords[0],)
            y_coords += (y_coords[0],)
            z_coords += (z_coords[0],)

            ax.plot(x_coords, y_coords, z_coords, color='b')

            control_point = np.mean(panel, axis=0)
            normal_vector = self.normal_vectors[i]

            ax.quiver(
                control_point[0], control_point[1], control_point[2],
                normal_vector[0], normal_vector[1], normal_vector[2],
                color='r', length=0.1, normalize=True
            )

        ax.set_title("Wing Geometry with Normal Vectors")
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Axis')
        ax.set_zlabel('Z Axis')

        plt.show()

    def compute_induced_velocity_vortex_segment(self, P0, P1, P3, gamma):

        """
        Computes the velocity induced by the bound vortex at a control point (eq. 7.38 a)

        Arguments:
            P0, P1, P3: Coordinates of the control point, and start/end points of the vortex line.
            gamma: unit circulation strength

        Returns: Induced velocity.
        """

        P0 = np.array(P0)
        P1 = np.array(P1)
        P3 = np.array(P3)

        r1 = P0 - P1
        r2 = P0 - P3
        r0 = P3 - P1

        r1_cross_r2 = np.cross(r1, r2)
        norm_r1_cross_r2_squared = np.linalg.norm(r1_cross_r2) ** 2

        if norm_r1_cross_r2_squared < 1e-8:
            return np.array([0.0, 0.0, 0.0])

        Fac1_AB = r1_cross_r2 / norm_r1_cross_r2_squared

        norm_r1 = np.linalg.norm(r1)
        norm_r2 = np.linalg.norm(r2)

        Fac2_AB = np.dot(r0, (r1 / norm_r1)) - np.dot(r0, (r2 / norm_r2))

        induced_velocity = (gamma / (4 * np.pi)) * Fac1_AB * Fac2_AB

        return induced_velocity

    def compute_inner_trailing_vortex_induced_velocity(self, gamma, P0, P1):

        """
        Computes the induced velocity from the inner trailing vortex segment (eq. 7.38 b).

        Arguments:
            P0, P1: Coordinates of the control point and vortex point.
            gamma: unit circulation strength

        Returns: Induced velocity.
        """

        x, y, z = P0
        x1, y1, z1 = P1

        dx = x - x1
        dy = y - y1
        dz = z - z1

        r_squared = dx ** 2 + dy ** 2 + dz ** 2
        r_xy_squared = (y1 - y) ** 2 + dz ** 2

        if r_squared < 1e-8 or r_xy_squared < 1e-8:
            return np.array([0.0, 0.0, 0.0])

        factor1 = ((dz) + (y1 - y)) / r_xy_squared
        factor2 = 1 + (dx / np.sqrt(r_squared))

        induced_velocity = (gamma / (4 * np.pi)) * factor1 * factor2

        return induced_velocity

    def compute_outer_trailing_vortex_induced_velocity(self, gamma, P0, P3):

        """
        Computes the induced velocity from the outer trailing vortex segment (eq. 7.38 c)

        Arguments:
            P0, P3: Coordinates of the control point and vortex point.
            gamma: unit circulation strength

        Returns: Induced velocity.
        """

        x, y, z = P0
        x2, y2, z2 = P3

        dx = x - x2
        dy = y - y2
        dz = z - z2

        r_squared = dx ** 2 + dy ** 2 + dz ** 2
        r_xy_squared = (y2 - y) ** 2 + dz ** 2

        if r_squared < 1e-8 or r_xy_squared < 1e-8:
            return np.array([0.0, 0.0, 0.0])

        factor1 = (dz + (y2 - y)) / r_xy_squared
        factor2 = 1 + (dx / np.sqrt(r_squared))

        induced_velocity = -(gamma / (4 * np.pi)) * factor1 * factor2

        return induced_velocity

    def compute_induced_velocities_at_control_points(self):

        """
        Calculates the total induced velocity at each control point due to all vortex segments, including the mirrored half-wing.

        Returns: Induced velocities in y and z direction, and downwash.
        """

        num_panels = len(self.control_points)
        induced_velocities_y = np.zeros((num_panels, num_panels))
        induced_velocities_z = np.zeros((num_panels, num_panels))
        downwash = []

        table_data = []

        for m, P0 in enumerate(self.control_points):
            total_downwash = 0

            for n in range(num_panels):
                P1 = self.vortex_points[n][0]
                P3 = self.vortex_points[n][1]

                gamma = 1

                velocity_bound_vortex_right = self.compute_induced_velocity_vortex_segment(P0, P1, P3, gamma)
                velocity_inner_trailing_vortex_right = self.compute_inner_trailing_vortex_induced_velocity(gamma, P0,
                                                                                                           P1)
                velocity_outer_trailing_vortex_right = self.compute_outer_trailing_vortex_induced_velocity(gamma, P0,
                                                                                                           P3)

                total_velocity_right = velocity_bound_vortex_right + velocity_inner_trailing_vortex_right + velocity_outer_trailing_vortex_right

                P1_mirrored = np.array([P3[0], -P3[1], P3[2]])
                P3_mirrored = np.array([P1[0], -P1[1], P1[2]])

                velocity_bound_vortex_left = self.compute_induced_velocity_vortex_segment(P0, P1_mirrored, P3_mirrored,
                                                                                          gamma)

                velocity_outer_trailing_vortex_left = self.compute_outer_trailing_vortex_induced_velocity(gamma, P0,
                                                                                                          P3_mirrored)

                velocity_inner_trailing_vortex_left = self.compute_inner_trailing_vortex_induced_velocity(gamma, P0,
                                                                                                          P1_mirrored)
                total_velocity_left = velocity_bound_vortex_left + velocity_inner_trailing_vortex_left + velocity_outer_trailing_vortex_left
                total_induced_velocity = total_velocity_right + total_velocity_left

                induced_velocities_y[m, n] = total_induced_velocity[1]
                induced_velocities_z[m, n] = total_induced_velocity[2]

                right_downwash = velocity_inner_trailing_vortex_right + velocity_outer_trailing_vortex_right
                left_downwash = velocity_inner_trailing_vortex_left + velocity_outer_trailing_vortex_left
                total_downwash += right_downwash + left_downwash

                xm, ym, _ = P0.tolist()
                x1n, y1n, _ = P1.tolist()
                x2n, y2n, _ = P3.tolist()
                table_data.append(f"{m + 1:<10}{xm:>10.3f}{ym:>10.3f}{x1n:>10.3f}{y1n:>10.3f}{x2n:>10.3f}{y2n:>10.3f}")

            downwash.append(total_downwash)

        downwash = np.array(downwash)

        self.text_field.delete(1.0, tk.END)
        self.text_field.insert(tk.END,
                               f"{'Panel':<10}{'xm':>10}{'ym':>10}{'x1n':>10}{'y1n':>10}{'x2n':>10}{'y2n':>10}\n")
        self.text_field.insert(tk.END, "-" * 70 + "\n")
        self.text_field.insert(tk.END, "\n".join(table_data) + "\n")

        return induced_velocities_y, induced_velocities_z, downwash

    def compute_circulation(self):

        """
        Assembles the AIC matrix and solves for circulation strengths (Gamma) over each panel.

        Returns: Circulation (Gamma) array for each panel.
        """

        induced_velocities_y, induced_velocities_z, downwash = self.compute_induced_velocities_at_control_points()
        num_panels = len(self.control_points)
        self.AIC_matrix = np.zeros((num_panels, num_panels))

        dihedral_angle = np.radians(float(self.geometry.dihedral_entry.get()))

        for m in range(num_panels):
            for n in range(num_panels):
                aic_value = -induced_velocities_y[m, n] * np.sin(dihedral_angle) + \
                            induced_velocities_z[m, n] * np.cos(dihedral_angle)
                self.AIC_matrix[m, n] = aic_value

        alpha = np.radians(float(self.freestream.aoa_entry.get()))
        U_inf = float(self.freestream.velocity_entry.get())
        rhs = -U_inf * np.sin(alpha) * np.cos(dihedral_angle)

        try:
            self.Gamma = np.linalg.solve(self.AIC_matrix, np.full(num_panels, rhs))
        except np.linalg.LinAlgError as e:
            print("Error in solving for circulation strengths:", e)
            self.Gamma = None

        self.text_field.delete(1.0, tk.END)
        self.text_field.insert(tk.END, f"Circulation (Gamma):\n{self.Gamma}\n")

    def plot_circulation_over_span(self):

        """
        Plots the chord-wise summed circulation (Gamma) distribution along the span of the wing.

        Returns: Spanwise circulation distribution.
        """

        if self.Gamma is None:
            self.text_field.delete(1.0, tk.END)
            self.text_field.insert(tk.END, "Please compute circulation first.\n")
            return

        num_spanwise_panels = int(self.geometry.n_panels_span_entry.get()) - 1
        num_chordwise_panels = int(self.geometry.n_panels_chord_entry.get())

        total_panels = num_spanwise_panels * num_chordwise_panels

        if len(self.Gamma) != total_panels:
            self.text_field.delete(1.0, tk.END)
            self.text_field.insert(tk.END, "Mismatch between Gamma size and panel grid dimensions.\n")
            return

        Gamma_reshaped = self.Gamma.reshape((num_spanwise_panels, num_chordwise_panels))

        self.Gamma_spanwise = np.sum(Gamma_reshaped, axis=1)

        span = float(self.geometry.span_entry.get()) / 2
        y_stations = np.linspace(0, span, num_spanwise_panels)

        self.text_field.delete(1.0, tk.END)
        self.text_field.insert(tk.END, f"Spanwise Circulation (Gamma):\n{self.Gamma_spanwise}\n")

        self.figure.clf()
        ax = self.figure.add_subplot(111)
        ax.plot(y_stations, self.Gamma_spanwise, '-', label='VLM')
        ax.set_xlabel('y (m)', fontsize=12)
        ax.set_ylabel('Circulation (m^2/s)', fontsize=12)
        ax.tick_params(axis='both', labelsize=10)
        #ax.set_title('Spanwise Circulation Distribution')
        ax.grid(True)
        ax.legend(fontsize=11)

        self.canvas.draw()

        if self.export_option.get():
            if self.export_folder:
                filepath = f"{self.export_folder}/circulation_distribution.png"
                self.figure.savefig(filepath)
                self.text_field.insert(tk.END, f"Circulation distribution plot exported to: {filepath}\n")
            else:
                self.text_field.insert(tk.END, "No folder selected for exporting plots. Please select a folder.\n")

        return self.Gamma_spanwise

    def calculate_total_lift(self, rho_inf, U_inf, q_inf, dy, S):

        """
        Calculates the total lift and lift coefficient from the computed circulation values according to Bertin (eq. 7.49 and 7.50 b).

        Arguments:
            rho_inf: Freestream density.
            U_inf: Freestream velocity.
            q_inf: Dynamic pressure.
            dy: Spanwise widths of each panel.
            S: Total wing surface area.

        Returns: Lift per panel, total lift, and lift coefficient (CL).
        """

        panel_lift = rho_inf * U_inf * self.Gamma * dy

        total_lift = 2 * rho_inf * U_inf * np.sum(self.Gamma * dy)

        CL = total_lift / (q_inf * S)

        return panel_lift, total_lift, CL

    def compute_panel_properties(self):

        """
        Computes the spanwise width and area of each panel

        Returns: Array of panel widths (dy), total surface area (S), and panel areas.
        """

        dy = []
        panel_areas = []

        for coords in self.panel_coordinates:
            y_a = coords[0][1]
            y_b = coords[3][1]

            delta_y_n = abs(y_b - y_a)
            dy.append(delta_y_n)

            P1 = np.array(coords[0])
            P2 = np.array(coords[1])
            P3 = np.array(coords[2])
            P4 = np.array(coords[3])

            triangle_1_area = 0.5 * np.linalg.norm(np.cross(P2 - P1, P3 - P1))
            triangle_2_area = 0.5 * np.linalg.norm(np.cross(P3 - P1, P4 - P1))

            panel_area = triangle_1_area + triangle_2_area
            panel_areas.append(panel_area)

        dy = np.array(dy)
        panel_areas = np.array(panel_areas)

        S = np.sum(panel_areas) * 2

        return dy, S, panel_areas

    def calculate_and_display_aerodynamics(self):

        """
        Computes aerodynamic properties, including lift, induced drag, and pressure distribution, and
        displays results in the GUI text field.
        """

        induced_velocities_y, induced_velocities_z, downwash = self.compute_induced_velocities_at_control_points()

        dy, S, panel_areas = self.compute_panel_properties()

        rho_inf = float(self.freestream.density_entry.get())
        U_inf = float(self.freestream.velocity_entry.get())

        q_inf = 0.5 * rho_inf * U_inf ** 2

        panel_lift, L, CL = self.calculate_total_lift(rho_inf, U_inf, q_inf, dy, S)

        panel_CL = panel_lift / (q_inf * panel_areas)

        P = panel_lift / panel_areas

        pressure_differences, delta_lift = self.calculate_pressure_differences(rho_inf, U_inf, dy, panel_areas)

        Di, Di_panel = self.calculate_induced_drag(rho_inf, downwash, dy)

        Di = Di * 2

        self.text_field.delete(1.0, tk.END)
        self.text_field.insert(tk.END, f"Panel Widths: {dy}\n")
        self.text_field.insert(tk.END, f"Total Surface Area: {S:.3f} m^2\n")
        self.text_field.insert(tk.END, f"Freestream Density: {rho_inf:.3f} kg/m^3\n")
        self.text_field.insert(tk.END, f"Freestream Velocity: {U_inf:.3f} m/s\n")
        self.text_field.insert(tk.END, f"Dynamic Pressure: {q_inf:.3f} Pa\n")
        self.text_field.insert(tk.END, f"Panel Lift Forces: {panel_lift} N\n")
        self.text_field.insert(tk.END, f"Panel Lift Coefficients: {panel_CL}\n ")
        self.text_field.insert(tk.END, f"Total Lift: {L:.3f} N\n")
        self.text_field.insert(tk.END, f"Wing Lift Coefficient: {CL:.4f}\n")
        self.text_field.insert(tk.END, f"Panel Pressure: {P} Pa\n")
        self.text_field.insert(tk.END, f"Pressure Differences across Panels: {pressure_differences} Pa\n")
        self.text_field.insert(tk.END, f"Induced Drag: {Di:.3f} N\n")
        self.text_field.insert(tk.END, f"Panel Induced Drag: {Di_panel} N\n")

        return panel_lift, L, CL, panel_CL, Di, Di_panel

    def plot_panel_CL_over_span(self):

        """
        Plots the chordwise-summed panel lift coefficient (CL) distribution along the span of the wing.
        """

        panel_lift, L, CL, panel_CL, Di, Di_panel = self.calculate_and_display_aerodynamics()

        num_spanwise_panels = int(self.geometry.n_panels_span_entry.get()) - 1
        num_chordwise_panels = int(self.geometry.n_panels_chord_entry.get())

        panel_CL_reshaped = panel_CL.reshape((num_spanwise_panels, num_chordwise_panels))
        self.CL_spanwise = np.sum(panel_CL_reshaped, axis=1)

        span = float(self.geometry.span_entry.get()) / 2
        y_stations = np.linspace(0, span, num_spanwise_panels)

        self.figure.clf()
        ax = self.figure.add_subplot(111)
        ax.plot(y_stations, self.CL_spanwise, '-', label='VLM')
        #ax.set_title('Spanwise Lift Coefficient Distribution')
        ax.set_xlabel('y (m)')
        ax.set_ylabel('Lift Coefficient')
        ax.grid(True)
        ax.legend()

        self.canvas.draw()

        if self.export_option.get():
            if self.export_folder:
                filepath = f"{self.export_folder}/cl_distribution.png"
                self.figure.savefig(filepath)
                self.text_field.insert(tk.END, f"CL distribution plot exported to: {filepath}\n")
            else:
                self.text_field.insert(tk.END, "No folder selected for exporting plots. Please select a folder.\n")

        return self.CL_spanwise

    def calculate_pressure_differences(self, rho_inf, U_inf, dy, panel_areas):

        """
        Calculates pressure differences across panels and lift contribution from circulation change along the span (eq. 12.26 from Katz and Plotkin)

        Arguments:
            rho_inf: Freestream density.
            U_inf: Freestream velocity.
            dy: Panel widths.
            panel_areas: Panel areas.

        Returns: Pressure differences and delta lift across panels.
        """

        delta_gamma = np.zeros_like(self.Gamma)
        delta_lift = np.zeros_like(self.Gamma)

        Q_inf = 0.5 * rho_inf * U_inf ** 2

        delta_gamma[0] = self.Gamma[0]
        delta_lift[0] = Q_inf * self.Gamma[0] * dy[0]

        for i in range(1, len(self.Gamma)):
            delta_gamma[i] = self.Gamma[i] - self.Gamma[i - 1]

            delta_lift[i] = Q_inf * delta_gamma[i] * dy[i]

        pressure_differences = delta_lift / panel_areas

        return pressure_differences, delta_lift

    def calculate_induced_drag(self, rho_inf, downwash, dy):

        """
        Calculates the induced drag based on the obtained circulation from the lift after equation 12.27 by Katz and Plotkin,

        Low Speed Aerodynamics

        Returns: Panel induced drag array, total induced drag.
        """

        num_panels = len(self.Gamma)
        panel_induced_drag = np.zeros(num_panels)

        for i in range(num_panels):
            if i == 0:
                delta_drag = -rho_inf * downwash[i] * self.Gamma[i] * dy[i]
            else:
                delta_gamma = self.Gamma[i] - self.Gamma[i - 1]
                delta_drag = -rho_inf * downwash[i] * delta_gamma * dy[i]

            panel_induced_drag[i] = delta_drag

        total_induced_drag = np.sum(panel_induced_drag)

        return total_induced_drag, panel_induced_drag

    def get_forces_bertin_vlm(self):

        """
        Collects the summed spanwise aerodynamic forces in a dictionary to pass to the FEM class
        """

        num_chordwise_panels = int(self.geometry.n_panels_chord_entry.get())
        span = float(self.geometry.span_entry.get()) / 2
        num_spanwise_panels = int(self.geometry.n_panels_span_entry.get()) - 1

        self.y_stations = np.linspace(0, span, num_spanwise_panels)

        panel_lift, L, CL, panel_CL, Di, Di_panel = self.calculate_and_display_aerodynamics()

        panel_lift_reshaped = panel_lift.reshape((num_spanwise_panels, num_chordwise_panels))

        panel_lift_spanwise = np.sum(panel_lift_reshaped, axis=1)

        panel_drag_reshaped = Di_panel.reshape((num_spanwise_panels, num_chordwise_panels))

        panel_drag_spanwise = np.sum(panel_drag_reshaped, axis=1)

        self.L_panel_tot = panel_lift_spanwise
        self.Di_panel_tot = panel_drag_spanwise

        forces = {
            'spanwise_lift_distribution': self.L_panel_tot,
            'spanwise_drag_distribution': self.Di_panel_tot,
        }

        self.text_field.delete(1.0, tk.END)

        self.text_field.insert(tk.END, "Spanwise Stations (y) and Force Distributions:\n")
        self.text_field.insert(tk.END, "-----------------------------------------------\n")
        self.text_field.insert(tk.END, f"{'y':>10}  |  {'Lift Distribution':>18}  |  {'Drag Distribution':>18}\n")

        for i in range(len(self.y_stations)):
            self.text_field.insert(tk.END,
                                   f"{self.y_stations[i]:>10.3f}  |  {forces['spanwise_lift_distribution'][i]:>18.3f}  |  {forces['spanwise_drag_distribution'][i]:>18.3f}\n")

        return forces, self.y_stations

    def toggle_export_folder(self):

        """
        Enables or disables the folder selection button based on the checkbox state.
        """

        if self.export_option.get():
            self.select_folder_button.config(state=tk.NORMAL)
        else:
            self.select_folder_button.config(state=tk.DISABLED)

    def select_export_folder(self):

        """
        Opens a dialog to select the folder where exported plots will be saved.
        """

        folder = filedialog.askdirectory(title="Select Folder for Exported Plots")
        if folder:
            self.export_folder = folder
            self.text_field.insert(tk.END, f"Export folder selected: {folder}\n")

    def get_results(self):

        """
        Returns the aerodynamic results computed by the VLM class.
        """

        return {
            "spanwise_positions": self.y_stations,
            "circulation_distribution": self.Gamma_spanwise,
            "lift_distribution": self.L_panel_tot,
            "cl_distribution": self.CL_spanwise,
            "Di_distribution": self.Di_panel_tot,
        }


