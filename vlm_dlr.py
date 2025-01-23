#####################################################################

# VLM implementation following DLR script

#####################################################################

import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

class VLMDLR(ttk.Frame):

    """
    This class implements a VLM based on a description by the DLR titled
    "An Implementation of the Vortex Lattice Method and the Doublet Lattice Method"
    """

    def __init__(self, parent, geometry, freestream, *args, **kwargs):

        """
        Initializes the VLM class, linking geometry and freestream objects and setting up the GUI for control.
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
        self.total_downwash = None

        button_frame = tk.Frame(self)
        button_frame.grid(row=0, column=0, padx=5, pady=5, sticky="nsw")

        self.import_geometry_button = ttk.Button(button_frame, text="Import Geometry Data", command=self.import_geometry_data)
        self.import_geometry_button.grid(row=0, column=0, padx=5, pady=5, sticky="w")

        self.import_freestream_button = ttk.Button(button_frame, text="Import Freestream Data", command=self.import_freestream_data)
        self.import_freestream_button.grid(row=1, column=0, padx=5, pady=5, sticky="w")

        self.compute_freestream_button = ttk.Button(button_frame, text="Compute Velocity Vector", command=self.compute_freestream_velocity)
        self.compute_freestream_button.grid(row=2, column=0, padx=5, pady=5, sticky="w")

        self.status_label = ttk.Label(button_frame, text="Status: Awaiting data import", font=("Arial", 10, 'italic'))
        self.status_label.grid(row=3, column=0, padx=5, pady=5, sticky="w")

        self.compute_normals_button = ttk.Button(button_frame, text="Compute Normal Vectors", command=self.visualize_normal_vectors)
        self.compute_normals_button.grid(row=4, column=0, padx=5, pady=5, sticky="w")

        self.compute_distances_button = ttk.Button(button_frame, text="Compute Distance Bound Vortex to Control Point", command=self.display_distances_and_dihedral)
        self.compute_distances_button.grid(row=5, column=0, padx=5, pady=5, sticky="w")

        self.output_label = ttk.Label(button_frame, text="Status: Awaiting data import", font=("Arial", 10, 'italic'))
        self.output_label.grid(row=6, column=0, padx=5, pady=5, sticky="w")

        self.compute_vortex_induced_velocities_button = ttk.Button(button_frame, text="Compute Induced Velocities from Bound Vortex", command=self.display_induced_velocities_vortex_line)
        self.compute_vortex_induced_velocities_button.grid(row=7, column=0, padx=5, pady=5, sticky="w")

        self.compute_inner_induced_velocities_button = ttk.Button(button_frame, text="Compute Induced Velocities from Inner Horseshoe Vortex", command=self.display_inner_horseshoe_induced_velocities)
        self.compute_inner_induced_velocities_button.grid(row=8, column=0, padx=5, pady=5, sticky="w")

        self.compute_outer_induced_velocities_button = ttk.Button(button_frame, text="Compute Induced Velocities from Outer Horseshoe Vortex", command=self.display_outer_horseshoe_induced_velocities)
        self.compute_outer_induced_velocities_button.grid(row=9, column=0, padx=5, pady=5, sticky="w")

        self.compute_cell_area_button = ttk.Button(button_frame, text="Compute Panel Area", command=self.display_panel_areas)
        self.compute_cell_area_button.grid(row=10, column=0, padx=5, pady=5, sticky="w")

        self.compute_cell_span_button = ttk.Button(button_frame, text="Compute Panel Width", command=self.display_spanwise_widths)
        self.compute_cell_span_button.grid(row=11, column=0, padx=5, pady=5, sticky="w")

        self.compute_pressure_aic_button = ttk.Button(button_frame, text="Compute Pressure Mapped AIC Matrix", command=self.assemble_pressure_aic)
        self.compute_pressure_aic_button.grid(row=12, column=0, padx=5, pady=5, sticky="w")

        self.compute_pressure_dist_button = ttk.Button(button_frame, text="Compute Pressure Distribution", command=self.compute_pressure_distribution)
        self.compute_pressure_dist_button.grid(row=13, column=0, padx=5, pady=5, sticky="w")

        self.plot_pressure_dist_button = ttk.Button(button_frame, text="Plot Pressure Distribution", command=self.plot_pressure_distribution)
        self.plot_pressure_dist_button.grid(row=14, column=0, padx=5, pady=5, sticky="w")

        self.compute_lift_button = ttk.Button(button_frame, text="Compute Lift",
                                                    command=self.compute_lift)
        self.compute_lift_button.grid(row=15, column=0, padx=5, pady=5, sticky="w")

        self.compute_induced_drag_button = ttk.Button(button_frame, text="Compute Induced Drag",
                                               command=self.compute_induced_drag)
        self.compute_induced_drag_button.grid(row=16, column=0, padx=5, pady=5, sticky="w")

        self.collect_forces_button = ttk.Button(button_frame, text="Collect Aerodynamic Forces",
                                                      command=self.get_forces_dlr_vlm)
        self.collect_forces_button.grid(row=17, column=0, padx=5, pady=5, sticky="w")

        text_frame = tk.Frame(self)
        text_frame.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

        self.scrollbar = tk.Scrollbar(text_frame)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        self.text_field = tk.Text(text_frame, height=40, width=80, yscrollcommand=self.scrollbar.set)
        self.text_field.pack(side=tk.LEFT, fill=tk.BOTH)

        self.scrollbar.config(command=self.text_field.yview)

        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

    def import_geometry_data(self):

        """
        Imports the geometry data from the Geometry class. This includes panel coordinates, vortex points,
        control points, and trailing vortex points.
        """

        self.panel_coordinates = None
        self.vortex_points = None
        self.control_points = None
        self.trailing_vortex_points = None

        wing_geometry_data = self.geometry.build_aerodynamic_mesh()

        if wing_geometry_data:

            self.panel_coordinates = wing_geometry_data.get('panel_coordinates', None)
            self.vortex_points = wing_geometry_data.get('vortex_points', None)
            self.control_points = wing_geometry_data.get('control_points', None)
            self.trailing_vortex_points = wing_geometry_data.get('trailing_vortex_points', None)

            # No mirroring right now

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
        self.output_label.config(text=f"Freestream Velocity Vector: {self.freestream_velocity}")
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

        self.text_field.delete(1.0, tk.END)
        self.text_field.insert(tk.END, "Normal Vectors and Their Magnitudes:\n")
        self.text_field.insert(tk.END, "-------------------------------------\n")

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

            norm_magnitude = np.linalg.norm(normal_vector)

            self.text_field.insert(
                tk.END, f"Panel {i + 1}: Normal Vector = {normal_vector}, Magnitude = {norm_magnitude:.4f}\n"
            )

            ax.quiver(
                control_point[0], control_point[1], control_point[2],
                normal_vector[0], normal_vector[1], normal_vector[2],
                color='r', length=0.01, normalize=True
            )

        ax.set_title("Wing Geometry with Normal Vectors")
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Axis')
        ax.set_zlabel('Z Axis')

        plt.show()

    def compute_panel_distances_and_dihedral(self, vortex_points, control_points):

        """
        Computes the distances between bound vortex segments and control points, including the dihedral angle.

        Arguments:
            vortex_points: List of vortex segment points.
            control_points: List of control points.

        Returns: Distances and dihedral angles for each panel as tuple
        """

        distances = []
        dihedral_angles = []

        dihedral_angle = float(self.geometry.dihedral_entry.get())

        for vortex_segment, control_point in zip(vortex_points, control_points):
            P1, P3 = vortex_segment

            P_mid = (P1 + P3) * 0.5

            r = control_point - P_mid
            distance = np.linalg.norm(r)
            distances.append(distance)

            dihedral_angles.append(dihedral_angle)

        return distances, dihedral_angles

    def display_distances_and_dihedral(self):

        """
        Displays the distances and dihedral for each panel in the GUI for verification
        """

        if self.vortex_points is None or self.control_points is None:
            self.output_label.config(text="Please import geometry data first.")
            return

        if self.normal_vectors is None:
            self.output_label.config(text="Please compute the normal vectors first.")
            return

        distances, dihedral_angles = self.compute_panel_distances_and_dihedral(self.vortex_points, self.control_points)

        self.panel_distances = distances
        self.panel_dihedrals = dihedral_angles

        return self.panel_distances, self.panel_dihedrals

    def compute_induced_velocity_vortex_line(self, P0, P1, P3):

        """
        Computes the velocity induced by the bound vortex at a control point (eq. 2.4 in the DLR script)

        Arguments:
            P0, P1, P3: Coordinates of the control point, and start/end points of the vortex line.

        Returns: Induced velocity vector at the control point.
        """

        r0 = np.array(P3) - np.array(P1)
        r1 = np.array(P0) - np.array(P1)
        r2 = np.array(P0) - np.array(P3)

        r1_cross_r2 = np.cross(r1, r2)
        r0_dot_r1 = np.dot(r0, r1)
        r0_dot_r2 = np.dot(r0, r2)

        r1_norm = np.linalg.norm(r1)
        r2_norm = np.linalg.norm(r2)
        r1_cross_r2_norm = np.linalg.norm(r1_cross_r2)

        epsilon = 1e-6
        if r1_norm < epsilon or r2_norm < epsilon or r1_cross_r2_norm < epsilon:
            return np.array([0.0, 0.0, 0.0])

        term1 = r1_cross_r2 / (4 * np.pi * r1_cross_r2_norm ** 2)
        term2 = (r0_dot_r1 / r1_norm) - (r0_dot_r2 / r2_norm)

        induced_velocity_bound_vortex = term1 * term2

        return induced_velocity_bound_vortex

    def display_induced_velocities_vortex_line(self):

        """
        Displays the bound vortex induced velocity for verification
        """

        if self.vortex_points is None or self.control_points is None:
            self.text_field.delete(1.0, tk.END)
            self.text_field.insert(tk.END, "Please import geometry data first.\n")
            return

        induced_velocities_bound_vortex = []

        for i, (vortex_segment, control_point) in enumerate(zip(self.vortex_points, self.control_points)):

            P1, P3 = vortex_segment
            P0 = control_point

            induced_velocity = self.compute_induced_velocity_vortex_line(np.array(P0), np.array(P1), np.array(P3))

            induced_velocities_bound_vortex.append(induced_velocity)

        self.text_field.delete(1.0, tk.END)

        for i, velocity in enumerate(induced_velocities_bound_vortex):
            result_text = f"Panel {i + 1} Induced Velocity: {velocity}\n"
            self.text_field.insert(tk.END, result_text)

        return induced_velocities_bound_vortex

    def compute_inner_horseshoe_induced_velocity(self, P0, P1, distance, dihedral_angle):

        """
        Computes the induced velocity from the inner trailing vortex segment (eq. 2.8 and 2.9 in the DLR script).

        Arguments:
            P0, P1: Coordinates of the control point and vortex point.
            distance: Distance from vortex to control point.
            dihedral_angle: Dihedral angle in radians.

        Returns: Induced velocity vector.
        """

        r1 = np.array(P0) - np.array(P1)
        r1_x = r1[0]
        r1_norm = np.linalg.norm(r1)

        cos_beta_1 = 1.0
        cos_beta_2 = -r1_x / r1_norm

        cos_difference = cos_beta_1 - cos_beta_2

        Dv_2 = (-np.sin(np.radians(dihedral_angle)) / (4 * np.pi * distance)) * cos_difference
        Dw_2 = (-np.cos(np.radians(dihedral_angle)) / (4 * np.pi * distance)) * cos_difference

        return [0, Dv_2, Dw_2]

    def display_inner_horseshoe_induced_velocities(self):

        """
        Displays the velocity induced by the inner trailing vortex in the GUI for verification
        """

        if self.panel_distances is None or self.control_points is None or self.vortex_points is None:
            self.text_field.delete(1.0, tk.END)
            self.text_field.insert(tk.END, "Please compute panel distances and normal vectors first.\n")
            return

        dihedral_angle = np.radians(float(self.geometry.dihedral_entry.get()))
        inner_horseshoe_velocities = []

        for i, (distance, control_point, vortex_segment) in enumerate(
                zip(self.panel_distances, self.control_points, self.vortex_points)):

            P1, P3 = vortex_segment
            P0 = control_point
            velocity = self.compute_inner_horseshoe_induced_velocity(P0, P1, distance, dihedral_angle)
            inner_horseshoe_velocities.append(velocity)

        self.text_field.delete(1.0, tk.END)
        self.text_field.insert(tk.END, f"Inner Horseshoe Induced Velocities: {inner_horseshoe_velocities}\n")

    def compute_outer_horseshoe_induced_velocity(self, P0, P3, distance, dihedral_angle):

        """
        Computes the induced velocity from the outer trailing vortex segment (eq. 2.8 and 2.9 in the DLR script modified)

        Arguments:
            P0, P3: Coordinates of the control point and vortex point.
            distance: Distance from vortex to control point.
            dihedral_angle: Dihedral angle in radians.

        Returns: Induced velocity vector.
        """

        r2 = np.array(P0) - np.array(P3)
        r2_x = r2[0]
        r2_norm = np.linalg.norm(r2)

        cos_beta_1 = r2_x / r2_norm
        cos_beta_2 = -1.0

        cos_difference = cos_beta_1 - cos_beta_2

        Dv_3 = (-np.sin(np.radians(dihedral_angle)) / (4 * np.pi * distance)) * cos_difference
        Dw_3 = (-np.cos(np.radians(dihedral_angle)) / (4 * np.pi * distance)) * cos_difference

        return [0, Dv_3, Dw_3]

    def display_outer_horseshoe_induced_velocities(self):

        """
        Displays the velocity induced by the outer trailing vortex in the GUI for verification
        """

        if self.panel_distances is None or self.control_points is None or self.vortex_points is None:
            self.text_field.delete(1.0, tk.END)
            self.text_field.insert(tk.END, "Please compute panel distances and normal vectors first.\n")
            return

        dihedral_angle = np.radians(float(self.geometry.dihedral_entry.get()))
        outer_horseshoe_velocities = []

        for i, (distance, control_point, vortex_segment) in enumerate(
                zip(self.panel_distances, self.control_points, self.vortex_points)):
            P1, P3 = vortex_segment
            P0 = control_point

            velocity = self.compute_outer_horseshoe_induced_velocity(P0, P3, distance, dihedral_angle)
            outer_horseshoe_velocities.append(velocity)

        self.text_field.delete(1.0, tk.END)
        self.text_field.insert(tk.END, f"Outer Horseshoe Induced Velocities: {outer_horseshoe_velocities}\n")

    def compute_panel_area(self, panel_coords):

        """
        Calculates the area of a panel given its coordinates.

        Arguments:
            panel_coords: List of coordinates defining the panel vertices.

        Returns: Area of the panel.
        """

        P1 = np.array(panel_coords[0])
        P2 = np.array(panel_coords[1])
        P3 = np.array(panel_coords[2])
        P4 = np.array(panel_coords[3])

        self.triangle_1_area = 0.5 * np.linalg.norm(np.cross(P2 - P1, P3 - P1))
        self.triangle_2_area = 0.5 * np.linalg.norm(np.cross(P3 - P1, P4 - P1))

        self.A = self.triangle_1_area + self.triangle_2_area

        return self.A

    def compute_spanwise_width(self, panel_coords):

        """
        Calculates the spanwise width of a panel.

        Arguments:
            panel_coords: List of coordinates defining the panel vertices.

        Returns: Spanwise width of the panel.
        """

        y_a = panel_coords[0][1]
        y_b = panel_coords[3][1]

        self.b_panel = np.abs(y_b - y_a)

        return self.b_panel

    def display_panel_areas(self):

        """
        Displays the panel areas in the GUI for verification
        """

        if self.panel_coordinates is None:
            self.text_field.delete(1.0, tk.END)
            self.text_field.insert(tk.END, "Please import geometry data first.\n")
            return

        self.text_field.delete(1.0, tk.END)
        self.text_field.insert(tk.END, "Panel Areas:\n")

        num_panels = min(10, len(self.panel_coordinates))
        for i in range(num_panels):
            area = self.compute_panel_area(self.panel_coordinates[i])
            self.text_field.insert(tk.END, f"Panel {i + 1} Area: {area:.4f} m²\n")

    def display_spanwise_widths(self):

        """
        Displays the panel areas in the GUI for verification
        """

        self.spanwise_widths = []

        if self.panel_coordinates is None:
            self.text_field.delete(1.0, tk.END)
            self.text_field.insert(tk.END, "Please import geometry data first.\n")
            return

        self.text_field.delete(1.0, tk.END)
        self.text_field.insert(tk.END, "Spanwise Widths:\n")

        num_panels = min(10, len(self.panel_coordinates))
        for i in range(num_panels):
            spanwise_width = self.compute_spanwise_width(self.panel_coordinates[i])
            self.text_field.insert(tk.END, f"Panel {i + 1} Spanwise Width: {spanwise_width:.4f} m\n")
            self.spanwise_widths.append(spanwise_width)

    '''
    def mirror_geometry(self):
    
        # Full wing geometry mirroring for complete AIC matrix - in this code, symmetry is assumed and only the half
        # wing is modelled

        mirrored_panel_coordinates = [[(x, -y, z) for (x, y, z) in panel] for panel in self.panel_coordinates]
        mirrored_vortex_points = [[(P1[0], -P1[1], P1[2]), (P3[0], -P3[1], P3[2])] for P1, P3 in self.vortex_points]
        mirrored_control_points = [(x, -y, z) for (x, y, z) in self.control_points]
        mirrored_normal_vectors = [(Nx, -Ny, -Nz) for (Nx, Ny, Nz) in self.normal_vectors]

        self.panel_coordinates += mirrored_panel_coordinates
        self.vortex_points += mirrored_vortex_points
        self.control_points += mirrored_control_points
        self.normal_vectors += mirrored_normal_vectors
        self.panel_distances *= 2
    '''

    def assemble_pressure_aic(self):

        """
        Assembles the pressure-mapped AIC matrix (eq. 2.10 - 2.13 in the DLR script)

        Returns: Pressure-mapped AIC matrix.
        """

        if self.panel_coordinates is None or self.vortex_points is None or self.control_points is None:
            self.text_field.delete(1.0, tk.END)
            self.text_field.insert(tk.END, "Please import geometry and vortex data first.\n")
            return

        num_panels = len(self.control_points)
        self.pressure_AIC_matrix = np.zeros((num_panels, num_panels))
        self.A_matrix = np.zeros((num_panels, num_panels))

        dihedral_angle = np.radians(float(self.geometry.dihedral_entry.get()))

        self.text_field.delete(1.0, tk.END)
        self.text_field.insert(tk.END,
                               'Panel i | Panel j |  Dv_total  |  Dw_total  |  D_weighted  | Panel Area | Spanwise Width\n')
        self.text_field.insert(tk.END, '--------------------------------------------------------------------------\n')

        for i in range(num_panels):

            P0 = self.control_points[i]
            normal_vector = self.normal_vectors[i]
            N_y = normal_vector[1]
            N_z = normal_vector[2]

            for j in range(num_panels):

                vortex_segment = self.vortex_points[j]
                P1, P3 = vortex_segment
                distance = self.panel_distances[j]

                bound_vortex_velocity = self.compute_induced_velocity_vortex_line(P0, P1, P3)
                inner_horseshoe_velocity = self.compute_inner_horseshoe_induced_velocity(P0, P1, distance,
                                                                                         dihedral_angle)
                outer_horseshoe_velocity = self.compute_outer_horseshoe_induced_velocity(P0, P3, distance,
                                                                                         dihedral_angle)

                Du_1, Dv_1, Dw_1 = bound_vortex_velocity
                Dv_2, Dw_2 = inner_horseshoe_velocity[1], inner_horseshoe_velocity[2]
                Dv_3, Dw_3 = outer_horseshoe_velocity[1], outer_horseshoe_velocity[2]

                Dv_total_starboard = Dv_1 + Dv_2 + Dv_3
                Dw_total_starboard = Dw_1 + Dw_2 + Dw_3

                P1_mirrored = np.array([P3[0], -P3[1], P3[2]])
                P3_mirrored = np.array([P1[0], -P1[1], P1[2]])

                bound_vortex_velocity_mirrored = self.compute_induced_velocity_vortex_line(P0, P1_mirrored, P3_mirrored)
                outer_horseshoe_velocity_mirrored = self.compute_outer_horseshoe_induced_velocity(
                    P0, P3_mirrored, distance, dihedral_angle)


                inner_horseshoe_velocity_mirrored = self.compute_inner_horseshoe_induced_velocity(
                        P0, P1_mirrored, distance, dihedral_angle)


                Du_1_mirrored, Dv_1_mirrored, Dw_1_mirrored = bound_vortex_velocity_mirrored
                Dv_2_mirrored, Dw_2_mirrored = inner_horseshoe_velocity_mirrored[1], inner_horseshoe_velocity_mirrored[
                    2]
                Dv_3_mirrored, Dw_3_mirrored = outer_horseshoe_velocity_mirrored[1], outer_horseshoe_velocity_mirrored[2]

                Dv_total_port = Dv_1_mirrored + Dv_2_mirrored + Dv_3_mirrored
                Dw_total_port = Dw_1_mirrored + Dw_2_mirrored + Dw_3_mirrored

                Dv_total = Dv_total_starboard + Dv_total_port
                Dw_total = Dw_total_starboard + Dw_total_port

                D_weighted = Dv_total * N_y + Dw_total * N_z

                panel_area = self.compute_panel_area(self.panel_coordinates[j])
                panel_spanwise_width = self.compute_spanwise_width(self.panel_coordinates[j])

                scaling_factor = panel_area / (2 * panel_spanwise_width)
                D_transformed = D_weighted * scaling_factor

                self.A_matrix[i, j] = D_transformed

                self.text_field.insert(
                    tk.END,
                    f'{i:^7} | {j:^7} | {Dv_total:^9.4f} | {Dw_total:^9.4f} | {D_weighted:^11.4f} | {panel_area:^10.4f} | {panel_spanwise_width:^14.4f}\n')

        try:
            self.pressure_AIC_matrix = -np.linalg.inv(self.A_matrix)
            condition_number = np.linalg.cond(self.pressure_AIC_matrix)
        except np.linalg.LinAlgError:
            self.text_field.insert(tk.END, "Matrix inversion failed (singular matrix).\n")
            return

        matrix_size = self.pressure_AIC_matrix.shape

        self.text_field.insert(tk.END, f"AIC Matrix with Pressure Size: {matrix_size}\n")
        self.text_field.insert(tk.END, f"AIC Matrix Condition Number: {condition_number}\n")
        self.text_field.insert(tk.END, f"AIC Matrix (First 5 Rows):\n{self.pressure_AIC_matrix[:5, :5]}\n")

        return self.pressure_AIC_matrix

    def compute_pressure_distribution(self):

        """
        Computes the pressure coefficient (ΔCp) distribution across panels based on the AIC matrix (eq. 2.3 in the DLR script).

        Returns: Pressure coefficient distribution.
        """

        if self.pressure_AIC_matrix is None:
            self.text_field.delete(1.0, tk.END)
            self.text_field.insert(tk.END, "Please ensure that the AIC matrix is computed first.\n")
            return

        self.total_downwash = []

        for i, (P0, vortex_segment) in enumerate(zip(self.control_points, self.vortex_points)):

            P1, P3 = vortex_segment
            distance = self.panel_distances[i % (len(self.control_points) // 2)]
            dihedral_angle = np.radians(float(self.geometry.dihedral_entry.get()))

            inner_horseshoe_velocity = self.compute_inner_horseshoe_induced_velocity(P0, P1, distance, dihedral_angle)
            outer_horseshoe_velocity = self.compute_outer_horseshoe_induced_velocity(P0, P3, distance, dihedral_angle)

            Dv_2, Dw_2 = inner_horseshoe_velocity[1], inner_horseshoe_velocity[2]
            Dv_3, Dw_3 = outer_horseshoe_velocity[1], outer_horseshoe_velocity[2]

            P1_mirrored = np.array([P3[0], -P3[1], P3[2]]) # mirrored vortex point
            P3_mirrored = np.array([P1[0], -P1[1], P1[2]]) # mirrored vortex point

            # Influence of left wing on downwash of right wing

            inner_horseshoe_velocity_mirrored = self.compute_inner_horseshoe_induced_velocity(
                P0, P1_mirrored, distance, dihedral_angle)
            outer_horseshoe_velocity_mirrored = self.compute_outer_horseshoe_induced_velocity(
                P0, P3_mirrored, distance, dihedral_angle)

            Dv_2_mirrored, Dw_2_mirrored = inner_horseshoe_velocity_mirrored[1], inner_horseshoe_velocity_mirrored[2]
            Dv_3_mirrored, Dw_3_mirrored = outer_horseshoe_velocity_mirrored[1], outer_horseshoe_velocity_mirrored[2]

            Dv_total = Dv_2 + Dv_3 + Dv_2_mirrored + Dv_3_mirrored
            Dw_total = Dw_2 + Dw_3 + Dw_2_mirrored + Dw_3_mirrored

            normal_vector = self.normal_vectors[i]
            N_y = normal_vector[1]
            N_z = normal_vector[2]

            total_downwash_value = Dw_total * N_z
            self.total_downwash.append(total_downwash_value)

        self.total_downwash = np.array(self.total_downwash)

        try:
            self.delta_Cp = np.dot(self.pressure_AIC_matrix, self.total_downwash)
        except np.linalg.LinAlgError as e:
            self.text_field.delete(1.0, tk.END)
            self.text_field.insert(tk.END, f"Error solving for ΔCp: {str(e)}\n")
            return None

        self.text_field.delete(1.0, tk.END)
        self.text_field.insert(tk.END, f"Pressure Coefficient Distribution (ΔCp):\n{self.delta_Cp}\n")

        return self.delta_Cp

    def plot_pressure_distribution(self):

        """
        Plots the computed pressure coefficient distribution as a color-coded 3D plot over the wing geometry.
        """

        if self.panel_coordinates is None or self.delta_Cp is None:
            self.text_field.delete(1.0, tk.END)
            self.text_field.insert(tk.END, "Please compute pressure distribution first.\n")
            return

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title("Pressure Coefficient Distribution")

        min_cp = np.min(self.delta_Cp[:len(self.panel_coordinates) // 2])
        max_cp = np.max(self.delta_Cp[:len(self.panel_coordinates) // 2])

        for i, panel_coords in enumerate(self.panel_coordinates):

            panel = np.array(panel_coords)
            verts = [panel]

            pressure = self.delta_Cp[i]

            color_value = plt.cm.jet((pressure - min_cp) / (max_cp - min_cp))

            poly = Poly3DCollection(verts, facecolors=color_value, edgecolor='k')
            ax.add_collection3d(poly)

        ax.set_xlabel("X Axis")
        ax.set_ylabel("Y Axis")
        ax.set_zlabel("Z Axis")

        mappable = plt.cm.ScalarMappable(cmap=plt.cm.jet, norm=plt.Normalize(vmin=min_cp, vmax=max_cp))
        mappable.set_array(self.delta_Cp[:len(self.panel_coordinates) // 2])  # Only for the original panels
        fig.colorbar(mappable, ax=ax, label="Pressure Coefficient (ΔCp)")

        plt.show()

    def compute_lift(self):

        """
        Calculates the panel lift, total lift and lift coefficient for the wing based on the pressure distribution.

        from Delta Cp = (p-p_\infty)/(0.5*\rho*V_\infty^2) -> pressure distribution * panel area

        Returns: Panel lift array, total lift, lift coefficient, and circulation.
        """

        if self.delta_Cp is None or self.panel_coordinates is None:
            self.text_field.delete(1.0, tk.END)
            self.text_field.insert(tk.END, "Please compute pressure distribution first.\n")
            return

        if self.velocity is None or self.density is None:
            self.text_field.delete(1.0, tk.END)
            self.text_field.insert(tk.END, "Please import freestream data first.\n")
            return

        self.q_inf = 0.5 * self.density * self.velocity ** 2

        self.L = 0.0
        self.S = 0.0

        Gamma = []
        panel_lift = []

        self.text_field.delete(1.0, tk.END)
        self.text_field.insert(tk.END, '\nPanel |   Cp   |  Panel Circulation | Panel Lift (N) |  Panel cl  |\n')
        self.text_field.insert(tk.END, '-----------------------------------------------\n')

        for i, panel in enumerate(self.panel_coordinates):

            v1 = np.array(panel[1]) - np.array(panel[0])
            v2 = np.array(panel[3]) - np.array(panel[0])
            panel_area = np.linalg.norm(np.cross(v1, v2))

            y_a = panel[0][1]
            y_b = panel[3][1]

            dy = abs(y_b - y_a)

            self.S += panel_area

            n_z = self.normal_vectors[i][2]

            L_panel = -self.delta_Cp[i] * self.q_inf * panel_area * n_z
            panel_lift.append(L_panel)

            self.L += L_panel

            Gamma_panel = L_panel / (self.density * self.velocity)

            Gamma.append(Gamma_panel)

            cl_panel = L_panel / (self.q_inf * panel_area)

            self.text_field.insert(tk.END, ' %3d  | %6.3f | %12.3f |   %12.3f  |  %7.4f   |\n' % (
                i, self.delta_Cp[i], Gamma_panel, L_panel, cl_panel))

        Gamma = np.array(Gamma)

        self.S = self.S * 2
        self.L = self.L * 2  # symmetry assumption
        self.CL = self.L / (self.q_inf * self.S)

        self.text_field.insert(tk.END, f"Total Lift: {self.L:.2f} N\n")
        self.text_field.insert(tk.END, 'Total Wing Area: %.3f m²\n' % self.S)
        self.text_field.insert(tk.END, 'Total Lift Coefficient (CL): %.4f\n' % self.CL)

        return panel_lift, self.L, self.CL, Gamma

    def compute_induced_drag(self):

        """
        Calculates the induced drag based on the obtained circulation from the lift after equation 12.27 by Katz and Plotkin,

        Low Speed Aerodynamics

        Returns: Panel induced drag array, total induced drag.
        """

        self.induced_drag = 0.0
        panel_induced_drag = []

        L_panel, L, CL, Gamma = self.compute_lift()

        for j, gamma_current in enumerate(Gamma):

            panel_spanwise_width = self.compute_spanwise_width(self.panel_coordinates[j])

            w_ind = self.total_downwash[j]

            if j == 0:
                induced_drag_panel = -self.density * w_ind * gamma_current * panel_spanwise_width
            else:
                gamma_previous = Gamma[j - 1]
                induced_drag_panel = -self.density * w_ind * (gamma_current - gamma_previous) * panel_spanwise_width

            panel_induced_drag.append(induced_drag_panel)
            self.induced_drag += induced_drag_panel

        self.induced_drag = self.induced_drag * 2

        self.text_field.delete(1.0, tk.END)
        self.text_field.insert(tk.END, f"Total Induced Drag: {self.induced_drag:.2f} N\n")

        return self.induced_drag, panel_induced_drag

    def get_forces_dlr_vlm(self):

        """
        Collects the summed spanwise aerodynamic forces in a dictionary to pass to the FEM class
        """

        panel_lift, total_lift, CL, _ = self.compute_lift()
        _, panel_induced_drag = self.compute_induced_drag()

        num_chordwise_panels = int(self.geometry.n_panels_chord_entry.get())
        span = float(self.geometry.span_entry.get()) / 2
        num_panels_spanwise = len(self.control_points) // num_chordwise_panels

        panel_lift_chordwise = np.zeros((num_chordwise_panels, num_panels_spanwise))
        panel_drag_chordwise = np.zeros((num_chordwise_panels, num_panels_spanwise))
        panel_pressure_chordwise = np.zeros((num_chordwise_panels, num_panels_spanwise))

        y_stations = np.linspace(0, span, num_panels_spanwise)

        panel_areas = np.zeros(num_chordwise_panels * num_panels_spanwise)

        for j in range(num_panels_spanwise):

            for i in range(num_chordwise_panels):

                idx = i * num_panels_spanwise + j
                panel_coords = self.panel_coordinates[idx]

                panel_area = self.compute_panel_area(panel_coords)
                panel_areas[idx] = panel_area

                panel_lift_chordwise[i, j] = panel_lift[idx]
                panel_drag_chordwise[i, j] = panel_induced_drag[idx]
                panel_pressure_chordwise[i, j] = panel_lift[idx] / panel_area

        L_panel_tot = np.sum(panel_lift_chordwise, axis=0)
        Di_panel_tot = np.sum(panel_drag_chordwise, axis=0)
        total_panel_pressure = np.sum(panel_pressure_chordwise, axis=0)

        forces = {
            'spanwise_lift_distribution': L_panel_tot,
            'spanwise_drag_distribution': Di_panel_tot,
            'spanwise_pressure_distribution': total_panel_pressure,
        }

        self.text_field.delete(1.0, tk.END)
        self.text_field.insert(tk.END, "Spanwise Stations (y) and Force Distributions:\n")
        self.text_field.insert(tk.END, "-----------------------------------------------\n")
        self.text_field.insert(tk.END, f"{'y':>10}  |  {'Lift Distribution':>18}  |  {'Drag Distribution':>18}\n")

        for i in range(len(y_stations)):
            self.text_field.insert(tk.END,
                                   f"{y_stations[i]:>10.3f}  |  {forces['spanwise_lift_distribution'][i]:>18.3f}  |  {forces['spanwise_drag_distribution'][i]:>18.3f}\n")

        return forces, y_stations

    def get_panel_forces_dlr_vlm(self):

        """
        Collects the panel aerodynamic forces in a dictionary to pass to the FEM class
        """

        panel_lift, total_lift, CL, _ = self.compute_lift()
        _, panel_induced_drag = self.compute_induced_drag()

        num_chordwise_panels = int(self.geometry.n_panels_chord_entry.get())
        span = float(self.geometry.span_entry.get()) / 2
        num_panels_spanwise = len(self.control_points) // num_chordwise_panels

        y_stations = np.linspace(0, span, num_panels_spanwise)

        forces = {
            'lift_distribution': panel_lift,
            'drag_distribution': panel_induced_drag,
        }

        return forces, y_stations
