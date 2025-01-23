#########################################################################################

# Aerodynamic geometry and mesh module

#########################################################################################

import tkinter as tk
from tkinter import ttk, filedialog
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class Panel:

    """
    Represents a single panel for an airfoil section, including geometric properties like endpoints,
    center point, length, and orientation.
    """

    def __init__(self, xa, ya, xb, yb):

        """
        Initializes a panel with endpoints and calculates its center, length, and orientation.
        """

        self.xa, self.ya = xa, ya
        self.xb, self.yb = xb, yb
        self.xc, self.yc = (xa + xb) / 2, (ya + yb) / 2
        self.length = np.sqrt((xb - xa) ** 2 + (yb - ya) ** 2)

        if xb - xa <= 0.0:
            self.beta = np.arccos((yb - ya) / self.length)
        else:
            self.beta = np.pi + np.arccos(-(yb - ya) / self.length)

        self.loc = 'upper' if self.beta <= np.pi else 'lower'


class Geometry(ttk.Frame):

    """
    Class to handle the user interface and calculations for wing geometry and airfoil data, including
    panel creation, 2D airfoil visualization, and 3D wing visualization.
    """

    def __init__(self, parent, *args, **kwargs):

        """
        Initializes the Geometry class, setting up input fields, buttons, and placeholders for airfoil
        and wing geometry visualizations.
        """

        super().__init__(parent, *args, **kwargs)

        self.airfoil_data = {}
        self.panels_data = {'root': None, 'tip': None}
        self.current_airfoil = 'root'

        ttk.Label(self, text="Root Airfoil:", font=("Arial", 12)).grid(row=0, column=0, padx=5, pady=5, sticky="w")
        root_file_button = ttk.Button(self, text="Select File", command=self.load_root_airfoil)
        root_file_button.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

        ttk.Label(self, text="Tip Airfoil:", font=("Arial", 12)).grid(row=1, column=0, padx=5, pady=5, sticky="w")
        tip_file_button = ttk.Button(self, text="Select File", command=self.load_tip_airfoil)
        tip_file_button.grid(row=1, column=1, padx=5, pady=5, sticky="ew")

        ttk.Label(self, text="Number of chordwise Panels:", font=("Arial", 12)).grid(row=2, column=0, padx=5, pady=5, sticky="w")
        self.n_panels_chord_entry = ttk.Entry(self)
        self.n_panels_chord_entry.grid(row=2, column=1, padx=5, pady=5, sticky="ew")

        ttk.Label(self, text="Number of spanwise Panels:", font=("Arial", 12)).grid(row=3, column=0, padx=5, pady=5, sticky="w")
        self.n_panels_span_entry = ttk.Entry(self)
        self.n_panels_span_entry.grid(row=3, column=1, padx=5, pady=5, sticky="ew")

        define_panels_button = ttk.Button(self, text="Build Panels", command=self.define_panels_and_visualize)
        define_panels_button.grid(row=4, column=0, columnspan=2, padx=5, pady=10, sticky="ew")

        show_root_button = ttk.Button(self, text="Show Root Airfoil", command=self.show_root_airfoil)
        show_root_button.grid(row=5, column=0, padx=5, pady=5, sticky="ew")

        show_tip_button = ttk.Button(self, text="Show Tip Airfoil", command=self.show_tip_airfoil)
        show_tip_button.grid(row=5, column=1, padx=5, pady=5, sticky="ew")

        ttk.Label(self, text="Wing Geometry Inputs", font=("Arial", 14, 'bold')).grid(row=0, column=2, padx=5, pady=5, sticky="w")

        ttk.Label(self, text="Span (m):", font=("Arial", 12)).grid(row=1, column=2, padx=5, pady=5)
        self.span_entry = ttk.Entry(self)
        self.span_entry.grid(row=1, column=3, padx=2, pady=5)

        ttk.Label(self, text="Root Chord Length (m):", font=("Arial", 12)).grid(row=2, column=2, padx=5, pady=5)
        self.root_chord_entry = ttk.Entry(self)
        self.root_chord_entry.grid(row=2, column=3, padx=2, pady=5)

        ttk.Label(self, text="Tip Chord Length (m):", font=("Arial", 12)).grid(row=3, column=2, padx=5, pady=5)
        self.tip_chord_entry = ttk.Entry(self)
        self.tip_chord_entry.grid(row=3, column=3, padx=2, pady=5)

        ttk.Label(self, text="Sweep Angle (degrees):", font=("Arial", 12)).grid(row=4, column=2, padx=5, pady=5)
        self.sweep_entry = ttk.Entry(self)
        self.sweep_entry.grid(row=4, column=3, padx=2, pady=5)

        ttk.Label(self, text="Dihedral Angle (degrees):", font=("Arial", 12)).grid(row=5, column=2, padx=5, pady=5)
        self.dihedral_entry = ttk.Entry(self)
        self.dihedral_entry.grid(row=5, column=3, padx=2, pady=5)

        build_wing_button = ttk.Button(self, text="Visualize Aerodynamic Mesh", command=self.build_aerodynamic_mesh)
        build_wing_button.grid(row=6, column=2, padx=5, pady=10, columnspan = 1)

        reset_button = ttk.Button(self, text="Reset", command=self.reset_plots)
        reset_button.grid(row=6, column=3, padx=5, pady=10, columnspan=1)

        self.airfoil_plot_frame = ttk.Frame(self)
        self.airfoil_plot_frame.grid(row=7, column=0, columnspan=2, padx=5, pady=5, sticky="nsew")
        self.airfoil_plot_frame.columnconfigure(0, weight=1)
        self.airfoil_plot_frame.rowconfigure(0, weight=1)

        self.init_empty_airfoil_plot()

        self.wing_plot_frame = ttk.Frame(self)
        self.wing_plot_frame.grid(row=7, column=2, columnspan=2, padx=5, pady=5, sticky="ns")
        self.wing_plot_frame.columnconfigure(0, weight=1)
        self.wing_plot_frame.rowconfigure(0, weight=1)

        self.init_empty_wing_plot()

    def reset_plots(self):

        """
        Resets the airfoil and wing plot frames to an empty state.
        """

        self.init_empty_airfoil_plot()

        self.init_empty_wing_plot()

        self.airfoil_data = {}
        self.panels_data = {'root': None, 'tip': None}
        self.current_airfoil = 'root'

    def init_empty_airfoil_plot(self):

        """
        Initializes an empty plot for the airfoil section.
        """

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.set_title("Airfoil Geometry")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_ylim(-0.2, 0.2)
        ax.grid()

        for widget in self.airfoil_plot_frame.winfo_children():
            widget.destroy()
        canvas = FigureCanvasTkAgg(fig, master=self.airfoil_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def init_empty_wing_plot(self):

        """
        Initializes an empty 3D plot for the wing geometry visualization.
        """

        fig = plt.figure(figsize=(13, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title("Wing Geometry Visualization")
        ax.set_xlabel("X Axis (m)")
        ax.set_ylabel("Y Axis (m)")
        ax.set_zlabel("Z Axis (m)")
        ax.grid(True)

        for widget in self.wing_plot_frame.winfo_children():
            widget.destroy()
        canvas = FigureCanvasTkAgg(fig, master=self.wing_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def load_airfoil_coordinates(self, filepath):

        """
        Loads airfoil coordinates from a .dat file.

        Arguments:
           filepath (str): Path to the .dat file containing airfoil coordinates.

        Returns:
           tuple: Arrays for x and y coordinates of the airfoil.

        Raises:
           FileNotFoundError: If the file does not exist.
           ValueError: If the file cannot be parsed correctly.
        """

        if not os.path.isfile(filepath):
            raise FileNotFoundError(f"Airfoil file not found: {filepath}")

        try:
            x, y = np.loadtxt(filepath, dtype=float, skiprows=1, unpack=True)
            return x, y
        except ValueError as e:
            print(f"Error loading airfoil coordinates: {e}")
            raise ValueError(f"Failed to load airfoil data from {filepath}. Please check the file format.")

    def load_root_airfoil(self):

        """
        Opens a file dialog to load the root airfoil coordinates.
        """

        file_path = filedialog.askopenfilename(defaultextension=".dat", filetypes=[("DAT files", "*.dat")])
        if file_path:
            x, y = self.load_airfoil_coordinates(file_path)
            self.airfoil_data['root'] = (x, y)
            self.visualize_airfoil_geometry()

    def load_tip_airfoil(self):

        """
        Opens a file dialog to load the tip airfoil coordinates.
        """

        file_path = filedialog.askopenfilename(defaultextension=".dat", filetypes=[("DAT files", "*.dat")])
        if file_path:
            x, y = self.load_airfoil_coordinates(file_path)
            self.airfoil_data['tip'] = (x, y)
            self.visualize_airfoil_geometry()

    '''
    def define_panels(self, x, y, N=40):
    
        """ Cosine airfoil shape spacing 
            
            Spaces the control points along the camber lines of the airfoil
        """
        
        R = (x.max() - x.min()) / 2.0  # radius
        x_center = (x.max() + x.min()) / 2.0  # center

        theta = np.linspace(0.0, 2.0 * np.pi, N + 1)
        x_circle = x_center + R * np.cos(theta)

        x_ends = np.copy(x_circle)
        y_ends = np.empty_like(x_ends)

        x, y = np.append(x, x[0]), np.append(y, y[0])

        I = 0
        for i in range(N):
            while I < len(x) - 1:
                if (x[I] <= x_ends[i] <= x[I + 1]) or (x[I + 1] <= x_ends[i] <= x[I]):
                    break
                else:
                    I += 1
            a = (y[I + 1] - y[I]) / (x[I + 1] - x[I])
            b = y[I + 1] - a * x[I + 1]
            y_ends[i] = a * x_ends[i] + b
        y_ends[N] = y_ends[0]

        # Create panels
        panels = np.array([Panel(x_ends[i], y_ends[i], x_ends[i + 1], y_ends[i + 1]) for i in range(N)])
        return panels
    '''

    '''
    def define_panels(self, x, y, N=40):
    
        """
        Chosine chord spacing 
        
        Uses cosine spacing to place the panel control points along the chord line of the airfoil
        
        """

        beta = np.linspace(0.0, np.pi, N + 1)
        x_ends = 0.5 * (1 - np.cos(beta))

        x_ends = x.min() + x_ends * (x.max() - x.min())

        y_interp = np.interp(x_ends, x, y)

        panels = []
        for i in range(N):
            panels.append(Panel(x_ends[i], y_interp[i], x_ends[i + 1], y_interp[i + 1]))

        return panels
    '''


    def define_panels(self, x, y, N):

        """
        Defines panel control points along an airfoil using linear spacing.

        Returns:
            List of `Panel` objects representing each panel along the airfoil.
        """

        x_ends = np.linspace(x.min(), x.max(), N + 1)

        y_interp = np.interp(x_ends, x, y)

        panels = []
        for i in range(N):
            panels.append(Panel(x_ends[i], y_interp[i], x_ends[i + 1], y_interp[i + 1]))

        return panels

    def define_panels_and_visualize(self):

        """
        Retrieves the number of panels from the entry, creates panels for the current airfoil,
        and visualizes the airfoil geometry.
        """

        try:
            N = int(self.n_panels_chord_entry.get())
            if self.current_airfoil in self.airfoil_data:
                x, y = self.airfoil_data[self.current_airfoil]
                self.panels_data[self.current_airfoil] = self.define_panels(x, y, N)
                self.visualize_airfoil_geometry()
            else:
                print(f"{self.current_airfoil.capitalize()} airfoil data not loaded yet.")
        except ValueError:
            print("Invalid number of panels.")

    def show_root_airfoil(self):

        """
        Displays the root airfoil plot
        """

        self.current_airfoil = 'root'
        self.visualize_airfoil_geometry()

    def show_tip_airfoil(self):

        """
        Displays the tip airfoil plot
        """

        self.current_airfoil = 'tip'
        self.visualize_airfoil_geometry()

    def visualize_airfoil_geometry(self):

        """
        Visualizes the airfoil geometry, including airfoil shape and the defined panels.
        """

        if self.current_airfoil not in self.airfoil_data:
            print(f"{self.current_airfoil.capitalize()} airfoil data not loaded yet.")
            return

        x, y = self.airfoil_data[self.current_airfoil]
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.grid()
        ax.set_xlabel('x', fontsize=16)
        ax.set_ylabel('y', fontsize=16)
        ax.plot(x, y, color='k', linestyle='-', linewidth=2, label=f"{self.current_airfoil.capitalize()} Airfoil")

        panels = self.panels_data.get(self.current_airfoil)
        if panels is not None:
            ax.plot(np.append([panel.xa for panel in panels], panels[0].xa),
                    np.append([panel.ya for panel in panels], panels[0].ya),
                    linestyle='-', linewidth=1, marker='o', markersize=6, color='#CD2305', label="Panels")

        ax.axis('scaled')
        ax.set_xlim(0, 1.1)
        ax.set_ylim(-0.25, 0.25)

        for widget in self.airfoil_plot_frame.winfo_children():
            widget.destroy()
        canvas = FigureCanvasTkAgg(fig, master=self.airfoil_plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def build_aerodynamic_mesh(self):

        """
        Builds and visualizes the 3D aerodynamic mesh based on user input for span, chords,
        sweep, and dihedral. Uses both root and tip airfoils to interpolate along the
        span and creates panels and vortex points for visualization.

        Returns:
            dict: A dictionary with geometry data, including panel coordinates,
                  vortex points, control points, and trailing vortex points.
        """

        try:
            span = float(self.span_entry.get()) / 2
            root_chord = float(self.root_chord_entry.get())
            tip_chord = float(self.tip_chord_entry.get())
            sweep = np.radians(float(self.sweep_entry.get()))
            dihedral = np.radians(float(self.dihedral_entry.get()))
            dihedral_deg = float(self.dihedral_entry.get())

            if 'root' not in self.airfoil_data or 'tip' not in self.airfoil_data:
                print("Please ensure both root and tip airfoil data are loaded.")
                return

            num_spanwise_panels = int(self.n_panels_span_entry.get())
            num_chordwise_panels = int(self.n_panels_chord_entry.get())

            root_x, root_y = self.airfoil_data['root']
            tip_x, tip_y = self.airfoil_data['tip']

            y_stations = np.linspace(0, span, num_spanwise_panels)
            chords = np.linspace(root_chord, tip_chord, num_spanwise_panels)

            fig = plt.figure(figsize=(13, 8))
            ax = fig.add_subplot(111, projection='3d')
            ax.set_title("Aerodynamic Mesh")
            ax.set_xlabel("X (m)")
            ax.set_ylabel("Y (m)")
            ax.set_zlabel("Z (m)")
            ax.grid(True)

            wing_panels = []
            vortex_points = []
            control_points = []
            panel_coordinates = []
            trailing_vortex_points = []

            for i, y_station in enumerate(y_stations):

                chord_length = chords[i]

                airfoil_x = (1 - y_station / span) * root_x + (y_station / span) * tip_x
                airfoil_y = (1 - y_station / span) * root_y + (y_station / span) * tip_y

                airfoil_x_scaled = airfoil_x * chord_length

                leading_edge_offset = (root_chord - chord_length) / 2
                airfoil_x_swept = airfoil_x_scaled + leading_edge_offset + y_station * np.tan(sweep)

                airfoil_z_dihedral = y_station * np.sin(dihedral)

                panels = self.define_panels(airfoil_x_swept, airfoil_y, num_chordwise_panels)

                if dihedral_deg != 0:
                    for panel in panels:
                        panel.ya += airfoil_z_dihedral
                        panel.yb += airfoil_z_dihedral
                else:
                    for panel in panels:
                        panel.ya = 0
                        panel.yb = 0

                wing_panels.append(panels)

            for i in range(len(wing_panels)):

                if i == len(wing_panels) - 1:

                    current_panels = wing_panels[i]
                    y_current = y_stations[i]

                    for j in range(len(current_panels)):
                        panel_a = current_panels[j]

                        P1_x = panel_a.xa + 0.25 * (panel_a.xb - panel_a.xa)
                        P1_y = y_current
                        P1_z = panel_a.ya + 0.25 * (panel_a.yb - panel_a.ya)
                        P1_z += y_current * np.sin(dihedral)

                        P3_x = panel_a.xb + 0.25 * (panel_a.xb - panel_a.xa)
                        P3_y = y_current
                        P3_z = panel_a.yb + 0.25 * (panel_a.yb - panel_a.ya)
                        P3_z += y_current * np.sin(dihedral)

                        P1 = np.array([P1_x, P1_y, P1_z])
                        P3 = np.array([P3_x, P3_y, P3_z])

                        P2_x = P1_x + 1e6
                        P2_y = P1_y
                        P2_z = P1_z

                        P4_x = P3_x + 1e6
                        P4_y = P3_y
                        P4_z = P3_z

                        P2 = np.array([P2_x, P2_y, P2_z])
                        P4 = np.array([P4_x, P4_y, P4_z])

                        vortex_points.append((P1, P3))
                        trailing_vortex_points.append((P2, P4))

                else:

                    current_panels = wing_panels[i]
                    y_current = y_stations[i]

                    next_panels = wing_panels[i + 1]
                    y_next = y_stations[i + 1]

                    y_mid = y_current + (y_next - y_current) / 2

                    for j in range(len(current_panels)):
                        panel_a = current_panels[j]
                        panel_b = next_panels[j]

                        x_a, y_a, z_a = panel_a.xa, y_current, panel_a.ya
                        x_b, y_b, z_b = panel_a.xb, y_current, panel_a.yb
                        x_c, y_c, z_c = panel_b.xb, y_next, panel_b.yb
                        x_d, y_d, z_d = panel_b.xa, y_next, panel_b.ya

                        z_a += y_current * np.sin(dihedral)
                        z_b += y_current * np.sin(dihedral)
                        z_c += y_next * np.sin(dihedral)
                        z_d += y_next * np.sin(dihedral)

                        P1_x = panel_a.xa + 0.25 * (panel_a.xb - panel_a.xa)
                        P1_y = y_current
                        P1_z = panel_a.ya + 0.25 * (panel_a.yb - panel_a.ya)

                        P3_x = panel_b.xa + 0.25 * (panel_b.xb - panel_b.xa)
                        P3_y = y_next
                        P3_z = panel_b.ya + 0.25 * (panel_b.yb - panel_b.ya)

                        P1_z += y_current * np.sin(dihedral)
                        P3_z += y_next * np.sin(dihedral)

                        P1 = np.array([P1_x, P1_y, P1_z])
                        P3 = np.array([P3_x, P3_y, P3_z])

                        x_mid = (panel_a.xa + panel_b.xa) / 2

                        P0_x = x_mid + 0.75 * (panel_a.xb - panel_a.xa)
                        P0_y = y_mid
                        P0_z = (panel_a.ya + 0.75 * (panel_a.yb - panel_a.ya)) + y_mid * np.sin(dihedral)
                        P0 = np.array([P0_x, P0_y, P0_z])

                        vortex_points.append((P1, P3))
                        control_points.append(P0)

                        panel_coordinates.append([(x_a, y_a, z_a), (x_b, y_b, z_b), (x_c, y_c, z_c), (x_d, y_d, z_d)])

                        P2_x = P1_x + 1e6
                        P2_y = P1_y
                        P2_z = P1_z

                        P4_x = P3_x + 1e6
                        P4_y = P3_y
                        P4_z = P3_z

                        P2 = np.array([P2_x, P2_y, P2_z])
                        P4 = np.array([P4_x, P4_y, P4_z])

                        trailing_vortex_points.append((P2, P4))

                        ax.plot([P1_x, P3_x], [P1_y, P3_y], [P1_z, P3_z], color='red',
                            label='Bound Vortex' if j == 0 and i == 0 else "")
                        ax.plot([x_a, x_b, x_c, x_d, x_a], [y_a, y_b, y_c, y_d, y_a], [z_a, z_b, z_c, z_d, z_a],
                            color='b')

            vortex_x, vortex_y, vortex_z = zip(*[v[0] for v in vortex_points])
            ax.scatter(vortex_x, vortex_y, vortex_z, color='k', s=10, label='Vortex Points', marker='o')

            control_x, control_y, control_z = zip(*control_points)
            ax.scatter(control_x, control_y, control_z, color='g', s=10, label='Control Points', marker='o')

            vortex_y_mirrored = [-y for y in vortex_y]
            control_y_mirrored = [-y for y in control_y]

            ax.scatter(vortex_x, vortex_y_mirrored, vortex_z, color='k', s=10, marker='o')
            ax.scatter(control_x, control_y_mirrored, control_z, color='g', s=10, marker='o')

            ax.legend()

            for line in ax.get_lines():
                x_data, y_data, z_data = line.get_data_3d()
                y_data_mirrored = [-y for y in y_data]
                ax.plot(x_data, y_data_mirrored, z_data, color='b')

                if line.get_color() == 'red':
                    ax.plot(x_data, [-y for y in y_data], z_data, color='red')

            for widget in self.wing_plot_frame.winfo_children():
                widget.destroy()
            canvas = FigureCanvasTkAgg(fig, master=self.wing_plot_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

            wing_geometry_data = {
                'panel_coordinates': panel_coordinates,
                'vortex_points': vortex_points,
                'control_points': control_points,
                'trailing_vortex_points': trailing_vortex_points
            }

            return wing_geometry_data

        except ValueError:
            print("Please enter valid numerical values for aerodynamic mesh.")

    def get_input_values(self):

        """
        Retrieves raw user input values for wing geometry.
        Returns:
            dict: Wing geometry data.
        """

        try:
            n_panels_chord = int(self.n_panels_chord_entry.get())
            n_panels_span = int(self.n_panels_span_entry.get())
            span = float(self.span_entry.get())
            root_chord = float(self.root_chord_entry.get())
            tip_chord = float(self.tip_chord_entry.get())
            sweep = float(self.sweep_entry.get())
            dihedral = float(self.dihedral_entry.get())

            return {
                'n_panels_chord': n_panels_chord,
                'n_panels_span': n_panels_span,
                'span': span,
                'root_chord': root_chord,
                'tip_chord': tip_chord,
                'sweep': sweep,
                'dihedral': dihedral
            }
        except ValueError as e:
            print(f"Error getting geometry input values: {e}")
            return None

    def set_input_values(self, input_values):

        """
        Sets the input fields based on provided values for wing geometry.
        """

        if 'n_panels_chord' in input_values:
            self.n_panels_chord_entry.delete(0, tk.END)
            self.n_panels_chord_entry.insert(0, input_values['n_panels_chord'])

        if 'n_panels_span' in input_values:
            self.n_panels_span_entry.delete(0, tk.END)
            self.n_panels_span_entry.insert(0, input_values['n_panels_span'])

        if 'span' in input_values:
            self.span_entry.delete(0, tk.END)
            self.span_entry.insert(0, input_values['span'])

        if 'root_chord' in input_values:
            self.root_chord_entry.delete(0, tk.END)
            self.root_chord_entry.insert(0, input_values['root_chord'])

        if 'tip_chord' in input_values:
            self.tip_chord_entry.delete(0, tk.END)
            self.tip_chord_entry.insert(0, input_values['tip_chord'])

        if 'sweep' in input_values:
            self.sweep_entry.delete(0, tk.END)
            self.sweep_entry.insert(0, input_values['sweep'])

        if 'dihedral' in input_values:
            self.dihedral_entry.delete(0, tk.END)
            self.dihedral_entry.insert(0, input_values['dihedral'])


