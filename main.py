######################################

# Main Page

######################################

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import json
import pandas as pd
from vlm_bertin import VLM
from geometry import Geometry
from freestream import Freestream
from fem import FEMAnalysisAlt
from clt import CLT

class HomePage(ttk.Frame):

    def __init__(self, parent, main_app, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.main_app = main_app

        self.pack(fill=tk.BOTH, expand=True)

        intro_text = (
            "This is a static aeroelastic model for lifting surface analysis in preliminary design stages.\n\n"
            "It features a VLM aerodynamic model and a spatial beam element FEM structural model.\n\n"
            "To begin, you can load an existing study or start a new one by configuring the inputs.\n\n"
        )
        label = ttk.Label(self, text=intro_text, font=("Arial", 14))
        label.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")

        save_button = ttk.Button(self, text="Save Study", command=self.save_input_fields)
        save_button.grid(row=1, column=0, padx=5, pady=5, sticky="ew")

        load_button = ttk.Button(self, text="Load Study", command=self.load_input_fields)
        load_button.grid(row=2, column=0, padx=5, pady=5, sticky="ew")

        export_button = ttk.Button(self, text="Export Results", command=self.export_results)
        export_button.grid(row=3, column=0, padx=5, pady=5, sticky="ew")

    def save_input_fields(self):

        file_path = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON files", "*.json")])
        if file_path:
            input_fields = {}

            if hasattr(self.main_app, 'freestream'):
                input_fields['freestream'] = self.main_app.freestream.get_input_values()
            if hasattr(self.main_app, 'geometry'):
                input_fields['geometry'] = self.main_app.geometry.get_input_values()
            if hasattr(self.main_app, 'fem'):
                input_fields['fem'] = self.main_app.fem.get_input_values()
            if hasattr(self.main_app, 'vlm') and hasattr(self.main_app.vlm, 'get_input_values'):
                input_fields['vlm'] = self.main_app.vlm.get_input_values()
            if hasattr(self.main_app, 'clt') and hasattr(self.main_app.clt, 'get_input_values'):
                input_fields['clt'] = self.main_app.clt.get_input_values()

            with open(file_path, 'w') as file:
                json.dump(input_fields, file)

    def load_input_fields(self):

        file_path = filedialog.askopenfilename(defaultextension=".json", filetypes=[("JSON files", "*.json")])
        if file_path:
            try:
                with open(file_path, 'r') as file:
                    input_fields = json.load(file)

                if 'freestream' in input_fields and hasattr(self.main_app, 'freestream'):
                    self.main_app.freestream.set_input_values(input_fields['freestream'])
                if 'geometry' in input_fields and hasattr(self.main_app, 'geometry'):
                    self.main_app.geometry.set_input_values(input_fields['geometry'])
                if 'fem' in input_fields and hasattr(self.main_app, 'fem'):
                    self.main_app.fem.set_input_values(input_fields['fem'])
                if 'vlm' in input_fields and hasattr(self.main_app, 'vlm') and hasattr(self.main_app.vlm,
                                                                                                 'set_input_values'):
                    self.main_app.clt.set_input_values(input_fields['clt'])
                if 'clt' in input_fields and hasattr(self.main_app, 'clt') and hasattr(self.main_app.clt,
                                                                                                 'set_input_values'):
                    self.main_app.clt.set_input_values(input_fields['clt'])

            except FileNotFoundError:
                print("No saved input fields found.")

    def export_results(self):

        """
        Exports aerodynamic and structural results to an Excel or text file,
        including displacements, twist, strains, and stresses.
        """

        file_path = filedialog.asksaveasfilename(
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx"), ("Text files", "*.txt")],
        )
        if not file_path:
            return

        vlm_results = self.main_app.vlm.get_results()
        fem_results = self.main_app.fem.get_results()

        spanwise_positions_vlm = vlm_results.get("spanwise_positions", [])
        circulation = vlm_results.get("circulation_distribution", [])
        lift = vlm_results.get("lift_distribution", [])
        cl = vlm_results.get("cl_distribution", [])
        Di = vlm_results.get("Di_distribution", [])

        spanwise_positions_fem = fem_results.get("spanwise_positions", [])
        displacements = fem_results.get("displacements", [])
        twist = fem_results.get("twist", [])
        strains = fem_results.get("strains", [])
        stresses = fem_results.get("stresses", [])

        if file_path.endswith(".xlsx"):
            with pd.ExcelWriter(file_path) as writer:
                aero_df = pd.DataFrame({
                    "Spanwise Position (m)": spanwise_positions_vlm,
                    "Circulation (m^2/s)": circulation,
                    "Lift (N)": lift,
                    "CL": cl,
                    "Induced Drag (N)": Di,
                })
                aero_df.to_excel(writer, sheet_name="Aerodynamic Results", index=False)

                struct_df = pd.DataFrame({
                    "Spanwise Position (m)": spanwise_positions_fem,
                    "X-Displacement (m)": displacements[0::6],
                    "Y-Displacement (m)": displacements[1::6],
                    "Z-Displacement (m)": displacements[2::6],
                    "Twist (deg)": twist,
                    "Strain ε_xx": [strain[0] for strain in strains],
                    "Strain ε_yy": [strain[1] for strain in strains],
                    "Strain γ_xy": [strain[3] for strain in strains],
                    "Stress σ_xx (MPa)": [stress[0] for stress in stresses],
                    "Stress σ_yy (MPa)": [stress[1] for stress in stresses],
                    "Shear Stress τ_xy (MPa)": [stress[2] for stress in stresses],
                })
                struct_df.to_excel(writer, sheet_name="Structural Results", index=False)

        elif file_path.endswith(".txt"):
            with open(file_path, "w") as file:
                file.write("Aerodynamic Results:\n")
                file.write("Spanwise Position (m), Circulation (m^2/s), Lift (N), CL, Induced Drag\n")
                for i in range(len(spanwise_positions_vlm)):
                    file.write(
                        f"{spanwise_positions_vlm[i]:.6f}, {circulation[i]:.6f}, {lift[i]:.6f}, {cl[i]:.6f}, {Di[i]:.6f}\n"
                    )

                file.write("\nStructural Results:\n")
                file.write("Spanwise Position (m), X-Displacement (m), Y-Displacement (m), Z-Displacement (m), "
                           "Twist (deg), Strain ε_xx, Strain ε_yy, Strain γ_xy, "
                           "Stress σ_xx (MPa), Stress σ_yy (MPa), Shear Stress τ_xy (MPa)\n")
                for i in range(len(spanwise_positions_fem)):
                    file.write(
                        f"{spanwise_positions_fem[i]:.6f}, "
                        f"{displacements[6 * i]:.6f}, {displacements[6 * i + 1]:.6f}, {displacements[6 * i + 2]:.6f}, "
                        f"{twist[i]:.6f}, {strains[i][0]:.6e}, {strains[i][1]:.6e}, {strains[i][3]:.6e}, "
                        f"{stresses[i][0]:.6e}, {stresses[i][1]:.6e}, {stresses[i][2]:.6e}\n"
                    )


class AeroelasticApplication(tk.Tk):

    def __init__(self):
        super().__init__()

        self.title("Aeroelastic Model")
        self.geometry("1000x800")

        style = ttk.Style(self)
        style.configure('TNotebook.Tab', font=('Arial', 14))

        self.notebook = ttk.Notebook(self, style='TNotebook')
        self.notebook.pack(fill=tk.BOTH, expand=True)

        home_page = HomePage(self.notebook, self)
        self.notebook.add(home_page, text="Home Page")

        self.freestream = Freestream(self.notebook)
        self.notebook.add(self.freestream, text="Freestream Conditions")

        self.geometry = Geometry(self.notebook)
        self.notebook.add(self.geometry, text="Wing Geometry")

        self.clt = CLT(self.notebook)
        self.notebook.add(self.clt, text="Classical Laminate Theory")

        self.vlm = VLM(self.notebook, self.geometry, self.freestream)
        self.notebook.add(self.vlm, text="Vortex Lattice Method")

        self.fem = FEMAnalysisAlt(self.notebook, self.freestream, self.geometry, self.vlm, self.clt)
        self.notebook.add(self.fem, text="Finite Element Method")

if __name__ == "__main__":
    app = AeroelasticApplication()
    app.mainloop()


