# NUC.REAC.PHY

This directory contains all the files and scripts developed to simulate a **Pressurized Water Reactor (PWR)** using *OpenMC*. The different subfolders correspond to specific stages and parameter studies carried out during the project.

> ⚠️ **Note:** The most important folder is **`Reactor_core/`**. It contains the core implementation and summarizes the work performed in all the other subfolders. For a clear and complete understanding of the project, it is recommended to focus on this directory.

---

## Folder structure

### **1. `Change_Enrichment/`**
This folder contains the code used to run *pin-cell simulations* with different **fuel types** and **fuel enrichments**. The main file is `Project3.py`.
It also includes:
- All `.xml` input files (geometry, materials, settings, tallies)
- Output files and results obtained for each enrichment case
- Data summarized and discussed in the main project report
Additional Python scripts are included for post-processing:
- `postprocess_energy.py`: generates the *Neutron energy spectrum in the fuel*
- `postprocess.py`: generates the *Neutron flux distribution / Fission site heatmap*

---

### **2. `Change_moderator/`**
Follows the same logic as `Change_Enrichment`, but here the **moderator material** is varied instead of the fuel. The main file is `Project2.py`. 


---

### **3. `cross_section_plot/`**
Contains the Python code and data needed to produce the **macroscopic cross-section plots** of **Boron Carbide (B₄C)**.  
Includes:
- `.txt` data files with energy and cross-section values  
- The corresponding plotting script and the resulting figure

---

### **4. `Light_water/`**
Holds the scripts used to simulate a *pin-cell* configuration with **light water (H₂O)** as moderator.  
This folder corresponds to the first stage of the project and contains the initial results presented in the report.

---

### **5. `Reactor_core/`**
This is the **main and most complete** part of the project. It integrates and summarizes the work done in the other folders.  
Contents include:
- **Images** of the reactor assemblies and the full reactor core  
- **Important Python scripts:**
  - `core.py`: **main script** — builds the core model using helper functions from `pins.py` and `utilities.py`
  - `assembly.py`: generates an assembly and the corresponding figure `assembly_xy.png`
  - `keff_vs_T.py`: produces the plot `keff_vs_T.pdf`, showing how the effective multiplication factor varies with temperature
- **XML files** with full definitions of materials, geometry, settings, tallies, and the depletion chain  
- **PDF results** containing key figures and output summaries

---

## Recommendations

- Focus on **`Reactor_core/`**, which contains the complete, functional, and well-organized implementation of the PWR simulation.  
- Other folders can be ignored.

---

