# TFG – DC-OPF-UC & Unit Commitment (Julia / JuMP)

This repository contains the implementation of several **power system optimization models** developed in **Julia** using **JuMP**, as part of a **Bachelor’s Thesis (TFG)**.

The project focuses on the **comparison of different modeling granularities** for power system operation, ranging from **Copper Plate Unit Commitment** models to **DC Optimal Power Flow with Unit Commitment** on multi-bus systems.

The objective is to analyze how modeling detail (network constraints, unit commitment formulation, and temporal constraints) affects system operation costs and dispatch decisions.

---

## 🔧 Technologies and Tools
<details>
<summary> 👉Click here to see more details</summary>
- **Julia** (v1.11 recommended)
- **JuMP** optimization framework
- **HiGHS**, **Gurobi**, **Ipopt** solvers
- CSV-based input data
- Multi-period optimization (24h horizon)
</details>


## 🚀 Run the Project Online (MyBinder)
<details>
<summary> 👉Click here to see more details</summary>
Click the badge below to launch the project directly in your browser using **MyBinder** (no local installation required):

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mjimenezdevega/TFG_OPF_UC/main?filepath=OPF-Multi-Period)

### 🔹 What is MyBinder?

MyBinder allows you to run this project **online in your browser**, without installing Julia or any dependencies locally.  
It automatically uses the `Project.toml` and `Manifest.toml` files included in this repository.

---

### 🔹 How to run the code on MyBinder
<br>

1️⃣ Click on the **MyBinder badge** above  
2️⃣ Wait until the environment is built (first launch may take a few minutes)  
3️⃣ A **JupyterLab interface** will open in your browser  
4️⃣ Open a **New -> Terminal** inside JupyterLab  
5️⃣ Navigate to the project folder: cd OPF-Multi-Period  
6️⃣ Run the main Julia script: julia --project=../binder main.jl  
7️⃣ Follow the interactive menu:  
    - Select study case  
    - Select OPF type (in this proyect you have to select DC-OPF for all the cases excep UC Copper that has ist own)  
    - Select solver → HiGHS is recommended for online runs  
    The optimization results will be displayed directly in the terminal.

⚠️ Notes  
• MyBinder sessions are temporary  
• Computation time may be slower than local execution  
• Commercial solvers (e.g. Gurobi) may not be available on MyBinder

</details>

## ▶️ How to Run the Code

<details>
<summary> 👉Click here to see more details</summary>

### 1️⃣ Open the integrated terminal

In **VS Code**, navigate to the project folder and **right-click** on: OPF-Multi-period

Then select: Open integrated terminal

Make sure the terminal is opened inside the `OPF-Multi-Period` directory.

---

### 2️⃣ Run the main Julia script

In the integrated terminal, execute: julia main.jl
---
</details>

## 📁 Case Studies Structure and Input Data

<details>
<summary> 👉Click here to see more details</summary>

This folder contains all the **case studies** used in the project, including input data files and system 
definitions for each simulation scenario.

Each subfolder inside `Cases/` represents an **independent study case**, with its own demand, generation,
and (if applicable) network data.

---

## 🔹 UC Copper 24h

**Copper Plate Unit Commitment (24-hour horizon)**

This case represents a **single-node (copper plate)** system, where network constraints are ignored and only **generation and demand balance** is considered.

### Files included:

#### 📄 `generatorData.csv`
Defines the **thermal generation technologies** aggregated by type.

| Column | Description |
|------|------------|
| `tech` | Generator technology name (e.g. CC, Nuclear, Biomasa) |
| `Ng_total` | Total number of units available |
| `Pub_MW` | Installed capacity (informative) |
| `Pmin_MW` | Minimum power per unit |
| `Pmax_MW` | Maximum power per unit |
| `Cm_EUR_MWh` | Marginal generation cost |
| `CnL_EUR_h` | No-load cost per unit |
| `C_start_EUR` | Startup cost per unit |
| `RR_MW_per_h` | Ramp rate |
| `Tst_h` | Startup time |
| `Tmut_h` | Minimum up time |
| `Tmdt_h` | Minimum down time |
| `Ng0_units` | Number of units online at hour 0 |
| `Pg0_total_MW` | Initial total power output |

---

#### 📄 `demandData.csv`
Hourly system demand.

| Column | Description |
|------|------------|
| `PD_MW` | System demand (MW) |

---

#### 📄 `solarData.csv`
Hourly solar generation availability.

| Column | Description |
|------|------------|
| `PSolar_MW` | Available solar power (MW) |

---

#### 📄 `windData.csv`
Hourly wind generation availability.

| Column | Description |
|------|------------|
| `PWind_MW` | Available wind power (MW) |

---

## 🔹 3 Nudos UC

**Three-bus system with Unit Commitment and DC-OPF**

This case introduces **network constraints**, modelling power flows between buses using a DC approximation.

### Files included:

#### 📄 `busData.csv`
Defines buses and their demand.

| Column | Description |
|------|------------|
| `bus_id` | Bus identifier |
| `PD_MW` | Demand at the bus |

---

#### 📄 `generatorData.csv`
Defines individual generators connected to buses.

| Column | Description |
|------|------------|
| `gen_id` | Generator identifier |
| `bus_id` | Connected bus |
| `Pmin_MW` | Minimum power |
| `Pmax_MW` | Maximum power |
| `Cm_EUR_MWh` | Marginal cost |
| `CnL_EUR_h` | No-load cost |
| `C_start_EUR` | Startup cost |
| `RR_MW_per_h` | Ramp rate |
| `Tst_h` | Startup time |
| `Tmut_h` | Minimum up time |
| `Tmdt_h` | Minimum down time |

---

#### 📄 `branchData.csv`
Transmission line data.

| Column | Description |
|------|------------|
| `from_bus` | Origin bus |
| `to_bus` | Destination bus |
| `x_pu` | Line reactance |
| `rate_MW` | Thermal limit |

---

## 🔹 3 Nudos Ampliación

**Three-bus system with transmission expansion planning**

This case extends the previous one by allowing **investment decisions** on transmission lines.

### Additional files:

#### 📄 `expansionData.csv`
Defines candidate transmission lines.

| Column | Description |
|------|------------|
| `from_bus` | Origin bus |
| `to_bus` | Destination bus |
| `cost_EUR` | Investment cost |
| `rate_MW` | Line capacity |

---

## 🔹 19 Nudos

**Multi-bus system with DC-OPF and Unit Commitment**

This is the **most detailed case**, used as the reference model for comparing the Copper Plate formulation.

### Files included:

- `busData.csv`
- `generatorData.csv`
- `branchData.csv`
- `solarData.csv`
- `windData.csv`

Each generator is connected to a specific bus, and all **network, ramping, startup, and minimum time constraints** 
are enforced.

---

## 📌 Notes

- All input files are in **CSV format** for transparency and reproducibility.
- Each case can be run independently via the interactive menu in `main.jl`.
- Results are automatically saved in the `Results/` folder with timestamps.

---

## 🎓 Purpose within this proyect

These case studies are designed to **compare different modeling granularities**, from simplified Copper Plate 
Unit Commitment to detailed DC-OPF-UC formulations, analyzing:

- Total system cost
- Commitment schedules
- Impact of network constraints
- Differences in marginal prices and dispatch

</details>
