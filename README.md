# TEnKF Lorenz-96 Experiments

Reproducible code for the Lorenz-96 data assimilation experiments in:

> *Tempered Ensemble Kalman Filter (TEnKF)
> for convective-scale radar data assimilation*

---

## Structure

```
.
├── compile.sh                  ← step 0: build Fortran modules
├── run_nature.py               ← step 1: generate nature run + observations
├── run_example.py              ← step 2: run the reference experiment
├── reproduce_example.ipynb     ← step 3: load results and plot diagnostics
│
├── fortran/
│   ├── model/
│   │   └── lorenzN.f90         ← Lorenz-96 single-scale model
│   ├── da/
│   │   ├── common_da_tools_1d.f90  ← TEnKF / LETKF DA engine
│   │   ├── common_letkf.f90        ← LETKF core
│   │   └── common_obs_lorenzN.f90  ← observation operator
│   └── common/
│       ├── common_tools.f90    ← shared math/utility routines
│       ├── common_mtx.f90      ← matrix utilities
│       └── netlib.f90          ← linear algebra (LAPACK subset)
│
├── codes_experiments/
│   ├── assimilation_letkf_module.py ← Python DA driver
│   ├── nature_module.py             ← nature run driver
│   ├── sensitivity_conf_default.py  ← default DA configuration
│   └── default_nature_conf.py       ← default nature run configuration
│
├── environment.yml             ← conda environment
└── data/                       ← created at runtime, not tracked by git
    ├── Nature/
    └── output/
```

---

## Quickstart

### 1. Environment

```bash
conda env create -f environment.yml
conda activate tenkf_l96
```

### 2. Compile Fortran modules

```bash
bash compile.sh
```

### 3. Generate the nature run

```bash
python run_nature.py
```

Produces `data/Nature/Paper_Nature_Freq4_Den1.0_Type3_ObsErr5.npz`.
To generate additional σ_obs values (0.3, 1, 25), add them to
`ObsErrorList` in `run_nature.py`.

### 4. Run the reference experiment

```bash
python run_example.py
```

Runs the (inflation × localization) sweep for the reference configuration:

| Parameter | Value |
|-----------|-------|
| σ_obs | 5 |
| N_ens | 20 |
| Obs density | 1.0 |
| α | 2 |
| N_temp | 2 |

To reproduce a different condition, edit the `USER PARAMETERS` block at the
top of `run_example.py`.

### 5. Visualize results

Open `reproduce_example.ipynb` in Jupyter. It loads the output from step 4
and produces RMSE/spread heatmaps and calibration diagnostics.

```bash
jupyter notebook reproduce_example.ipynb
```

---

## Key method

The **TEnKF** splits the Bayesian update into `N_temp` pseudo-time steps
controlled by a tempering schedule parameterized by `α`. At `N_temp = 1`
the method reduces to the standard LETKF. The tempering schedule is computed
in `codes_experiments/assimilation_letkf_module.get_temp_steps()` and the
Fortran DA code is in `fortran/da/common_da_tools_1d.f90`.
