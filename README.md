# Pupil GLM

MATLAB analysis of pupil dilation (PD) during an emotional memory (EM) task,
comparing three groups: younger adults (YAs), older adults (OAs), and patients
with mild cognitive impairment (MCI).

The main script (`pupil_GLM_permutations.m`) runs a two-level GLM:

1. **First level (within subject):** each subject's single-trial pupil timeseries
   is regressed on task variables (emotionality, recognition, certainty,
   testing period, and two interactions), yielding one beta timeseries per
   regressor per subject.
2. **Second level (across subjects):** the per-subject beta timeseries are
   regressed on group indicators (YAs/OAs/MCI) plus, where applicable, an
   interindividual behavioral covariate. Significance is assessed with
   sign-flipping / label-shuffling permutation tests at both levels.

## Setup

Open `EM_pupil_GLM.m` and set one variable at the top:

```matlab
projroot = pwd;   % or an explicit path to the project folder
```

Every other path is built relative to `projroot`, so nothing else needs editing.

## Expected folder layout

```
projroot/
├── EM_pupil_GLM.m
├── eye_glm.m                          # not included – see below
├── ols.m                              # not included – see below
├── raacampbell-shadedErrorBar-19cf3fe/   # external dependency (see below)
├── fieldtrip/                            # external dependency (see below)
└── EM_task/
    ├── EM_cleaneddata_test/           # cleaned pupil data – not included
    │   ├── r<ID>_eye_rv200_1.mat
    │   └── r<ID>_eye_rv200_2.mat
    ├── behav_em/                      # encoding-session behaviour – not included
    │   └── <ID>/behav/EM/<ID>_emomem.mat
    ├── files_for_GLM/                 # behavioural covariates – not included
    │   ├── RTs_Emo_minus_Neu_T1_for_GLM.mat
    │   ├── RTs_Emo_minus_Neu_T2_for_GLM.mat
    │   ├── dprime_Emo_minus_Neu_T1_for_GLM.mat
    │   └── dprime_Emo_minus_Neu_T2_for_GLM.mat
    └── GLM_outputs/                   # figures are written here
```

## Not included in this repository

The following are required to run the script but are **not distributed here**:

- **Data files** — the cleaned/raw pupil `.mat` files
  (`EM_cleaneddata_test/`, `behav_em/`) and the behavioral covariate files
  (`files_for_GLM/`). These are not shared for data-protection reasons.
- **Analysis helpers** — `eye_glm.m` (first-level GLM) and `ols.m`
  (second-level OLS with contrasts). Place them somewhere on the MATLAB path
  (e.g. in `projroot`).

## External dependencies

- **FieldTrip** — used here for `ft_resampledata`.
  https://www.fieldtriptoolbox.org/
- **shadedErrorBar** (Rob Campbell) — used for the shaded SEM plots.
  https://github.com/raacampbell/shadedErrorBar

Install these separately.
