#!/bin/bash
# Shared settings for the cycles pipeline. Every numbered stage sources this
# file, so paths and the session list only ever get edited here.

# ---- people and places ----
lab_root=/orcd/data/ldlewis/001
work_root=${lab_root}/users/jpatel/cycles                    # everything we produce goes under here
raw_root=${work_root}/nifty                                  # converted scans; stage 01 finds them with or without per-session folders
bids_root=${work_root}/bids
deriv_root=${work_root}/derivatives
code_dir=${work_root}/code                                   # scripts + sessions.txt live here
log_dir=${work_root}/logs

# ---- subject / sessions ----
# One participant scanned many times. Session names match the raw folders
# (cycles006 becomes ses-006). Add a line per new session.
subject=sub-01
sessions_file=${code_dir}/sessions.txt   # one raw session name per line, e.g. cycles006

# ---- scan settings (checked against the sidecars in stage 01) ----
# stage 08 reads the real TR from the sidecar; these are only for checks.
expected_tr=0.45
n_slices=48
mb_factor=8

# ---- tool env ----
# pinned to what module avail reports on ORCD as of Sep 2026
load_fsl()        { module use /orcd/compute/bcs/001/modulefiles; module load fsl/6.0.7; }
load_freesurfer() { module use /orcd/compute/bcs/001/modulefiles; module load freesurfer/8.0.0; export SUBJECTS_DIR=${deriv_root}/freesurfer; }
load_ants()       { module use /orcd/compute/bcs/001/modulefiles; module load ants/2.1.0-3.8bed08; }
load_matlab()     { module use /orcd/compute/bcs/001/modulefiles; module load matlab/matlab-2023b; }
spm_path=${lab_root}/software/spm12

# ---- where each stage writes ----
topup_dir=${deriv_root}/topup
spm_dir=${deriv_root}/spm_preproc
fs_dir=${deriv_root}/freesurfer
ants_dir=${deriv_root}/ants_template
reg_dir=${deriv_root}/func_anat_reg
atlas_dir=${deriv_root}/atlas
qc_dir=${deriv_root}/qc
ts_dir=${deriv_root}/timeseries
conn_dir=${deriv_root}/connectivity

# turn the SLURM array index into a session name
session_from_array() {
    sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$sessions_file"
}

mkdir -p "$log_dir"
