#!/bin/bash
#SBATCH -J recon_all
#SBATCH -t 24:00:00
#SBATCH --mem=16G
#SBATCH -c 8
#SBATCH -o logs/04_fs_%a.out
#SBATCH --array=0-0

# Stage 04: FreeSurfer segmentation of each session's T1. Stage 08 cuts its
# white matter and CSF masks from the labels this produces. FreeSurfer was
# chosen over SPM for the job because its labels are sharper and steadier
# across sessions.

set -euo pipefail
source "$(dirname "$0")/config.sh"
load_freesurfer

ses_raw=$(session_from_array); ses=ses-${ses_raw#cycles}
t1=${bids_root}/${subject}/${ses}/anat/${subject}_${ses}_T1w.nii.gz
fs_id=${subject}_${ses}
mkdir -p "$SUBJECTS_DIR"

if [[ -d ${SUBJECTS_DIR}/${fs_id} ]]; then
    echo "recon dir exists for $fs_id; resuming with -make all"
    recon-all -s "$fs_id" -make all -parallel -openmp 8
else
    recon-all -s "$fs_id" -i "$t1" -all -parallel -openmp 8
fi

# later stages read nifti, so convert the two volumes they need
mri_convert "${SUBJECTS_DIR}/${fs_id}/mri/aseg.mgz"  "${SUBJECTS_DIR}/${fs_id}/mri/aseg.nii.gz"
mri_convert "${SUBJECTS_DIR}/${fs_id}/mri/brain.mgz" "${SUBJECTS_DIR}/${fs_id}/mri/brain.nii.gz"

echo "stage 04 done: $fs_id"
