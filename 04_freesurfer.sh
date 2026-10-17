#!/bin/bash
#SBATCH -J recon_all
#SBATCH -t 24:00:00
#SBATCH --mem=32G
#SBATCH -c 8
#SBATCH -o logs/04_fs_%a.out
#SBATCH --array=0-0

# Stage 04: FreeSurfer segmentation of each session's T1. Stage 08 cuts its
# white matter and CSF masks from the labels this produces. FreeSurfer was
# chosen over SPM for the job because its labels are sharper and steadier
# across sessions, and a segmentation that wobbles session to session could
# masquerade as a cycle effect.
#
# Not done here: extra bias correction before recon-all. Its built-in
# correction has handled these T1s fine; revisit if QC shows patchy white
# matter in any session.

set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/config.sh"
load_freesurfer

ses_raw=$(session_from_array); ses=ses-${ses_raw#cycles}
t1=${bids_root}/${subject}/${ses}/anat/${subject}_${ses}_T1w.nii.gz
fs_id=${subject}_${ses}
mkdir -p "$SUBJECTS_DIR"

# FS 8 dropped -make all, so a partial dir cannot be resumed; delete and restart.
# recon-all.done appears only after a successful -all, so it marks a finished session.
done_flag=${SUBJECTS_DIR}/${fs_id}/scripts/recon-all.done
if [[ -d ${SUBJECTS_DIR}/${fs_id} && ! -e $done_flag ]]; then
    echo "partial recon dir for $fs_id; deleting and restarting"
    rm -rf "${SUBJECTS_DIR:?}/${fs_id:?}"
fi
if [[ -e $done_flag ]]; then
    echo "recon already done for $fs_id; skipping recon-all"
else
    recon-all -s "$fs_id" -i "$t1" -all -parallel -openmp 8
fi

# later stages read nifti, so convert the two volumes they need
mri_convert "${SUBJECTS_DIR}/${fs_id}/mri/aseg.mgz"  "${SUBJECTS_DIR}/${fs_id}/mri/aseg.nii.gz"
mri_convert "${SUBJECTS_DIR}/${fs_id}/mri/brain.mgz" "${SUBJECTS_DIR}/${fs_id}/mri/brain.nii.gz"

echo "stage 04 done: $fs_id"
