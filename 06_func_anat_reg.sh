#!/bin/bash
#SBATCH -J func_reg
#SBATCH -t 02:00:00
#SBATCH --mem=8G
#SBATCH -o logs/06_reg_%a.out
#SBATCH --array=0-0

# Stage 06: line each run's mean image up with that session's T1 using
# bbregister, which snaps the fit to the gray/white boundary from the
# FreeSurfer run and holds up better than plain intensity matching at this
# resolution. The result is saved in a format ANTs can chain in stage 06b.
# The functional data itself never gets resampled: the atlas comes to it.

set -euo pipefail
source "$(dirname "$0")/config.sh"
load_freesurfer
load_ants

ses_raw=$(session_from_array); ses=ses-${ses_raw#cycles}
fs_id=${subject}_${ses}
out=${reg_dir}/${subject}/${ses}
mkdir -p "$out"

for run in 1 2; do
    # align using the run's mean image, computed before smoothing
    mean=${spm_dir}/${subject}/${ses}/mean_run-${run}.nii.gz
    if [[ ! -f $mean ]]; then
        load_fsl
        fslmaths "${spm_dir}/${subject}/${ses}/rrun-${run}_bold.nii.gz" -Tmean "$mean"
    fi

    bbregister --s "$fs_id" --mov "$mean" --reg "${out}/run-${run}_bbr.dat" \
               --init-coreg --bold --lta "${out}/run-${run}_bbr.lta"

    # save the alignment as a plain text matrix for stage 06b to chain
    lta_convert --inlta "${out}/run-${run}_bbr.lta" \
                --outitk "${out}/run-${run}_func_to_anat.txt"

    # bbregister scores its own fit; above about 0.8 the alignment probably
    # failed, so print the score where stage 07 can collect it
    cost=$(cat "${out}/run-${run}_bbr.dat.mincost" | awk '{print $1}')
    echo "$ses run-$run bbr mincost = $cost"
done

echo "stage 06 done: $ses"
