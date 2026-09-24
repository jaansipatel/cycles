#!/bin/bash
#SBATCH -J spm_preproc
#SBATCH -t 06:00:00
#SBATCH --mem=24G
#SBATCH -c 4
#SBATCH -o logs/03_spm_%a.out
#SBATCH --array=0-0

# Stage 03: SPM motion correction plus 4 mm smoothing on the topup-corrected runs.
#
# No unwarp: stage 02 already fixed the distortion, so this step should only need
# to track head motion, which stage 08 reuses as noise regressors.

set -euo pipefail
source "$(dirname "$0")/config.sh"
load_matlab

ses_raw=$(session_from_array); ses=ses-${ses_raw#cycles}
in_dir=${topup_dir}/${subject}/${ses}
out=${spm_dir}/${subject}/${ses}
mkdir -p "$out"

# SPM edits files where they sit and wants them unzipped, so work on copies
for run in 1 2; do
    src=${in_dir}/${subject}_${ses}_dir-AP_task-rest_run-${run}_bold_topup.nii.gz
    tgt=${out}/run-${run}_bold.nii
    [[ -f $tgt ]] || { cp "$src" "${tgt}.gz"; gunzip "${tgt}.gz"; }
done

matlab -nodisplay -nosplash <<MAT
addpath('${spm_path}');
spm('defaults','fmri'); spm_jobman('initcfg');

r1 = spm_select('expand', '${out}/run-1_bold.nii');
r2 = spm_select('expand', '${out}/run-2_bold.nii');

% motion correct both runs together, write corrected copies plus a mean
mb.spm.spatial.realign.estwrite.data = {cellstr(r1); cellstr(r2)};
mb.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
mb.spm.spatial.realign.estwrite.eoptions.sep = 4;
mb.spm.spatial.realign.estwrite.eoptions.rtm = 0;    % align to the first image, not the mean: faster and just as good here
mb.spm.spatial.realign.estwrite.eoptions.interp = 2;
mb.spm.spatial.realign.estwrite.roptions.which = [2 1];
mb.spm.spatial.realign.estwrite.roptions.interp = 4;
spm_jobman('run', {mb}); clear mb;

% smooth at 4 mm, same as the reference paper
for run = 1:2
    src = spm_select('expand', sprintf('${out}/rrun-%d_bold.nii', run));
    mb.spm.spatial.smooth.data = cellstr(src);
    mb.spm.spatial.smooth.fwhm = [4 4 4];
    mb.spm.spatial.smooth.prefix = 's';
    spm_jobman('run', {mb}); clear mb;
end
exit;
MAT

# motion estimates per run, renamed so stage 08 can find them
mv "${out}/rp_run-1_bold.txt" "${out}/motion_run-1.txt"
mv "${out}/rp_run-2_bold.txt" "${out}/motion_run-2.txt"

for run in 1 2; do
    gzip -f "${out}/srrun-${run}_bold.nii"
    gzip -f "${out}/rrun-${run}_bold.nii"
done
gzip -f "${out}/meanrun-1_bold.nii" 2>/dev/null || gzip -f "${out}"/mean*.nii

echo "stage 03 done: $ses"
