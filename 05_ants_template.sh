#!/bin/bash
#SBATCH -J ants_template
#SBATCH -t 48:00:00
#SBATCH --mem=32G
#SBATCH -c 16
#SBATCH -o logs/05_ants.out

# Stage 05: average all the session T1s into one reference brain for this
# participant, then work out how each session's T1 maps onto that average,
# and how the standard MNI brain maps onto it. With one person scanned many
# times, their own average is a better common space than a population brain.
#
# Files left behind for stage 06b to chain together:
#   MNI -> template:      mni_to_tpl_1Warp.nii.gz + mni_to_tpl_0GenericAffine.mat
#   session T1 -> tpl:    ${ses}_to_tpl_1Warp / _0GenericAffine  (+ inverses)

set -euo pipefail
source "$(dirname "$0")/config.sh"
load_ants

mkdir -p "${ants_dir}"; cd "${ants_dir}"

# use the skull-stripped T1s so skull and neck don't drag the average around
t1_list=()
while read -r ses_raw; do
    ses=ses-${ses_raw#cycles}
    brain=${fs_dir}/${subject}_${ses}/mri/brain.nii.gz
    [[ -f $brain ]] || { echo "missing $brain, run stage 04 first"; exit 1; }
    cp -n "$brain" "t1_${ses}.nii.gz"
    t1_list+=("t1_${ses}.nii.gz")
done < "$sessions_file"

if [[ ! -f tpl_template0.nii.gz ]]; then
    antsMultivariateTemplateConstruction2.sh -d 3 -o tpl_ \
        -i 4 -g 0.2 -c 2 -j 16 -k 1 -t SyN -m CC -n 0 "${t1_list[@]}"
fi
tpl=tpl_template0.nii.gz

# map each session's T1 onto the average
while read -r ses_raw; do
    ses=ses-${ses_raw#cycles}
    [[ -f ${ses}_to_tpl_1Warp.nii.gz ]] && continue
    antsRegistrationSyN.sh -d 3 -f "$tpl" -m "t1_${ses}.nii.gz" \
        -o "${ses}_to_tpl_" -n 16 -t s
done < "$sessions_file"

# map the standard brain onto the average too; the atlas ships in that space
mni=${FSLDIR:-/orcd/software/fsl}/data/standard/MNI152_T1_1mm_brain.nii.gz
if [[ ! -f mni_to_tpl_1Warp.nii.gz ]]; then
    antsRegistrationSyN.sh -d 3 -f "$tpl" -m "$mni" -o mni_to_tpl_ -n 16 -t s
fi

echo "stage 05 done: template + $(wc -l < "$sessions_file") session registrations"
