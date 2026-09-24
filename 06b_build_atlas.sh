#!/bin/bash
#SBATCH -J atlas_115
#SBATCH -t 02:00:00
#SBATCH --mem=8G
#SBATCH -o logs/06b_atlas_%a.out
#SBATCH --array=0-0

# Stage 06b: build the 115-region atlas once (100 Schaefer cortical regions
# plus 15 subcortical ones from Harvard-Oxford), write the labels table that
# stage 08 treats as the source of truth, then carry the atlas into each
# run's own space through the chain of alignments:
#   MNI -> subject average -> session T1 -> functional
# The atlas is what moves; the functional data never gets resampled.

set -euo pipefail
source "$(dirname "$0")/config.sh"
load_fsl
load_ants

mkdir -p "${atlas_dir}"
schaefer=${atlas_dir}/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz
schaefer_lut=${atlas_dir}/Schaefer2018_100Parcels_7Networks_order.txt
combined=${atlas_dir}/atlas115_mni.nii.gz
labels=${atlas_dir}/atlas115_labels.tsv

# ---- build the combined atlas once ----
if [[ ! -f $combined ]]; then
    base=https://raw.githubusercontent.com/ThomasYeoLab/CBIG/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations
    [[ -f $schaefer ]]     || curl -fsSL -o "$schaefer"     "$base/MNI/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii.gz"
    [[ -f $schaefer_lut ]] || curl -fsSL -o "$schaefer_lut" "$base/MNI/freeview_lut/Schaefer2018_100Parcels_7Networks_order.txt"
    ho=${FSLDIR}/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr25-1mm.nii.gz

    python3 - "$schaefer" "$schaefer_lut" "$ho" "$combined" "$labels" <<'PY'
import nibabel as nib, numpy as np, sys
sch_f, lut_f, ho_f, out_f, tsv_f = sys.argv[1:6]

sch = nib.load(sch_f); atlas = np.asarray(sch.dataobj).astype(np.int16)
ho  = np.asarray(nib.load(ho_f).dataobj).astype(np.int16)

# the Harvard-Oxford labels we keep: 7 structures per hemisphere plus the
# brainstem, numbered 101 to 115. That atlas also labels white matter,
# cortex and ventricles, which we skip.
keep = [(10,"L","Thalamus"), (11,"L","Caudate"), (12,"L","Putamen"), (13,"L","Pallidum"),
        (9,"L","Hippocampus"), (18,"L","Amygdala"), (26,"L","Accumbens"),
        (49,"R","Thalamus"), (50,"R","Caudate"), (51,"R","Putamen"), (52,"R","Pallidum"),
        (48,"R","Hippocampus"), (54,"R","Amygdala"), (58,"R","Accumbens"),
        (8,"B","Brainstem")]
rows = []
for i, line in enumerate(open(lut_f)):
    p = line.split()
    if len(p) < 2: continue
    idx, name = int(p[0]), p[1]                      # e.g. 7Networks_LH_Vis_1
    _, hemi, net = name.split("_")[0], name.split("_")[1], name.split("_")[2]
    rows.append((idx, name, net, "L" if hemi == "LH" else "R"))

for j, (ho_idx, hemi, name) in enumerate(keep):
    node = 101 + j
    mask = (ho == ho_idx) & (atlas == 0)             # where the two atlases overlap, cortex wins
    atlas[mask] = node
    rows.append((node, f"Subcortical_{hemi}_{name}", "Subcortical", hemi))

assert len(rows) == 115 and sorted(r[0] for r in rows) == list(range(1, 116))
nib.save(nib.Nifti1Image(atlas, sch.affine, sch.header), out_f)
with open(tsv_f, "w") as f:
    f.write("node_id\tname\tnetwork\themisphere\n")
    for r in sorted(rows): f.write("\t".join(map(str, r)) + "\n")
print("atlas115 written:", int((atlas > 0).sum()), "labeled voxels")
PY
fi

# ---- carry the atlas into this session's runs ----
ses_raw=$(session_from_array); ses=ses-${ses_raw#cycles}
out=${atlas_dir}/${subject}/${ses}
mkdir -p "$out"

for run in 1 2; do
    ref=${spm_dir}/${subject}/${ses}/mean_run-${run}.nii.gz
    bbr=${reg_dir}/${subject}/${ses}/run-${run}_func_to_anat.txt
    antsApplyTransforms -d 3 -i "$combined" -r "$ref" \
        -o "${out}/atlas115_run-${run}_native.nii.gz" -n GenericLabel \
        -t ["$bbr",1] \
        -t ["${ants_dir}/${ses}_to_tpl_0GenericAffine.mat",1] \
        -t "${ants_dir}/${ses}_to_tpl_1InverseWarp.nii.gz" \
        -t "${ants_dir}/mni_to_tpl_1Warp.nii.gz" \
        -t "${ants_dir}/mni_to_tpl_0GenericAffine.mat"

    # bring the FreeSurfer segmentation along too. It already lives in this
    # session's T1 space, so only the last step of the chain applies. Stage
    # 08 cuts its white matter and CSF masks from this file.
    antsApplyTransforms -d 3 -i "${fs_dir}/${subject}_${ses}/mri/aseg.nii.gz" \
        -r "$ref" -o "${out}/aseg_run-${run}_native.nii.gz" -n GenericLabel \
        -t ["$bbr",1]

    n=$(fslstats "${out}/atlas115_run-${run}_native.nii.gz" -R | awk '{print $2}')
    [[ ${n%.*} -eq 115 ]] || echo "WARNING: max label ${n} != 115 in run-${run}, nodes lost in warp"
done

echo "stage 06b done: $ses"
