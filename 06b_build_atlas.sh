#!/bin/bash
#SBATCH -J atlas_115
#SBATCH -t 02:00:00
#SBATCH --mem=8G
#SBATCH -o logs/06b_atlas_%a.out
#SBATCH --array=0-0

# Stage 06b: build the 115-region atlas once (100 Schaefer cortical regions
# plus 15 subcortical ones from Harvard-Oxford), then carry the labels into
# each run's own space through the chain of alignments:
#   MNI -> subject average -> session T1 -> functional

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
# cortex and ventricles, which we skip. The numbers are the voxel values
# in the sub-maxprob image, which run 1 to 21 in the order of the
# HarvardOxford-Subcortical.xml entries. An earlier version of this list
# carried FreeSurfer aseg numbers by mistake, which grabbed the wrong
# structures below 22 and nothing at all above it.
keep = [(4,"L","Thalamus"),  (5,"L","Caudate"),  (6,"L","Putamen"),  (7,"L","Pallidum"),
        (9,"L","Hippocampus"), (10,"L","Amygdala"), (11,"L","Accumbens"),
        (15,"R","Thalamus"), (16,"R","Caudate"), (17,"R","Putamen"), (18,"R","Pallidum"),
        (19,"R","Hippocampus"), (20,"R","Amygdala"), (21,"R","Accumbens"),
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
    assert mask.sum() > 0, f"node {node} ({hemi} {name}, HO value {ho_idx}) got no voxels"
    print(f"node {node} {hemi} {name}: {int(mask.sum())} voxels")
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

    # bring the FreeSurfer segmentation along too
    antsApplyTransforms -d 3 -i "${fs_dir}/${subject}_${ses}/mri/aseg.nii.gz" \
        -r "$ref" -o "${out}/aseg_run-${run}_native.nii.gz" -n GenericLabel \
        -t ["$bbr",1]

    # the old check here only looked at the maximum label, which stays 115
    # as long as the brainstem survives, so missing nodes slipped past it
    missing=$(fslstats "${out}/atlas115_run-${run}_native.nii.gz" -H 116 0 116 \
              | awk 'NR > 1 && $1 == 0 {printf "%d ", NR - 1}')
    [[ -z $missing ]] || echo "WARNING: run-${run} nodes with no voxels after warp: ${missing}"
done

echo "stage 06b done: $ses"
