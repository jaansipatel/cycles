#!/usr/bin/env python3
#SBATCH -J extract_ts
#SBATCH -t 02:00:00
#SBATCH --mem=24G
#SBATCH -o logs/08_ts_%a.out
#SBATCH --array=0-0

# Stage 08: turn each preprocessed run into clean signals for the 115 atlas
# regions, following the reference code step for step: put the run on a
# common intensity scale, pull one signal per region, regress out motion and
# the white matter and CSF signals, then keep only the slow frequency band
# the analysis uses.
# Writes one table per run (volumes x 115) plus a json recording every
# setting. Columns follow the labels table from stage 06b exactly.

import json, os, subprocess, glob
import numpy as np
import nibabel as nib
from scipy.ndimage import binary_erosion
import pritschet_core as pc

# ---- knobs ----
wavelet_scales = [3, 4, 5, 6]   # the settled band, 0.017 to 0.278 Hz at this
                                # TR. The solver picks 3-7 for the paper's
                                # exact band; scale 7 was dropped on purpose,
                                # it sits mostly below 0.01 Hz and costs a lot
                                # of usable signal on runs this short
extract_method = "eig"          # "eig" matches the reference; "mean" is the plain fallback
erode_iters = 1                 # shave one voxel off the mask edges to keep gray matter out
min_node_vox = 5                # regions smaller than this give a shaky signal
wm_labels = (2, 41)
csf_labels = (4, 43, 14, 15)    # the ventricles only, not CSF at the brain surface


def cfg(*names):
    here = os.path.dirname(os.path.abspath(__file__))
    out = subprocess.check_output(
        ["bash", "-c", f"source {here}/config.sh >/dev/null 2>&1; " +
         " ".join(f'echo ${n};' for n in names)], text=True).splitlines()
    return out if len(names) > 1 else out[0]


def roi_mean(data2d_getter, mask):
    return data2d_getter(mask).mean(axis=1)


def main():
    (subject, sessions_file, bids_root, spm_dir, atlas_dir, ts_dir) = cfg(
        "subject", "sessions_file", "bids_root", "spm_dir", "atlas_dir", "ts_dir")
    idx = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
    ses_raw = open(sessions_file).read().split()[idx]
    ses = "ses-" + ses_raw.removeprefix("cycles")
    out = os.path.join(ts_dir, subject, ses)
    os.makedirs(out, exist_ok=True)

    labels_tsv = os.path.join(atlas_dir, "atlas115_labels.tsv")
    node_ids = [int(l.split("\t")[0]) for l in open(labels_tsv).read().splitlines()[1:]]
    assert node_ids == list(range(1, 116)), "labels TSV is the source of truth and it disagrees"

    for run in (1, 2):
        sesdir = os.path.join(spm_dir, subject, ses)
        adir = os.path.join(atlas_dir, subject, ses)

        tr = json.load(open(os.path.join(
            bids_root, subject, ses, "func",
            f"{subject}_{ses}_dir-AP_task-rest_run-{run}_bold.json")))["RepetitionTime"]

        smooth = np.asarray(nib.load(os.path.join(sesdir, f"srrun-{run}_bold.nii.gz")).dataobj)
        raw    = np.asarray(nib.load(os.path.join(sesdir, f"rrun-{run}_bold.nii.gz")).dataobj)
        atlas  = np.asarray(nib.load(os.path.join(adir, f"atlas115_run-{run}_native.nii.gz")).dataobj).astype(int)
        aseg   = np.asarray(nib.load(os.path.join(adir, f"aseg_run-{run}_native.nii.gz")).dataobj).astype(int)
        n_vols = smooth.shape[3]

        # one scale factor per run, applied to both copies so region and
        # noise signals sit on the same footing
        brain = aseg > 0
        _, factor = pc.grand_median_scale(raw, np.broadcast_to(brain[..., None], raw.shape))
        smooth = smooth * factor
        raw = raw * factor

        # one signal per atlas region, from the smoothed copy
        ts = np.zeros((n_vols, 115))
        thin = []
        for node in range(1, 116):
            m = atlas == node
            if m.sum() < min_node_vox:
                thin.append(node)
                ts[:, node - 1] = np.nan if m.sum() == 0 else smooth[m].T.mean(axis=1)
                continue
            block = smooth[m].T                      # time x voxels
            ts[:, node - 1] = pc.first_eigenvariate(block) if extract_method == "eig" \
                else block.mean(axis=1)
        if thin:
            print(f"WARNING {ses} run-{run}: nodes under {min_node_vox} voxels: {thin}")

        # noise signals come from the unsmoothed copy: smoothing smears real
        # brain signal into the white matter and CSF, and it would then get
        # regressed back out of the data
        wm = binary_erosion(np.isin(aseg, wm_labels), iterations=erode_iters)
        csf = binary_erosion(np.isin(aseg, csf_labels), iterations=erode_iters)
        if not csf.any():
            csf = np.isin(aseg, csf_labels)          # small ventricles can vanish when shaved; use the full mask
            print(f"NOTE {ses} run-{run}: csf erosion emptied the mask, using uneroded")
        rp = np.loadtxt(os.path.join(sesdir, f"motion_run-{run}.txt"))
        nuis = np.column_stack([pc.friston24_ar1(rp),
                                raw[wm].T.mean(axis=1),
                                raw[csf].T.mean(axis=1)])

        resid = pc.nuisance_regress(ts, nuis)
        clean = pc.modwt_band(resid, wavelet_scales)

        base = os.path.join(out, f"{subject}_{ses}_run-{run}")
        np.savetxt(base + "_ts115.tsv", clean, delimiter="\t", fmt="%.6f")
        json.dump({
            "tr": tr, "n_vols": n_vols, "scaling": "grand_median_1000",
            "scale_factor": float(factor), "extract": extract_method,
            "wavelet": "db8", "modwt_scales": wavelet_scales,
            "band_hz": [pc.scale_band_hz(max(wavelet_scales), tr)[0],
                        pc.scale_band_hz(min(wavelet_scales), tr)[1]],
            "solver_check": pc.solve_scales(tr, 0.01, 0.17),
            "nuisance": "friston24_ar1 + wm + csf, detrend projector both sides",
            "wm_vox": int(wm.sum()), "csf_vox": int(csf.sum()),
            "thin_nodes": thin, "labels": os.path.relpath(labels_tsv, out),
        }, open(base + "_ts115.json", "w"), indent=2)
        print(f"{ses} run-{run}: wrote {clean.shape} timeseries")

    print(f"stage 08 done: {ses}")


if __name__ == "__main__":
    main()
