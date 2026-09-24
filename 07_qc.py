#!/usr/bin/env python3
#SBATCH -J qc
#SBATCH -t 01:00:00
#SBATCH --mem=8G
#SBATCH -o logs/07_qc_%a.out
#SBATCH --array=0-0

# Stage 07: quality checks per session. Numbers first (head motion, signal
# to noise, alignment score, atlas coverage), then pictures for eyeballing
# the alignment. Everything lands in one summary table row per run so
# sessions can be compared across the cycle: a quality measure that itself
# tracks cycle phase would contaminate the analysis, so catch it early.

import json, os, subprocess, sys
import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def cfg(*names):
    # config.sh stays the single source of truth; read it through bash
    here = os.path.dirname(os.path.abspath(__file__))
    out = subprocess.check_output(
        ["bash", "-c", f"source {here}/config.sh >/dev/null 2>&1; " +
         " ".join(f'echo ${n};' for n in names)], text=True).splitlines()
    return out if len(names) > 1 else out[0]


def framewise_displacement(rp, radius=50.0):
    # motion per volume; rotations become mm of movement at the head surface
    d = np.diff(rp, axis=0)
    d[:, 3:] *= radius
    return np.abs(d).sum(axis=1)


def main():
    (subject, sessions_file, spm_dir, atlas_dir, reg_dir, qc_dir) = cfg(
        "subject", "sessions_file", "spm_dir", "atlas_dir", "reg_dir", "qc_dir")
    idx = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
    ses_raw = open(sessions_file).read().split()[idx]
    ses = "ses-" + ses_raw.removeprefix("cycles")
    out = os.path.join(qc_dir, subject, ses)
    os.makedirs(out, exist_ok=True)

    summary_path = os.path.join(qc_dir, "qc_summary.tsv")
    if not os.path.exists(summary_path):
        with open(summary_path, "w") as f:
            f.write("session\trun\tmean_fd\tmax_fd\tpct_fd_over_0.3\tmedian_tsnr\t"
                    "bbr_cost\tn_nodes_present\tmin_node_vox\n")

    for run in (1, 2):
        sesdir = os.path.join(spm_dir, subject, ses)
        rp = np.loadtxt(os.path.join(sesdir, f"motion_run-{run}.txt"))
        fd = framewise_displacement(rp)

        func = nib.load(os.path.join(sesdir, f"rrun-{run}_bold.nii.gz"))
        y = np.asarray(func.dataobj)
        mu, sd = y.mean(axis=3), y.std(axis=3)
        brain = mu > np.percentile(mu, 75)          # rough brain mask, fine for this number
        tsnr = np.median(mu[brain] / np.maximum(sd[brain], 1e-6))

        atlas = np.asarray(nib.load(os.path.join(
            atlas_dir, subject, ses, f"atlas115_run-{run}_native.nii.gz")).dataobj).astype(int)
        counts = np.bincount(atlas.ravel(), minlength=116)[1:116]
        cost = float(open(os.path.join(
            reg_dir, subject, ses, f"run-{run}_bbr.dat.mincost")).read().split()[0])

        with open(summary_path, "a") as f:
            f.write(f"{ses}\t{run}\t{fd.mean():.4f}\t{fd.max():.4f}\t"
                    f"{(fd > 0.3).mean()*100:.2f}\t{tsnr:.1f}\t{cost:.4f}\t"
                    f"{(counts > 0).sum()}\t{counts.min()}\n")

        # motion over time
        fig, ax = plt.subplots(figsize=(9, 2.2))
        ax.plot(fd, lw=0.6); ax.axhline(0.3, color="r", lw=0.6)
        ax.set(xlabel="volume", ylabel="FD (mm)", title=f"{ses} run-{run}")
        for s in ("top", "right"): ax.spines[s].set_visible(False)
        fig.tight_layout(); fig.savefig(os.path.join(out, f"fd_run-{run}.png"), dpi=120)
        plt.close(fig)

        # atlas outlines drawn over the mean image, three slices
        fig, axes = plt.subplots(1, 3, figsize=(9, 3.2))
        edges = (np.diff(atlas, axis=0, prepend=0) != 0) | \
                (np.diff(atlas, axis=1, prepend=0) != 0)
        for ax, z in zip(axes, np.percentile(np.arange(atlas.shape[2]), [30, 50, 70]).astype(int)):
            ax.imshow(mu[:, :, z].T, cmap="gray", origin="lower")
            ax.contour(edges[:, :, z].T, levels=[0.5], colors="r", linewidths=0.4)
            ax.axis("off")
        fig.suptitle(f"{ses} run-{run}: atlas on mean EPI")
        fig.tight_layout(); fig.savefig(os.path.join(out, f"atlas_overlay_run-{run}.png"), dpi=120)
        plt.close(fig)

        zero = np.where(counts == 0)[0] + 1
        if zero.size:
            print(f"WARNING {ses} run-{run}: nodes with zero voxels: {zero.tolist()}")

    print(f"stage 07 done: {ses}")


if __name__ == "__main__":
    main()
