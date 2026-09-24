#!/usr/bin/env python3
#SBATCH -J connectivity
#SBATCH -t 00:20:00
#SBATCH --mem=4G
#SBATCH -o logs/09_conn_%a.out
#SBATCH --array=0-0

# Stage 09: connectivity per run, one 115 x 115 matrix of correlations
# between region signals, transformed so values compare fairly across runs.

import json, os, subprocess
import numpy as np


def cfg(*names):
    here = os.path.dirname(os.path.abspath(__file__))
    out = subprocess.check_output(
        ["bash", "-c", f"source {here}/config.sh >/dev/null 2>&1; " +
         " ".join(f'echo ${n};' for n in names)], text=True).splitlines()
    return out if len(names) > 1 else out[0]


def fisher_z(r):
    # clip so perfect correlations don't blow up to infinity
    return np.arctanh(np.clip(r, -0.999999, 0.999999))


def main():
    subject, sessions_file, ts_dir, conn_dir = cfg(
        "subject", "sessions_file", "ts_dir", "conn_dir")
    idx = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
    ses_raw = open(sessions_file).read().split()[idx]
    ses = "ses-" + ses_raw.removeprefix("cycles")
    out = os.path.join(conn_dir, subject, ses)
    os.makedirs(out, exist_ok=True)

    zs = []
    for run in (1, 2):
        ts = np.loadtxt(os.path.join(
            ts_dir, subject, ses, f"{subject}_{ses}_run-{run}_ts115.tsv"))
        bad = np.isnan(ts).any(axis=0)
        if bad.any():
            print(f"WARNING {ses} run-{run}: nan nodes {np.where(bad)[0] + 1}, edges set to nan")
        z = fisher_z(np.corrcoef(ts, rowvar=False))
        np.fill_diagonal(z, 0.0)
        z[bad, :] = np.nan; z[:, bad] = np.nan
        np.savetxt(os.path.join(out, f"{subject}_{ses}_run-{run}_fc115_z.tsv"),
                   z, delimiter="\t", fmt="%.6f")
        zs.append(z)

    np.savetxt(os.path.join(out, f"{subject}_{ses}_mean_fc115_z.tsv"),
               np.nanmean(zs, axis=0), delimiter="\t", fmt="%.6f")
    print(f"stage 09 done: {ses}")


if __name__ == "__main__":
    main()
