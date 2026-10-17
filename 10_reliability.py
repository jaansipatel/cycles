#!/orcd/data/ldlewis/001/users/jpatel/software/cycles_venv/bin/python3
#SBATCH -J reliability
#SBATCH -t 00:15:00
#SBATCH --mem=4G
#SBATCH -o logs/10_rel.out

# Stage 10: calculate connectome-level reliability over stage 09's FC matrices.
# This script checks how well the two same-day runs agree and how similar each day's
# mean connectome is to every other day's.

import json, os, subprocess, sys
import numpy as np

def main():
    # config.sh stays the single source of truth; read it through bash
    here = os.environ.get("SLURM_SUBMIT_DIR", os.path.dirname(os.path.abspath(__file__)))
    out = subprocess.check_output(
        ["bash", "-c", f"source {here}/config.sh >/dev/null 2>&1; "
         "echo $subject; echo $conn_dir; echo $qc_dir; echo $sessions_file"],
        text=True).split()
    subject, conn_dir, qc_dir, sessions_file = out

    sessions = ["ses-" + s.removeprefix("cycles")
                for s in open(sessions_file).read().split()]

    iu = None            # upper-triangle index, set from the first matrix
    within = {}
    means = {}
    for ses in sessions:
        d = os.path.join(conn_dir, subject, ses)
        r1 = np.loadtxt(os.path.join(d, f"{subject}_{ses}_run-1_fc115_z.tsv"))
        r2 = np.loadtxt(os.path.join(d, f"{subject}_{ses}_run-2_fc115_z.tsv"))
        if iu is None:
            iu = np.triu_indices(r1.shape[0], k=1)
        within[ses] = np.corrcoef(r1[iu], r2[iu])[0, 1]
        means[ses] = np.loadtxt(os.path.join(d, f"{subject}_{ses}_mean_fc115_z.tsv"))[iu]

    with open(os.path.join(qc_dir, "reliability.tsv"), "w") as f:
        f.write("session\twithin_session_r\n")
        for ses in sessions:
            f.write(f"{ses}\t{within[ses]:.4f}\n")
            print(f"{ses}: within-session r = {within[ses]:.3f}")

    # between-session similarity of mean connectomes, full square matrix
    with open(os.path.join(qc_dir, "similarity_between_sessions.tsv"), "w") as f:
        f.write("session\t" + "\t".join(sessions) + "\n")
        for a in sessions:
            row = [f"{np.corrcoef(means[a], means[b])[0, 1]:.4f}" for b in sessions]
            f.write(a + "\t" + "\t".join(row) + "\n")
            print(a, " ".join(row))

    print(f"stage 10 done: {len(sessions)} sessions")

if __name__ == "__main__":
    main()
