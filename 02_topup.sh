#!/bin/bash
#SBATCH -J topup
#SBATCH -t 04:00:00
#SBATCH --mem=16G
#SBATCH -c 4
#SBATCH -o logs/02_topup_%a.out
#SBATCH --array=0-0

# Stage 02: Topup runs before motion correction so that step sees clean images.

set -euo pipefail
source "$(dirname "$0")/config.sh"
load_fsl

ses_raw=$(session_from_array); ses=ses-${ses_raw#cycles}
func=${bids_root}/${subject}/${ses}/func
fmap=${bids_root}/${subject}/${ses}/fmap
out=${topup_dir}/${subject}/${ses}
mkdir -p "$out"; cd "$out"

# each rest run gets its own correction, from the PA scanned right before
# it (stage 01 already sorted out which is which), so a distortion that
# drifts over the session can't push run 1's fix onto run 2
for run in 1 2; do
    ap_src=${func}/${subject}_${ses}_dir-AP_task-rest_run-${run}_bold.nii.gz
    pa_src=${fmap}/${subject}_${ses}_dir-PA_acq-rest_run-${run}_epi.nii.gz

    # 5 volumes per direction is enough and keeps topup quick
    fslroi "$ap_src" ap5_r${run} 0 5
    fslroi "$pa_src" pa5_r${run} 0 5
    fslmerge -t appa_r${run} ap5_r${run} pa5_r${run}

    # acqparams rows follow the merge order: AP volumes first, then PA
    trt=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['TotalReadoutTime'])" \
          "${ap_src%.nii.gz}.json")
    { for i in 1 2 3 4 5; do echo "0 -1 0 $trt"; done
      for i in 1 2 3 4 5; do echo "0  1 0 $trt"; done; } > acqparams_r${run}.txt

    # the standard config halves the image early on to save time, which needs
    # even dimensions; this data is 92 x 92 x 48, so that is fine
    topup --imain=appa_r${run} --datain=acqparams_r${run}.txt --config=b02b0.cnf \
          --out=topup_r${run} --fout=fieldmap_hz_r${run} --iout=appa_unwarped_r${run} -v

    for suf in bold sbref; do
        img=${func}/${subject}_${ses}_dir-AP_task-rest_run-${run}_${suf}.nii.gz
        applytopup --imain="$img" --datain=acqparams_r${run}.txt --inindex=1 \
                   --topup=topup_r${run} --method=jac \
                   --out=${out}/$(basename "${img%.nii.gz}")_topup
    done

    # rough plausibility check: a correction this strong would probably mean the
    # pairing or a sign in acqparams is wrong
    p98=$(fslstats fieldmap_hz_r${run} -a -P 98)
    echo "run-${run} fieldmap |Hz| p98 = $p98"
    python3 -c "import sys; sys.exit(0 if float('$p98') < 150 else 1)" \
        || echo "WARNING: run-${run} field magnitude suspicious, check pairing/sign"
done

# the two corrections come from the same head minutes apart, so they
# should nearly match; a big gap means one pairing is off
load_fsl 2>/dev/null || true
fslmaths fieldmap_hz_r1 -sub fieldmap_hz_r2 -abs fieldmap_diff
echo "field |r1 - r2| p98 = $(fslstats fieldmap_diff -P 98)"

echo "stage 02 done: $ses"
