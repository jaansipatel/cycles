#!/bin/bash
#SBATCH -J bids_org
#SBATCH -t 00:30:00
#SBATCH --mem=4G
#SBATCH -o logs/01_bids_%a.out
#SBATCH --array=0-0            # raise this as sessions are added

# Stage 01: sort one session of converted scans into BIDS folders and check
# the header fields the rest of the pipeline relies on. Files are copied,
# never moved, so the raw folder stays untouched.

set -euo pipefail
source "$(dirname "$0")/config.sh"

ses_raw=$(session_from_array)                 # e.g. cycles006
ses=ses-${ses_raw#cycles}                     # e.g. ses-006
# the raw folder has moved before, so try the likely layouts one by one
src=""
for cand in "${raw_root}/${ses_raw}/mri/nifty" "${raw_root}/${ses_raw}" "${raw_root}"; do
    if compgen -G "${cand}/run_*_T1_MEMPRAGE*.nii*" > /dev/null; then src=$cand; break; fi
done
[[ -n $src ]] || { echo "no run_* files found under ${raw_root} for ${ses_raw}"; exit 1; }
echo "raw source: $src"
dst=${bids_root}/${subject}/${ses}
mkdir -p "$dst"/{anat,func,fmap}

# copy an image and its json under a new BIDS name
place() {
    local pattern=$1 out=$2 sub=$3
    local hits=( "$src"/${pattern}.nii* )
    [[ -e ${hits[0]} ]] || { echo "MISSING: $pattern in $src"; return 1; }
    # rescans show up as extra copies; keep the last one acquired
    local nii=${hits[-1]}
    local base=${nii%.nii*}
    cp "$nii" "$dst/$sub/${subject}_${ses}_${out}.nii${nii##*.nii}"
    cp "${base}.json" "$dst/$sub/${subject}_${ses}_${out}.json"
}

# two T1 copies per session; keep the later one
place "run_*_T1_MEMPRAGE_1.0mm_p2_RMS"                          T1w                              anat

# the two rest runs and their reference volumes
place "run_*_SMS_CMRR_2.5mm_S8pe4_slice48_TR450_RS_run1"        dir-AP_task-rest_run-1_bold      func
place "run_*_SMS_CMRR_2.5mm_S8pe4_slice48_TR450_RS_run1_SBRef"  dir-AP_task-rest_run-1_sbref     func
place "run_*_SMS_CMRR_2.5mm_S8pe4_slice48_TR450_RS_run2"        dir-AP_task-rest_run-2_bold      func
place "run_*_SMS_CMRR_2.5mm_S8pe4_slice48_TR450_RS_run2_SBRef"  dir-AP_task-rest_run-2_sbref     func

# PA runs: the scanner reuses names, so two different PA runs can both be
# called RS_run1. Names get ignored here; each rest run is paired with the
# PA scanned right before it, going by series number.
python3 - "$src" "$dst" "$subject" "$ses" <<'PY'
import json, sys, glob, os, shutil
src, dst, sub, ses = sys.argv[1:5]

def series(j): return json.load(open(j))["SeriesNumber"]

pa = sorted((f for f in glob.glob(f"{src}/run_*_PA_SMS_CMRR_2.5mm_*_RS_run*.json")
             if "SBRef" not in f), key=series)
assert pa, f"no PA rest runs in {src}"
for run in (1, 2):
    ap_sn = series(f"{dst}/func/{sub}_{ses}_dir-AP_task-rest_run-{run}_bold.json")
    before = [f for f in pa if series(f) < ap_sn]
    assert before, f"no PA series precedes AP series {ap_sn}"
    pick = before[-1][:-5]                # strip .json
    out = f"{dst}/fmap/{sub}_{ses}_dir-PA_acq-rest_run-{run}_epi"
    nii = glob.glob(pick + ".nii*")[0]
    shutil.copy(nii, out + ".nii" + nii.split(".nii")[1])
    j = json.load(open(pick + ".json"))
    j["IntendedFor"] = os.path.relpath(glob.glob(f"{dst}/func/*run-{run}_bold.nii*")[0], os.path.dirname(dst))
    json.dump(j, open(out + ".json", "w"), indent=2)
    print(f"run-{run}: AP series {ap_sn} <- PA series {series(pick + '.json')} ({os.path.basename(pick)})")
PY

# the face task and breath hold runs stay in raw for now; this analysis
# only needs rest

# make sure the header fields later stages depend on are what we expect
python3 - "$dst" "$subject" "$ses" "$expected_tr" <<'PY'
import json, sys, glob, os
dst, sub, ses, tr_expect = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])

bad = []
for f in glob.glob(f"{dst}/func/*_bold.json"):
    s = json.load(open(f))
    if abs(s["RepetitionTime"] - tr_expect) > 1e-6: bad.append((f, "TR", s["RepetitionTime"]))
    if s.get("PhaseEncodingDirection") != "j-":     bad.append((f, "PE", s.get("PhaseEncodingDirection")))
    if "TotalReadoutTime" not in s:                 bad.append((f, "no TotalReadoutTime", None))
for epi in glob.glob(f"{dst}/fmap/*_epi.json"):
    if json.load(open(epi)).get("PhaseEncodingDirection") != "j":
        bad.append((epi, "PE", "expected j"))
for b in bad: print("SIDECAR CHECK FAILED:", *b)
sys.exit(1 if bad else 0)
PY

# only needs writing once
desc=${bids_root}/dataset_description.json
[[ -f $desc ]] || cat > "$desc" <<EOF
{"Name": "Lewis Lab cycles dense-sampling", "BIDSVersion": "1.8.0"}
EOF

echo "stage 01 done: $ses"
