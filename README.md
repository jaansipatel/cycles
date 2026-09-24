# cycles

Preprocessing and functional connectivity for a dense-sampling menstrual cycle neuroimaging study.

The design closely follows the methodology of Pritschet, Santander et al. (2020, *NeuroImage*).

## Study design

- One subject, many sessions (`sub-01`, `ses-006` and onward)
- Two resting-state runs per session, 1000 volumes each
- CMRR multiband EPI: MB8, 48 slices, 2.5 mm isotropic, TR 450 ms, matrix 92 x 92
- A 5-volume reversed phase encode (PA) run acquired immediately before each
  rest run, for distortion correction
- Two MEMPRAGE T1 series per session

Slice timing, the wavelet band, and the parcellation size all follow from the TR, the run length, and the single-subject structure.

## Requirements

Built for the ORCD/OpenMind cluster with SLURM. Module versions are pinned in `config.sh`:

| Tool | Version |
|------|---------|
| FSL | 6.0.7 |
| FreeSurfer | 8.0.0 |
| ANTs | 2.1.0-3.8bed08 |
| MATLAB | 2023b (with SPM12) |

Python needs numpy, scipy, nibabel, pywavelets, and matplotlib.

## Stages

Each stage is one file, run in order. Bash stages carry SLURM headers and run under `sbatch`.

| Stage | File | What it does |
|-------|------|--------------|
| 01 | `01_bids_organize.sh` | Sorts converted scans into BIDS, pairs each rest run with its PA scan, checks header fields |
| 02 | `02_topup.sh` | Distortion correction, estimated separately per run |
| 03 | `03_spm_preproc.sh` | SPM motion correction and 4 mm smoothing |
| 04 | `04_freesurfer.sh` | recon-all on the session T1 |
| 05 | `05_ants_template.sh` | Builds a subject-specific anatomical average and registers sessions and MNI to it |
| 06 | `06_func_anat_reg.sh` | bbregister, mean functional to session T1 |
| 06b | `06b_build_atlas.sh` | Builds the 115-region atlas, carries it and the FreeSurfer segmentation into native functional space |
| 07 | `07_qc.py` | Motion, tSNR, alignment score, atlas coverage, overlay images |
| 08 | `08_extract_timeseries.py` | Region signals, nuisance regression, wavelet band-limiting |
| 09 | `09_connectivity.py` | 115 x 115 correlation matrices, Fisher z |

Shared code is in `config.sh` (all paths and settings) and `pritschet_core.py` (the signal processing functions).

## Running it

```bash
echo cycles006 > sessions.txt          # one raw session name per line
sbatch 01_bids_organize.sh             # adjust --array as sessions accrue
```

Stage 05 is not an array job, since the anatomical average needs every session at once. This can be rerun when new sessions land, then rerun 06b and everything after.

For a first session, it is best to run the stages interactively on a compute node rather than through `sbatch`, so failures are visible as they happen.


## Methodological decisions

**No slice timing correction.** At a 450 ms TR the largest between-slice offset is a small fraction of the frequencies analyzed here. 

**Distortion correction per run, not per session.** Each rest run is corrected using the PA scan acquired immediately before it.

**FreeSurfer, not SPM, for the noise masks.** Labels are sharper at 2.5 mm, and they do not wobble session to session.

**115 regions.** 100 Schaefer cortical parcels (7-network ordering) plus 15 subcortical regions from Harvard-Oxford. `atlas115_labels.tsv` indicates choices on node count, ordering, and network assignment.


### Changes from reference papers

1. The frequency band is isolated with a wavelet transform, not a bandpass filter
2. Intensity scaling divides by the median, not the mean
3. Detrending is applied to the nuisance regressors as well as the data
4. Motion regressors take the form `[R, R(t-1), R^2, R(t-1)^2]` with no intercept


## References

Pritschet, L., Santander, T., et al. (2020). Functional reorganization of brain
networks across the human menstrual cycle. *NeuroImage*, 220, 117091.

Taylor, C. M., Pritschet, L., et al. (2020). Progesterone shapes medial temporal
lobe volume across the human menstrual cycle. *NeuroImage*, 220, 117125.

Reference implementation: `tsantander/PritschetSantander2020_NI_Hormones`
