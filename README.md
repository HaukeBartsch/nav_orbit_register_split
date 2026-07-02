# Load / register / split / export

Task: Load a set of 2 CT images with their corresponding masks. Split the masks after registration and export some statistics.

Using pandas, SimpleITK and SimpleITK-SimpleElastix.

```bash
conda activate nav_orbit_register_split
./analyze.py \
  -i1 data/NAV_ORBIT/images/1.3.6.1.4.1.45037.004171102544140049240402150654051107991535320.nii.gz \
  -i2 data/NAV_ORBIT/images/1.3.6.1.4.1.45037.417744054018056302213130162265012103225857091.nii.gz \
  -m1 data/NAV_ORBIT/labels/1.3.6.1.4.1.45037.004171102544140049240402150654051107991535320.nii.gz \
  -m2 data/NAV_ORBIT/labels/1.3.6.1.4.1.45037.417744054018056302213130162265012103225857091.nii.gz \
  --output /tmp/NAV_ORBIT_test/
cat /tmp/volumes.json
```

As output the following information is generated (/tmp/volumes.json).

```json
{
  "final_mutual_information": 0.671422,
  "full_masks_dice": 0.8309363178873083,
  "volume_change_ratio": {
    "Left eye / medial": 0.3948,
    "Right eye / medial": 0.9031,
    "Left eye / lateral": 0.7542,
    "Right eye / lateral": 0.8061
  },
  "moved": {
    "Left eye / medial": 7.89731547915303,
    "Right eye / medial": 4.151688998790885,
    "Left eye / lateral": 17.828592184085956,
    "Right eye / lateral": 17.480299698230382
  },
  "fixed": {
    "Left eye / medial": 8.385168201593249,
    "Right eye / medial": 7.072876303773375,
    "Left eye / lateral": 20.90119778497074,
    "Right eye / lateral": 20.97036979765703
  }
}
```

## Three-stage data processing pipeline

- **Intensity-based affine registration** — bone mask extracted from pre-op CT (≥200 HU), dilated 7 voxels, used to initialize and constrain Elastix affine transform (12 DOF) registering post-op CT to pre-op CT space
- **Mask transformation** — post-op mask re-sampled via nearest-neighbor interpolation using the computed transform
- **Eye splitting & volume analysis** — hierarchical split (medial/lateral, then anterior/posterior) into 4 labeled regions, volumes computed in cm³

## Generated Files

| File |	What it is |
|------|-------------|
|ct_fixed.nii.gz	| Pre-op CT (reference) |
|ct_moved_resampled.nii.gz	| Post-op CT warped to pre-op space |
|mask_moved_resampled.nii.gz	| Post-op mask registered to pre-op space |
|split_mask_fixed_resampled.nii.gz	| Pre-op mask split into 4 regions |
|split_mask_moved_resampled.nii.gz	| Post-op mask split into 4 regions |
|volumes.json	| Volume per region, per image with QC measures |
|&lt;NAV_ORBIT_test&gt;.csv | Volume per region, per image with QC measures (csv format) |

## How to verify visually / numerically

- A larger mutual information value indicates a better fit. A value close to 0 would mean that both images are still misaligned.
- The dice coefficient between the aligned masks (before splitting) should be closer to 1 for a good fit. A value of 1 indicates that both mask fit perfectly, which should not happen as both images are pre/post surgery. A value of 0 mean there is no overlap between the orbit masks after registration.
- The volume change ratios should be closer to 0 for the sides without surgery.
- Check that all four regions have non-zero volumes in both `fixed` and `moved` (unless the surgery intentionally removed tissue).
- Check that the total volume (sum of all four regions) is in a plausible range for orbital volumes (typically 25–45 cm³ total per eye in adults).
- **Failure indicators:**
  - All volumes are 0 — the mask was empty or the splitting logic failed.
  - One region has an implausibly large volume (e.g., > 50 cm³) — the split may have assigned too many voxels to one region.
  - The `fixed` and `moved` volumes are identical — the registration transform may be the identity (no transformation applied).
