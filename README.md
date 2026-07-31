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
  "final_mutual_information": 1.24027,
  "full_masks_dice": 0.9029907246512155,
  "affine_scale": {
    "x": 1,
    "y": 1,
    "z": 1
  },
  "volume_change_ratio": {
    "Left eye / medial": 1.718,
    "Right eye / medial": 1.7642,
    "Left eye / lateral": 1.0277,
    "Right eye / lateral": 1.0083
  },
  "moved": {
    "Left eye / medial": 6.094768261645225,
    "Right eye / medial": 5.3841309189096425,
    "Left eye / lateral": 21.377972060092763,
    "Right eye / lateral": 20.881426106104357
  },
  "fixed": {
    "Left eye / medial": 3.547581782437404,
    "Right eye / medial": 3.0518850604551395,
    "Left eye / lateral": 20.8024475295336,
    "Right eye / lateral": 20.708862162457287
  },
  "fixed_image_name": "NAV-ORBIT_0XX",
  "moving_image_name": "NAV-ORBIT_0XX",
  "FixedSeriesInstanceUID": "1.2.752.2771",
  "MovingSeriesInstanceUID": "1.2.752.3771",
  "FixedSeriesDescription": "2",
  "MovingSeriesDescription": "1",
  "FixedEvent": "EventName:1",
  "MovingEvent": "EventName:2",
  "FixedStudyDate": "20180000",
  "MovingStudyDate": "20250000",
  "DayDifference": 2352
}
```

## Three-stage data processing pipeline

- **Intensity-based rigid registration** — bone mask extracted from pre-op CT (≥200 HU), dilated 7 voxels, used to initialize and constrain Elastix rigid transform (6 DOF: 3 translation + 3 rotation) registering post-op CT to pre-op CT space, transformation was initialized using a rigid post-op mask to pre-op mask registration
- **Mask transformation** — post-op mask re-sampled via nearest-neighbor interpolation using the computed transform
- **Eye splitting & volume analysis** — hierarchical split (medial/lateral, then anterior/posterior) into 4 labeled regions, volumes computed in cm³

```mermaid
flowchart TD
    subgraph Inputs["Inputs"]
        CT1["CT fixed (pre-op)"]
        CT2["CT moving (post-op)"]
        M1["Mask fixed (pre-op)"]
        M2["Mask moving (post-op)"]
    end

    subgraph Stage1["Stage 1: Mask-based rigid registration"]
        G1["Geometry check: mask1 ↔ ct1"]
        G2["Geometry check: mask2 ↔ ct2"]
        R1["Rigid registration: mask2 → mask1"]
        T1["EulerTransform, 6 DOF\nGeometricalCenter init\n500 iterations, 2000 samples"]
        TX1["Initial transform\nmask_initial_transform.txt"]
        MR1["Registered mask\nmask2OnMask1.nii.gz"]
    end

    subgraph Stage2["Stage 2: Intensity-based affine CT registration"]
        BM1["Bone mask (fixed/pre-op)\nthreshold 200–1000 HU\nsmooth σ=0.5"]
        BM1D["Bone mask dilated ×4"]
        BM2["Bone mask (moving/post-op)\nthreshold 300–600 HU"]
        BM2D["Bone mask dilated ×4"]
        R2["Rigid registration: ct2 → ct1\nAdvancedMattesMutualInformation\n6 DOF, 7 resolutions\n2000 samples, ASGD"]
        CR2["Registered CT\nct_moved_resampled.nii.gz"]
    end

    subgraph Stage3["Stage 3: Mask transformation & volume analysis"]
        TM["Transform mask2 via Transformix\nnearest-neighbor interpolation"]
        MM["Registered mask\nmask_moved_resampled.nii.gz"]
        QC["Dice coefficient\nmask1 ↔ registered mask2"]
        BB["Bounding box of\nregistered mask\ncompute center_x/y/z, split1, split2"]

        SPLIT_M["Split registered mask (label 1–4)\n• z > center_z → label 2 (medial)\n• z < split1 → label 3 (lateral)\n• z > split2 → label 4 (lateral)"]

        SPLIT_F["Split fixed mask (same bounds)\n• z > center_z → label 2\n• z < split1 → label 3\n• z > split2 → label 4"]

        STATS["LabelShapeStatistics\nvolume per label (mm³ → cm³)"]

        RATIO["Volume change ratio\nmoved / fixed per region"]

        OUT1["split_mask_moved_resampled.nii.gz"]
        OUT2["split_mask_fixed_resampled.nii.gz"]
        VJ["volumes.json\nmutualInformation, dice,\nvolumes, ratios, metadata"]
        CSV["<folder>.csv"]
    end

    CT1 --> G1
    CT2 --> R2
    M1 --> G1
    M2 --> R1
    G1 --> G2
    G2 --> R1
    R1 --> T1
    T1 --> TX1
    R2 --> CR2
    BM1 --> BM1D
    BM2 --> BM2D
    TX1 --> R2
    BM1D -. fixed mask .-> R2
    BM2D -. moving mask .-> R2
    CR2 --> TM
    TM --> MM
    MM --> QC
    MM --> BB
    BB --> SPLIT_M
    M1 --> SPLIT_F
    SPLIT_M --> OUT1
    SPLIT_F --> OUT2
    SPLIT_M --> STATS
    SPLIT_F --> STATS
    STATS --> RATIO
    RATIO --> VJ
    OUT1 --> VJ
    OUT2 --> VJ
    QC --> VJ
    VJ --> CSV
```

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
- The dice coefficient between the aligned masks (before splitting) should be close to 1 for a good fit. A value of 1 indicates that both mask fit perfectly, which should not happen as both images are pre/post surgery. A value of 0 mean there is no overlap between the orbit masks after registration.
- The volume change ratios should be closer to 1 for the sides without surgery. Larger than one indicates that volume after surgery is larger by that factor compared to before surgery.
- Check that all four regions have non-zero volumes in both `fixed` and `moved` (unless the surgery intentionally removed tissue).
- Processing might have failed at an initial stage, check for a log file that indicates the reason for the failure, check the raw data using a tool like 3D Slicer or ITK Snap.
- Check that the total volume (sum of all four regions) is in a plausible range for orbital volumes (typically 25–45 cm³ total per eye in adults).
- **Failure indicators:**
  - All volumes are 0 — the mask was empty or the splitting logic failed.
  - One region has an implausibly large volume (e.g., > 50 cm³) — the split may have assigned too many voxels to one region.
  - The `fixed` and `moved` volumes are identical — the registration transform may be the identity (no transformation applied).
