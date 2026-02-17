# Load / register / split / export

Load a set of 2 CT images with their corresponding masks. Split the masks after registration and export some statistics.

Using pandas, SimpleITK and SimpleITK-SimpleElastix.

```bash
./analyze.py \
  -i1 data/NAV_ORBIT/images/1.3.6.1.4.1.45037.004171102544140049240402150654051107991535320.nii.gz \
  -i2 data/NAV_ORBIT/images/1.3.6.1.4.1.45037.417744054018056302213130162265012103225857091.nii.gz \
  -m1 data/NAV_ORBIT/labels/1.3.6.1.4.1.45037.004171102544140049240402150654051107991535320.nii.gz \
  -m2 data/NAV_ORBIT/labels/1.3.6.1.4.1.45037.417744054018056302213130162265012103225857091.nii.gz \
  --output /tmp/
```

As output the following information is generated (/tmp/nav_orbit_volumes.json).

```json
{
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