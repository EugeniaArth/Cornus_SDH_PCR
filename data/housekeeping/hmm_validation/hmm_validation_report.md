# Comprehensive HMM Validation Report (corrected profiles)

Profiles: H3=PF00125, eIF5A=PF01287, TEF1/EF1a=PF00009, UBQ=PF00240

> NOTE: eIF5A previously used PF01331 (mRNA_cap_enzyme, WRONG). Correct profile PF01287 (eIF-5a) is used here.

## Summary

| Family | Profile | Total hits | Validated |
|--------|---------|-----------|-----------|
| H3 | PF00125 | 45 | 6 |
| eIF5A | PF01287 | 3 | 3 |
| TEF1 / EF1a | PF00009 | 46 | TEF1:0, EF1a:9 |
| Ubiquitin | PF00240 | 54 | polyUb:12 |

## eIF5A candidates (PF01287)

- `NODE_39166_length_1416_cov_4264.255304_g17868_i0` (score 101.2): Eukaryotic translation initiation factor 5A-2
- `NODE_50859_length_1062_cov_1528.604146_g3934_i1` (score 98.6): Eukaryotic translation initiation factor 5A-4
- `NODE_53930_length_980_cov_1766.262084_g13316_i2` (score 100.8): Eukaryotic translation initiation factor 5A

## TEF1 candidates

None found by HMM+NR classification (Cornus lacks a distinct TEF1).

## H3 candidates (canonical)

- `NODE_39897_length_1392_cov_4023.766195_g12629_i2`: Histone H3.3
- `NODE_53397_length_994_cov_1811.162963_g15731_i1`: Histone H3.3
- `NODE_60621_length_815_cov_632.157963_g30276_i0`: Histone H3.2
- `NODE_60817_length_810_cov_532.544021_g30423_i0`: Histone H3.2
- `NODE_63066_length_759_cov_701.432394_g32145_i0`: Histone H3.2
- `NODE_69948_length_625_cov_309.848958_g38046_i0`: Histone H3.2

## Polyubiquitin candidates

- `NODE_27696_length_1817_cov_1648.611991_g12754_i0` n_dom=5 cov=0.94 spacing=[76, 76, 76, 76]
- `NODE_28715_length_1778_cov_1424.970503_g13204_i0` n_dom=5 cov=0.94 spacing=[76, 76, 76, 76]
- `NODE_30540_length_1710_cov_16320.173390_g14002_i0` n_dom=5 cov=0.94 spacing=[76, 76, 76, 76]
- `NODE_31256_length_1685_cov_142.812347_g14326_i0` n_dom=5 cov=0.78 spacing=[74, 77, 73, 73]
- `NODE_31412_length_1680_cov_12.937462_g14395_i0` n_dom=4 cov=0.72 spacing=[74, 76, 74]
- `NODE_32737_length_1631_cov_360.464602_g14995_i0` n_dom=6 cov=0.87 spacing=[76, 76, 76, 76, 76]
- `NODE_34465_length_1571_cov_4544.060447_g11248_i2` n_dom=5 cov=0.84 spacing=[76, 76, 76, 76]
- `NODE_36068_length_1517_cov_697.857629_g16432_i0` n_dom=2 cov=0.93 spacing=[76]
- `NODE_38635_length_1433_cov_2316.104046_g17618_i0` n_dom=4 cov=0.94 spacing=[76, 76, 76]
- `NODE_40207_length_1381_cov_1448.836336_g14875_i2` n_dom=5 cov=0.8 spacing=[76, 76, 76, 76]
- `NODE_45641_length_1215_cov_633.473413_g16925_i1` n_dom=2 cov=0.92 spacing=[76]
- `NODE_63731_length_745_cov_1210.834770_g32674_i0` n_dom=2 cov=0.85 spacing=[76]
