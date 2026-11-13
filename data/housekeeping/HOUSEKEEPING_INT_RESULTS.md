# Combined BLAST / Keyword / HMM Validation Results

## Methodology
1. **BLASTn** — Agrostemma githago Table S7 sequences vs Cornus transcriptome (e-value 1e-5)
2. **Keyword search** — product-name patterns in GFF3 annotation
3. **HMM validation** — hmmscan with Pfam profiles (gathering thresholds):
   - PF00125 (Histone) -> H3 (filtered to canonical H3.2/H3.3 by NR product)
   - PF01287 (eIF-5a) -> eIF5A (corrected; replaces wrongly-used PF01331/mRNA_cap_enzyme)
   - PF00009 (GTP_EFTU) -> EF-1a / TEF1 / organellar EF-Tu (classified by NR product)
   - PF00240 (ubiquitin) -> ubiquitin fold; polyubiquitin by domain architecture

## Summary Table

| Family | BLAST hits | Keyword hits | Intersection | HMM validated | Notes |
|--------|-----------|-------------|-------------|---------------|-------|
| ACT (Actin) | 2 | 155 | 2 | 0 |  |
| BTUB (Beta-tubulin) | 3 | 32 | 2 | 0 |  |
| EF1a (Elongation factor 1-alpha) | 4 | 26 | 4 | 0 | EF1a: 0 |
| GAPDH (GAPDH cytosolic) | 2 | 7 | 2 | 0 |  |
| H3 (Histone H3) | 0 | 11 | 0 | 6 | canonical H3: 6 (45 raw PF00125 hits) |
| eIF5A (eIF5A) | 0 | 3 | 0 | 3 | eIF5A (PF01287): 3 |
| TEF1 (Translation elongation factor 1) | 0 | 12 | 0 | 0 | TEF1: 0 |
| UBQ (Ubiquitin / polyubiquitin) | 0 | 472 | 0 | 54 | polyUb: 12/54 |

## ACT — Actin

**BLAST hits (2):**
  - `NODE_19648_length_2181_cov_78.058630_g9247_i0` (id 83.05%, bitscore 1044.0)
  - `NODE_33831_length_1593_cov_726.832902_g10573_i1` (id 82.58%, bitscore 1014.0)
**Keyword hits (155):**
  - `NODE_101333_length_322_cov_3.260073_g68952_i0`
  - `NODE_105337_length_305_cov_1.542969_g72956_i0`
  - `NODE_10642_length_2810_cov_13.667150_g3306_i1`
  - `NODE_11020_length_2776_cov_226.300330_g2855_i2`
  - `NODE_11441_length_2739_cov_22.293680_g5595_i0`
  - ... plus 150 more
**Intersection:** 2
**HMM validated:** 0

## BTUB — Beta-tubulin

**BLAST hits (3):**
  - `NODE_23938_length_1973_cov_540.610187_g9818_i1` (id 81.13%, bitscore 1059.0)
  - `NODE_24156_length_1963_cov_264.541275_g9334_i1` (id 80.14%, bitscore 968.0)
  - `NODE_903_length_5512_cov_73.338093_g441_i0` (id 76.28%, bitscore 259.0)
**Keyword hits (32):**
  - `NODE_13743_length_2555_cov_85.273743_g6152_i1`
  - `NODE_18997_length_2216_cov_2235.678819_g2331_i1`
  - `NODE_1932_length_4591_cov_42.699251_g949_i0`
  - `NODE_21001_length_2114_cov_130.399516_g6197_i2`
  - `NODE_21309_length_2100_cov_55.309605_g9061_i1`
  - ... plus 27 more
**Intersection:** 2
**HMM validated:** 0

## EF1a — Elongation factor 1-alpha

**BLAST hits (4):**
  - `NODE_26996_length_1845_cov_2241.225501_g6372_i1` (id 87.59%, bitscore 457.0)
  - `NODE_26674_length_1858_cov_391.308458_g12320_i0` (id 86.58%, bitscore 435.0)
  - `NODE_27894_length_1810_cov_321.861442_g5609_i1` (id 86.79%, bitscore 431.0)
  - `NODE_5275_length_3536_cov_5178.433611_g2630_i0` (id 86.31%, bitscore 422.0)
**Keyword hits (26):**
  - `NODE_12049_length_2687_cov_48.741092_g5856_i0`
  - `NODE_15213_length_2448_cov_84.235098_g7282_i0`
  - `NODE_15669_length_2418_cov_11.002955_g7480_i0`
  - `NODE_16779_length_2346_cov_112.272529_g6682_i1`
  - `NODE_18779_length_2228_cov_406.409821_g8226_i1`
  - ... plus 21 more
**Intersection:** 4
**HMM validated:** 0

## GAPDH — GAPDH cytosolic

**BLAST hits (2):**
  - `NODE_10643_length_2809_cov_958.134058_g2667_i2` (id 83.40%, bitscore 909.0)
  - `NODE_19071_length_2212_cov_4499.645400_g3603_i3` (id 81.89%, bitscore 832.0)
**Keyword hits (7):**
  - `NODE_10643_length_2809_cov_958.134058_g2667_i2`
  - `NODE_19071_length_2212_cov_4499.645400_g3603_i3`
  - `NODE_22975_length_2019_cov_1623.357868_g846_i3`
  - `NODE_25416_length_1909_cov_574.066667_g10164_i1`
  - `NODE_25657_length_1899_cov_3954.542162_g4634_i5`
  - ... plus 2 more
**Intersection:** 2
**HMM validated:** 0

## H3 — Histone H3

**BLAST hits:** 0 (none at e-value 1e-5)
**Keyword hits (11):**
  - `NODE_29197_length_1761_cov_15.182827_g13396_i0`
  - `NODE_33730_length_1597_cov_2765.946382_g15430_i0`
  - `NODE_39897_length_1392_cov_4023.766195_g12629_i2`
  - `NODE_44095_length_1262_cov_118.439406_g20262_i0`
  - `NODE_50008_length_1086_cov_31.182257_g23380_i0`
  - ... plus 6 more
**HMM validated (6):**
  - `NODE_39897_length_1392_cov_4023.766195_g12629_i2` score=132.0 E=2e-42  _Histone H3.3_
  - `NODE_53397_length_994_cov_1811.162963_g15731_i1` score=130.2 E=6e-42  _Histone H3.3_
  - `NODE_60621_length_815_cov_632.157963_g30276_i0` score=131.6 E=2e-42  _Histone H3.2_
  - `NODE_60817_length_810_cov_532.544021_g30423_i0` score=132.3 E=1e-42  _Histone H3.2_
  - `NODE_63066_length_759_cov_701.432394_g32145_i0` score=131.4 E=3e-42  _Histone H3.2_
  - `NODE_69948_length_625_cov_309.848958_g38046_i0` score=132.3 E=1e-42  _Histone H3.2_

## eIF5A — eIF5A

**BLAST hits:** 0 (none at e-value 1e-5)
**Keyword hits (3):**
  - `NODE_39166_length_1416_cov_4264.255304_g17868_i0`
  - `NODE_50859_length_1062_cov_1528.604146_g3934_i1`
  - `NODE_53930_length_980_cov_1766.262084_g13316_i2`
**HMM validated (3):**
  - `NODE_39166_length_1416_cov_4264.255304_g17868_i0` score=101.2 E=6e-33  _Eukaryotic translation initiation factor 5A-2_
  - `NODE_50859_length_1062_cov_1528.604146_g3934_i1` score=98.6 E=4e-32  _Eukaryotic translation initiation factor 5A-4_
  - `NODE_53930_length_980_cov_1766.262084_g13316_i2` score=100.8 E=8e-33  _Eukaryotic translation initiation factor 5A_

## TEF1 — Translation elongation factor 1

**BLAST hits:** 0 (none at e-value 1e-5)
**Keyword hits (12):**
  - `NODE_12049_length_2687_cov_48.741092_g5856_i0`
  - `NODE_15669_length_2418_cov_11.002955_g7480_i0`
  - `NODE_20746_length_2127_cov_9.090953_g9714_i0`
  - `NODE_21746_length_2078_cov_137.641695_g10158_i0`
  - `NODE_25561_length_1903_cov_54.217907_g11828_i0`
  - ... plus 7 more
**HMM validated:** 0

## UBQ — Ubiquitin / polyubiquitin

**BLAST hits:** 0 (none at e-value 1e-5)
**Keyword hits (472):**
  - `NODE_10197_length_2853_cov_37.823823_g4113_i2`
  - `NODE_10237_length_2848_cov_41.178635_g5009_i0`
  - `NODE_103905_length_311_cov_1.450382_g71524_i0`
  - `NODE_10503_length_2823_cov_60.649243_g5132_i0`
  - `NODE_10526_length_2821_cov_32.762987_g1471_i1`
  - ... plus 467 more
**HMM validated (54):**
  - `NODE_11021_length_2776_cov_189.479648_g5383_i0` n_dom=2 cov=0.25 kind=singleUb/fusion  _Phosphatidylinositol 4-kinase gamma 3_
  - `NODE_11122_length_2768_cov_183.803972_g5432_i0` n_dom=1 cov=0.22 kind=singleUb/fusion  _BAG family molecular chaperone regulator 1_
  - `NODE_11690_length_2716_cov_152.215598_g5697_i0` n_dom=2 cov=0.25 kind=singleUb/fusion  _Phosphatidylinositol 4-kinase gamma 4_
  - `NODE_12462_length_2654_cov_155.001919_g6042_i0` n_dom=2 cov=0.25 kind=singleUb/fusion  _Phosphatidylinositol 4-kinase gamma 4_
  - `NODE_12893_length_2621_cov_230.645412_g6242_i0` n_dom=1 cov=0.13 kind=singleUb/fusion  _Ubiquitin domain-containing protein DSK2b_
  - `NODE_15893_length_2402_cov_341.591160_g571_i1` n_dom=1 cov=0.13 kind=singleUb/fusion  _Ubiquitin domain-containing protein DSK2b_
  - `NODE_18713_length_2231_cov_361.979835_g5205_i1` n_dom=1 cov=0.18 kind=singleUb/fusion  _Ubiquitin receptor RAD23c_
  - `NODE_20592_length_2134_cov_155.982254_g8014_i2` n_dom=1 cov=0.21 kind=singleUb/fusion  _BAG family molecular chaperone regulator 1_
  - `NODE_21389_length_2096_cov_25.707377_g10011_i0` n_dom=1 cov=0.25 kind=singleUb/fusion  _BAG family molecular chaperone regulator 1_
  - `NODE_21441_length_2093_cov_229.857632_g10030_i0` n_dom=1 cov=0.16 kind=singleUb/fusion  _XP_059660312.1 uncharacterized protein LOC132306772 [Cornus florida]_
  - `NODE_21623_length_2084_cov_479.932187_g6561_i2` n_dom=1 cov=0.14 kind=singleUb/fusion  _Ubiquitin carboxyl-terminal hydrolase 7_
  - `NODE_24946_length_1930_cov_48.686869_g9674_i1` n_dom=1 cov=0.19 kind=singleUb/fusion  _Ubiquitin receptor RAD23d_
  - `NODE_27696_length_1817_cov_1648.611991_g12754_i0` n_dom=5 cov=0.94 kind=polyubiquitin  _Polyubiquitin_
  - `NODE_28144_length_1801_cov_570.082763_g10574_i2` n_dom=1 cov=0.17 kind=singleUb/fusion  _Ubiquitin receptor RAD23c_
  - `NODE_28715_length_1778_cov_1424.970503_g13204_i0` n_dom=5 cov=0.94 kind=polyubiquitin  _Polyubiquitin_
  - `NODE_29762_length_1739_cov_63.603550_g13220_i1` n_dom=1 cov=0.19 kind=singleUb/fusion  _LRR repeats and ubiquitin-like domain-containing protein At2g30105_
  - `NODE_30268_length_1720_cov_270.873728_g10269_i4` n_dom=1 cov=0.18 kind=singleUb/fusion  _Protein DNA-DAMAGE INDUCIBLE 1_
  - `NODE_30540_length_1710_cov_16320.173390_g14002_i0` n_dom=5 cov=0.94 kind=polyubiquitin  _Polyubiquitin_
  - `NODE_31191_length_1687_cov_289.449328_g7660_i1` n_dom=1 cov=0.2 kind=singleUb/fusion  _Ubiquitin receptor RAD23b_
  - `NODE_31256_length_1685_cov_142.812347_g14326_i0` n_dom=5 cov=0.78 kind=polyubiquitin  _Polyubiquitin_
  - `NODE_31412_length_1680_cov_12.937462_g14395_i0` n_dom=4 cov=0.72 kind=polyubiquitin  _Polyubiquitin (Fragment)_
  - `NODE_32025_length_1657_cov_276.667289_g1516_i3` n_dom=1 cov=0.13 kind=singleUb/fusion  _Ubiquitin domain-containing protein DSK2b_
  - `NODE_32737_length_1631_cov_360.464602_g14995_i0` n_dom=6 cov=0.87 kind=polyubiquitin  _Polyubiquitin_
  - `NODE_33910_length_1590_cov_353.295263_g8863_i2` n_dom=1 cov=0.74 kind=singleUb/fusion  _Small ubiquitin-related modifier 1_
  - `NODE_34465_length_1571_cov_4544.060447_g11248_i2` n_dom=5 cov=0.84 kind=polyubiquitin  _Polyubiquitin_
  - `NODE_36068_length_1517_cov_697.857629_g16432_i0` n_dom=2 cov=0.93 kind=polyubiquitin  _Ubiquitin-NEDD8-like protein RUB1_
  - `NODE_38234_length_1446_cov_28.372226_g17446_i0` n_dom=1 cov=0.3 kind=singleUb/fusion  _BAG family molecular chaperone regulator 4_
  - `NODE_38635_length_1433_cov_2316.104046_g17618_i0` n_dom=4 cov=0.94 kind=polyubiquitin  _Polyubiquitin (Fragment)_
  - `NODE_39742_length_1398_cov_355.910304_g14213_i1` n_dom=1 cov=0.22 kind=singleUb/fusion  _Ubiquitin-like domain-containing CTD phosphatase_
  - `NODE_40207_length_1381_cov_1448.836336_g14875_i2` n_dom=5 cov=0.8 kind=polyubiquitin  _Polyubiquitin_
  - `NODE_43628_length_1275_cov_293.643556_g20030_i0` n_dom=1 cov=0.27 kind=singleUb/fusion  _BAG family molecular chaperone regulator 2_
  - `NODE_44600_length_1247_cov_636.883973_g20510_i0` n_dom=1 cov=0.73 kind=singleUb/fusion  _Small ubiquitin-related modifier 2_
  - `NODE_45641_length_1215_cov_633.473413_g16925_i1` n_dom=2 cov=0.92 kind=polyubiquitin  _Ubiquitin-NEDD8-like protein RUB1_
  - `NODE_49946_length_1087_cov_259.983622_g23207_i1` n_dom=1 cov=0.57 kind=singleUb/fusion  _Small ubiquitin-related modifier 2_
  - `NODE_50400_length_1075_cov_81.334308_g21388_i1` n_dom=1 cov=0.33 kind=singleUb/fusion  _XP_059643647.1 uncharacterized protein LOC132285473 [Cornus florida]_
  - `NODE_5443_length_3504_cov_200.533430_g398_i1` n_dom=1 cov=0.07 kind=singleUb/fusion  _Ubiquitin-like domain-containing protein CIP73_
  - `NODE_5464_length_3501_cov_127.212051_g2726_i0` n_dom=1 cov=0.07 kind=singleUb/fusion  _Protein STABILIZED1_
  - `NODE_548_length_6152_cov_48.633131_g97_i1` n_dom=1 cov=0.04 kind=singleUb/fusion  _Transcription initiation factor TFIID subunit 1_
  - `NODE_54988_length_953_cov_204.388274_g26356_i0` n_dom=1 cov=0.71 kind=singleUb/fusion  _XP_059635974.1 polyubiquitin 9 isoform X1 [Cornus florida]_
  - `NODE_55138_length_949_cov_782.843333_g23128_i1` n_dom=1 cov=0.51 kind=singleUb/fusion  _Ubiquitin-ribosomal protein eL40 fusion protein_
  - `NODE_55233_length_948_cov_1.812013_g26517_i0` n_dom=2 cov=0.57 kind=multiUb  _Polyubiquitin (Fragment)_
  - `NODE_56345_length_918_cov_2214.281933_g27255_i0` n_dom=1 cov=0.46 kind=singleUb/fusion  _Ubiquitin-ribosomal protein eS31 fusion protein_
  - `NODE_57837_length_881_cov_1402.335337_g28305_i0` n_dom=1 cov=0.56 kind=singleUb/fusion  _Ubiquitin-ribosomal protein eL40 fusion protein_
  - `NODE_58673_length_860_cov_85.849568_g28873_i0` n_dom=1 cov=0.44 kind=singleUb/fusion  _BAG family molecular chaperone regulator 2_
  - `NODE_59521_length_839_cov_1492.755696_g29464_i0` n_dom=1 cov=0.46 kind=singleUb/fusion  _Ubiquitin-ribosomal protein eS31 fusion protein_
  - `NODE_62534_length_771_cov_945.271468_g1749_i1` n_dom=1 cov=0.73 kind=singleUb/fusion  _Small ubiquitin-related modifier 2_
  - `NODE_62712_length_767_cov_3047.445682_g31861_i0` n_dom=1 cov=0.46 kind=singleUb/fusion  _Ubiquitin-ribosomal protein eS31 fusion protein_
  - `NODE_63731_length_745_cov_1210.834770_g32674_i0` n_dom=2 cov=0.85 kind=polyubiquitin  _Polyubiquitin (Fragment)_
  - `NODE_63832_length_743_cov_2833.184438_g30821_i1` n_dom=1 cov=0.56 kind=singleUb/fusion  _Ubiquitin-ribosomal protein eL40 fusion protein_
  - `NODE_7710_length_3150_cov_178.721703_g3815_i0` n_dom=1 cov=0.08 kind=singleUb/fusion  _E3 ubiquitin-protein ligase UPL5_
  - `NODE_7841_length_3131_cov_252.178456_g3878_i0` n_dom=1 cov=0.1 kind=singleUb/fusion  _Probable splicing factor 3A subunit 1_
  - `NODE_83233_length_446_cov_24.098237_g50867_i0` n_dom=1 cov=0.75 kind=singleUb/fusion  _XP_059664474.1 uncharacterized protein LOC132310308 [Cornus florida]_
  - `NODE_8745_length_3013_cov_195.032389_g3288_i2` n_dom=1 cov=0.1 kind=singleUb/fusion  _Ubiquitin-like domain-containing protein CIP73_
  - `NODE_9659_length_2906_cov_21.974799_g4722_i0` n_dom=2 cov=0.25 kind=singleUb/fusion  _Phosphatidylinositol 4-kinase gamma 4_

## Profile corrections

TEF1 is absent in Cornus sericea (0 of 46 PF00009 hits classify as TEF1).
