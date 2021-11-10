# Preparation of the maize Hmp 3.2.1 VCF files for GWAS analysis

## TODO: Complete section
* Upload the filtered Hmp 3.2.1 Hapmap file to Cyverse
* Upload the imputed file that can be used with rMVP
* Finish documenting the steps for generating the file


## After downloading the individual VCF chromsome files, TASSEL5 was used for:
* Reading the VCF files
* Homozygous Genotype (converts heterozygous SNPs to NA) 
* Filtering SNPs based on minimum count of 27 (10% of lines)
* Filtering SNPs based on minimum frequency of 0.1
* Filtering SNPs based on maximum frequency of 0.9
* Saving the individual filtered chromosome files as Hapmap files

## The Hapmap file was generated using the following command line code:
```
# chr00_header.txt is the first line of any Hapmap file.
head -n 1 chr01_filtered_27_01_09_indels.hmp.txt > chr00_header.txt # Should work, not tested
tail -q -n +2 chr00_header.txt chr01_filtered_27_01_09_indels.hmp.txt chr02_filtered_27_01_09_indels.hmp.txt chr03_filtered_27_01_09_indels.hmp.txt chr04_filtered_27_01_09_indels.hmp.txt chr05_filtered_27_01_09_indels.hmp.txt chr06_filtered_27_01_09_indels.hmp.txt chr07_filtered_27_01_09_indels.hmp.txt chr08_filtered_27_01_09_indels.hmp.txt chr09_filtered_27_01_09_indels.hmp.txt chr10_filtered_27_01_09_indels.hmp.txt > combined_filtered_27_01_09_indels.hmp.txt 
```

## Download source
```
Maize 282 association panel genotypes (7x, AGPv4 coordinates):
iplant/home/shared/panzea/hapmap3/hmp321/unimputed/282_libs_2015/uplifted_APGv4
hmp321_agpv4_chr1.vcf
hmp321_agpv4_chr...vcf
hmp321_agpv4_chr10.vcf
```

## Note from Panzea about these VCF files
```
"This directory contains VCF files with genotypes of the 277 maize lines from the "282" set, computed on HapMap 3.2.1 sites using the full depth obtained from sequencing of "Cornell_2015" libraries. The sequencing was done in 2 stages: 2x at Cornell, and then 4x at Novogene. Partial taxa name swap at Novogene initially led to incorrectly merged read depth and therefore incorrect genotypes.

The incorrect files c*_282_onHmp321.vcf.gz, previously stored in this directory have been replaced by their correct versions, c*_282_corrected_onHmp321.vcf.gz, on May 16 2017, after the name swap problem had been resolved.”

```

## Personal note
```
There are no c*_282_corrected_onHmp321.vcf.gz files in any of the hmp321 folders. On top of that there are no agpv3-aligned VCF files, only the agpv4 uplifted SNPs. The agpv4 uplifted SNPs were added/modified on June 2017 so I am assuming that they are based on the corrected hmp321 VCF files. 

Since the SNP name and coordinates are different, I assume the SNP name is based on agpv3 coordinates and that the actual coordinates are the uplifted agpv4.
```

## MD5 List
```
1bcd51a726a6aa5f64e7c703b42cb017  hmp321_agpv4_chr1.vcf.gz
e905af48900aaa70202f1709a0bfc87c  hmp321_agpv4_chr2.vcf.gz
415b282b181bee9bb7cb94455b167ef2  hmp321_agpv4_chr3.vcf.gz
80f84079d8c63068ab006228d890a645  hmp321_agpv4_chr4.vcf.gz
b58557fc827832faad5d8ac7f7bda036  hmp321_agpv4_chr5.vcf.gz
f929a880154ca411908202ce57f57041  hmp321_agpv4_chr6.vcf.gz
aa3e845417d6089ce426599c6cdb212c  hmp321_agpv4_chr7.vcf.gz
d12130af484672542b6307dec93708bb  hmp321_agpv4_chr8.vcf.gz
80f1b396bc501b6e3720801e66fea0f4  hmp321_agpv4_chr9.vcf.gz
26263a369717154d9091041168f90fe4  hmp321_agpv4_chr10.vcf.gz
```

## Taxa List
```
282set_33-16
282set_38-11Goodman-Buckler
282set_4226
282set_4722
282set_A188
282set_A214NGoodman-Buckler
282set_A239
282set_A441-5
282set_A554
282set_A556
282set_A6
282set_A619
282set_A632
282set_A634
282set_A635
282set_A641
282set_A654
282set_A659
282set_A661
282set_A679
282set_A680
282set_A682
282set_Ab28A
282set_B10
282set_B103
282set_B104
282set_B105
282set_B109
282set_B14A
282set_B164
282set_B2
282set_B37
282set_B46
282set_B52
282set_B57
282set_B64
282set_B68
282set_B73
282set_B73Htrhm
282set_B75
282set_B76
282set_B77
282set_B79
282set_B84
282set_B97
282set_C103
282set_C123
282set_C49A
282set_CH701-30
282set_CH9
282set_CI187-2
282set_CI21E
282set_CI28AGoodman-Buckler
282set_CI31A
282set_CI3A
282set_CI64
282set_CI66
282set_CM105
282set_CM174
282set_CM37
282set_CM7
282set_CML10
282set_CML103
282set_CML108
282set_CML11
282set_CML14
282set_CML154Q
282set_CML157Q
282set_CML158Q
282set_CML218
282set_CML220
282set_CML228
282set_CML238
282set_CML247
282set_CML254
282set_CML258
282set_CML261
282set_CML264
282set_CML277
282set_CML281
282set_CML287
282set_CML311
282set_CML314
282set_CML321
282set_CML322
282set_CML323
282set_CML328
282set_CML331
282set_CML332
282set_CML333
282set_CML341
282set_CML38
282set_CML45
282set_CML5
282set_CML52
282set_CML61
282set_CML69
282set_CML77
282set_CML91
282set_CML92
282set_CMV3
282set_CO106
282set_CO125
282set_CO255
282set_Ci7Goodman-Buckler
282set_Ci90CGoodman-Buckler
282set_Ci91BGoodman-Buckler
282set_D940Y
282set_DE-2
282set_DE1
282set_DE811
282set_E2558W
282set_EP1
282set_F2834T
282set_F44
282set_F6
282set_F7
282set_GA209
282set_GT112
282set_H105W
282set_H49
282set_H84
282set_H91
282set_H95
282set_H99
282set_HP301
282set_Hi27Goodman-Buckler
282set_I137TN
282set_I205
282set_I29
282set_IA2132Goodman-Buckler
282set_IDS28
282set_IDS69
282set_IDS91
282set_ILLHy
282set_Ia5125
282set_Il101
282set_Il14H
282set_K148
282set_K4
282set_K55
282set_K64
282set_KI3
282set_KY226
282set_KY228
282set_Ki11
282set_Ki14
282set_Ki2021
282set_Ki21
282set_Ki43
282set_Ki44
282set_Ky21
282set_L317
282set_L578
282set_M14
282set_M162W
282set_M37W
282set_MEF156-55-2
282set_MO1W
282set_MS1334
282set_MS153
282set_Mo17
282set_Mo18W
282set_Mo24W
282set_Mo44
282set_Mo45
282set_Mo46
282set_Mo47
282set_MoG
282set_Mp339
282set_Ms71
282set_Mt42
282set_N192
282set_N28Ht
282set_N6
282set_N7A
282set_NC222
282set_NC230
282set_NC232
282set_NC236
282set_NC238
282set_NC250
282set_NC258
282set_NC260
282set_NC262
282set_NC264
282set_NC290A
282set_NC294
282set_NC296
282set_NC296A
282set_NC298
282set_NC300
282set_NC302
282set_NC304
282set_NC306
282set_NC310
282set_NC314
282set_NC318
282set_NC320
282set_NC324
282set_NC326
282set_NC328
282set_NC33
282set_NC336
282set_NC338
282set_NC340
282set_NC342
282set_NC344
282set_NC346
282set_NC348
282set_NC350
282set_NC352
282set_NC354
282set_NC356
282set_NC358
282set_NC360
282set_NC362
282set_NC364
282set_NC366
282set_NC368
282set_ND246
282set_OH7B
282set_Oh40B
282set_Oh43
282set_Oh43E
282set_Oh603
282set_Os420
282set_P39Goodman-Buckler
282set_Pa762
282set_Pa875
282set_Pa880
282set_Pa91
282set_R109B
282set_R168
282set_R177
282set_R229
282set_R4
282set_SA24
282set_SC213R
282set_SC357
282set_SC55
282set_SD40
282set_SD44
282set_Sg1533
282set_Sg18
282set_T232
282set_T234
282set_T8
282set_Tx303
282set_Tx601
282set_Tzi10
282set_Tzi11
282set_Tzi16
282set_Tzi18
282set_Tzi25
282set_Tzi8
282set_Tzi9
282set_U267Y
282set_VA102
282set_Va14
282set_Va17
282set_Va22
282set_Va26
282set_Va35
282set_Va59
282set_Va85
282set_Va99
282set_VaW6
282set_W117Ht
282set_W153R
282set_W182B
282set_W22
282set_W64A
282set_WD
282set_Wf9
282set_Yu796
282set_i1677a
```