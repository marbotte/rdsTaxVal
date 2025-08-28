# rdsTaxVal


<!-- badges: start -->

[![R-CMD-check](https://github.com/marbotte/rdsTaxVal/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/marbotte/rdsTaxVal/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

``` r
require(knitr)
```

    Loading required package: knitr

## Installing the package

You may install the `rdsTaxVal` package with the following command:

``` r
devtools::install_github("marbotte/rdsTaxVal")
```

## Loading the package

``` r
library(rdsTaxVal)
```

## Reading the files

``` r
rdsBST<-readRDS("../data_google/4.rdsProyectos/dataRedBSTCol.rds")
rdsTDF<-readRDS("../data_google/4.rdsProyectos/dataTDF.rds")
#rdsTDF$taxonomy$family<-as.character(rdsTDF$taxonomy$family)
```

## Creating the oneTable objects

``` r
taxoBST <- new_taxo_oneTab(rdsBST$taxonomy)
taxoTDF <- new_taxo_oneTab(rdsTDF$taxonomy, currentFormat = "oneTable", taxonRanks_names = c(family="family", genus="genus",species="sp_epithet"), taxonRanks_epithetized = c("sp_epithet"))

taxoBST[arrayInd(which(as.matrix(taxoBST) %in% c("","N/D","NA","N/A")), .dim=dim(taxoBST), useNames = T)] <- NA
taxoTDF[arrayInd(which(as.matrix(taxoTDF) %in% c("","N/D","NA","N/A")), .dim=dim(taxoTDF), useNames = T)] <- NA
```

## Manage spaces in the taxonomy table

``` r
suggested_spaces_BST<-checkSpace(taxoBST)
kable(suggested_spaces_BST)
```

|  | id_suggest | row | ref_plot | ref_code | specificEpithet | suggest_specificEpithet |
|:---|---:|---:|:---|:---|:---|:---|
| 1298 | 1 | 1298 | HondaInicial | Chlomang | mangense# | mangense |
| 1301 | 2 | 1301 | HondaInicial | Cupalati | latifolia# | latifolia |
| 1302 | 3 | 1302 | HondaInicial | Cupalati | latifolia# | latifolia |
| 1308 | 4 | 1308 | HondaInicial | Rondpube | pubescens# | pubescens |
| 1404 | 5 | 1404 | HondaTemprano | Cupalati | latifolia# | latifolia |
| 1406 | 6 | 1406 | HondaTemprano | Eugeunif | uniflora# | uniflora |

``` r
taxoBST<-correct(taxoBST,suggested_spaces_BST)
```

``` r
suggested_spaces_TDF<-checkSpace(taxoTDF)
kable(suggested_spaces_TDF)
```

| id_suggest | row | ref_plot | ref_code |
|-----------:|----:|:---------|:---------|

``` r
taxoTDF<-correct(taxoTDF,suggested_spaces_TDF)
table(taxoRanks(taxoTDF))
```


     higher  family   genus species 
          0       0       0     891 

## Manage undetermined taxa

``` r
suggested_undetermined_BST<-checkUndetermited(taxoBST)
```

    Warning in checkUndetermited(taxoBST): Some taxa have an identifier of type
    "sp. n" in more than one column see rows: 2470

    Warning in checkUndetermited(taxoBST): Some taxa have an identifier of type
    "sp." or "spp." in more than one column see rows:2461 2471 2475 2477

``` r
kable(head(suggested_undetermined_BST))
```

|  | id_suggest | row | plot | code | family | genus | specificEpithet | infraspecificEpithet | suggest_family | suggest_genus | suggest_specificEpithet | suggest_infraspecificEpithet | suggest_identificationQualifier | suggest_verbatimTaxonRank |
|:---|---:|---:|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| 10 | 1 | 10 | AltoSanJorgeInicial | Cordsp2 | Boraginaceae | Cordia | sp2 | NA | Boraginaceae | Cordia | NA | NA | NA | sp2 |
| 11 | 2 | 11 | AltoSanJorgeInicial | Cupasp2 | Sapindaceae | Cupania | sp2 | NA | Sapindaceae | Cupania | NA | NA | NA | sp2 |
| 13 | 3 | 13 | AltoSanJorgeInicial | Ficusp | Moraceae | Ficus | sp | NA | Moraceae | Ficus | NA | NA | NA | sp |
| 19 | 4 | 19 | AltoSanJorgeInicial | Machsp5 | Fabaceae | Machaerium | sp5 | NA | Fabaceae | Machaerium | NA | NA | NA | sp5 |
| 22 | 5 | 22 | AltoSanJorgeInicial | Morfsp14 | Asteraceae | Morpho | sp14 | NA | Asteraceae | NA | NA | NA | NA | sp14 |
| 23 | 6 | 23 | AltoSanJorgeInicial | Morfsp19 | Indet | Morpho | sp19 | NA | NA | NA | NA | NA | NA | sp19 |

``` r
taxoBST<-correct(taxoBST,suggested_undetermined_BST)
```

``` r
suggested_undetermined_TDF<-checkUndetermited(taxoTDF)
kable(head(suggested_undetermined_TDF))
```

|  | id_suggest | row | plot | code | family | genus | sp_epithet | suggest_family | suggest_genus | suggest_sp_epithet | suggest_verbatimTaxonRank |
|:---|---:|---:|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| 2 | 1 | 2 | CardonalLoma | Ampesp1 | Ulmaceae | Ampelocera | sp1 | Ulmaceae | Ampelocera | NA | sp1 |
| 3 | 2 | 3 | CardonalLoma | Ampesp1 | Ulmaceae | Ampelocera | sp1 | Ulmaceae | Ampelocera | NA | sp1 |
| 11 | 3 | 11 | CardonalLoma | Buncsp | Malpighiaceae | Bunchosia | sp | Malpighiaceae | Bunchosia | NA | sp |
| 17 | 4 | 17 | CardonalLoma | Casesp1 | Salicaceae | Casearia | sp1 | Salicaceae | Casearia | NA | sp1 |
| 23 | 5 | 23 | CardonalLoma | Coccsp1 | Polygonaceae | Coccoloba | sp1 | Polygonaceae | Coccoloba | NA | sp1 |
| 27 | 6 | 27 | CardonalLoma | Crotsp | Euphorbiaceae | Croton | sp | Euphorbiaceae | Croton | NA | sp |

``` r
taxoTDF<-correct(taxoTDF,suggested_undetermined_TDF)
table(taxoRanks(taxoTDF))
```


     higher  family   genus species 
          7      10     221     653 

## Structural problems

### Are genera always in the same family

``` r
suggested_unicity_familyBST<-checkUnicityRankSup(taxoBST, rank = "genus", superior="family")
kable(head(suggested_unicity_familyBST))
```

|  | id_suggest | row | plot | code | genus | family | suggest_family |
|:---|---:|---:|:---|:---|:---|:---|:---|
| 195 | 1 | 195 | ArmeroIntermedio | Anacexce | Anacardium | Acardiaceae | Anacardiaceae |
| 198 | 2 | 198 | ArmeroIntermedio | Cocccoro | Coccoloba | Polygoceae | Polygonaceae |
| 201 | 3 | 201 | ArmeroIntermedio | Forsspic | Forsteronia | Apocyceae | Apocynaceae |
| 205 | 4 | 205 | ArmeroIntermedio | Oxanespi | Oxandra | Annoceae | Annonaceae |
| 206 | 5 | 206 | ArmeroIntermedio | Oxanespi | Oxandra | Annoceae | Annonaceae |
| 207 | 6 | 207 | ArmeroIntermedio | Oxanespi | Oxandra | Annoceae | Annonaceae |

``` r
#correctUndetermined(taxoBST,suggested_unicity_familyBST)
taxoBST<-correct(taxoBST, suggested_unicity_familyBST)
```

``` r
suggested_unicity_familyTDF<-checkUnicityRankSup(taxoTDF, rank="genus",superior="family")
kable(head(suggested_unicity_familyTDF))
```

| id_suggest | row | plot | code | genus | family | suggest_family |
|-----------:|----:|:-----|:-----|:------|:-------|:---------------|

``` r
taxoTDF<- correct(taxoTDF,suggested_unicity_familyTDF)
table(taxoRanks(taxoTDF))
```


     higher  family   genus species 
          7      10     221     653 

### Are taxonomic codes always associated with the same information

``` r
suggested_taxoCodeBST <- checkUnicityCodetax(taxoBST,"takeFirst")
```

    Warning in checkUnicityCodetax(taxoBST, "takeFirst"): These codes corresponds to different information, but we were unable to find the most represented set of variables:
    Bignsp3 Cappsp2 Cnidsp Eugemont Ficuobtu indet1 Lundsp Machinun Malvsp Morfsp15 Morfsp17 Morfsp3 Morfsp6 Morfsp7 Morfsp8 Morpsp1 Morpsp4 Morpsp6 Morpsp7 Morpsp8 Morpsp9 Nectline Paulglob Pipehisp Stylturb Termamaz Verbsp
     Note the THE FIRST SET OF VARIABLES ASSOCIATED WITH THE CODE WILL BE ARBITRARILY USED AS SUGGESTED CORRECTION FOR THE OTHER ONES

``` r
kable(head(showUnicityCodetax(suggested_taxoCodeBST,"suggested")))
```

|  | id_suggest | row | plot | code | family | genus | specificEpithet | infraspecificEpithet | suggest_family | suggest_genus | suggest_specificEpithet | suggest_infraspecificEpithet |
|:---|---:|---:|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| 4597 | 1 | 4597 | SanFelipeIntermedio | Anacexce | Acardiaceae | Acardium | excelsum | NA | Anacardiaceae | Anacardium | excelsum | NA |
| 4598 | 2 | 4598 | SanFelipeIntermedio | Annosp | Annoceae | Anno | NA | NA | Annonaceae | Annona | NA | NA |
| 1297 | 3 | 1297 | HondaInicial | Attabuty | Arecaceae | Attalea | butyraceae | NA | Arecaceae | Attalea | butyracea | NA |
| 1315 | 4 | 1315 | HondaIntermedio | Attabuty | Arecaceae | Attalea | butyraceae | NA | Arecaceae | Attalea | butyracea | NA |
| 1361 | 5 | 1361 | HondaTardio | Attabuty | Arecaceae | Attalea | butyraceae | NA | Arecaceae | Attalea | butyracea | NA |
| 1400 | 6 | 1400 | HondaTemprano | Attabuty | Arecaceae | Attalea | butyraceae | NA | Arecaceae | Attalea | butyracea | NA |

``` r
taxoBST<-correct(taxoBST, suggested_taxoCodeBST)
```

``` r
suggested_taxoCodeTDF <- checkUnicityCodetax(taxoTDF,"takeFirst")
kable(head(showUnicityCodetax(suggested_taxoCodeTDF,"suggested")))
```

``` r
taxoTDF<-correct(taxoTDF, suggested_taxoCodeTDF)
table(taxoRanks(taxoTDF))
```


     higher  family   genus species 
          7      10     221     653 

## Using the gbif backbone to get suggested corrections

``` r
suggested_species_gbif_BST<-checkGbif(taxoBST,rankCheck="species")
```

    Searching for 689 taxa in the GBIF Backbone
    ...

    done

    done
    Analysing GBIF Backbone information

    610 taxa are found without any modification needed

    13 taxa are found with suggested orthographic changes

    54 taxa are suggested synonyms

    16 taxa are found with suggested higher rank changes

    1 taxa were not found

``` r
kable(head(suggested_species_gbif_BST$suggested))
```

|  | row | code | plot | family | genus | specificEpithet | type | suggest_family | suggest_genus | suggest_specificEpithet | gbifid |
|:---|---:|:---|:---|:---|:---|:---|:---|:---|:---|:---|---:|
| 8 | 8 | Cochviti | AltoSanJorgeInicial | Bixaceae | Cochlospermum | vitifolium | exactMatch_changeHigherRanks | Cochlospermaceae | Cochlospermum | vitifolium | 2874865 |
| 9 | 9 | Cordbico | AltoSanJorgeInicial | Boraginaceae | Cordia | bicolor | exactMatch_changeHigherRanks | Cordiaceae | Cordia | bicolor | 5660135 |
| 18 | 28 | Myrcfall | AltoSanJorgeInicial | Myrtaceae | Myrcia | fallax | exactMatch_synonym | Myrtaceae | Myrcia | splendens | 3174691 |
| 21 | 32 | Schemoro | AltoSanJorgeInicial | Araliaceae | Schefflera | morototoni | exactMatch_synonym | Araliaceae | Didymopanax | morototoni | 3038523 |
| 23 | 35 | Tremmicr | AltoSanJorgeInicial | Cannabaceae | Trema | micrantha | fuzzyMatch | Cannabaceae | Trema | micranthum | 2984501 |
| 24 | 36 | Tremmicr | AltoSanJorgeInicial | Cannabaceae | Trema | micrantha | fuzzyMatch | Cannabaceae | Trema | micranthum | 2984501 |

``` r
taxoBST<-correct(taxoBST,suggested_species_gbif_BST)
```

``` r
suggested_genus_gbif_BST<-checkGbif(taxoBST,rankCheck="genus")
```

    Searching for 178 taxa in the GBIF Backbone
    ...

    done

    done
    Analysing GBIF Backbone information

    166 taxa are found without any modification needed

    3 taxa are found with suggested orthographic changes

    4 taxa are suggested synonyms

    2 taxa are found with suggested higher rank changes

    4 taxa were not found

``` r
kable(head(suggested_genus_gbif_BST$suggested))
```

|  | row | code | plot | family | genus | type | suggest_family | suggest_genus | gbifid |
|:---|---:|:---|:---|:---|:---|:---|:---|:---|---:|
| 1 | 10 | Cordsp2 | AltoSanJorgeInicial | Boraginaceae | Cordia | exactMatch_changeHigherRanks | Cordiaceae | Cordia | 2900865 |
| 55 | 455 | Cordsp1 | BelloHorizonteInicial | Boraginaceae | Cordia | exactMatch_changeHigherRanks | Cordiaceae | Cordia | 2900865 |
| 56 | 456 | Cordsp2 | BelloHorizonteInicial | Boraginaceae | Cordia | exactMatch_changeHigherRanks | Cordiaceae | Cordia | 2900865 |
| 58 | 462 | Amphsp | BeltranP10 | Bignoniaceae | Amphillophium | fuzzyMatch | Bignoniaceae | Amphilophium | 3172588 |
| 108 | 820 | Ampesp | CambaoP18 | Ulmaceae | Ampelocera | exactMatch_changeHigherRanks | Cannabaceae | Ampelocera | 7303603 |
| 109 | 829 | Cordsp | CambaoP18 | Boraginaceae | Cordia | exactMatch_changeHigherRanks | Cordiaceae | Cordia | 2900865 |

``` r
taxoBST<-correct(taxoBST,suggested_genus_gbif_BST)
```

``` r
suggested_family_gbif_BST<-checkGbif(taxoBST,rankCheck="family")
```

    Searching for 22 taxa in the GBIF Backbone
    ...

    done

    done
    Analysing GBIF Backbone information

    22 taxa are found without any modification needed

    0 taxa are found with suggested orthographic changes

    0 taxa are suggested synonyms

    0 taxa are found with suggested higher rank changes

    0 taxa were not found

``` r
kable(head(suggested_family_gbif_BST$suggested))
```

| row | code | plot | family | type | suggest_family | gbifid |
|----:|:-----|:-----|:-------|:-----|:---------------|-------:|

``` r
taxoBST<-correct(taxoBST,suggested_family_gbif_BST)
```

``` r
suggested_species_gbif_TDF<-checkGbif(taxoTDF,rankCheck="species")
```

    Searching for 292 taxa in the GBIF Backbone
    ...

    done

    done
    Analysing GBIF Backbone information

    259 taxa are found without any modification needed

    3 taxa are found with suggested orthographic changes

    23 taxa are suggested synonyms

    10 taxa are found with suggested higher rank changes

    1 taxa were not found

``` r
kable(head(suggested_species_gbif_TDF$suggested))
```

|  | row | code | plot | family | genus | sp_epithet | type | suggest_family | suggest_genus | suggest_sp_epithet | gbifid |
|:---|---:|:---|:---|:---|:---|:---|:---|:---|:---|:---|---:|
| 7 | 9 | Bauhhyme | CardonalLoma | Fabaceae | Bauhinia | hymenaeifolia | exactMatch_synonym | Fabaceae | Schnella | hymenaeifolia | 2953134 |
| 20 | 25 | Cordgera | CardonalLoma | Boraginaceae | Cordia | gerascanthus | exactMatch_changeHigherRanks | Cordiaceae | Cordia | gerascanthus | 5341264 |
| 44 | 62 | Seguamer | CardonalLoma | Petiveriaceae | Seguieria | americana | exactMatch_changeHigherRanks | Phytolaccaceae | Seguieria | americana | 4192680 |
| 57 | 77 | Zizistry | CardonalLoma | Rhamnaceae | Ziziphus | strychnifolia | exactMatch_synonym | Rhamnaceae | Sarcomphalus | strychnifolius | 8759836 |
| 65 | 88 | Bauhhyme | CardonalPlana | Fabaceae | Bauhinia | hymenaeifolia | exactMatch_synonym | Fabaceae | Schnella | hymenaeifolia | 2953134 |
| 78 | 102 | Cordgera | CardonalPlana | Boraginaceae | Cordia | gerascanthus | exactMatch_changeHigherRanks | Cordiaceae | Cordia | gerascanthus | 5341264 |

``` r
taxoTDF<-correct(taxoTDF,suggested_species_gbif_TDF)
```

``` r
suggested_genus_gbif_TDF<-checkGbif(taxoTDF,rankCheck="genus")
```

    Searching for 98 taxa in the GBIF Backbone
    ...

    done

    done
    Analysing GBIF Backbone information

    94 taxa are found without any modification needed

    1 taxa are found with suggested orthographic changes

    2 taxa are suggested synonyms

    2 taxa are found with suggested higher rank changes

    0 taxa were not found

``` r
kable(head(suggested_genus_gbif_TDF$suggested))
```

|  | row | code | plot | family | genus | type | suggest_family | suggest_genus | gbifid |
|:---|---:|:---|:---|:---|:---|:---|:---|:---|---:|
| 1 | 2 | Ampesp1 | CardonalLoma | Ulmaceae | Ampelocera | exactMatch_changeHigherRanks | Cannabaceae | Ampelocera | 7303603 |
| 2 | 3 | Ampesp1 | CardonalLoma | Ulmaceae | Ampelocera | exactMatch_changeHigherRanks | Cannabaceae | Ampelocera | 7303603 |
| 21 | 83 | Ampesp1 | CardonalPlana | Ulmaceae | Ampelocera | exactMatch_changeHigherRanks | Cannabaceae | Ampelocera | 7303603 |
| 22 | 84 | Ampesp1 | CardonalPlana | Ulmaceae | Ampelocera | exactMatch_changeHigherRanks | Cannabaceae | Ampelocera | 7303603 |
| 24 | 103 | Cordsp | CardonalPlana | Boraginaceae | Cordia | exactMatch_changeHigherRanks | Cordiaceae | Cordia | 2900865 |
| 25 | 104 | Cordsp | CardonalPlana | Boraginaceae | Cordia | exactMatch_changeHigherRanks | Cordiaceae | Cordia | 2900865 |

``` r
taxoTDF<-correct(taxoTDF,suggested_genus_gbif_TDF)
```

``` r
suggested_family_gbif_TDF<-checkGbif(taxoTDF,rankCheck="family")
```

    Searching for 3 taxa in the GBIF Backbone
    ...

    done

    done
    Analysing GBIF Backbone information

    3 taxa are found without any modification needed

    0 taxa are found with suggested orthographic changes

    0 taxa are suggested synonyms

    0 taxa are found with suggested higher rank changes

    0 taxa were not found

``` r
kable(head(suggested_family_gbif_TDF$suggested))
```

| row | code | plot | family | type | suggest_family | gbifid |
|----:|:-----|:-----|:-------|:-----|:---------------|-------:|

``` r
taxoTDF<-correct(taxoTDF,suggested_family_gbif_TDF)
```
