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
if(!require(devtools)){install.packages("devtools")}
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
taxoBST_initial<-taxoBST
taxoTDF_initial<-taxoTDF
```

## Full diagnostic

``` r
suggestedBST<-fullTaxonomicDiagnostic(taxoBST)
```

    cleaning space characters

    misplaced qualifiers for undetermined taxa

    Warning in checkUndetermited(taxo = structure(list(plot =
    c("AltoSanJorgeInicial", : Some taxa have an identifier of type "sp. n" in more
    than one column see rows: 2470

    Warning in checkUndetermited(taxo = structure(list(plot =
    c("AltoSanJorgeInicial", : Some taxa have an identifier of type "sp." or "spp."
    in more than one column see rows:2461 2471 2475 2477

    unicity of genus in family

    checking for unicity of taxonomic information associated with taxonomic code

    Warning in checkUnicityCodetax(taxo = structure(list(plot = c("AltoSanJorgeInicial", : These codes corresponds to different information, but we were unable to find the most represented set of variables:
    Bignsp3 Cappsp2 Cnidsp Eugemont Ficuobtu indet1 Lundsp Machinun Malvsp Morfsp15 Morfsp17 Morfsp3 Morfsp6 Morfsp7 Morfsp8 Morpsp1 Morpsp4 Morpsp6 Morpsp7 Morpsp8 Morpsp9 Nectline Paulglob Pipehisp Stylturb Termamaz Verbsp
     Note that NO CORRECTION WILL BE SUGGESTED FOR THESE CASES

    Comparing species information with gbif backbone

    Searching for 697 taxa in the GBIF Backbone
    ...

    done

    Analysing GBIF Backbone information

    617 taxa are found without any modification needed

    13 taxa are found with suggested orthographic changes

    55 taxa are suggested synonyms

    16 taxa are found with suggested higher rank changes

    1 taxa were not found

    Comparing genus information with gbif backbone

    Searching for 184 taxa in the GBIF Backbone
    ...

    done

    Analysing GBIF Backbone information

    171 taxa are found without any modification needed

    3 taxa are found with suggested orthographic changes

    4 taxa are suggested synonyms

    2 taxa are found with suggested higher rank changes

    5 taxa were not found

    Comparing family information with gbif backbone

    Searching for 27 taxa in the GBIF Backbone
    ...

    done

    Analysing GBIF Backbone information

    27 taxa are found without any modification needed

    0 taxa are found with suggested orthographic changes

    0 taxa are suggested synonyms

    0 taxa are found with suggested higher rank changes

    0 taxa were not found

``` r
kable(head(suggestedBST$suggested))
```

|  | id_suggest | row | suggestDescription | code | plot | family | genus | specificEpithet | infraspecificEpithet | comments | suggest_family | suggest_genus | suggest_specificEpithet | suggest_infraspecificEpithet | suggest_identificationQualifier | suggest_verbatimTaxonRank | gbifid |
|:---|---:|---:|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| 8 | 1 | 8 | Comparing species information with gbif backbone (exactMatch_changeHigherRanks) | Cochviti | AltoSanJorgeInicial | Bixaceae | Cochlospermum | vitifolium | NA | Cochlospermum vitifolium | Cochlospermaceae | Cochlospermum | vitifolium | NA | NA | NA | 2874865 |
| 9 | 2 | 9 | Comparing species information with gbif backbone (exactMatch_changeHigherRanks) | Cordbico | AltoSanJorgeInicial | Boraginaceae | Cordia | bicolor | NA | Cordia bicolor | Cordiaceae | Cordia | bicolor | NA | NA | NA | 5660135 |
| 10 | 3 | 10 | misplaced qualifiers for undetermined taxa -\> Comparing genus information with gbif backbone (exactMatch_changeHigherRanks) | Cordsp2 | AltoSanJorgeInicial | Boraginaceae | Cordia | sp2 | NA | Cordia sp.2 | Cordiaceae | Cordia | NA | NA | NA | sp2 | 2900865 |
| 11 | 4 | 11 | misplaced qualifiers for undetermined taxa | Cupasp2 | AltoSanJorgeInicial | Sapindaceae | Cupania | sp2 | NA | Cupania sp.2 | Sapindaceae | Cupania | NA | NA | NA | sp2 | NA |
| 13 | 5 | 13 | misplaced qualifiers for undetermined taxa | Ficusp | AltoSanJorgeInicial | Moraceae | Ficus | sp | NA | Ficus sp. | Moraceae | Ficus | NA | NA | NA | sp | NA |
| 19 | 6 | 19 | misplaced qualifiers for undetermined taxa | Machsp5 | AltoSanJorgeInicial | Fabaceae | Machaerium | sp5 | NA | Machaerium sp.5 | Fabaceae | Machaerium | NA | NA | NA | sp5 | NA |

``` r
suggestedTDF<-fullTaxonomicDiagnostic(taxoTDF)
```

    cleaning space characters

    misplaced qualifiers for undetermined taxa

    unicity of genus in family

    checking for unicity of taxonomic information associated with taxonomic code

    Comparing species information with gbif backbone

    Searching for 292 taxa in the GBIF Backbone
    ...

    done

    Analysing GBIF Backbone information

    259 taxa are found without any modification needed

    3 taxa are found with suggested orthographic changes

    23 taxa are suggested synonyms

    10 taxa are found with suggested higher rank changes

    1 taxa were not found

    Comparing genus information with gbif backbone

    Searching for 98 taxa in the GBIF Backbone
    ...

    done

    Analysing GBIF Backbone information

    94 taxa are found without any modification needed

    1 taxa are found with suggested orthographic changes

    2 taxa are suggested synonyms

    2 taxa are found with suggested higher rank changes

    0 taxa were not found

    Comparing family information with gbif backbone

    Searching for 3 taxa in the GBIF Backbone
    ...

    done

    Analysing GBIF Backbone information

    3 taxa are found without any modification needed

    0 taxa are found with suggested orthographic changes

    0 taxa are suggested synonyms

    0 taxa are found with suggested higher rank changes

    0 taxa were not found

``` r
kable(head(suggestedTDF$suggested))
```

|  | id_suggest | row | suggestDescription | code | plot | family | genus | sp_epithet | suggest_family | suggest_genus | suggest_sp_epithet | suggest_verbatimTaxonRank | gbifid |
|:---|---:|---:|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| 2 | 1 | 2 | misplaced qualifiers for undetermined taxa -\> Comparing genus information with gbif backbone (exactMatch_changeHigherRanks) | Ampesp1 | CardonalLoma | Ulmaceae | Ampelocera | sp1 | Cannabaceae | Ampelocera | NA | sp1 | 7303603 |
| 3 | 2 | 3 | misplaced qualifiers for undetermined taxa -\> Comparing genus information with gbif backbone (exactMatch_changeHigherRanks) | Ampesp1 | CardonalLoma | Ulmaceae | Ampelocera | sp1 | Cannabaceae | Ampelocera | NA | sp1 | 7303603 |
| 9 | 3 | 9 | Comparing species information with gbif backbone (exactMatch_synonym) | Bauhhyme | CardonalLoma | Fabaceae | Bauhinia | hymenaeifolia | Fabaceae | Schnella | hymenaeifolia | NA | 2953134 |
| 11 | 4 | 11 | misplaced qualifiers for undetermined taxa | Buncsp | CardonalLoma | Malpighiaceae | Bunchosia | sp | Malpighiaceae | Bunchosia | NA | sp | NA |
| 17 | 5 | 17 | misplaced qualifiers for undetermined taxa | Casesp1 | CardonalLoma | Salicaceae | Casearia | sp1 | Salicaceae | Casearia | NA | sp1 | NA |
| 23 | 6 | 23 | misplaced qualifiers for undetermined taxa | Coccsp1 | CardonalLoma | Polygonaceae | Coccoloba | sp1 | Polygonaceae | Coccoloba | NA | sp1 | NA |

### Exporting the full Diagnostic in Excel

``` r
fileDiagnosticBST<- "../suggestedBST.xlsx"
exportXlDiagnostic(suggestedBST,fileDiagnosticBST,overwrite = T)
```

    Writing sheets: suggested failed_step_4_unicityCodeTax failed_step_5_gbif_species failed_step_5_gbif_genus failed_step_5_gbif_family
    into file:/home/marius/Travail/traitementDonnees/2024_parcelas_permanentes/suggestedBST.xlsx

``` r
fileDiagnosticTDF<-"../suggestedTDF.xlsx"
exportXlDiagnostic(suggestedTDF,fileDiagnosticTDF,overwrite = T)
```

    Writing sheets: suggested failed_step_5_gbif_species failed_step_5_gbif_genus failed_step_5_gbif_family
    into file:/home/marius/Travail/traitementDonnees/2024_parcelas_permanentes/suggestedTDF.xlsx

Here you might want to:

- check on the excel file which correction you do not want to apply
  (delete the rows in the sheet “suggested”)
- load back the file with the function `importXlDiagnostic`
- use the function `correct` to apply the corrections to the R object
  directly

## Going into details

### Manage spaces in the taxonomy table

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

### Manage undetermined taxa

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

### Structural problems

#### Are genera always in the same family

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
##correctUndetermined(taxoBST,suggested_unicity_familyBST)
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

#### Are taxonomic codes always associated with the same information

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

| row | plot | code | family | genus | sp_epithet | suggest_family | suggest_genus | suggest_sp_epithet |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|

``` r
taxoTDF<-correct(taxoTDF, suggested_taxoCodeTDF)
table(taxoRanks(taxoTDF))
```


     higher  family   genus species 
          7      10     221     653 

### Using the gbif backbone to get suggested corrections

``` r
suggested_species_gbif_BST<-checkGbif(taxoBST,rankCheck="species")
```

    Searching for 689 taxa in the GBIF Backbone
    ...

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

## Other useful functions

### Information at a particular rank

The `getRank` function returns the name of taxa at a particular rank (NA
if the lower rank with information is higher than the one queried).

``` r
table(getRank(taxoBST, rank="family"))
```


         Acanthaceae      Achariaceae  Achatocarpaceae    Amaranthaceae 
                  22               23               13                5 
       Anacardiaceae       Annonaceae      Apocynaceae     Aptandraceae 
                 241              115              161                1 
          Araliaceae        Arecaceae Aristolochiaceae       Asteraceae 
                  41               54               13               30 
         Basellaceae     Bignoniaceae         Bixaceae     Bromeliaceae 
                   3              171                2                3 
         Burseraceae        Cactaceae      Cannabaceae      Capparaceae 
                 108               12               57              116 
          Caricaceae     Celastraceae Chrysobalanaceae       Clusiaceae 
                   7               43               24               42 
    Cochlospermaceae     Combretaceae      Connaraceae   Convolvulaceae 
                  17               33               10               10 
          Cordiaceae    Cyclanthaceae     Dilleniaceae        Ebenaceae 
                  98                1               14                5 
         Ehretiaceae   Elaeocarpaceae  Erythropalaceae  Erythroxylaceae 
                   1                1                1               50 
       Euphorbiaceae         Fabaceae    Hernandiaceae     Hypericaceae 
                 215              943               24                2 
     Lacistemataceae        Lamiaceae        Lauraceae    Lecythidaceae 
                   9               31              106               30 
         Loganiaceae       Lythraceae    Malpighiaceae        Malvaceae 
                   2                1              117              240 
      Marcgraviaceae  Melastomataceae        Meliaceae   Menispermaceae 
                   1               27              191                2 
         Monimiaceae         Moraceae    Muntingiaceae    Myristicaceae 
                   1              275                4                7 
           Myrtaceae    Nyctaginaceae        Ochnaceae         Oleaceae 
                 313              116                2                1 
      Passifloraceae         Peraceae   Phyllanthaceae   Phytolaccaceae 
                   4                3               27               17 
       Picramniaceae       Piperaceae          Poaceae     Polygalaceae 
                  11               56                2                5 
        Polygonaceae      Primulaceae       Rhamnaceae        Rubiaceae 
                 136               59               20              389 
            Rutaceae       Salicaceae      Sapindaceae       Sapotaceae 
                 226              274              171               77 
       Schoepfiaceae    Simaroubaceae     Siparunaceae      Smilacaceae 
                   1                2                4                2 
          Solanaceae    Stemonuraceae    Thymelaeaceae     Trigoniaceae 
                  29                2                3                1 
            Ulmaceae       Urticaceae      Verbenaceae      Viburnaceae 
                   3               72               87                1 
           Violaceae         Vitaceae     Vochysiaceae      Ximeniaceae 
                  19                9                1                5 
      Zygophyllaceae 
                   5 

### Lower rank for each row

It is often useful to know the lowest rank that might apply for each
register:

``` r
taxoRanks(taxoBST)[1:100]
```

      [1] species species species species species species species species species
     [10] genus   genus   species genus   species species species species species
     [19] genus   species species higher  higher  higher  higher  higher  higher 
     [28] species species genus   species species genus   species species species
     [37] species species species species species species species species species
     [46] species genus   species species species species species species species
     [55] species species higher  higher  higher  higher  genus   family  species
     [64] genus   genus   species species species species species species species
     [73] species species species genus   species genus   species species species
     [82] species species species species species species species species species
     [91] genus   species species genus   species species species species genus  
    [100] species
    Levels: higher < family < genus < species < infraspecies

Note that no analysis is done here, which means that if there is wrong
information in a column represented a rank, it will not be discarded.
See the following example where we count the ranks before and after
automatic corrections:

``` r
table(taxoRanks(taxoBST_initial))
```

    Warning in taxoRanks(taxoBST_initial): Some of the taxa have missing higher
    information (see rows 62 149 150 151 154 362 371 430 558 580 590 634 636 676
    683 684 706 771 799 1688 1745 1746 1755 1788 1816 1825 1893 2470 2638 2730 2744
    2766 2814 2854 2855 2888 2902 2928 2929 2952 3056 4206 4487 4526 4771 4777 4820
    4821 4822 4845 4866 5601 5602 5658 5674 6066


          higher       family        genus      species infraspecies 
             113            0            4         5952            5 

``` r
table(taxoRanks(taxoBST))
```


          higher       family        genus      species infraspecies 
             143          127          918         4886            0 

### Extracting an epithet from a complete name

The `extractEpithet` function allows us to extract an epithet:

``` r
extractEpithet("Espeletia grandiflora")
```

    [1] "grandiflora"

``` r
extractEpithet("Espeletia grandiflora var. boyacana")
```

    [1] "boyacana"

While the latest example works well , you should prefer the form without
qualifiers:

``` r
extractEpithet("Espeletia grandiflora boyacana")
```

    [1] "boyacana"

Used without any other argument, the function will basically extract the
last word of a string, however, if you give the higher ranks, it will
check whether they actually correspond:

``` r
extractEpithet("Espeletia grandiflora", higherTaxon = "Espeletia")
```

    [1] "grandiflora"

``` r
extractEpithet("Magnolia grandiflora", higherTaxon = "Magnolia")
```

    [1] "grandiflora"

``` r
extractEpithet("Espeletia grandiflora", higherTaxon = "Magnolia")
```

    Warning in extractEpithet("Espeletia grandiflora", higherTaxon = "Magnolia"): While extracting the epithet we found discrepancies between the taxon names and higher taxon names 
    Problematic cases:
    taxon: Espeletia grandiflora ; epithet: grandiflora ; higher taxon: Magnolia

    [1] "grandiflora"

Of course this function may be applied on more than one name:

``` r
speciesTableFromGbif<-suggested_species_gbif_BST$analysedGbif[1:30]
(spNames<-sapply(speciesTableFromGbif,function(x)x$canonicalname))
```

              Attalea butyracea          Bellucia pentamera 
            "Attalea butyracea"        "Bellucia pentamera" 
                  Bixa orellana         Casearia sylvestris 
                "Bixa orellana"       "Casearia sylvestris" 
               Cecropia peltata    Cochlospermum vitifolium 
             "Cecropia peltata"  "Cochlospermum vitifolium" 
                 Cordia bicolor        Dendropanax arboreus 
               "Cordia bicolor"      "Dendropanax arboreus" 
               Genipa americana             Guarea guidonia 
             "Genipa americana"           "Guarea guidonia" 
              Guazuma ulmifolia       Machaerium biovulatum 
            "Guazuma ulmifolia"     "Machaerium biovulatum" 
             Miconia spicellata               Myrcia fallax 
           "Miconia spicellata"          "Myrcia splendens" 
             Ochroma pyramidale      Phyllanthus acuminatus 
           "Ochroma pyramidale"    "Phyllanthus acuminatus" 
          Schefflera morototoni             Spondias mombin 
       "Didymopanax morototoni"           "Spondias mombin" 
                Trema micrantha           Xylopia aromatica 
             "Trema micranthum"         "Xylopia aromatica" 
          Acalypha diversifolia   Adenocalymma aspericarpum 
        "Acalypha diversifolia" "Adenocalymma aspericarpum" 
                Aegiphila laeta         Aristolochia maxima 
              "Aegiphila laeta"       "Aristolochia maxima" 
           Astronium graveolens       Bignonia diversifolia 
         "Astronium graveolens"     "Bignonia diversifolia" 
            Carludovica palmata            Clarisia biflora 
          "Carludovica palmata"          "Clarisia biflora" 
              Cupania americana              Ficus insipida 
            "Cupania americana"            "Ficus insipida" 

``` r
(gnNames<-sapply(speciesTableFromGbif,function(x)
  {hr<-x$higherRanks
  return(hr$canonicalname[hr$rank=="genus"])}))
```

            Attalea butyracea        Bellucia pentamera             Bixa orellana 
                    "Attalea"                "Bellucia"                    "Bixa" 
          Casearia sylvestris          Cecropia peltata  Cochlospermum vitifolium 
                   "Casearia"                "Cecropia"           "Cochlospermum" 
               Cordia bicolor      Dendropanax arboreus          Genipa americana 
                     "Cordia"             "Dendropanax"                  "Genipa" 
              Guarea guidonia         Guazuma ulmifolia     Machaerium biovulatum 
                     "Guarea"                 "Guazuma"              "Machaerium" 
           Miconia spicellata             Myrcia fallax        Ochroma pyramidale 
                    "Miconia"                  "Myrcia"                 "Ochroma" 
       Phyllanthus acuminatus     Schefflera morototoni           Spondias mombin 
                "Phyllanthus"             "Didymopanax"                "Spondias" 
              Trema micrantha         Xylopia aromatica     Acalypha diversifolia 
                      "Trema"                 "Xylopia"                "Acalypha" 
    Adenocalymma aspericarpum           Aegiphila laeta       Aristolochia maxima 
               "Adenocalymma"               "Aegiphila"            "Aristolochia" 
         Astronium graveolens     Bignonia diversifolia       Carludovica palmata 
                  "Astronium"                "Bignonia"             "Carludovica" 
             Clarisia biflora         Cupania americana            Ficus insipida 
                   "Clarisia"                 "Cupania"                   "Ficus" 

``` r
extractEpithet(spNames,gnNames)
```

            Attalea butyracea        Bellucia pentamera             Bixa orellana 
                  "butyracea"               "pentamera"                "orellana" 
          Casearia sylvestris          Cecropia peltata  Cochlospermum vitifolium 
                 "sylvestris"                 "peltata"              "vitifolium" 
               Cordia bicolor      Dendropanax arboreus          Genipa americana 
                    "bicolor"                "arboreus"               "americana" 
              Guarea guidonia         Guazuma ulmifolia     Machaerium biovulatum 
                   "guidonia"               "ulmifolia"              "biovulatum" 
           Miconia spicellata             Myrcia fallax        Ochroma pyramidale 
                 "spicellata"               "splendens"              "pyramidale" 
       Phyllanthus acuminatus     Schefflera morototoni           Spondias mombin 
                 "acuminatus"              "morototoni"                  "mombin" 
              Trema micrantha         Xylopia aromatica     Acalypha diversifolia 
                 "micranthum"               "aromatica"            "diversifolia" 
    Adenocalymma aspericarpum           Aegiphila laeta       Aristolochia maxima 
               "aspericarpum"                   "laeta"                  "maxima" 
         Astronium graveolens     Bignonia diversifolia       Carludovica palmata 
                 "graveolens"            "diversifolia"                 "palmata" 
             Clarisia biflora         Cupania americana            Ficus insipida 
                    "biflora"               "americana"                "insipida" 

### Getting clean taxonomy/taxonomical classification from gbif and adding higher ranks to the table

If the Gbif check has been done, through the `gbifCheck` function

Whenever you have done a search in Gbif, you may extract a complete list
of taxa and their higher taxa, from the GBIF backbone:

``` r
(tdfFamilies<-taxoTDF[taxoRanks(taxoTDF)=="family",])
```

            code        family genus leaf_phenology life_form collected.code  ind
    47   Morfsp1     Myrtaceae  <NA>      evergreen     shrub        RG-1881  376
    123 Morfsp49     Rubiaceae  <NA>      deciduous     liana           <NA> <NA>
    234 Morfsp17     Rubiaceae  <NA>      evergreen      tree           <NA> <NA>
    369 Morfsp11     Rubiaceae  <NA>      evergreen     liana           <NA> <NA>
    535 Morfsp47     Rubiaceae  <NA>      evergreen      tree           <NA> <NA>
    778  Morfsp2     Myrtaceae  <NA>      evergreen      tree        RG-2372   47
    779  Morfsp2     Myrtaceae  <NA>      evergreen      tree        RG-2372 <NA>
    780  Morfsp2     Myrtaceae  <NA>      evergreen      tree        RG-2388 <NA>
    781  Morfsp3     Myrtaceae  <NA>      evergreen      tree        RG-2458 <NA>
    857  Morfsp5 Amaranthaceae  <NA>      evergreen      forb           <NA>   T4
        sp_epithet          plot verbatimTaxonRank
    47        <NA>  CardonalLoma               sp1
    123       <NA> CardonalPlana              sp49
    234       <NA>     Colorados              sp17
    369       <NA>        Jabiru              sp11
    535       <NA>        Tambor              sp47
    778       <NA>       Tuparro               sp2
    779       <NA>       Tuparro               sp2
    780       <NA>       Tuparro               sp2
    781       <NA>       Tuparro               sp3
    857       <NA>       Vinculo               sp5

``` r
kable(tabClassifTdfFamilies<-extractCompleteTaxo(suggested_family_gbif_BST$analysedGbif))
```

|     | rank    | canonicalname  |  gbifid | parent_gbifid |
|:----|:--------|:---------------|--------:|--------------:|
| 1   | kingdom | Plantae        |       6 |            NA |
| 2   | phylum  | Tracheophyta   | 7707728 |             6 |
| 3   | class   | Magnoliopsida  |     220 |       7707728 |
| 4   | order   | Laurales       |     407 |           220 |
| 5   | family  | Lauraceae      |    6688 |           407 |
| 9   | order   | Malvales       |     941 |           220 |
| 10  | family  | Malvaceae      |    6685 |           941 |
| 14  | order   | Rosales        |     691 |           220 |
| 15  | family  | Ulmaceae       |    2382 |           691 |
| 19  | order   | Lamiales       |     408 |           220 |
| 20  | family  | Bignoniaceae   |    6655 |           408 |
| 24  | order   | Gentianales    |     412 |           220 |
| 25  | family  | Rubiaceae      |    8798 |           412 |
| 29  | order   | Fabales        |    1370 |           220 |
| 30  | family  | Fabaceae       |    5386 |          1370 |
| 34  | order   | Sapindales     |     933 |           220 |
| 35  | family  | Sapindaceae    |    6657 |           933 |
| 40  | family  | Apocynaceae    |    6701 |           412 |
| 44  | order   | Asterales      |     414 |           220 |
| 45  | family  | Asteraceae     |    3065 |           414 |
| 49  | order   | Ericales       |    1353 |           220 |
| 50  | family  | Sapotaceae     |    8802 |          1353 |
| 55  | family  | Urticaceae     |    6639 |           691 |
| 59  | order   | Myrtales       |     690 |           220 |
| 60  | family  | Myrtaceae      |    5014 |           690 |
| 64  | order   | Brassicales    | 7225535 |           220 |
| 65  | family  | Capparaceae    |    3111 |       7225535 |
| 70  | family  | Meliaceae      |    2397 |           933 |
| 75  | family  | Rhamnaceae     |    2407 |           691 |
| 79  | order   | Malpighiales   |    1414 |           220 |
| 80  | family  | Violaceae      |    6631 |          1414 |
| 84  | order   | Celastrales    |     950 |           220 |
| 85  | family  | Celastraceae   |    6715 |           950 |
| 90  | family  | Verbenaceae    |    6689 |           408 |
| 95  | family  | Primulaceae    |    6674 |          1353 |
| 100 | family  | Malpighiaceae  |    6676 |          1414 |
| 105 | family  | Moraceae       |    6640 |           691 |
| 109 | order   | Caryophyllales |     422 |           220 |
| 110 | family  | Amaranthaceae  |    3064 |           422 |

Note that since the `fullTaxonomicDiagnostic` function (when applied
with the checks=“gbif” option) has an analysedGbif component, function
`extractCompleteTaxo` may be applied on it too!

Then it might be useful to extract a taxonomical table on a column-based
format:

``` r
tabTaxoFromRank(tabClassifTdfFamilies, rank = "order")
```

    $names
            kingdom   phylum         class           order           
    407     "Plantae" "Tracheophyta" "Magnoliopsida" "Laurales"      
    941     "Plantae" "Tracheophyta" "Magnoliopsida" "Malvales"      
    691     "Plantae" "Tracheophyta" "Magnoliopsida" "Rosales"       
    408     "Plantae" "Tracheophyta" "Magnoliopsida" "Lamiales"      
    412     "Plantae" "Tracheophyta" "Magnoliopsida" "Gentianales"   
    1370    "Plantae" "Tracheophyta" "Magnoliopsida" "Fabales"       
    933     "Plantae" "Tracheophyta" "Magnoliopsida" "Sapindales"    
    414     "Plantae" "Tracheophyta" "Magnoliopsida" "Asterales"     
    1353    "Plantae" "Tracheophyta" "Magnoliopsida" "Ericales"      
    690     "Plantae" "Tracheophyta" "Magnoliopsida" "Myrtales"      
    7225535 "Plantae" "Tracheophyta" "Magnoliopsida" "Brassicales"   
    1414    "Plantae" "Tracheophyta" "Magnoliopsida" "Malpighiales"  
    950     "Plantae" "Tracheophyta" "Magnoliopsida" "Celastrales"   
    422     "Plantae" "Tracheophyta" "Magnoliopsida" "Caryophyllales"

    $gbifid
            kingdom  phylum class   order
    407           6 7707728   220     407
    941           6 7707728   220     941
    691           6 7707728   220     691
    408           6 7707728   220     408
    412           6 7707728   220     412
    1370          6 7707728   220    1370
    933           6 7707728   220     933
    414           6 7707728   220     414
    1353          6 7707728   220    1353
    690           6 7707728   220     690
    7225535       6 7707728   220 7225535
    1414          6 7707728   220    1414
    950           6 7707728   220     950
    422           6 7707728   220     422

And finally we can merge with the register table, to obtain:

``` r
res<-addHigherRanks(tdfFamilies,suggested_family_gbif_BST$analysedGbif,ranks = c("phylum","class","order"),mergeOn="family")
kable(res)
```

|  | code | family | genus | leaf_phenology | life_form | collected.code | ind | sp_epithet | plot | verbatimTaxonRank | phylum | class | order |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| 47 | Morfsp1 | Myrtaceae | NA | evergreen | shrub | RG-1881 | 376 | NA | CardonalLoma | sp1 | Tracheophyta | Magnoliopsida | Myrtales |
| 123 | Morfsp49 | Rubiaceae | NA | deciduous | liana | NA | NA | NA | CardonalPlana | sp49 | Tracheophyta | Magnoliopsida | Gentianales |
| 234 | Morfsp17 | Rubiaceae | NA | evergreen | tree | NA | NA | NA | Colorados | sp17 | Tracheophyta | Magnoliopsida | Gentianales |
| 369 | Morfsp11 | Rubiaceae | NA | evergreen | liana | NA | NA | NA | Jabiru | sp11 | Tracheophyta | Magnoliopsida | Gentianales |
| 535 | Morfsp47 | Rubiaceae | NA | evergreen | tree | NA | NA | NA | Tambor | sp47 | Tracheophyta | Magnoliopsida | Gentianales |
| 778 | Morfsp2 | Myrtaceae | NA | evergreen | tree | RG-2372 | 47 | NA | Tuparro | sp2 | Tracheophyta | Magnoliopsida | Myrtales |
| 779 | Morfsp2 | Myrtaceae | NA | evergreen | tree | RG-2372 | NA | NA | Tuparro | sp2 | Tracheophyta | Magnoliopsida | Myrtales |
| 780 | Morfsp2 | Myrtaceae | NA | evergreen | tree | RG-2388 | NA | NA | Tuparro | sp2 | Tracheophyta | Magnoliopsida | Myrtales |
| 781 | Morfsp3 | Myrtaceae | NA | evergreen | tree | RG-2458 | NA | NA | Tuparro | sp3 | Tracheophyta | Magnoliopsida | Myrtales |
| 857 | Morfsp5 | Amaranthaceae | NA | evergreen | forb | NA | T4 | NA | Vinculo | sp5 | Tracheophyta | Magnoliopsida | Caryophyllales |
