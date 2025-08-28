# rdsTaxVal


``` r
require(knitr)
```

    Loading required package: knitr

## sourcing the function (development phase)

``` r
lapply(dir("R"),function(x)source(paste("R",x,sep="/")))
```

    [[1]]
    [[1]]$value
    function () 
    {
        print("Hello, world!")
    }

    [[1]]$visible
    [1] FALSE


    [[2]]
    [[2]]$value
    function (obj) 
    {
        NA
    }

    [[2]]$visible
    [1] FALSE

## Reading the files

``` r
rdsBST<-readRDS("../data_google/4.rdsProyectos/dataRedBSTCol.rds")
rdsTDF<-readRDS("../data_google/4.rdsProyectos/dataTDF.rds")
rdsTDF$taxonomy$family<-as.character(rdsTDF$taxonomy$family)
```

## Creating the oneTable objects

``` r
taxoBST <- new_taxo_oneTab(rdsBST$taxonomy)
taxoTDF <- new_taxo_oneTab(rdsTDF$taxonomy, "oneTable", taxonRanks_names = c(family="family", genus="genus",species="sp_epithet"), taxonRanks_epithetized = c("sp_epithet"))

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
```

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
```

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
```

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
```
