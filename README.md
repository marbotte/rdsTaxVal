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
```

## Creating the oneTable objects

``` r
taxoBST <- new_taxo_oneTab(rdsBST$taxonomy)
taxoTDF <- new_taxo_oneTab(rdsTDF$taxonomy, "oneTable", epithet="sp_epithet")
```

## Manage spaces in the taxonomy table

``` r
suggested_spaces_BST<-checkSpace(taxoBST)
kable(suggested_spaces_BST)
```

| id_suggest | row | col | type | field | code | plot | problem | show_space | suggested |
|---:|---:|---:|:---|:---|:---|:---|:---|:---|:---|
| 1 | 1298 | 8 | end | specificEpithet | Chlomang | HondaInicial | mangense | mangense# | mangense |
| 2 | 1301 | 8 | end | specificEpithet | Cupalati | HondaInicial | latifolia | latifolia# | latifolia |
| 3 | 1302 | 8 | end | specificEpithet | Cupalati | HondaInicial | latifolia | latifolia# | latifolia |
| 4 | 1308 | 8 | end | specificEpithet | Rondpube | HondaInicial | pubescens | pubescens# | pubescens |
| 5 | 1404 | 8 | end | specificEpithet | Cupalati | HondaTemprano | latifolia | latifolia# | latifolia |
| 6 | 1406 | 8 | end | specificEpithet | Eugeunif | HondaTemprano | uniflora | uniflora# | uniflora |

``` r
taxoBST<-correctSpace(taxoBST,suggested_spaces_BST)
```

``` r
suggested_spaces_TDF<-checkSpace(taxoTDF)
kable(suggested_spaces_TDF)
```

| id_suggest | row | col | type | field | code | plot | problem | show_space | suggested |
|-----------:|----:|----:|:-----|:------|:-----|:-----|:--------|:-----------|:----------|

``` r
taxoTDF<-correctSpace(taxoTDF,suggested_spaces_TDF)
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
taxoBST<-correctUndetermined(taxoBST,suggested_undetermined_BST)
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
taxoTDF<-correctUndetermined(taxoTDF,suggested_undetermined_TDF)
```

## Structural problems

### Are genera always in the same family

``` r
suggested_unicity_familyBST<-unicityGnInFam(taxoBST)
kable(head(suggested_unicity_familyBST))
```

| id_suggest | row | plot | code | genus | family | suggest_family |
|---:|---:|:---|:---|:---|:---|:---|
| 1 | 195 | ArmeroIntermedio | Anacexce | Anacardium | Acardiaceae | Anacardiaceae |
| 2 | 198 | ArmeroIntermedio | Cocccoro | Coccoloba | Polygoceae | Polygonaceae |
| 3 | 201 | ArmeroIntermedio | Forsspic | Forsteronia | Apocyceae | Apocynaceae |
| 4 | 205 | ArmeroIntermedio | Oxanespi | Oxandra | Annoceae | Annonaceae |
| 5 | 206 | ArmeroIntermedio | Oxanespi | Oxandra | Annoceae | Annonaceae |
| 6 | 207 | ArmeroIntermedio | Oxanespi | Oxandra | Annoceae | Annonaceae |

``` r
#correctUndetermined(taxoBST,suggested_unicity_familyBST)
```

``` r
suggested_unicity_familyTDF<-unicityGnInFam(taxoTDF)
kable(head(suggested_unicity_familyTDF))
```

| id_suggest | row | plot | code | genus | family | suggest_family |
|-----------:|----:|:-----|:-----|:------|:-------|:---------------|
