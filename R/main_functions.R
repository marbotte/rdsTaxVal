#! Put the taxonomy from the RDS in one table with attributes
#!
#! Technically this creates an object with a new class inherited from data.frame
#! The attributes of this class will ensure that we can always understand the
#! structure of the object, and come back to the initial format
#!
#! @param obj the base object with the taxonomy, may be a list with taxonomy
#! data.frame for each plot, or a single data.frame
#! @param currentFormat describe the format of obj: "listPlot" for a list of
#! plots and "oneTable" for a single
#! @param family name of the family column
#! @param genus name of the genus column
#! @param epithet name of the specific epithet column
#! @param infraspecies name of the infraspecific epithet column
#! @param taxoCode name of the column code corresponding to code of the taxon
#! @param plot name of the plot (name of the permanent plot) column
#! @param cf_aff name of the column containing the cf. or aff. expression
#! @param sp_specif name of the column containg information about the morphospecies


# obj=rdsBST$taxonomy
# currentFormat="listPlot"
# family="family"
# genus="genus"
# epithet="specificEpithet"
# infraspecies="infraspecificEpithet"
# taxoCode="code"
# plot="plot"
# cf_aff="identificationQualifier"
# sp_specif="verbatimTaxonRank"

new_taxo_oneTab <- function(obj,currentFormat=c("listPlot","oneTable"),family="family", genus="genus", epithet="specificEpithet", infraspecies="infraspecificEpithet", taxoCode="code", plot="plot", cf_aff="identificationQualifier", sp_specif="verbatimTaxonRank")
{
  currentFormat <- match.arg(currentFormat)
  if(currentFormat == "listPlot")
  {
    cn<-lapply(obj,colnames)
    if(any(match(cn,cn) != 1))  {stop("Object from the list do not have the same colnames")}
    if(plot%in%cn[[1]])
    {
      oneTab <- Reduce(rbind,obj)
    }else{
      oneTab<-data.frame(plot=rep(names(obj),sapply(obj,nrow)),
                         Reduce(rbind,obj))
      colnames(oneTab)[1]<-plot
    }
  }
  if(currentFormat == "oneTable")
  {
    oneTab <- obj
  }
  # Checking whether minimal information is in the table
  stopifnot(family %in% colnames(oneTab))
  stopifnot(genus %in% colnames(oneTab))
  stopifnot(epithet %in% colnames(oneTab))
  stopifnot(taxoCode %in% colnames(oneTab))
  stopifnot(plot %in% colnames(oneTab))
  # Give the object its attributes
  oneTab<-structure(oneTab,
            origin=currentFormat,
            family=family,
            genus=genus,
            epithet=epithet,
            taxoCode=taxoCode,
            plot=plot,
            infraspecies=ifelse(infraspecies %in% colnames(oneTab),infraspecies,NA),
            cf_aff=ifelse(cf_aff %in% colnames(oneTab),cf_aff,NA),
            sp_specif=ifelse(sp_specif %in% colnames(oneTab), sp_specif,NA))
  class(oneTab) <- c("taxo_oneTab","data.frame")
  return(oneTab)
}



checkSpace <- function(taxo,show_space="#")
{
  stopifnot(is(taxo,"taxo_oneTab"))
  mat<-as.matrix(taxo[na.omit(unlist(attributes(taxo)[c("family","genus","epithet","taxoCode","plot","infraspecies","cf_aff","sp_specif")]))])
  m_col<-match(colnames(mat),colnames(taxo))
  w_space_start <- which(`dim<-`(grepl("^[[:space:]]+",mat),dim(mat)),arr.ind=T)
  w_space_end <- which(`dim<-`(grepl("[[:space:]]+$",mat),dim(mat)),arr.ind=T)
  w_space_multi <- which(`dim<-`(grepl("[[:space:]]{2,}",mat),dim(mat)),arr.ind=T)
  numCases<-nrow(w_space_start)+nrow(w_space_end)+nrow(w_space_multi)
  tabReport<-data.frame(
    id_suggest=seq(1,numCases,length.out=numCases),
    Reduce(rbind,list(w_space_start,w_space_end,w_space_multi)),
    type=c(rep("start",nrow(w_space_start)),rep("end",nrow(w_space_end)),rep("multi",nrow(w_space_multi))),
    field=colnames(mat)[c(w_space_start[,"col"],w_space_end[,"col"],w_space_multi[,"col"])],
    taxoCode=mat[c(w_space_start[,"row"],w_space_end[,"row"],w_space_multi[,"row"]),attr(taxo,"taxoCode")],
    plot=mat[c(w_space_start[,"row"],w_space_end[,"row"],w_space_multi[,"row"]),attr(taxo,"plot")],
    problem=c(mat[w_space_start],mat[w_space_end],mat[w_space_multi]),
    show_space=gsub("[[:space:]]",show_space,c(mat[w_space_start],mat[w_space_end],mat[w_space_multi])),
    suggested=gsub("[[:space:]]{2,}"," ",gsub("[[:space:]]+$","",gsub("^[[:space:]]+","",c(mat[w_space_start],mat[w_space_end],mat[w_space_multi])))),
    row.names=NULL)
  tabReport$col<-m_col[tabReport$col]
  colnames(tabReport)[c(6,7)]<-c(attr(taxo,"taxoCode"),attr(taxo,"plot"))
  return(tabReport)
}

correctSpace<-function(taxo,suggested)
{
  if(nrow(suggested)>0){
    taxo[as.matrix(suggested[c("row","col")])]<-suggested$suggested
  }
  return(taxo)
}

checkUndetermited <- function(taxo, correct = F)
# Note: when there are cases of Indet or Morpho, species code should be checked for relevant information
{
  stopifnot(is(taxo,"taxo_oneTab"))
  mat<-as.matrix(taxo[na.omit(unlist(attributes(taxo)[c("family","genus","epithet","infraspecies")]))])
  #if(is.na(attr(taxo,"infraspecies"))) {mat<-cbind(mat, infraspecificEpithet = rep(NA,nrow(mat)))}
  w_Indet <- which(`dim<-`(grepl("^[Ii]ndet",mat),dim(mat)),arr.ind=T)
  w_empty <- which(`dim<-`(grepl("^$",mat),dim(mat)),arr.ind=T)
  w_cf <- which(`dim<-`(grepl("^cf\\. ?(.*)$",mat),dim(mat)),arr.ind=T)
  w_aff <- which(`dim<-`(grepl("^aff\\. ?(.*)$",mat),dim(mat)),arr.ind=T)
  w_sp <- which(`dim<-`(grepl("^spp?\\.?$",mat),dim(mat)),arr.ind=T)
  w_spNb <- which(`dim<-`(grepl("^sp?\\.?[0-9]{1,3}\\.?$",mat),dim(mat)),arr.ind=T)
  w_morfo <- which(`dim<-`(grepl("^[Mm]or?(f)|(ph)o?\\.?$",mat),dim(mat)),arr.ind=T)

  matPb <- matrix(data = NA, nrow = nrow(mat), ncol = ncol(mat))
  matPb[w_Indet] <- "indet"
  #matPb[w_empty] <- "empty"
  matPb[w_cf] <- "cf"
  matPb[w_aff] <- "aff"
  matPb[w_sp] <- "sp"
  matPb[w_spNb] <- "spNb"
  matPb[w_morfo] <- "morpho"

  mat2<-mat
  mat2[w_Indet] <- NA
  mat2[w_empty] <- NA
  mat2[w_cf] <- NA
  mat2[w_aff] <- NA
  mat2[w_sp] <- NA
  mat2[w_spNb] <- NA
  mat2[w_morfo] <- NA
  notNALevels <- apply(apply(mat2,1,function(x)!is.na(x)),2,which,simplify = F)
  NALevels <- apply(apply(mat2,1,is.na),2,which,simplify = F)
  maxInfoLevel <- sapply(notNALevels,function(x)ifelse(length(x)==0,0,max(x)))+1
  minNotInfoLevel <- sapply(NALevels,function(x,y)ifelse(length(x)==0,y,min(x)),y=ncol(mat)+1)
  if(any(minNotInfoLevel != maxInfoLevel)) {stop("Seems there is a taxon with low level information and no high-level:\nlook at cases",paste(which(minNotInfoLevel != maxInfoLevel),collapse=" "))}
  #Now we have all the information we need about potential "undetermination" error in the taxonomic table, we may try and find suggested solutions, case by case
  suggested_cf_aff<-character(length=nrow(taxo))
  suggested_sp_specif <- character(length=nrow(taxo))
  # Case 1: we have a spNb (ejemplo sp1)
  case1 <- unique(w_spNb[,"row"])
  colCase1 <- tapply(w_spNb[,"col"], w_spNb[,"row"],function(x)x,simplify=F)
  if(any(sapply(colCase1,length)!=1)){warning("Some taxa have an identifier of type \"sp. n\" in more than one column see rows: ",paste(names(colCase1)[which(sapply(colCase1,length)!=1)],collapse=" "))}
  colCase1<-sapply(colCase1,function(x)x[1])[match(as.integer(names(colCase1)),case1)]
  suggested_sp_specif[case1]<-gsub("^sp?\\.?([0-9]{1,3})\\.?$","sp\\1",mat[cbind(row=case1,col=colCase1)])
  # Case 2: we have a sp. or spp. case
  case2 <- unique(w_sp[,"row"])
  if(any(case2%in%case1)){stop("Some taxa have an identifier of type \"sp.\" or \"spp.\" and an identifier of type \"sp. n\"",paste(case2[case2%in%case1],collapse = " "))}
  colCase2<- tapply(w_sp[,"col"], w_sp[,"row"],function(x)x,simplify=F)
  if(any(sapply(colCase2,length)!=1)){warning("Some taxa have an identifier of type \"sp.\" or \"spp.\" in more than one column see rows:",paste(names(colCase2)[which(sapply(colCase2,length)!=1)],collapse=" "))}
  colCase2<-sapply(colCase2,function(x)x[1])[match(as.integer(names(colCase2)),case2)]
  suggested_sp_specif[case2]<-gsub("^(spp?)\\.?$","\\1",mat[cbind(row=case2,col=colCase2)])
  # Manejo de los cf. and aff. should be possible to be true in more than one level
  cf_affMat<-rep("",length(mat))
  dim(cf_affMat)<-dim(mat)
  cf_affMat[w_cf]<-gsub("^cf\\. ?(.*)$","cf. \\1",mat[w_cf])
  cf_affMat[w_aff] <- gsub("^aff\\. ?(.*)$","aff. \\1",mat[w_aff])
  suggested_cf_aff<-trimws(apply(cf_affMat,1,paste,sep=" ",collapse=" "))

  if(!is.na(attr(taxo,"cf_aff")))
  {
    raw_cf_aff<-taxo[attr(taxo,"cf_aff")]
    raw_cf_aff[grepl("^[[:space:]]*$",raw_cf_aff)]<-NA
    suggested_cf_aff[suggested_cf_aff==""]<-raw_cf_aff[suggested_cf_aff==""]
  }else{
    suggested_cf_aff[suggested_cf_aff==""]<-NA
  }

  if(!is.na(attr(taxo,"sp_specif")))
  {
    raw_sp_specif<-taxo[attr(taxo,"sp_specif")]
    raw_sp_specif[grepl("^[[:space:]]*$",raw_sp_specif)]<-NA
    suggested_sp_specif[suggested_sp_specif==""]<-raw_sp_specif[suggested_sp_specif==""]
  }else{
    suggested_sp_specif[suggested_sp_specif==""]<-NA
  }

  anyPb<-sort(unique(c(w_Indet[,"row"], w_empty[,"row"], w_cf[,"row"], w_aff[,"row"], w_sp[,"row"], w_spNb[,"row"], w_morfo[,"row"])))
  tabRaw<-taxo[na.omit(unlist(attributes(taxo)[c("plot","taxoCode","family","genus","epithet","infraspecies","cf_aff","sp_specif")]))]
  tabSuggest<-as.data.frame(mat2)
  if(any(!is.na(suggested_cf_aff)))
  {
    tabSuggest[[ifelse(is.na(attr(taxo,"cf_aff")),"identificationQualifier",attr(taxo,"cf_aff"))]] <- suggested_cf_aff
  }
  if(any(!is.na(suggested_sp_specif)))
  {
    tabSuggest[[ifelse(is.na(attr(taxo,"sp_specif")),"verbatimTaxonRank",attr(taxo,"sp_specif"))]] <- suggested_sp_specif
  }
  colnames(tabSuggest)<-paste0("suggest_",colnames(tabSuggest))
  return(data.frame(id_suggest=seq(1,length(anyPb),length.out=length(anyPb)),row=anyPb,tabRaw[anyPb,],tabSuggest[anyPb,]))
}

correctUndetermined <- function(taxo, suggested)
{
  stopifnot(is(taxo,"taxo_oneTab"))
  cn <- colnames(suggested)[grep("suggest_",colnames(suggested))]
  resCol<-gsub("suggest_","",cn)
  colToAdd<-resCol[!resCol%in%colnames(taxo)]
  if(length(colToAdd)>0){
    taxo[colToAdd]<-NA
    if(any(colToAdd=="identificationQualifier")&is.na(attr(taxo,"cf_aff")))
    {attr(taxo,"cf_aff")<-"identificationQualifier"}
    if(any(colToAdd=="verbatimTaxonRank"&is.na(attr(taxo,"sp_specif"))))
    {attr(taxo,"sp_specif")<-"verbatimTaxonRank"}
    }
  w_suggested<-as.matrix(cbind(row=as.numeric(row(suggested[,cn])),col=as.numeric(col(suggested[,cn]))))
  w_taxo <- cbind(row=suggested$row[w_suggested[,"row"]],col=match(resCol,colnames(taxo))[w_suggested[,"col"]])
  taxo[w_taxo]<-suggested[,cn][w_suggested]
  return(taxo)
}

unicityCodetax <- function(taxo)
{
  stopifnot(is(taxo,"taxo_oneTab"))
  mat<-as.matrix(taxo[na.omit(unlist(attributes(taxo)[c("family","genus","epithet","infraspecies","cf_aff","sp_specif")]))])
  corres<-match(split(mat,row(mat)),split(mat,row(mat)))
  corresByCode<-tapply(corres,taxo[[attr(taxo,"taxoCode")]],unique,simplify=F)
  ln_corresByCode<-sapply(corresByCode,length)
  if(any(ln_corresByCode!=1)){
    lapply(corresByCode[ln_corresByCode>1],function(x,m){m[x,]},m=taxo)
    }

}

unicityGnInFam <- function(taxo,simplified=F)
{
  stopifnot(is(taxo,"taxo_oneTab"))
  gnFamTab<-data.frame(row=1:nrow(taxo),taxo[c(attr(taxo,"plot"),attr(taxo,"taxoCode"),attr(taxo,"family"),attr(taxo,"genus"))])
  gnFamTab<-gnFamTab[!is.na(gnFamTab[,attr(taxo,"genus")]),]
  gnInFam<-tapply(as.character(gnFamTab[[attr(taxo,"family")]]),as.character(gnFamTab[[attr(taxo,"genus")]]),function(x)sort(table(x),decreasing = T),simplify=F)
  if(any(sapply(gnInFam,length)!=1))
  {
    listProblems<-lapply(gnInFam[sapply(gnInFam,length)!=1],function(x)names(x[2:length(x)]))
    listSolutions<-lapply(gnInFam[sapply(gnInFam,length)!=1],function(x)names(x[1]))
    suggested_generic<-data.frame(genus=rep(names(listProblems),sapply(listProblems,length)),
                         family=unlist(listProblems),
                         suggested=rep(unlist(listSolutions),sapply(listProblems,length)))
    matGnFam<-as.matrix(gnFamTab[c(attr(taxo,"family"),attr(taxo,"genus"))])
    matSuggested_GnFam<-as.matrix(suggested_generic[c("family","genus")])
    m <- match(split(matGnFam,row(matGnFam)),split(matSuggested_GnFam,row(matSuggested_GnFam)))
    suggested <- data.frame(
      id_suggest=seq(1,sum(!is.na(m)),length.out=sum(!is.na(m))),
      row=gnFamTab$row[!is.na(m)],
      plot=gnFamTab[[attr(taxo,"plot")]][!is.na(m)],
      taxoCode=gnFamTab[[attr(taxo,"taxoCode")]][!is.na(m)],
      genus=suggested_generic$genus[m[!is.na(m)]],
      family=suggested_generic$family[m[!is.na(m)]],
      suggest_family=suggested_generic$suggested[m[!is.na(m)]]
    )
  }else{
    suggested_generic<-data.frame(genus=NULL,family=NULL,suggested=NULL)
    suggested<-data.frame(id_suggest=numeric(0),row=numeric(0),plot=character(0),taxoCode=character(0),genus=character(0),family=character(0),suggest_family=character(0))
  }
  colnames(suggested)[colnames(suggested)=="genus"]<-attr(taxo,"genus")
  colnames(suggested)[colnames(suggested)=="family"]<-attr(taxo,"family")
  colnames(suggested)[colnames(suggested)=="plot"]<-attr(taxo,"plot")
  colnames(suggested)[colnames(suggested)=="taxoCode"]<-attr(taxo,"taxoCode")
  colnames(suggested)[colnames(suggested)=="suggest_family"]<-paste("suggest",attr(taxo,"family"),sep="_")
  if(simplified)
  {return(suggested_generic)}
  return(suggested)
}



familtySearchGbif<-function(obj)
{NA}

genusSearchGbif <- function(obj)
{NA}

spSearchGbif <- function(obj)
{NA}
