
#' Put the taxonomy from the RDS in one table with attributes
#'
#' Technically this creates an object with a new class inherited from data.frame
#' The attributes of this class will ensure that we can always understand the
#' structure of the object, and come back to the initial format
#'
#' @param obj the base object with the taxonomy, may be a list with taxonomy
#' data.frame for each plot, or a single data.frame
#' @param currentFormat describe the format of obj: "listPlot" for a list of
#' plots and "oneTable" for a single data.frame containing all the taxonomy
#' @param taxoCode name of the column code corresponding to code of the taxon
#' @param taxonRanks_names named character vector (names:rank as defined in taxize, content name of the column)
#' @param taxonRanks_epithetized names of the columns (as in taxonRanks_names) which correspond to epithets
#' @param morphoQualifiers named character vector, the names are the types of qualifiers, the contents are the names of the object fields. For now it only work with cf_aff and sp_specif
#' @param comments name of the column containing the comments
#' @param plot name of the plot (name of the permanent plot) column
#'
#' @export
#'

new_taxo_oneTab <- function(obj,currentFormat=c("listPlot","oneTable"), taxonRanks_names = c(family = "family", genus = "genus", species = "specificEpithet", infraspecies = "infraspecificEpithet"), taxonRanks_epithetized=c("specificEpithet", "infraspecificEpithet"), taxoCode="code", plot="plot", morphoQualifiers=c(cf_aff = "identificationQualifier", sp_specif = "verbatimTaxonRank"), comments="comments")
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
  # clean function parameters in function of object colnames
  taxonRanks_names <- taxonRanks_names[taxonRanks_names %in% colnames(oneTab)]
  taxonRanks_epithetized <- taxonRanks_epithetized [taxonRanks_epithetized %in% taxonRanks_names]
  morphoQualifiers <- morphoQualifiers [morphoQualifiers %in% colnames(oneTab)]
  # Checking whether minimal information is in the table
  #stopifnot(family %in% colnames(oneTab))
  stopifnot(length(taxonRanks_names) >= 1)
  stopifnot(!all(taxonRanks_names %in% taxonRanks_epithetized))
  stopifnot(taxonRanks_epithetized %in% taxonRanks_names)
  #stopifnot(genus %in% colnames(oneTab))
  #stopifnot(species_epithet %in% colnames(oneTab))
  #stopifnot(taxoCode %in% colnames(oneTab))
  #stopifnot(plot %in% colnames(oneTab))
  # Create taxonRanks object
  separ<-strsplit(taxize::rank_ref$ranks,",")
  taxize_rank_ref<-data.frame(rankid=rep(as.integer(taxize::rank_ref$rankid),sapply(separ,length)),
                              rankname=unlist(separ))
  taxize_rank_ref<-taxize_rank_ref[1:which(taxize_rank_ref$rankname=="form"),]
  stopifnot(names(taxonRanks_names) %in% taxize_rank_ref$rankname)
  taxonRanks<-data.frame(
    rank=names(taxonRanks_names),
    column=taxonRanks_names,
    taxize_rankid=taxize_rank_ref$rankid[match(names(taxonRanks_names),taxize_rank_ref$rankname)]
  )
  taxonRanks$epithetized <- taxonRanks$column %in% taxonRanks_epithetized
  taxonRanks<-taxonRanks[order(taxonRanks$taxize_rankid),]
  rownames(taxonRanks)<-NULL
  if(!comments %in% colnames(oneTab)){comments<-NA}
  if(!plot %in% colnames(oneTab)){plot<-NA}
  if(!taxoCode %in% colnames(oneTab)){taxoCode<-NA}
  # Give the object its attributes
  oneTab<-structure(oneTab,
            origin=currentFormat,
            taxonRanks=taxonRanks,
            taxoCode=taxoCode,
            plot=plot,
            morphoQualifiers=morphoQualifiers,
            comments=comments)
  class(oneTab) <- c("taxo_oneTab","data.frame")
  return(oneTab)
}

#' Extract information from a taxonomic table
#'
#' @param taxo taxonomic table (of class taxo_oneTab)
#' @param parts parts of the taxonomic information you want to extract (from the following list: "taxonRanks","taxoCode","plot","morphoQualifiers","comments")
#' @param onlyRanks ranks (as defined in `taxize` package) to be extracted when "taxonRanks" in the `parts` if null all the ranks are extracted
#' @param onlyQualifiers qualifiers (cf_aff and/or sp_specif) to be extracted when "taxonRanks" in the `parts` if null all the ranks are extracted
#' @export
extract <-function(taxo,parts=c("taxonRanks","taxoCode","plot","morphoQualifiers","comments"),onlyRanks=NULL,onlyQualifiers=NULL)
{
  stopifnot(methods::is(taxo,"taxo_oneTab"))
  parts=match.arg(parts,choices=c("plot","taxoCode","taxonRanks","morphoQualifiers","comments"),several.ok = T)
  colToGet<-data.frame(gp=character(),cn=character())
  if("plot" %in% parts)
  {
    colToGet<-rbind(colToGet,data.frame(gp="plot",cn=attr(taxo,"plot")))
  }
  if("comments" %in% parts)
  {
    colToGet<-rbind(colToGet,data.frame(gp="comments",cn=attr(taxo,"comments")))
  }
  if("taxoCode" %in% parts)
  {
    colToGet<-rbind(colToGet,data.frame(gp="taxoCode",cn=attr(taxo,"taxoCode")))
  }
  if("taxonRanks" %in% parts)
  {
    if(!is.null(onlyRanks))
    {
      stopifnot(onlyRanks %in% attr(taxo,"taxonRanks")$rank)
      colToGet<-rbind(colToGet,data.frame(gp="taxonRanks",cn=attr(taxo,"taxonRanks")$column[match(onlyRanks,attr(taxo,"taxonRanks")$rank)]))
    }else{
      colToGet<-rbind(colToGet,data.frame(gp="taxonRanks",cn=attr(taxo,"taxonRanks")$column))
    }
  }
  if("morphoQualifier" %in% parts)
  {
    if(!is.null(onlyQualifiers))
    {
      stopifnot(onlyQualifiers %in% names(attr(taxo,"morphoQualifiers")))
      colToGet<-rbind(colToGet,data.frame(gp="taxonRanks",cn=attr(taxo,"morphoQualifiers")[match(onlyQualifiers,names(attr(taxo,"morphoQualifiers")))]))
    }else{
      colToGet<-rbind(colToGet,data.frame(gp="taxonRanks",cn=attr(taxo,"morphoQualifiers")))
    }
  }
  colToGet<-colToGet[order(match(colToGet$gp,parts)),]
  colToGet<-stats::na.omit(colToGet)
  return(taxo[colToGet$cn])
}

#' Get taxon information at a particular taxonomic rank from a taxonomic table
#'
#' @param taxo taxonomic table (of class taxo_oneTab)
#' @param rank TODO: document
#' @export
getRank <- function(taxo, rank)
{
  stopifnot(methods::is(taxo,"taxo_oneTab"))
  ATTR_TR <- attr(taxo, "taxonRanks")
  stopifnot(rank %in% ATTR_TR$rank)
  epi <- ATTR_TR[ATTR_TR$rank==rank,"epithetized"]
  if(!epi)
  {return(as.character(taxo[,ATTR_TR$column[ATTR_TR$rank==rank],drop=T]))}
  if(rank=="species") {otherRanks<-"genus"} else {otherRanks<-c("genus", "species")}
  tab <- extract(taxo,parts="taxonRanks",onlyRanks = c(otherRanks,rank))
  res <- do.call(paste, tab)
  res[unique(which(is.na(tab),arr.ind=T)[,"row"])]<-NA
  return(as.character(res))
}

#' Correct a taxonomic table according to a table of suggested corrections
#'
#' @param taxo taxonomic table (of class taxo_oneTab)
#' @param suggested TODO: document
#' @export
correct<-function(taxo, suggested)
{
  stopifnot(methods::is(taxo,"taxo_oneTab"))
  if(!methods::is(suggested,"data.frame"))
  {
    if("suggested" %in% names(suggested)) { suggested <- suggested$suggested }
  }
  if(!nrow(suggested)) { return(taxo) }
  stopifnot("row" %in% colnames(suggested) & any(grepl("suggest",colnames(suggested))))
  if("col" %in% colnames(suggested)) { type= "cell" } else { type <- "row" }
  if(type=="row")
  {
    colsToChange<-gsub("suggest_","",colnames(suggested)[grepl("suggest_",colnames(suggested))])
    colsToAdd<-colsToChange[!colsToChange%in%colnames(taxo)]
    separ<-strsplit(taxize::rank_ref$ranks,",")
    taxize_rank_ref<-data.frame(rankid=rep(as.integer(taxize::rank_ref$rankid),sapply(separ,length)),
                                rankname=unlist(separ))
    taxize_rank_ref <- taxize_rank_ref[1:which(taxize_rank_ref$rankname=="form"),]
    match_mat <- cbind(match(colsToAdd, taxize_rank_ref$rankname), match(colsToAdd, gsub(" ", "_", taxize_rank_ref$rankname)), match(colsToAdd, gsub("[[:space:]]","_",paste0(taxize_rank_ref$rankname,"epithet",sep="_"))))
    taxizeMatch_m <- apply(match_mat, 1, function(x)
      {
        res=stats::na.omit(x)
        if(length(res)==0){res<-NA}
        return(res)
      })
    if(any(!is.na(taxizeMatch_m)))
    {
      add<-data.frame(rank = taxize_rank_ref$rankname[stats::na.omit(taxizeMatch_m)],
                      column = colsToAdd[!is.na(taxizeMatch_m)],
                      taxize_rankid = taxize_rank_ref$rankid[stats::na.omit(taxizeMatch_m)],
                      epithetized=grepl("epithet",colsToAdd[!is.na(taxizeMatch_m)]))
      stopifnot(!add$taxize_rankid %in% attr(taxo,"taxonRanks")$taxize_rankid)
      attr(taxo,"taxonRanks") <- rbind(attr(taxo,"taxonRanks"), add)
      attr(taxo, "taxonRanks") <- attr(taxo, "taxonRanks")[order( attr(taxo, "taxonRanks")$taxize_rankid),]
    }
    if(any(colsToAdd=="identificationQualifier"))
    {
      stopifnot(!"cf_aff" %in% names(attr(taxo,"morphoQualifiers")))
      attr(taxo,"morphoQualifiers")<-c(attr(taxo,"morphoQualifiers"),cf_aff="identificationQualifier")
    }
    if(any(colsToAdd=="verbatimTaxonRanks"))
    {
      stopifnot(!"sp_specif" %in% names(attr(taxo,"morphoQualifiers")))
      attr(taxo,"morphoQualifiers")<-c(attr(taxo,"morphoQualifiers"),sp_specif="verbatimTaxonRanks")
    }
    if(any(colsToAdd == "comments"))
    {
      stopifnot(is.na(attr(taxo, "comments")))
      attr(taxo, "comments") <- "comments"
    }
    if(any(colsToAdd == "plot"))
    {
      stopifnot(is.na(attr(taxo, "plot")))
      attr(taxo, "plot") <- "plot"
    }
    if(any(colsToAdd == "taxoCode"))
    {
      stopifnot(is.na(attr(taxo, "taxoCode")))
      attr(taxo, "taxoCode") <- "taxoCode"
    }
    taxo[colsToAdd]<-NA
    taxo[suggested$row,colsToChange]<-suggested[,paste("suggest",colsToChange,sep="_")]
    return(taxo)
  }
}

#' Check for irregular spaces in a taxonomic table and suggest corrections
#'
#' @param taxo taxonomic table (of class taxo_oneTab)
#' @param parts TODO: document
#' @param show_ref TODO: document
#' @param show_space TODO: document
#' @export
checkSpace <- function(taxo, parts=c("plot","taxoCode","taxonRanks","morphoQualifiers"), show_ref = c(attr(taxo,"plot"),attr(taxo,"taxoCode")),show_space="#")
{
  stopifnot(methods::is(taxo,"taxo_oneTab"))

  mat<-as.matrix(extract(taxo,parts))
  m_col<-match(colnames(mat),colnames(taxo))
  w_space_start <- which(`dim<-`(grepl("^[[:space:]]+",mat),dim(mat)),arr.ind=T)
  w_space_end <- which(`dim<-`(grepl("[[:space:]]+$",mat),dim(mat)),arr.ind=T)
  w_space_multi <- which(`dim<-`(grepl("[[:space:]]{2,}",mat),dim(mat)),arr.ind=T)
  numCases<-nrow(w_space_start)+nrow(w_space_end)+nrow(w_space_multi)
  concernedCol<-unique(c(w_space_start[,"col"],w_space_end[,"col"],w_space_multi[,"col"]))
  concernedRow<-unique(c(w_space_start[,"row"],w_space_end[,"row"],w_space_multi[,"row"]))
  types<-data.frame(Reduce(rbind,list(w_space_start,w_space_end,w_space_multi)),
  type=c(rep("start",nrow(w_space_start)),rep("end",nrow(w_space_end)),rep("multi",nrow(w_space_multi))))
  stopifnot(show_ref %in% colnames(taxo))
  references <- taxo[concernedRow,show_ref]
  colnames(references)<-paste("ref",show_ref,sep="_")
  suggested<-as.data.frame(gsub("[[:space:]]{2,}"," ",gsub("[[:space:]]+$","",gsub("^[[:space:]]+","",mat[concernedRow,concernedCol,drop=F]))))
  if(length(concernedRow))
  {colnames(suggested)<-paste("suggest",colnames(suggested),sep="_")}
  res<-data.frame(id_suggest=seq(1,length(concernedRow),length.out=length(concernedRow)),
             row=concernedRow,
             references
             )
  if(length(concernedRow))
  {
    res<-cbind(res,
              as.data.frame(gsub("[[:space:]]",show_space,mat[concernedRow,concernedCol,drop=F])),
              suggested
               )
  }
  return(res)
}

#' Check for bad formatted undetermination qualifiers and suggest corrections
#'
#' @param taxo taxonomic table (of class taxo_oneTab)
#' @export
checkUndetermited <- function(taxo)
# Note: when there are cases of Indet or Morpho, species code should be checked for relevant information
{
  stopifnot(methods::is(taxo,"taxo_oneTab"))
  mat<-as.matrix(extract(taxo,"taxonRanks"))
  #if(is.na(attr(taxo,"infraspecies_epithet"))) {mat<-cbind(mat, infraspecificEpithet = rep(NA,nrow(mat)))}
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

  if("cf_aff" %in% names(attr(taxo,"morphoQualifiers")))
  {
    raw_cf_aff<-taxo[attr(taxo,"morphoQualifiers")["cf_aff"]]
    raw_cf_aff[grepl("^[[:space:]]*$",raw_cf_aff)]<-NA
    suggested_cf_aff[suggested_cf_aff==""]<-raw_cf_aff[suggested_cf_aff==""]
  }else{
    suggested_cf_aff[suggested_cf_aff==""]<-NA
  }

  if("sp_specif" %in% names(attr(taxo,"morphoQualifiers")))
  {
    raw_sp_specif<-taxo[attr(taxo,"morphoQualifiers")["sp_specif"]]
    raw_sp_specif[grepl("^[[:space:]]*$",raw_sp_specif)]<-NA
    suggested_sp_specif[suggested_sp_specif==""]<-raw_sp_specif[suggested_sp_specif==""]
  }else{
    suggested_sp_specif[suggested_sp_specif==""]<-NA
  }

  anyPb<-sort(unique(c(w_Indet[,"row"], w_empty[,"row"], w_cf[,"row"], w_aff[,"row"], w_sp[,"row"], w_spNb[,"row"], w_morfo[,"row"])))
  tabRaw<-extract(taxo,c("plot","taxoCode","taxonRanks","morphoQualifiers"))
  tabSuggest<-as.data.frame(mat2)
  if(any(!is.na(suggested_cf_aff)))
  {
    name_cf_aff<-ifelse("cf_aff" %in% names(attr(taxo,"morphoQualifiers")),attr(taxo,"morphoQualifiers")["cf_aff"],"identificationQualifier")
    tabSuggest[[name_cf_aff]] <- suggested_cf_aff
  }
  if(any(!is.na(suggested_sp_specif)))
  {
    name_sp_specif<-ifelse("sp_specif" %in% names(attr(taxo,"morphoQualifiers")),attr(taxo,"morphoQualifiers")["sp_specif"],"verbatimTaxonRank")
    tabSuggest[[name_sp_specif]] <- suggested_sp_specif
  }
  colnames(tabSuggest)<-paste0("suggest_",colnames(tabSuggest))
  return(data.frame(id_suggest=seq(1,length(anyPb),length.out=length(anyPb)),row=anyPb,tabRaw[anyPb,],tabSuggest[anyPb,]))
}


#' Check whether taxa of a particular rank are consistently always in the same higher rank and suggest corrections
#'
#' @param taxo taxonomic table (of class taxo_oneTab)
#' @param rank TODO: document
#' @param superior TODO: document
#' @export
checkUnicityRankSup <- function(taxo,rank="genus",superior="family")
{
  stopifnot(methods::is(taxo,"taxo_oneTab"))
  ATTR_TR<-attr(taxo,"taxonRanks")
  stopifnot(c(rank,superior) %in% ATTR_TR$rank)
  stopifnot(ATTR_TR$taxize_rankid[ATTR_TR$rank==rank]>ATTR_TR$taxize_rankid[ATTR_TR$rank==superior])
  TAB<-data.frame(row=1:nrow(taxo),superior=getRank(taxo,superior), rank=getRank(taxo,rank))
  TAB<-TAB[!is.na(TAB$rank),]
  majority<-tapply(TAB$superior,TAB$rank,function(x)names(sort(table(x),decreasing=T))[1])
  m<-match(TAB$rank,names(majority))
  rowConcerned<-TAB$row[TAB$superior!=majority[m]]
  suggest<-data.frame(suggest=majority[m][TAB$superior!=majority[m]])
  colnames(suggest)<-paste("suggest",ATTR_TR$column[ATTR_TR$rank==superior],sep="_")
  return(data.frame(id_suggest=seq(1,length(rowConcerned),length.out=length(rowConcerned)),
                  row=rowConcerned,
                  extract(taxo,parts=c("plot","taxoCode","taxonRanks"),onlyRanks=c(rank,superior))[rowConcerned,],
                  suggest
    ))
}


#' Check whether taxonomic codes are alway associated with the same taxonomic information and suggests corrections
#'
#' @param taxo taxonomic table (of class taxo_oneTab)
#' @param noMajority TODO: document
#' @export
checkUnicityCodetax <- function(taxo, noMajority=c("skip","takeFirst","stop"))
{
  stopifnot(methods::is(taxo,"taxo_oneTab"))
  noMajority<-match.arg(noMajority)
  mat<-as.matrix(extract(taxo,c("taxonRanks","morphoQualifiers")))
  corres<-match(split(mat,row(mat)),split(mat,row(mat)))
  corresByCode<-tapply(corres,taxo[[attr(taxo,"taxoCode")]],unique,simplify=F)
  ln_corresByCode<-sapply(corresByCode,length)
  variablesToShow<-colnames(extract(taxo,c("plot","taxoCode","taxonRanks","morphoQualifiers")))
  varSuggest<-colnames(extract(taxo,c("taxonRanks","morphoQualifiers")))
  if(any(ln_corresByCode!=1)){
    nbCasesByCode <- tapply(corres,extract(taxo,"taxoCode")[,1,drop=T], table)
    orderedCases <- lapply(nbCasesByCode[sapply(nbCasesByCode, length) > 1], sort, decreasing=T)
    problems<-lapply(names(orderedCases),function(n,t,c,o,m){
      res<-data.frame(
        row=which(t[,attr(t,"taxoCode")]==n),
        t[t[,attr(t,"taxoCode")]==n,c]
      )
      m2<-as.matrix(res[colnames(m)])
      matchMat<-match(split(m2,row(m2)),split(m,row(m)))
      res$suggestedCase<-matchMat==as.integer(names(o[[n]])[1])
      return(res)
      },
      t=taxo,
      c=variablesToShow,
      o=orderedCases,
      m=mat)
    names(problems)<-names(orderedCases)
    noMaj<-sapply(orderedCases,function(x)x[1]==x[2])
    if(any(noMaj) & noMajority=="skip")
    {
      warning("These codes corresponds to different information, but we were unable to find the most represented set of variables:\n", paste(names(orderedCases)[noMaj],collapse=" "),"\n Note that NO CORRECTION WILL BE SUGGESTED FOR THESE CASES\n")
      problems<-mapply(function(x,n){if(n)x$suggestedCase<-NA;return(x)},problems,noMaj,SIMPLIFY=F)
    }
    if(any(noMaj) & noMajority=="takeFirst")
    {
      warning("These codes corresponds to different information, but we were unable to find the most represented set of variables:\n", paste(names(orderedCases)[noMaj],collapse=" "),"\n Note the THE FIRST SET OF VARIABLES ASSOCIATED WITH THE CODE WILL BE ARBITRARILY USED AS SUGGESTED CORRECTION FOR THE OTHER ONES\n")
    }
    if(any(noMaj) & noMajority=="stop")
    {
      stop("These codes corresponds to different information, but we were unable to find the most represented set of variables:\n", paste(names(orderedCases)[noMaj],collapse=" "))
    }
    suggested<-Reduce(rbind,lapply(problems[!(noMaj&noMajority=="skip")],function(pb,vs){
      tab1<-pb[!pb$suggestedCase,,drop=F]
      row2<-unique(pb[pb$suggestedCase,vs,drop=F])
      stopifnot(nrow(row2)==1)
      colnames(row2)<-paste("suggest",colnames(row2),sep="_")
      return(cbind(tab1,row2))
    },vs=varSuggest))
    suggested<-data.frame(id_suggest=seq(1,nrow(suggested),length.out=nrow(suggested)),suggested[names(suggested)!="suggestedCase"])
  return(list(problems=problems,suggested=suggested))
  }
  cn<-c("row",variablesToShow,paste("suggest",varSuggest,sep="_"))
  names(cn)<-cn
  return(list(problems=NULL,suggested=as.data.frame(lapply(cn,function(x){return(NULL)}))))
}

#' Format and show results from the analysis of unicity of information associated with taxonomic codes
#'
#' @param resUnicityCodeTax TODO: document
#' @param type TODO: document
#' @param code TODO: document
#' @export
showUnicityCodetax<-function(resUnicityCodeTax,type=c("code","problems","suggested"),code=NA)
{
  type<-match.arg(type)
  if(type=="code")
  {
    return(names(resUnicityCodeTax$problems))
  }
  if(type=="problems")
  {
    if(is.na(code)){return(resUnicityCodeTax$problems)}else{
        if(length(code)==1){return(resUnicityCodeTax$problems[[code]])}else{return(resUnicityCodeTax$problems[code])
        }}
  }
  if(type=="suggested")
  {
    if(is.na(code)){return(resUnicityCodeTax$suggested)}else{return(resUnicityCodeTax$suggested[resUnicityCodeTax$suggested[,3] %in% code])
    }}
}



#' Analyse a table obtained from a gbif backbone search
#'
#' @param searched TODO: document
#' @param tabGbif TODO: document
#' @param rank TODO: document
#' @param obligatory TODO: document
#' @param expected TODO: document
#' @param returnGbifRes TODO: document
#' @export
analyseGbifTable <- function(searched, tabGbif,rank, obligatory, expected, returnGbifRes=T)
{
  separ<-strsplit(taxize::rank_ref$ranks,",")
  taxize_rank_ref<-data.frame(rankid=rep(as.integer(taxize::rank_ref$rankid),sapply(separ,length)),
                              rankname=unlist(separ))
  taxize_rank_ref<-taxize_rank_ref[1:which(taxize_rank_ref$rankname=="form"),]
  orderedRanks=taxize_rank_ref$rankname
  if(sum(sapply(obligatory,length)))
  {
    stopifnot(names(obligatory) %in% names(tabGbif) & methods::is(obligatory,"character"))
    matObliGbif<-as.matrix(tabGbif[names(obligatory)])
    matchOblig<-!is.na(match(split(matObliGbif,row(matObliGbif)),list(unname(obligatory))))
  }
  rankOk <- (tabGbif$rank %in% rank)
  # case: perfect match without any discrepancy nor problematic taxonomic information
  minimumRequired<-which(matchOblig & rankOk)
  if(!length(minimumRequired))
  {
    res<-list(type = "Failed",  canonicalname = NA,  authorship = NA, gbifid = NA, initialRank = NA, finalRank = NA, synonym = data.frame(gbifid = numeric(), scientificname = character()), modifiedHR = NULL, higherRanks = data.frame(rank = character(),  canonicalname = character(),  gbifid = numeric()), parent_gbifid = numeric())
    if(returnGbifRes){res$gbifTab <- tabGbif}
    return(res)
  }
  tabFiltered <- tabGbif[minimumRequired,,drop=F]
  exact<-(searched == tabFiltered$canonicalname)
  if(length(expected)&&sum(!is.na(expected)))
  {
    stopifnot(names(expected) %in% names(tabGbif) & methods::is(expected,"character"))
    matExpectGbif<-as.matrix(tabFiltered[names(expected)])
    matchExpect<-!is.na(match(split(matExpectGbif,row(matExpectGbif)),list(unname(expected))))
  }else{matchExpect<-rep(T,nrow(tabFiltered))}
  chosen<-order(exact,matchExpect,decreasing=T)[1]
  rowChosen<-tabFiltered[chosen,,drop=F]
  # Determining category
  exactChosen <- exact[chosen]
  matchExpectChosen <- matchExpect[chosen]
  statusChosen<-c(ACCEPTED="ok","SYNONYM"="synonym","HOMOTYPIC_SYNONYM"="synonym")[rowChosen$status]
  if(is.na(statusChosen))
  {
    res<-list(type = "Failed",  canonicalname = NA,  authorship = NA, gbifid = NA, initialRank = NA, finalRank = NA, synonym = data.frame(gbifid = numeric(), scientificname = character()), modifiedHR = NULL, higherRanks = data.frame(rank = character(),  canonicalname = character(),  gbifid = numeric()), parent_gbifid = numeric())
    if(returnGbifRes){res$gbifTab <- tabGbif}
    return(res)
  }
  type_var<-c(ifelse(exactChosen,"exactMatch","fuzzyMatch"), if(!statusChosen=="ok") {statusChosen} else {NULL}, if(matchExpectChosen) {NULL} else {"changeHigherRanks"})
  addrRank<-which(orderedRanks == rowChosen$rank)
  potentialHR<-orderedRanks[1:(addrRank-1)]
  colHR<-potentialHR[potentialHR %in% names(rowChosen)]
  colHR<-colHR[!is.na(rowChosen[colHR])]
  colHRkey <- paste0(colHR,"key")
  stopifnot(all(colHRkey %in% colnames(rowChosen)))
  canonicalname<-rowChosen[,rowChosen$rank]
  gbifid<-rowChosen[,paste0(rowChosen$rank,"key")]
  extractableAuthorship<-grepl(paste0("^",canonicalname," ?(.+)$"), rowChosen$scientificname)
  authorship <- ifelse(extractableAuthorship,gsub(paste0("^",canonicalname," ?(.+)$"),"\\1",rowChosen$scientificname),NA)
  if(statusChosen=="synonym")
  {
    synonym=data.frame(gbifid=rowChosen$usagekey,scientificname=rowChosen$scientificname)
  }else{
    synonym=data.frame(gbifid=numeric(),scientificname=character())
  }
  if(!matchExpectChosen)
  {
    modifiedHR<-names(expected)
  }else{
    modifiedHR<-NULL
  }
  higherRanks <- data.frame(rank=colHR, canonicalname= as.character(rowChosen[colHR]), gbifid=as.numeric(rowChosen[colHRkey]), row.names=NULL)
  higherRanks$parent_gbifid <- c(NA,higherRanks$gbifid[1:(nrow(higherRanks )-1)])
  res=list(
    type=paste(type_var,collapse="_"),
    canonicalname=canonicalname,
    authorship=authorship,
    gbifid=gbifid,
    initialRank=rank,
    finalRank=rowChosen$rank,
    synonym=synonym,
    modifiedHR=modifiedHR,
    higherRanks=higherRanks)
  if(returnGbifRes){res$gbifTab <- tabGbif}
  return(res)
}

#' Get suggested taxonomic information from an analysed table from a search in the Gbif backbone
#'
#' @param analysedGbif TODO: document
#' @param ranks TODO: document
#' @export
getSuggestsFromResGbif <- function(analysedGbif,ranks=c("family","genus","species","subspecies","variety"))
{
  separ<-strsplit(taxize::rank_ref$ranks,",")
  taxize_rank_ref<-data.frame(rankid=rep(as.integer(taxize::rank_ref$rankid),sapply(separ,length)),
                              rankname=unlist(separ))
  taxize_rank_ref<-taxize_rank_ref[1:which(taxize_rank_ref$rankname=="form"),]
  orderedRanks=taxize_rank_ref$rankname
  stopifnot(ranks %in% orderedRanks)
  if(!is.na(analysedGbif$finalRank) && !analysedGbif$finalRank %in% ranks)
  {warning("The rank of the taxon found in the GBIF backbone (",analysedGbif$finalRank,") is not in the extracted ranks")}
  res<-sapply(ranks,function(analysedGbif)NA,USE.NAMES = T)
  if(analysedGbif$finalRank %in% names(res)) {res[analysedGbif$finalRank]<-analysedGbif$canonicalname}
  inHR<-match(ranks,analysedGbif$higherRanks$rank)
  res[!is.na(inHR)]<-analysedGbif$higherRanks$canonicalname[stats::na.omit(inHR)]
  return(data.frame(type=analysedGbif$type,as.data.frame(as.list(res)),gbifid=analysedGbif$gbifid))
}

#' Get suggested taxa from a list of analysed tables from Gbif backbone searches
#'
#' @param listAnalysedGbif TODO: document
#' @param ranks TODO: document
#' @param exclude TODO: document
#' @export
getSuggestsFromListGbif <- function(listAnalysedGbif,ranks=c("family","genus","species","subspecies","variety"), exclude=c("exactMatch"))
{
  suggested<-Reduce(rbind,lapply(listAnalysedGbif,getSuggestsFromResGbif,ranks=ranks))
  rownames(suggested)<-names(listAnalysedGbif)
  suggested<-suggested[!suggested$type%in%exclude,]
  return(suggested)
}

#' Extract epithets from full name taxa
#'
#' @param taxon TODO: document
#' @param higherTaxon TODO: document
#' @export
extractEpithet <- function(taxon,higherTaxon=NULL)
{
  epithet <- sapply(strsplit(taxon,"[[:space:]]"),function(x)x[length(x)])
  if(!is.null(higherTaxon))
  {
    w_problem <- which(!is.na(taxon) & paste(higherTaxon,epithet,sep=" ")!=taxon)
    if(length(w_problem)){
      warning("While extracting the epithet we found discrepancies between the taxon names and higher taxon names \nProblematic cases:\n",
            paste("taxon:",taxon[w_problem],"; epithet:",epithet[w_problem],"; higher taxon:",higherTaxon[w_problem],collapse = "\n"))
    }
  }
  return(epithet)
}


#' Extract the taxonomic ranks of taxa from a taxonomic table
#'
#' @param taxo taxonomic table (of class taxo_oneTab)
#' @export
taxoRanks<-function(taxo)
{
  stopifnot(methods::is(taxo,"taxo_oneTab"))
  ATTR_TR<-attr(taxo,"taxonRanks")
  stopifnot(order(ATTR_TR$taxize_rankid)==1:nrow(ATTR_TR))
  mat_taxRanks<-as.matrix(taxo[ATTR_TR$column])
  lowerInfo<-apply(mat_taxRanks,1,function(x,n)
  {
    if(all(is.na(x))){return(1)}
    max(which(!is.na(x)))+1
  })
  higherMissing<-apply(mat_taxRanks,1,function(x,n)
  {
    if(!any(is.na(x))){return(n+1)}
    return(min(which(is.na(x))))
  },n=nrow(ATTR_TR))
  if(any(lowerInfo!=higherMissing))
  {warning("Some of the taxa have missing higher information (see rows ",paste0(which(lowerInfo!=higherMissing),collapse=" "))}
  LEV<-c("higher",ATTR_TR$rank)
  return(factor(LEV[lowerInfo],levels=LEV, ordered = T))
}

#' Determine which higher ranks can be found from a taxonomic table
#'
#' @param taxo taxonomic table (of class taxo_oneTab)
#' @param rank TODO: document
#' @param excludeEpithComponents TODO: document
#'
#' @export
higherRanks <- function(taxo, rank, excludeEpithComponents=T)
{
  ATTR_TR<-attr(taxo,"taxonRanks")
  w<-which(ATTR_TR$rank==rank)
  if(w==1){return(character(0))}
  HR<-ATTR_TR$rank[1:(w-1)]
  if(excludeEpithComponents & ATTR_TR[w,"epithetized"])
  {
    toSupp<-c("genus","species")
    HR<-HR[!HR %in% toSupp]
  }
  return(HR)
}



#' Check in the gbif backbone whether taxa from a particular rank in a table are erroneous and suggests corrections
#'
#' @param taxo taxonomic table (of class taxo_oneTab)
#' @param rankCheck TODO: document
#' @param higherRankCheck TODO: document
#' @param messagesGbif TODO: document
#' @param messagesOther TODO: document
#' @param filters TODO: document
#' @param excludeFromSuggests TODO: document
#' @param epithetizeSuggests TODO: document
#' @param returnGbifRes TODO: document
#' @param returnAnalysedGbif TODO: document
#' @param returnFailed TODO: document
#' @param returnExactMatches TODO: document
#' @export
checkGbif<-function(taxo,rankCheck,higherRankCheck=higherRanks(taxo,rankCheck),messagesGbif=F,messagesOther=T,filters=c(kingdom="Plantae"), excludeFromSuggests=c("exactMatch","Failed"), epithetizeSuggests=T, returnGbifRes=F, returnAnalysedGbif = !returnGbifRes, returnFailed="Failed" %in% excludeFromSuggests, returnExactMatches=F)
{
  #Checking conditions for application
  stopifnot(methods::is(taxo,"taxo_oneTab"))
  ATTR_TR<-attr(taxo,"taxonRanks")
  stopifnot(length(rankCheck==1) & rankCheck %in% ATTR_TR$rank)
  stopifnot(higherRankCheck %in% ATTR_TR$rank)
  stopifnot(ATTR_TR$taxize_rankid[ATTR_TR$rank %in% higherRankCheck] < ATTR_TR$taxize_rankid[ATTR_TR$rank == rankCheck])
  # preparing the filter for rows corresponding to the rank
  rankOk <- taxoRanks(taxo)==rankCheck
  stopifnot(sum(rankOk)>0)

  # preparing the taxa to search
  tax <- getRank(taxo[rankOk,,drop=F], rankCheck)
  toSearch <- unique(tax)
  if(messagesOther)
  {message("Searching for ",length(toSearch)," taxa in the GBIF Backbone\n...")}
  resGbifTax<-taxize::get_gbifid_(toSearch,messages=messagesGbif)
  if(messagesOther)
  {message("done")}

  # preparing the data for analysing gbif search results (higher ranks etc)
  HR_tab<-as.data.frame(lapply(higherRankCheck,getRank,taxo=taxo[rankOk,,drop=F]))
  if(!length(higherRankCheck)){HR_tab<-data.frame(matrix(nrow=sum(rankOk),ncol=0))}
  colnames(HR_tab)<-higherRankCheck
  tabRanks<-data.frame(cbind(HR_tab,tax))
  unTabRanks<-unique(tabRanks)
  m<-match(unTabRanks$tax,toSearch)
  searchedTax<-toSearch[m]
  resGbifTax<-resGbifTax[m]
  HR<-apply(unTabRanks,1,function(x,cn){res<-x[cn!="tax"]; names(res)<-cn[cn!="tax"] ; return(res)},cn=colnames(unTabRanks),simplify = F)

  # analysing gbif results
  if(messagesOther)
  {message("done\nAnalysing GBIF Backbone information")}
  gbifAnalysed<-mapply(analyseGbifTable,searched=searchedTax, tabGbif=resGbifTax,expected=HR, MoreArgs=list(rank=rankCheck, obligatory= filters),SIMPLIFY = F)
  diagnostic<-table(sapply(gbifAnalysed,function(x)x$type))
  if(messagesOther)
  {
    message(diagnostic["exactMatch"]," taxa are found without any modification needed")
    message(sum(diagnostic[grep("fuzzy",names(diagnostic))])," taxa are found with suggested orthographic changes")
    message(sum(diagnostic[grep("synonym",names(diagnostic))])," taxa are suggested synonyms")
    message(sum(diagnostic[grep("changeHigherRanks",names(diagnostic))])," taxa are found with suggested higher rank changes")
    message(sum(diagnostic[grep("Failed",names(diagnostic))])," taxa were not found")
  }

  # getting suggested taxa (epithetized if needed)
  rankSuggests<-c(higherRanks(taxo,rankCheck,excludeEpithComponents = F),rankCheck)
  suggestedTax<-getSuggestsFromListGbif(gbifAnalysed,ranks=rankSuggests,exclude=NULL)
  if(epithetizeSuggests)
  {
    epithetizable<-c("species","infraspecies","subspecies","variety","form")
    ATTR_epithetizable<-ATTR_TR[stats::na.omit(match(epithetizable,ATTR_TR$rank)),]
    supp<-ATTR_epithetizable$rank[!ATTR_epithetizable$epithetized]
    toEpithetize<-epithetizable[!epithetizable %in% supp & epithetizable %in% rankSuggests]
    if(length(toEpithetize))
    {
      higherRanks_epi<-rep("species",length(toEpithetize))
      higherRanks_epi[toEpithetize=="species"]<-"genus"
      epithetized_suggest <- as.data.frame(mapply(extractEpithet,as.list(suggestedTax[,toEpithetize,drop=F]),suggestedTax[,higherRanks_epi,drop=F],SIMPLIFY = F))
      suggestedTax[toEpithetize]<-epithetized_suggest
    }
  }else{toEpithetize=character()}
  # define colnames for ranks
  tabDecision<-data.frame(cn=rankSuggests,
             inATTR=rankSuggests%in%ATTR_TR$rank,
             epiHere=rankSuggests%in%toEpithetize,
             epiATTR=NA,
             nameATTR=NA,
             naturalName=ifelse(rankSuggests%in%toEpithetize, paste(rankSuggests,"epithet",sep="_"), rankSuggests),
             finalName=NA
             )
  m<-match(tabDecision$cn[tabDecision$inATTR],ATTR_TR$rank)
  tabDecision$nameATTR[tabDecision$inATTR]<-ATTR_TR[m,"column"]
  tabDecision$epiATTR[tabDecision$inATTR]<-ATTR_TR[m,"epithetized"]
  tabDecision$finalName[is.na(tabDecision$nameATTR)|(tabDecision$epiHere!=tabDecision$epiATTR)]<-tabDecision$naturalName[is.na(tabDecision$nameATTR)|(tabDecision$epiHere!=tabDecision$epiATTR)]
  tabDecision$finalName[!is.na(tabDecision$nameATTR)&(tabDecision$epiHere==tabDecision$epiATTR)]<-tabDecision$nameATTR[!is.na(tabDecision$nameATTR)&(tabDecision$epiHere==tabDecision$epiATTR)]
  stopifnot(!is.na(tabDecision$finalName))
  tabDecision$finalName<-paste("suggest",tabDecision$finalName,sep="_")
  mcn<-match(colnames(suggestedTax),tabDecision$cn)
  colnames(suggestedTax)[!is.na(mcn)]<-tabDecision$finalName[stats::na.omit(mcn)]

  # Getting what rows are concerned
  suggestedAllRows<-data.frame(row=(1:nrow(taxo))[rankOk],
                               extract(taxo,c("taxoCode","plot","taxonRanks"),onlyRanks = rankSuggests)[rankOk,],
    suggestedTax[match(split(as.matrix(tabRanks),row(as.matrix(tabRanks))), split(as.matrix(unTabRanks),row(as.matrix(unTabRanks)))),]
  ,row.names=NULL)
  suggested<-suggestedAllRows[!suggestedAllRows$type%in%excludeFromSuggests,]
  res<-list()
  if(returnGbifRes)
  {
    res$gbifRes <- resGbifTax
  }
  if(returnAnalysedGbif)
  {
    res$analysedGbif <- gbifAnalysed
  }
  if(returnFailed)
  {
    res$failed <- suggestedAllRows[suggestedAllRows$type=="Failed",]
  }
  if(returnExactMatches)
  {
    res$exactMatches <- suggestedAllRows[suggestedAllRows$type=="exactMatch",]
  }

  res$rowRankOk<-which(rankOk)
  res$suggested<-suggested
  return(res)
}
