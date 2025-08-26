
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
#! @param infraspecies_epithet name of the infraspecific epithet column
#! @param taxoCode name of the column code corresponding to code of the taxon
#! @param plot name of the plot (name of the permanent plot) column
#! @param cf_aff name of the column containing the cf. or aff. expression
#! @param sp_specif name of the column containg information about the morphospecies


# obj=rdsBST$taxonomy
# currentFormat="listPlot"
# family="family"
# genus="genus"
# species_epithet="specificEpithet"
# infraspecies_epithet="infraspecificEpithet"
# taxoCode="code"
# plot="plot"
# cf_aff="identificationQualifier"
# sp_specif="verbatimTaxonRank"

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
  taxonRanks<-taxonRanks[order(taxonRanks$taxize_rankid)]
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

extract <-function(taxo,parts=c("taxonRanks","taxoCode","plot","morphoQualifiers","comments"),onlyRanks=NULL,onlyQualifiers=NULL)
{
  stopifnot(is(taxo,"taxo_oneTab"))
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
  colToGet<-na.omit(colToGet)
  return(taxo[colToGet$cn])
}


checkSpace <- function(taxo, parts=c("plot","taxoCode","taxonRanks","morphoQualifiers"), show_ref = c(attr(taxo,"plot"),attr(taxo,"taxoCode")),show_space="#")
{
  stopifnot(is(taxo,"taxo_oneTab"))

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
  colnames(suggested)<-paste("suggest",colnames(suggested),sep="_")
  return(data.frame(id_suggest=seq(1,length(concernedRow),length.out=length(concernedRow)),
             row=concernedRow,
             references,
             as.data.frame(gsub("[[:space:]]",show_space,mat[concernedRow,concernedCol,drop=F])),
             suggested
             ))
}

# correctSpace<-function(taxo,suggested)
# {
#   if(nrow(suggested)>0){
#     taxo[as.matrix(suggested[c("row","col")])]<-suggested$suggested
#   }
#   return(taxo)
# }

checkUndetermited <- function(taxo, correct = F)
# Note: when there are cases of Indet or Morpho, species code should be checked for relevant information
{
  stopifnot(is(taxo,"taxo_oneTab"))
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

# correctUndetermined <- function(taxo, suggested)
# {
#   stopifnot(is(taxo,"taxo_oneTab"))
#   cn <- colnames(suggested)[grep("suggest_",colnames(suggested))]
#   resCol<-gsub("suggest_","",cn)
#   colToAdd<-resCol[!resCol%in%colnames(taxo)]
#   if(length(colToAdd)>0){
#     taxo[colToAdd]<-NA
#     if(any(colToAdd=="identificationQualifier")&is.na(attr(taxo,"cf_aff")))
#     {attr(taxo,"cf_aff")<-"identificationQualifier"}
#     if(any(colToAdd=="verbatimTaxonRank"&is.na(attr(taxo,"sp_specif"))))
#     {attr(taxo,"sp_specif")<-"verbatimTaxonRank"}
#     }
#   w_suggested<-as.matrix(cbind(row=as.numeric(row(suggested[,cn])),col=as.numeric(col(suggested[,cn]))))
#   w_taxo <- cbind(row=suggested$row[w_suggested[,"row"]],col=match(resCol,colnames(taxo))[w_suggested[,"col"]])
#   taxo[w_taxo]<-suggested[,cn][w_suggested]
#   return(taxo)
# }


checkUnicityGnInFam <- function(taxo,simplified=F)
{
  stopifnot(is(taxo,"taxo_oneTab"))
  stopifnot(c("genus","family") %in% attr(taxo, "taxonRanks")$rank)
  gnFamTab<-data.frame(row=1:nrow(taxo),extract(taxo,c("plot","taxoCode","taxonRanks"),onlyRanks = c("family","genus")))
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

correctUnicityGnInFam <- function(taxo, suggested)
{
  if(nrow(suggested)>0)
  {
    taxo[suggested$row,attr(taxo,"family")]<-suggested[,paste("suggest",attr(taxo,"family"),sep="_")]
  }
  return(taxo)
}

checkUnicityCodetax <- function(taxo, noMajority=c("skip","takeFirst","stop"))
{
  stopifnot(is(taxo,"taxo_oneTab"))
  noMajority<-match.arg(noMajority)
  mat<-as.matrix(taxo[na.omit(unlist(attributes(taxo)[c("family","genus","species_epithet","infraspecies_epithet","cf_aff","sp_specif")]))])
  corres<-match(split(mat,row(mat)),split(mat,row(mat)))
  corresByCode<-tapply(corres,taxo[[attr(taxo,"taxoCode")]],unique,simplify=F)
  ln_corresByCode<-sapply(corresByCode,length)
  if(any(ln_corresByCode!=1)){
    nbCasesByCode <- tapply(corres,taxo[[attr(taxo,"taxoCode")]], table)
    orderedCases <- lapply(nbCasesByCode[sapply(nbCasesByCode, length) > 1], sort, decreasing=T)
    variablesToShow<-na.omit(unlist(attributes(taxo)[c("plot","taxoCode","family","genus","species_epithet","infraspecies_epithet","cf_aff","sp_specif")]))
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
    varSuggest<-na.omit(unlist(attributes(taxo)[c("family","genus","species_epithet","infraspecies_epithet","cf_aff","sp_specif")]))
    suggested<-Reduce(rbind,lapply(problems[!(noMaj&noMajority=="skip")],function(pb,vs){
      tab1<-pb[!pb$suggestedCase,,drop=F]
      row2<-unique(pb[pb$suggestedCase,vs,drop=F])
      stopifnot(nrow(row2)==1)
      colnames(row2)<-paste("suggest",colnames(row2),sep="_")
      return(cbind(tab1,row2))
    },vs=varSuggest))
  return(list(problems=problems,suggested=suggested))
  }
  cn<-c("row",na.omit(unlist(attributes(taxo)[c("plot","taxoCode","family","genus","species_epithet","infraspecies_epithet","cf_aff","sp_specif")])),paste("suggest",na.omit(unlist(attributes(taxo)[c("family","genus","species_epithet","infraspecies_epithet","cf_aff","sp_specif")])),sep="_"))
  names(cn)<-cn
  return(list(problems=NULL,suggested=as.data.frame(lapply(cn,function(x){return(NULL)}))))
}

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

correct<-function(taxo, suggested)
{
  stopifnot(is(taxo,"taxo_oneTab"))
  if(!is(suggested,"data.frame"))
  {
    if("suggested" %in% names(suggested)) { suggested <- suggested$suggested }
  }
  stopifnot("row" %in% colnames(suggested) & any(grepl("suggest",colnames(suggested))))
  if("col" %in% colnames(suggested)) { type= "cell" } else { type <- "row" }
  if(type=="row")
  {
    colsToChange<-gsub("suggest_","",colnames(suggested)[grepl("suggest_",colnames(suggested))])
    colsToAdd<-colsToChange[!colsToChange%in%colnames(taxo)]
    if("infraspecificEpithet" %in% colsToAdd & is.na(attr(taxo,"infraspecies_epithet")))
    {
      taxo$infraspecificEpithet<-NA
      attr(taxo,"infraspecies_epithet")<-"infraspecificEpithet"
    }
    if(any(colsToAdd=="identificationQualifier")&is.na(attr(taxo,"cf_aff")))
    {
      taxo$identificationQualifier<-NA
      attr(taxo,"cf_aff")<-"identificationQualifier"
    }
    if(any(colsToAdd=="verbatimTaxonRank"&is.na(attr(taxo,"sp_specif"))))
    {
      taxo$verbatimTaxonRank<-NA
      attr(taxo,"sp_specif")<-"verbatimTaxonRank"
    }
    taxo[suggested$row,colsToChange]<-suggested[,paste("suggest",colsToChange,sep="_")]
    return(taxo)
  }
}

# n<-which(speciesToSearch == "Termilia amazonica")
# n<-which(speciesToSearch == "Aristolochia rigens")
# n<-which(speciesToSearch == "Protium guianense")
# n<-which(speciesToSearch == "Phanera guianensis")
# searched <- speciesToSearch[n]
# tabGbif <- resGbifSpe[[n]]
# rank <- "species"
# obligatory <- c(kingdom = "Plantae")
# #obligatory <- c(kingdom = "Plantae", phylum = "Tracheophyta")
# expected <- c(family = familySpecies[n])
# reference <- searched
# toCompare <- tabGbif$canonicalname

# naSuppressed <- function(reference, toCompare)
# {
#   sepChar <- strsplit(reference,"")[[1]]
#   sepCharToCompare <- strsplit(toCompare,"")
#   sapply(sepCharToCompare, function(x, y)
#     {
#       if(length(x)-2 != length(y)) {return(F)}
#       comparison <- y==x
#       if(!comparison[1]) {return(F)}
#       return(
#         sum(!comparison)==2 & (which(!comparison)[2]==(which(!comparison)[1] -1)) & grepl("^[Nn]a$",paste(y[!comparison],collapse=""))
#       )
#     }
#     , y=sepChar)
# }


# searched=speciesToSearch[i]
# tabGbif = resGbifSpe[[i]]
# expected=c(family=familySpecies[i])
# rank="species"
# obligatory = c(kingdom="Plantae")
# returnGbifRes=T
# orderedRanks=c("domain","superkingdom","kingdom","subkingdom","infrakingdom","superphylum","phylum","division","subphylum","subdivision","infradivision","superclass","class","subclass","infraclass","subterclass","parvclass","megacohort","supercohort","cohort","subcohort","infracohort","superorder","order","suborder","infraorder","parvorder","superfamily","family","subfamily","supertribe","tribe","subtribe","genus","subgenus","section","subsection","species group","series","species subgroup","species","infraspecies","subspecies","forma specialis","variety","varietas","subvariety","race","stirp","form","forma","morph","subform")



analyseGbifTable <- function(searched, tabGbif,rank, obligatory, expected, returnGbifRes=T,
                             orderedRanks=c("domain","superkingdom","kingdom","subkingdom","infrakingdom","superphylum","phylum","division","subphylum","subdivision","infradivision","superclass","class","subclass","infraclass","subterclass","parvclass","megacohort","supercohort","cohort","subcohort","infracohort","superorder","order","suborder","infraorder","parvorder","superfamily","family","subfamily","supertribe","tribe","subtribe","genus","subgenus","section","subsection","species group","series","species subgroup","species","infraspecies","subspecies","forma specialis","variety","varietas","subvariety","race","stirp","form","forma","morph","subform")
                             )
{
  if(sum(sapply(obligatory,length)))
  {
    stopifnot(names(obligatory) %in% names(tabGbif) & is(obligatory,"character"))
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
  if(sum(sapply(expected,length)))
  {
    stopifnot(names(expected) %in% names(tabGbif) & is(expected,"character"))
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

getSuggestsFromResGbif <- function(analysedGbif,ranks=c("family","genus","species","subspecies","variety"))
{
  orderedRanks=c("domain","superkingdom","kingdom","subkingdom","infrakingdom","superphylum","phylum","division","subphylum","subdivision","infradivision","superclass","class","subclass","infraclass","subterclass","parvclass","megacohort","supercohort","cohort","subcohort","infracohort","superorder","order","suborder","infraorder","parvorder","superfamily","family","subfamily","supertribe","tribe","subtribe","genus","subgenus","section","subsection","species group","series","species subgroup","species","infraspecies","subspecies","forma specialis","variety","varietas","subvariety","race","stirp","form","forma","morph","subform")
  stopifnot(ranks %in% orderedRanks)
  if(!is.na(analysedGbif$finalRank) && !analysedGbif$finalRank %in% ranks)
  {warning("The rank of the taxon found in the GBIF backbone (",analysedGbif$finalRank,") is not in the extracted ranks")}
  res<-sapply(ranks,function(analysedGbif)NA,USE.NAMES = T)
  if(analysedGbif$finalRank %in% names(res)) {res[analysedGbif$finalRank]<-analysedGbif$canonicalname}
  inHR<-match(ranks,analysedGbif$higherRanks$rank)
  res[!is.na(inHR)]<-analysedGbif$higherRanks$canonicalname[na.omit(inHR)]
  return(data.frame(type=analysedGbif$type,as.data.frame(as.list(res)),gbifid=analysedGbif$gbifid))
}

getSuggestsFromListGbif <- function(listAnalysedGbif,ranks=c("family","genus","species","subspecies","variety"), exclude=c("exactMatch"))
{
  suggested<-Reduce(rbind,lapply(listAnalysedGbif,getSuggestsFromResGbif,ranks=ranks))
  rownames(suggested)<-names(listAnalysedGbif)
  suggested<-suggested[!suggested$type%in%exclude,]
  return(suggested)
}

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

# taxo<- taxoTDF
# messagesGbif=F
# messageOther=T
# kingdom="Plantae"

analyseTaxoRanksInColAttrib <- function(taxo)
{
  isTaxo<-is(taxo_oneTab)
  orderedRanks <- c("domain","superkingdom","kingdom","subkingdom","infrakingdom","superphylum","phylum","division","subphylum","subdivision","infradivision","superclass","class","subclass","infraclass","subterclass","parvclass","megacohort","supercohort","cohort","subcohort","infracohort","superorder","order","suborder","infraorder","parvorder","superfamily","family","subfamily","supertribe","tribe","subtribe","genus","subgenus","section","subsection","species group","series","species subgroup","species","infraspecies","subspecies","forma specialis","variety","varietas","subvariety","race","stirp","form","forma","morph","subform")
  epithetized <- which(orderedRanks=="species"):length(orderedRanks)
}

taxoRanks <- function(taxo)
{
  orderedRanks <- c("domain","superkingdom","kingdom","subkingdom","infrakingdom","superphylum","phylum","division","subphylum","subdivision","infradivision","superclass","class","subclass","infraclass","subterclass","parvclass","megacohort","supercohort","cohort","subcohort","infracohort","superorder","order","suborder","infraorder","parvorder","superfamily","family","subfamily","supertribe","tribe","subtribe","genus","subgenus","section","subsection","species group","series","species subgroup","species","infraspecies","subspecies","forma specialis","variety","varietas","subvariety","race","stirp","form","forma","morph","subform")
  epithetized <- which(orderedRanks=="species"):length(orderedRanks)
  cnRanks <- cnRanks_epithetized <- orderedRanks
  cnRanks[orderedRanks=="class"]<-cnRanks_epithetized[orderedRanks=="class"]<-"taxonomic_class"
  cnRanks[orderedRanks=="order"]<-cnRanks_epithetized[orderedRanks=="order"]<-"taxonomic_order"
  cnRanks_epithetized[epithetized] <- paste(cnRanks_epithetized[epithetized],"epithet",sep="_")
  notNAattributes<-names(attributes(taxo))[!sapply(attributes(taxo),function(x)all(is.na(x)))]
  attributed_ranks <-which(cnRanks %in% notNAattributes | cnRanks_epithetized %in% notNAattributes)
  stopifnot(order(attributed_ranks)==1:length(attributed_ranks))
  attributes_ranks<-mapply(function(x,y)ifelse(is.na(x),y,x),match(cnRanks[attributed_ranks],names(attributes(taxo))),match(cnRanks_epithetized[attributed_ranks],names(attributes(taxo))))
  numcol_attributed_ranks<-match(unlist(attributes(taxo)[attributes_ranks]),colnames(taxo))
  minAttributedRank<-orderedRanks[attributed_ranks][1]
  maxAttributedRank<-orderedRanks[attributed_ranks][length(orderedRanks[attributed_ranks])]
  unattributed_ranks <-which(cnRanks %in% colnames(taxo) | cnRanks_epithetized %in% colnames(taxo) | orderedRanks %in% colnames(orderedRanks))
  unattributed_ranks <- unattributed_ranks[!unattributed_ranks %in% attributed_ranks]
  if(length(unattributed_ranks))
  {
    warning("We detected that you might have taxonomic information which may not be accounted for in the taxonomic object, ranks: ", paste0(orderedRanks[unattributed_ranks], collapse=" "))
  }
  lowerInfo<-apply(taxo[numcol_attributed_ranks],1,function(x,n)
  {
    if(all(is.na(x))){return(1)}
    max(which(!is.na(x)))+1
  })
  higherMissing<-apply(taxo[numcol_attributed_ranks],1,function(x,n)
    {
      if(!any(is.na(x))){return(n+1)}
      return(min(which(is.na(x))))
    },n=length(numcol_attributed_ranks))
  if(any(lowerInfo!=higherMissing))
  {warning("Some of the taxa have missing higher information (see rows ",paste0(which(lowerInfo!=higherMissing),collapse=" "))}
  LEV<-c("higher",orderedRanks[attributed_ranks])
  return(factor(LEV[lowerInfo],levels=LEV, ordered = T))
}

searchGbif<-function(taxo,messagesGbif=F,messageOther=T,kingdom="Plantae",excludeFromSuggests=c("exactMatch","Failed"))
{
  stopifnot(is(taxo,"taxo_oneTab"))
  mat<-as.matrix(taxo[na.omit(unlist(attributes(taxo)[c("family","genus","species_epithet","infraspecies_epithet","cf_aff","sp_specif")]))])
  m_mat<-match(split(mat,row(mat)),split(mat,row(mat)))
  corresCode<-tapply(m_mat,taxo[[attr(taxo,"taxoCode")]],unique,simplify = F)
  codePb<-names(corresCode)[sapply(corresCode,length)>1]
  if(length(codePb)>0)
  {
    warning("The taxa code do not always correspond to the same taxonomic information, we urge you to use the `checkUnicityCodetax` function to address these problems before querying GBIF\nProblematic codes:",paste(codePb,collapse=" "))
  }
  corresGenusFamily<-tapply(taxo[[attr(taxo,"family")]],taxo[[attr(taxo,"genus")]],unique)
  genusFamPb<-names(corresGenusFamily[sapply(corresGenusFamily,length)>1])
  if(length(genusFamPb)>0)
  {
    warning("The genera do not always have the same family, we urge you to use the `checkUnicityGnInFam` function to address these problems before querying GBIF\nProblematic genera:",paste(genusFamPb,collapse=" "))
  }
  ranks<-taxoRanks(taxo)
  if(any(ranks=="infraspecies_epithet")){
    tabInfraspecies<-unique(taxo[ranks=="infraspecies_epithet",unlist(attributes(taxo)[c("family","genus","species_epithet","infraspecies_epithet")])])
    infraspeciesToSearch<-apply(tabInfraspecies[2:4],1,paste)
    familiesInfraspecies<-tabInfraspecies[[1]]
    #todo:managing the search and error management of infraspecies
  }
  if(any(ranks=="species")){
    tabSpecies<-unique(taxo[ranks=="species",unlist(attributes(taxo)[c("family","genus","species_epithet")])])
    matSpecies<-as.matrix(tabSpecies)
    mSpe<-match(split(mat[,colnames(matSpecies)],row(mat[,colnames(matSpecies)])),split(matSpecies,row(matSpecies)))
    mSpe[ranks!="species"]<-NA
    speciesToSearch<-paste(tabSpecies[,2],tabSpecies[,3])
    familySpecies<-tabSpecies[,1]
    if(messageOther)
    {message("Searching for ",length(speciesToSearch)," species in the GBIF Backbone\n...")}
    resGbifSpe<-taxize::get_gbifid_(speciesToSearch,messages=messagesGbif)
    if(messageOther)
    {message("done\nAnalysing GBIF Backbone information")}
    gbifAnalysed<-mapply(analyseGbifTable,searched=speciesToSearch, tabGbif=resGbifSpe,expected= lapply(as.character(familySpecies),function(x)c(family=x)), MoreArgs=list(rank="species", obligatory=c(kingdom=kingdom)),SIMPLIFY = F)
    diagnostic<-table(sapply(gbifAnalysed,function(x)x$type))
    if(messageOther)
    {
      message(diagnostic["exactMatch"]," species are found without any modification needed")
      message(sum(diagnostic[grep("fuzzy",names(diagnostic))])," species are found with suggested orthographic changes")
      message(sum(diagnostic[grep("synonym",names(diagnostic))])," species are suggested synonyms")
      message(sum(diagnostic[grep("changeHigherRanks",names(diagnostic))])," species are found with suggested higher rank changes")
      message(sum(diagnostic[grep("Failed",names(diagnostic))])," species were not found")
    }
    suggestedSpecies<-getSuggestsFromListGbif(gbifAnalysed,ranks=c("family","genus","species"),exclude=NULL)
    suggestedSpecies[attr(taxo,"species_epithet")]<-extractEpithet(suggestedSpecies$species,suggestedSpecies$genus)
    concernedRowsTaxo<-which(!is.na(mSpe) & !suggestedSpecies$type[mSpe] %in% excludeFromSuggests)
    # TODO: working on the warnings depending on excludeFromSuggest
    # Note: here to find the correspondance between suggestedSpecies and the taxo row concerned you may want to use: data.frame(row_taxo=concernedRowsTaxo, row_suggestedSpecies=mSpe[concernedRowsTaxo])
    data.frame(row=concernedRowsTaxo, )
  }

  familyToSearch<-na.omit(unique(taxo[[attr(taxo,"family")]]))
  resGbifFam<-taxize::get_gbifid_(familyToSearch,messages=messagesGbif)
  perfectMatchOK<-lapply(resGbifFam,function(x){NA})

}



genusSearchGbif <- function(obj)
{NA}

spSearchGbif <- function(obj)
{NA}
