#function to annotate if pass variant
annotate_pass_rf_filter <- function(model, data, rescue=F){
  data %>% 
    mutate(gdc_mask=!(pon=='')) %>%
    mutate(is_indel=!(nchar(ref)==nchar(alt))) %>%
    mutate(in_dbsnp=!is.na(dbsnp)) %>%
    mutate(is_rmsk=!(rmsk_hit=='')) %>%
    replace_na(list(oxog=0,mqrs=0,mqs=0,sor=0,str_ru='')) -> rf_ready
  #
  p <- predict(model, rf_ready)
  p[is.na(p)] <- FALSE 
  p <- as.logical(p) & (rf_ready$impact %in% c('HIGH','MODERATE'))
  if(rescue){
    p <- p | with(rf_ready,
      (
       tlod>15 & nlod>20 & flod>16 & sor<3 &
       xfet<1e-4 & adt_fwd>=5 & adt_rev>=5 & adn<=2 &
       impact %in% c('HIGH','MODERATE')
       )
    )
  }
  stopifnot(length(p)==nrow(data))
  return(as.logical(p))
}

#' Import grouping from vault metadata
#'
#' @param vault a TPO vailt object
#'
#' @return id the "id" from the groupings table
#'
#' @export
import_meta <- function(vault){
  db <- parse_db_creds()
  con <- dbConnect(PostgreSQL(),
    host=db['host'],
    port=db['port'],
    dbname = db['db'],
    user=db['user'],
    password=db['password']
  )
  #get the existing groupings from the database
  dbgrps <- tbl(con, 'groupings') %>% select(-c(time,id)) %>% collect() %>% as.data.frame()

  grp <- with(vault$meta,data.frame(
                patient=case,
                cohort=ifelse(is.na(cohort),'Default',cohort),
                alignid_t=ifelse(is.null(tumor.dna.sample),'',tumor.dna.sample),
                alignid_n=ifelse(is.null(normal.dna.sample),'',normal.dna.sample),
                alignid_tr=ifelse(is.null(tumor.rna.sample),'',tumor.rna.sample),
                alignid_nr=ifelse(is.null(normal.rna.sample),'',normal.rna.sample),
                stringsAsFactors=F
  ))
  #remove the patient id from the RNA name if present
  if(grepl(grp$patient,grp$alignid_tr)){
    #grp$alignid_tr <- gsub(sprintf('%s[-_//.]',grp$patient),'',grp$alignid_tr)
  }
  put <- dbWriteTable(con,'groupings',grp,append=T,row.names=F)

  
  #get the new id from the db
  if(!is.na(grp$alignid_t)){
    gid <- tbl(con, 'groupings') %>% filter(
        alignid_t==!!grp$alignid_t & 
        (alignid_n==!!grp$alignid_n | is.na(alignid_n)) & 
        (alignid_tr==!!grp$alignid_tr | is.na(alignid_nr)) &
        (alignid_nr==!!grp$alignid_nr | is.na(alignid_nr))
        ) %>% 
      select(id) %>% collect() %>% .$id
  }else{
    gid <- tbl(con, 'groupings') %>% filter(
        (alignid_n==!!grp$alignid_n | is.na(alignid_n)) & 
        alignid_tr==!!grp$alignid_tr  &
        (alignid_nr==!!grp$alignid_nr | is.na(alignid_nr))
        ) %>% 
      select(id) %>% collect() %>% .$id
  }
  stopifnot(class(gid)=='integer' & length(gid)==1)

  meta <- with(vault$meta,data.frame(
                id=gid,
                assembly=ifelse(exists('assembly'),assembly,''),
                dna_cap=ifelse(exists('anno.targets'),basename(anno.targets),''),
                sex=sex,
                purity=purity,
                ploidy=ploidy
          ))
  dbWriteTable(con,'meta',meta,append=T,row.names=F)

  #clean up
  dbDisconnect(con)

  #return the ID from the groupings table
  return(gid)
}



#' Import variants from a TPO vault
#'
#' @param vault a TPO vault object
#' @param tp type (somatic or germline)
#' @param gid The id from the groupings table for the analysis
#'
#' @export
import_variants <- function(vault,tp,gid){
  db <- parse_db_creds()
  con <- dbConnect(PostgreSQL(),
    host=db['host'],
    port=db['port'],
    dbname = db['db'],
    user=db['user'],
    password=db['password']
  )

  dbcols <- colnames(tbl(con,tp))
  dbcols <- dbcols[dbcols!='iddb']

  #load the variants
  if(tp=='somatic'){
    vars <- vault$tables$somatic %>% select(-id)
  }else if(tp=='germline'){
    vars <- vault$tables$germline %>% select(-id)
  }
  colnames(vars) <- tolower(colnames(vars))
  colnames(vars) <- gsub('\\.','_',colnames(vars))

  vars <- filter(vars, (adt>2 | adn>2)) # clean some odd calls out

  # add the runid
  vars$runid <- gid

  #unlist rmsk and pon
  foo <- sapply(vars$rmsk_hit, paste, collapse=',')
  vars$rmsk_hit <- foo
  foo <- sapply(vars$pon, paste, collapse=',')
  vars$pon <- foo

  #set common TERT mutations to high impact
  vars$impact[vars$gene_name=='TERT' & vars$cosmic_cnt>500] <- 'HIGH'


  #sanity check columns
  #stopifnot(all(colnames(ds) %in% dbcols))
  w <- dbcols[!dbcols %in% colnames(vars)]
  if(length(w)>0){
    for(a in w){
      vars[,a] <- NA
    }
    warning(sprintf("Missing columns in vault, using NA:%s",paste0(w,collapse=',')))
  }
  if(!all(colnames(vars) %in% dbcols)){
    warning(sprintf("New columns in vault, not in DB:%s",paste0(setdiff(colnames(vars), dbcols),collapse=',')))
    vars <- vars %>% select(dbcols[-1])
  }


  #if(tp=='somatic'){
  ##TMB/MSI
  #  print("writing TMB")
  #  #tmb_to_add <- mioncoseq2::tmb(union, cap)
  #  tmb_to_add <- c(0,0)

  #  tmb_df <- data.frame(
  #    runid=gid,
  #    mutations=as.integer(tmb_to_add[1]),
  #    tmb=tmb_to_add[2],
  #    mask=F,
  #    msi_score=vault$meta$derived$msi_score,
  #    msi_call=vault$meta$derived$msi_suspected>.25 #change this
  #  )

  #  dbWriteTable(con,'tmb',tmb_df,append=T,row.names=F)
  #}

  # insert
  vars <- as.data.frame(vars)
  print(sprintf(
    "writing %d variants to database %s(%s) %s",
    nrow(vars),db[1],db[3],tp
  ))
  dbWriteTable(con,tp,vars,append=T,row.names=F)

  #clean up
  dbDisconnect(con)
}


#' Import structural variants from a TPO vault
#'
#' @param vault a TPO vault object
#' @param gid The id from the groupings table for the analysis
#'
#' @export
import_structural <- function(vault,gid){
  db <- parse_db_creds()
  con <- dbConnect(PostgreSQL(),
    host=db['host'],
    port=db['port'],
    dbname = db['db'],
    user=db['user'],
    password=db['password']
  )

  dbcols <- colnames(tbl(con,'structural'))
  struc <- vault$tables$structural %>% select(-id) %>%
          mutate(runid=gid) %>% as.data.frame()

  colnames(struc) <- gsub("\\.","_",tolower(colnames(struc)))

  #sanity check
  stopifnot(all(colnames(dbcols %in% colnames(struc))))

  dbWriteTable(con, 'structural', struc, append=T, row.names=F)
  #clean up
  dbDisconnect(con)

}

#' Import CNV segments from vault to database
#'
#' @param vault a TPO vault object
#' @param gid the grouping id in the database
#' @param type somatic ('som') or germline ('ger') analysis 
#'
#' @export
import_cnv <- function(vault,gid,type){

  db <- parse_db_creds()
  con <- dbConnect(PostgreSQL(),
    host=db['host'],
    port=db['port'],
    dbname = db['db'],
    user=db['user'],
    password=db['password']
  )

  dbcols <- colnames(tbl(con,'cnv'))

  #merge the genes into the segments
  segs <- vault$maps$segment %>% group_by(var_id) %>% 
    summarize(genes=sprintf("{%s}",paste(gene_id, collapse=","))) %>%
    ungroup() %>%
    merge(vault$tables$segment, by='var_id', all=T) %>%
    select(seg_start=start,seg_end=end,everything()) %>%
    mutate(runid=!!gid, tp=!!type)

  stopifnot(nrow(segs)==nrow(vault$tables$segment))
  colnames(segs) <- tolower(colnames(segs))
  segs <- segs[,dbcols]
  

  ft <- list("runid"="bigint",'tp'='text',"seg"="int","chr"="text",
            "seg_start"="int","seg_end"="int","c"="int","k"="int",
            "sc"="real","lr"="real","tl"="real","al"="real","d"="real",
            "anom"="real","mse"="real","nlr"="real","naf"="real","genes"="text[]"
            )
  dbWriteTable(con, 'cnv', segs, append=T, row.names=F, field.types=ft)

  #
  dbDisconnect(con)
}


#' Import gene expression from vault to database
#'
#' @param vault a TPO vault object
#' @param gid the grouping id in the database
#' @param type 'tumor' or 'normal'
#'
#' @export
import_gxp <- function(vault,gid,type){
  db <- parse_db_creds()
  con <- dbConnect(PostgreSQL(),
    host=db['host'],
    port=db['port'],
    dbname = db['db'],
    user=db['user'],
    password=db['password']
  )

  dbgenes <- tbl(con,'gxp_genes') %>% select(gene_id) %>% collect()
  alignid <- tbl(con,'groupings') %>% filter(id==!!gid) %>% collect() %>% .$alignid_tr
  stopifnot(length(alignid)==1)

  gxp <- vault$tables$gene.expression %>% 
            mutate(runid=gid,tp=type,alignid=alignid) %>%
            arrange(gene_id)

  #sanity check
  stopifnot(all(gxp$gene_id==dbgenes$gene_id))

  #
  gxp %>% group_by(tp,runid,alignid)  %>%
            tidyr::replace_na(list('tcount'=-1,'tcpm'=-1,'trpkm'=-1)) %>%
            summarize(
            count=sprintf("{%s}",paste(tcount, collapse=",")),
            cpm=sprintf("{%s}",paste(round(tcpm,2), collapse=",")),
            rpkm=sprintf("{%s}",paste(round(trpkm,2), collapse=",")),
            .groups='keep') %>%
            ungroup() -> to_add

  # TODO: leave comment on how it works
  ft <- list('id'='integer','tp'='text','alignid'='text','count'='integer[]',
            'cpm'='real[]','rpkm'='real[]')
  dbWriteTable(con, 'gxp', to_add, append=T, row.names=F, field.types=ft)

  #refresh the view
  dbExecute(con,'refresh materialized view gxp_pctile')
  #clean up
  dbDisconnect(con)
}

#' Import gene fusions from vault to database
#'
#' @param vault a TPO vault object
#' @param gid the grouping id in the database
#' @param type 'tumor' or 'normal'
#'
#' @export
import_fus <- function(vault,gid,type){
  db <- parse_db_creds()
  con <- dbConnect(PostgreSQL(),
    host=db['host'],
    port=db['port'],
    dbname = db['db'],
    user=db['user'],
    password=db['password']
  )

  dbcols <- colnames(tbl(con,'fusions'))
  fus <- vault$tables$fusion %>% 
          select(-c(sid,id)) %>% 
          mutate(runid=gid) %>%
          as.data.frame()

  # TODO: check if this is necessary
  fus$ctg.seq=unlist(sapply(fus$ctg.seq,
    function(x){
      if(length(x)>0){as.character(x)[[1]]}else{NA}
    }
  ))
  colnames(fus) <- gsub("\\.","_",tolower(colnames(fus)))

  #sanity check
  stopifnot(all(colnames(dbcols %in% colnames(fus))))

  dbWriteTable(con, 'fusions', fus, append=T, row.names=F)
  #clean up
  dbDisconnect(con)
}
