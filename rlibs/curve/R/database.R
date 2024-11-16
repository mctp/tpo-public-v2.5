#' Copy new groupings from a sheet to the database
#'
#' @export
db_put_groupings <- function(sheet,db=Sys.getenv('CURVE_DB')){
  if(class(db)=="character"){
    db <- strsplit(db,':')[[1]]
  }
  con <- dbConnect(PostgreSQL(),
    host=db[1],
    port=db[2],
    dbname = db[3],
    user=db[4],
    password=db[5]
  )

  #get the existing groupings from the database
  dbgrps <- tbl(con, 'groupings') %>% select(-c(time,id)) %>% collect() %>% as.data.frame()

  #load the groupings sheet
  sheet <- read.csv(sheet,stringsAsFactors=F) %>%
              mutate(varid=paste(patient,alignid_t,alignid_n,sep='.'))

  if(nrow(dbgrps)==0){
    #initialize
    dbWriteTable(con,'groupings',sheet,append=T,row.names=F)
  }else{
    #append
    to_do <- setdiff(sheet,dbgrps)
    if(nrow(to_do)>0){
        dbWriteTable(con,'groupings',to_do,append=T,row.names=F)
    }else{
    #nothing to do
        print('no new groupings found')
    }
  }
       
  #clean up
  dbDisconnect(con)
}

#' Parse DB credentials from an environment variable
#'
#' @export
parse_db_creds <- function(string=NULL){
  if(is.null(string)){
    string <- Sys.getenv('CURVE_DB')
  }
  db <- strsplit(string,":")[[1]]
  stopifnot(length(db)==5)
  stopifnot(!any(db==""))
  names(db) <- c('host','port','db','user','password')
  return(db)
}

#' Initialize the gxp_genes table using a GTF file
#' 
#' @export
init_gxp_genes <- function(gtf='/tpo/refs/grch38/ensembl/grch38.97.clean.gtf'){
  suppressPackageStartupMessages(require(dplyr))
  db <- parse_db_creds()
  con <- dbConnect(PostgreSQL(),
    host=db['host'],
    port=db['port'],
    dbname = db['db'],
    user=db['user'],
    password=db['password']
  )
  existing <- tbl(con, 'gxp_genes') %>% tally() %>% collect() %>% .$n
  if(existing==0){
    print(sprintf("Initializing gxp_genes table using %s", gtf))
    genes_gr <- robustImport(gtf,getGobj("hg38",NULL,F)$seqi,feature.type='gene')
    to_insert <- data.frame(gene_id=genes_gr$gene_id,gene_name=genes_gr$gene_name) %>% arrange(gene_id)
    to_insert$idx <- seq(1:nrow(to_insert))
    dbWriteTable(con, 'gxp_genes', to_insert,append=T, row.names=F)
  }
}


#### DB Retrieval Functions ######

#' Get the Patient metadata
#' 
#' @export
db_get_meta <- function(runid){
    meta <- suppressWarnings(
              tbl(db_pool, 'meta') %>% 
                filter(id %in% !!runid) %>% 
                collect() 
              )
    if(!is.na(meta$user_meta)){
      user_meta <- as.data.frame(fromJSON(meta$user_meta),stringsAsFactors=F)
      meta <- cbind(meta,user_meta)
    }
    meta <- meta %>%
              select(-user_meta) %>%
              mutate(purity=sprintf('%.0f%%',purity*100)) %>%
              mutate(ploidy=round(ploidy,1)) %>%
              t() %>% as.data.frame() %>%
              tibble::rownames_to_column('V0') %>%
              filter(!V0 %in% c('time','id','dna_cap','assembly')) %>%
              mutate(V0=tools::toTitleCase(V0))

  return(meta)   
}

#' Get the variants with recur,gxp and CN from the database
#' 
#' @export
db_get_variants <- function(runid,tp){
    fun <- ifelse(tp=='somatic','get_int_vars','get_int_ger')
    q <- sprintf("select * from %s(%d)",fun,as.integer(runid))
    res <- pool::dbGetQuery(db_pool,q) %>%
                mutate(rpkm=round(rpkm,2)) %>%
                arrange(desc(abs(cn-2))) %>% # keep the dup with the most interesting CN
                distinct(iddb, .keep_all=T)
    res$recur[res$recur==0] <- 1 # everything has a recur of 1
    return(res)
}


#' Get the gene expression with copynumber from the database
#' 
#' @export
db_get_gxp <- function(runid){
    data <- dbGetQuery(db_pool,
      sprintf("select gene_id,gene_name,rpkm,c from get_gxp(%s)",runid)
    )  %>%
        arrange(desc(abs(c-2))) %>% # keep the dup with the most interesting CN
        distinct(gene_id, .keep_all=T) %>%
        arrange(desc(rpkm))
    return(data)
}

#' Get the CNV segments from the database
#' 
#' @export
db_get_cnv <- function(runid){
    d <- tbl(db_pool, 'cnv') %>% 
      filter(runid %in% !!runid) %>% 
      collect()
    gl <- sapply(gsub('\\{|\\}|"','',d$genes),strsplit,split=',')
    names(gl) <- d$seg
    d$genes <- gl
    #
    return(d)
}

#' Get the SV from the database
#' 
#' @export
db_get_sv <- function(runid){
    sv <- tbl(db_pool, 'structural') %>% 
      filter(runid %in% !!runid) %>% 
      collect() %>%
      arrange(desc(aft), desc(qual))
    #
    return(sv)
}

#' Get the gene fusions from the database
#' 
#' @export
db_get_fus <- function(runid){
    data <- dbGetQuery(db_pool,
      sprintf("select * from get_fus(%s)",runid)
    ) %>%
        arrange(desc(sum_jnc+sum_bpt))
    return(data)
}

#' Get the list of what data is sideloaded for a patient
#' 
#' @export
db_get_sld_list <- function(runid){
    data <- dbGetQuery(db_pool,
      sprintf("select name from sideload_data where id in (%s)",runid)
    )
    return(data$name)
}

#' Get a sideloaded dataset for a patient
#' 
#' @export
db_get_sld <- function(runid,nm){
    q <- sprintf("select * from get_sideload_data('%s',%s);",nm,runid)
    data <- dbGetQuery(db_pool,q) %>% arrange(desc(value))
    data$value[data$value<0] <- NA
    colnames(data) <- c('id','gene',nm)
    return(data)
}

#' Get the distribution of a sideloaded feature
#' 
#' @export
db_get_sld_dist <- function(id,nm){
    idx <- dbGetQuery(db_pool,
      sprintf("select idx from(select unnest(id) as id, generate_subscripts(id,1) as idx from sideload_meta where name='%s') foo where id='%s'",nm,id)
    )$idx
    if(length(id)==1){
      dist <- dbGetQuery(db_pool,
        sprintf("select groupings.id,patient,cohort,value[%d] from sideload_data  left join groupings on groupings.id= sideload_data.id where name='Tumor Proteomics';",idx)
      )
    }else{
      dist <- data.frame()
    }
    return(dist)
}
