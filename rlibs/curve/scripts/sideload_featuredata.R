
option_list = list(
    optparse::make_option(c("-f", "--featuredata"), type="character",
                          help="csv containing feature data to add to database"),
    optparse::make_option(c("-g", "--gtf"), type="character",
                          help="path to gtf gene model"),
    optparse::make_option(c("-n", "--name"), type="character",
                          help="name of the feature for the database"),
    optparse::make_option(c("-d", "--database"), type="character",
                          default=NULL,
                          help="database credentials (will use CURVE_DB if missing)")
)
parser = optparse::OptionParser(
  "Rscript sideload_featuredata.R [options]",
  description=c("Import feature-level data (i.e. proteomics) to a CURVE database. The csv must contain a column called gene_id which maps to the gtf features. All other column names must map to the database alignid_t column (the tumor DNA library)."),
  epilogue=c(
      "Michigan Center for Translational Pathology (c) 2022\n"),
  option_list=option_list,
  add_help_option=T
  )

args <- optparse::parse_args(parser, positional_arguments=FALSE)

library(dplyr)
library(dbplyr)
library(jsonlite)
library(curve)
library(DBI)

# get DB creds (if needed)
if(!is.null(args$database)){
  db <- parse_db_creds(args$database)
}else{
  db <- parse_db_creds()
}
con <- DBI::dbConnect(RPostgreSQL::PostgreSQL(),
  host=db['host'],
  port=db['port'],
  dbname = db['db'],
  user=db['user'],
  password=db['password']
)

# Get the existing groupings and sideload metadata
grps <- tbl(con, 'groupings') %>% collect()
slm <- tbl(con, 'sideload_meta') %>% filter(name %in% !!args$name) %>% select(name) %>% collect()



#to_load <- read.csv('/data/cptac/cptac_pt_metadata.csv', stringsAsFactors=F) %>% filter(case %in% grps$patient)
to_load <- read.csv(args$featuredata)  %>% arrange(gene_id)

rnm <- rep(-1,ncol(to_load))[-1]
names(rnm) <- colnames(to_load)[-1]

w <- which(colnames(to_load) %in% c('gene_id', grps$alignid_t))
to_load <- to_load[,w] %>% tidyr::replace_na(as.list(rnm))


genes <- cnvex::robustImport(args$gtf,cnvex::getGobj("hg38",NULL,F)$seqi,feature.type='gene')
if(any(!to_load$gene_id %in% genes$gene_id)){
  w <- which(!to_load$gene_id %in% genes$gene_id)
  warning(sprintf('removing %d features not in the gtf',length(w)))
  to_load <- to_load[-w,]
}


#check that there is some data to add
stopifnot('gene_id' %in% colnames(to_load))
stopifnot(nrow(to_load)>0)
stopifnot(ncol(to_load)>1)

#if metadata has not been set up, add it
if(nrow(slm)==0){
  genes_to_add <- genes %>% 
    as.data.frame() %>%
    select(gene_name,gene_id) %>% 
    filter(gene_id %in% to_load$gene_id) %>% 
    arrange(gene_id)

  to_write <- data.frame(
    name=args$name,
    id=sprintf('{%s}',paste(genes_to_add$gene_id,collapse=',')),
    gene=sprintf('{%s}',paste(genes_to_add$gene_name,collapse=','))
  ) 
  ft <- list('name'='text','id'='text[]','gene'='text[]')
  DBI::dbWriteTable(con, 'sideload_meta', to_write, append=T, row.names=F, field.types=ft)
}

#data should be in the same order as the database
dbgenes <- dbGetQuery(con,sprintf("select unnest(id) as gene_id from sideload_meta where name='%s';",args$name))
stopifnot(all(to_load$gene_id==dbgenes$gene_id))

for(a in 2:ncol(to_load)){
  id <- grps %>% filter(alignid_t %in% !!colnames(to_load)[a]) %>% .$id
  to_add <- data.frame(
    name=args$name, 
    id=id,
    value=sprintf('{%s}',paste(round(to_load[,a],2),collapse=','))
  )
  existing <- tbl(con, 'sideload_data') %>% filter(id==!!id & name==!!args$name) %>% collect() %>% nrow()
  if(existing==0){
    print(sprintf('adding %s record for %s',args$name,colnames(to_load)[a]))
    ft <- list('name'='text','id'='text','value'='numeric[]')
    DBI::dbWriteTable(con, 'sideload_data', to_add, append=T, row.names=F, field.types=ft)
  }else{print(sprintf('skipping %s do to existing record',colnames(to_load)[a]))}
}


#clean up 
DBI::dbDisconnect(con)


