
option_list = list(
    optparse::make_option(c("-m", "--metadata"), type="character",
                          help="csv containing metadata to add to database"),
    optparse::make_option(c("-d", "--database"), type="character",
                          default=NULL,
                          help="database credentials (will use CURVE_DB if missing)")
)
parser = optparse::OptionParser(
  "Rscript sideload_metadata.R [options]",
  description=c("Import the metadata in a csv to the Curve database. The csv must contain a column called 'case' to match to the database 'patient' field. If a column 'alignid_t' exists in the tble, it will be used to match the tumor library in the database as well as case (if metadata is sample based rather than case based). Additional columns will be converted to JSON and stored in the meta table of the database for display. Missing or NA values are removed before storage. Current values will be overwritten in the database if they exist."),
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

# Get the existing groupings and metadata
grps <- tbl(con, 'groupings') %>% collect()
meta <- tbl(con, 'meta') %>% collect()


#to_load <- read.csv('/data/cptac/cptac_pt_metadata.csv', stringsAsFactors=F) %>% filter(case %in% grps$patient)
to_load <- read.csv(args$metadata) %>% filter(case %in% grps$patient) %>% select(-X)
if('alignid_t' %in% colnames(to_load)){
  to_load <- to_load %>% filter(alignid_t %in% grps$alignid_t)
}

stopifnot(nrow(to_load)>0)

for(a in 1:nrow(to_load)){
  #
  print(sprintf('Attempting to update metadata for %s',to_load$case[a]))
  #clean up the metadata and remove missing
  tmp <- as.list(to_load[a,])
  w <- which(tmp %in% c('NA','','case') | is.na(tmp) | is.null(tmp))
  if(length(w)>0){tmp <- tmp[-w]}
  names(tmp) <- gsub('\\.',' ',names(tmp))
  names(tmp) <- gsub('  ',' ',names(tmp))
  #make the JSON to write
  str <- toJSON(tmp,auto_unbox=T)

  #assemble the query
  if('alignid_t' %in% names(tmp)){
    id <- grps %>% filter(patient==to_load$case[a] & alignid_t==to_load$alignid_t[a]) %>% .$id
  }else{
    id <- grps %>% filter(patient==to_load$case[a]) %>% .$id
  }
  if(length(id)==0){
    print(sprintf('%s was not found in the database, skipping',to_load$case[a]))
  }else{
    q <- sprintf("update meta set user_meta='%s' where id=%s",str,id[1])
    res <- DBI::dbExecute(con,q)
    if(res!=length(id)){
      print('# of updated rows does not match expectation')
      print(str)
    }
    if(res>0){
      print(sprintf('updated %d rows for %s (id:%s):',res,to_load$case[a],paste0(id,collapse=',')))
      print(
        toJSON(fromJSON(str),auto_unbox=T,pretty=T)
      )
    }
  }

}

#clean up 
DBI::dbDisconnect(con)


