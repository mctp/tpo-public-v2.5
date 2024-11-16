library(shiny)
library(shinyBS)
library(dplyr)
library(curve)

color1 <- "#989C97" # umich ash
color2 <- "#00274C" # umich blue
color3 <- "#575294" # umich violet

## get hostname
hn <- system("hostname", intern=TRUE)

## get params from CURVE_DB environment variable
db <- parse_db_creds()
db_pool <- pool::dbPool(
  drv = RPostgreSQL::PostgreSQL(),
  host=db["host"],
  port=db["port"],
  dbname = db["db"],
  user=db["user"],
  password=db["password"]
)

## load global variables
gxp_ntile <- tbl(db_pool, "gxp_pctile") %>% select(-fooidx) %>% collect()
gxp_to_id <- tbl(db_pool, "gxp_genes") %>% select(-idx) %>% collect()
n_recur <- tbl(db_pool, "somatic_recur") %>% select(total) %>% head(1) %>% collect() %>% .$total
groupings <- tbl(db_pool, "groupings") %>% collect()

## define genes of interest
source(system.file("extdata","goi.R", package="curve"))

## data model mapping colums to display
dict <- read.csv(system.file("extdata", "data_model.csv", package="curve"))

## print debug info
print(sprintf("Package version: %s", packageVersion('curve')))
print(sprintf("cnvex version: %s", packageVersion('cnvex')))
print(sprintf("Database: %s",db['db']))
print(sprintf("Database Location: %s",db['host']))
print(sprintf("Hostname: %s",hn))
