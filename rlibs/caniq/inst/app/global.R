library(devtools)
load_all("../../../cargo")

library(shiny)
library(gridlayout)
library(bslib)
library(DT)

library(ggplot2)


pool <- pool::dbPool(
    drv = duckdb::duckdb(),
    dbdir = "/work/tmp/test-small33.db",
    read_only = TRUE
)

onStop(function() {
    print("Disconnecting DuckDB...")
    # con <- pool::poolCheckout(pool)
    # duckdb::dbDisconnect(con)
    # pool::poolReturn(con)
    pool::poolClose(pool)
    print("... disconnected.")

})
