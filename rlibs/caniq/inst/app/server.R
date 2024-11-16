

server <- function(input, output) {

  output$somatic_table <- renderDT({
    sel.sample <- input$somatic_select_sample
    got.sample <- sel.sample != ""
    sel.symbol <- input$somatic_select_symbol
    got.symbol <- sel.symbol != ""
    sel.pass <- input$somatic_filter_pass
    sel.impact <- input$somatic_filter_impact

    ## somatic table
    query_somatic <- dplyr::tbl(pool, "tables_somatic") |>
    ## flter by pass
    dplyr::filter( if (sel.pass) filter == "PASS" else TRUE) |>
    ## flter by impact
    dplyr::filter( if (sel.impact) IMPACT %in% c("MODERATE", "HIGH") else TRUE) |>
    ## filter by sample
    dplyr::filter( if (got.sample) id %in% sel.sample else TRUE) |>
    ## filter by bene
    dplyr::filter( if (got.symbol) SYMBOL %in% sel.symbol else TRUE)

    ## expression
    query_expression <- dplyr::tbl(pool, "tables_gene_expression") |>
    ## filter by sample
    dplyr::filter( if (got.sample) id %in% sel.sample else TRUE) |>
    ## filter by bene
    dplyr::filter( if (got.symbol) SYMBOL %in% sel.symbol else TRUE)

    ## join
    query_full <- query_somatic |> dplyr::inner_join(query_expression, by=c("id"="id", "GENE"="gene_id", "SYMBOL"="SYMBOL"))
    ## formatting
    query_full |>
    dplyr::select(id, SYMBOL, HGVSp, HGVSc, AFT, AFN, tlod, nlod, rpkm) |>
    head(n=input$somatic_limit_maxrow) |> dplyr::collect()

  }, rownames=FALSE)
}
