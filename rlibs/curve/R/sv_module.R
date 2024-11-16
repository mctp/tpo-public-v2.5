#'@export
svUI <- function(id){
  ns <- NS(id)
  tagList(
  uiOutput(ns("checkbox")),
  withSpinner(DT::dataTableOutput(ns('sv_table'),width='95%')),
  downloadButton(ns('download_sv_seg'),'Download'),
  )
}


#'@export
sv <- function(input, output, session, runid){

  output$checkbox <- renderUI({
    ns <- session$ns
    #
      checkboxGroupInput(ns("include"), "Include:",
        c(
          "Pass Filter"="pass"
          ),
        inline=T,selected=c('pass')
      )
  })

  output$sv_table <- DT::renderDataTable({
    out <- db_get_sv(runid()) %>% arrange(desc(adt),desc(aft))
    if('pass' %in% input$include){
      out <- out[out$triage=='PASS',]
    }
    out <- format_sv(out,dict)
    DT::datatable(out$table, 
      rownames=F,
      selection='single',
      filter='bottom',
      options=list(
        pageLength=22,
        autoWidth=T,
        dom='itlp'#,
        #columnDefs=list(list(visible=FALSE, targets=c(out$hide-1)))
      )
    ) %>%
      DT::formatPercentage(out$renderpct,2)
  })

}



