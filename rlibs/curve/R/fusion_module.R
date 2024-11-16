#'@export
fusionUI <- function(id){
  ns <- NS(id)
  tagList(
  uiOutput(ns("checkbox")),
  withSpinner(DT::dataTableOutput(ns('fusion_table'),width='95%')),
  downloadButton(ns('download_fusions'),'Download'),
  )
}


#'@export
fusion <- function(input, output, session, runid){
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


  FUSIONS <- reactive({
    out <- db_get_fus(runid())
    if('pass' %in% input$include){
      out <- out[out$triage=='PASS',]
    }
    out <- format_fus(out,dict)
    return(out$table)
  })
  output$fusion_table <- DT::renderDataTable({
    out <- FUSIONS()
    DT::datatable(out, 
      rownames=F,
      selection='single',
      filter='bottom',
      options=list(
        pageLength=22,
        autoWidth=T,
        dom='itlp'#,
        #columnDefs=list(list(visible=FALSE, targets=c(out$hide-1)))
      )
    )# %>%
      #DT::formatPercentage(out$renderpct,2)
  })

  output$download_fusions <- downloadHandler(
    filename=function(){
        paste0(sprintf('%s_fusions.csv',runid()))
      },
    content=function(f){write.csv(FUSIONS(),file=f,quote=F,row.names=F)}
  )

}



