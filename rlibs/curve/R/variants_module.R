#' @export
variantsUI <- function(id){
  ns <- NS(id)
  #select filter
  tagList(
    #make the filter boxes a little easier to read
    tags$head(tags$style(HTML("
      .form-control {font-size: 10px; padding: 1px 2px 1px 2px; height: 25px}
    "))),
    uiOutput(ns('tmb_text')),
    uiOutput(ns('purity_warning')),
    uiOutput(ns("checkbox")),
    actionButton(ns('view_info_modal'),'View'),
    downloadButton(ns('download_all'),'Download All'),
    downloadButton(ns('download_hgvs'),'Download HGVS'),
    withSpinner(
      DT::dataTableOutput(ns("variants"),width='95%')
    ),
    shinyBS::bsModal(id=ns("info_modal"),title="Variant Info",trigger=ns("view_info_modal"),size='large',
      withSpinner(uiOutput(ns("modal_contents")))
    ) 
 )
}



#' @export
variants <- function(input,output,session,runid,tp){

  VAR <- reactive({
    validate(need(!is.null(runid()), 'variants needs runid'))
    res <- db_get_variants(runid(),tp) %>% arrange(desc(cosmic_cnt),desc(tlod))
    validate(need(nrow(res)>0,"no variants found"))
    if(tp=='somatic'){
      hits <- GER_HITS()
    }else if(tp=='germline'){
      print('using ger hits')
      hits <- SOM_HITS()
    }
    res <- filter_variants(res,input$include,hits)
    return(res)
  })

  SOM_HITS <- reactive({
    validate(need(!is.null(runid()), 'variants needs runid'))
    res <- tbl(db_pool,'somatic') %>% filter(runid==!!runid()) %>% filter(triage=='pass') %>% select(gene_id) %>% collect()
    return(res$gene_id)

  })
  GER_HITS <- reactive({
    validate(need(!is.null(runid()), 'variants needs runid'))
    res <- tbl(db_pool,'germline') %>% filter(runid==!!runid()) %>% filter(triage=='pass') %>% select(gene_id) %>% collect()
    return(res$gene_id)

  })

  output$checkbox <- renderUI({
    ns <- session$ns
    #
    if(tp=='germline'){
      checkboxGroupInput(ns("include"), "Include:",
        c(
          "Pass Filter"="pass",
          "High/Med Impact Only"="imp",
          "Germline Panel Only"="ger",
          "Two Hits Only"="two"
          ),
        inline=T,selected=c('pass','imp','ger')
      )
    }else if(tp=='somatic'){
      checkboxGroupInput(ns("include"), "Include:",
        c(
          "Pass Filter"="pass",
          "High/Med Impact Only"="imp",
          "Two Hits Only"="two"
          ),
        inline=T,selected=c('pass','imp')
      )
    }
  })

  #the variants table
  var_table_opts <- list(
    pageLength=15,
    dom='itlp',
    scrollX=TRUE,
    autoWidth=TRUE
  )
  output$variants <- DT::renderDataTable({
#    if(tp=='somatic'){
#      hits <- GER_HITS()
#    }else if(tp=='germline'){
#      print('using ger hits')
#      hits <- SOM_HITS()
#    }
    #out <- filter_variants(VAR(),input$include,hits)
    out <- VAR()
    validate(need(nrow(out)>0, "no variants found"))
    out <- format_variants(out,dict,tp)
    #
    DT::datatable(out$table, 
      rownames=F, 
      escape=F,
      selection='single',
      filter=list(position='bottom',clear=FALSE),
      options=list(
        pageLength=15,
        dom='itlp',
        scrollX=TRUE,
        autoWidth=TRUE
        #columnDefs=cd
      )
    ) %>%
      DT::formatPercentage(out$renderpct,2)
#      DT::formatStyle('captured',target='row',fontWeight=DT::styleEqual(TRUE,'bold')) %>%
#      DT::formatStyle('Warnings',target='row',
#        color=DT::styleEqual('Possible Phasing Issue','grey')) %>%
      #DT::formatPercentage(renderpct,2)
  })


# #contents of the variant modal 
 fv <- reactive({
    w <- input$variants_rows_selected
    vars <- VAR()[w,]
    out <- format_variants(vars[1,],dict, tp,detail=TRUE)$table
    return(data.frame(Metric=colnames(out),Value=as.character(out[1,])))
 })

  output$modal_contents <- renderUI({
    ns <- session$ns
    validate(need(!is.null(input$variants_rows_selected),"Please Select a Variant"))
    tagList(
      DT::dataTableOutput(ns("ext_links2")),
      downloadButton(ns('download_some'),'Download Variant')
    )
  })
#
  output$ext_links2 <- DT::renderDataTable(
      {fv()},
      escape=F, 
      rownames=F,
      selection='none',
      options=list(dom='t', pageLength=50, autoWidth=TRUE, scrollY=TRUE)
 )


  #download buttons
  output$download_all <- downloadHandler(
    filename=function(){paste0(runid(),'_all_variants.csv')},
    content=function(f){write.csv(VAR(), f, row.names=F)}
  )
  output$download_some <- downloadHandler(
    filename=function(){
      paste0(
          runid(), '-',
          gsub('/','-',VAR()[input$variants_rows_selected,'id']),
          '.csv'
      )
    },
    content=function(f){
      write.csv(VAR()[input$variants_rows_selected,], f, row.names=F)
    }
  )
  
  output$download_hgvs <- downloadHandler(
    filename=function(){paste0(runid(),'_hgvs.csv')},
    content=function(f){
      out <- VAR() %>% 
        mutate(Variant=paste0(transcript,':',hgvsc,',',hgvsp)) %>%
        mutate(loc=paste0(chr,":",pos)) %>%
        mutate(aft=sprintf("%.1f%%",100*aft)) %>%
        mutate(afn=sprintf("%.1f%%",100*afn))
        out <- out %>% select(gene_name,loc,Variant,aft,afn) %>%
                        S4Vectors::rename(
                          'gene_name'='Gene', 'loc'='Location (GRCh38)', 
                          'aft'='Tumor Allelic Fraction', 'afn'='Normal Allelic Fraction'
                        )
    write.csv(out,file=f,row.names=F)
    }
  )

}


#' Filter variants according to checkboxes
#' @export
filter_variants <- function(vars,sel,hits=NULL){
  if('pass' %in% sel){
    vars <- vars[vars$triage=='pass',]
  }
  if('imp' %in% sel){
    vars <- vars[vars$impact %in% c('HIGH','MODERATE'),]
  }
  if('ger' %in% sel){
    vars <- vars[vars$gene_name %in% curve_goi$germline,]
  }
  if('two' %in% sel){
    w <- which(
      vars$gene_id %in% hits | vars$cn<2 | vars$k==0
    )
    vars <- vars[w,]
  }
return(vars)
}
