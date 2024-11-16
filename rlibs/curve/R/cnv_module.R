#'@export
cnvUI <- function(id){
  ns <- NS(id)
  tagList(
    withSpinner(plotOutput(ns('cnv_plot'))),
    downloadButton(ns('download_cnv_img'),'Save'),
    hr(),
    tabsetPanel(type='pills',
      tabPanel("Gene Table",
        checkboxInput(ns("goi_only"), "Notable Events Only", TRUE),
        withSpinner(DT::dataTableOutput(ns('gene_table'),width='70%')),
        downloadButton(ns('download_cnv_genes'),'Download')
      ),
      tabPanel("Segment Table",
        checkboxInput(ns("notable_only"), "Notable Events Only", TRUE),
        tags$br(),
        withSpinner(DT::dataTableOutput(ns('seg_table'),width='95%')),
        downloadButton(ns('download_cnv_seg'),'Download')
      )
    )
  )
}


#'@export
cnv <- function(input, output, session, runid){

  CNV_GENES <- reactive({
    validate(need(!is.null(runid()), 'gene table needs runid'))
    res <- db_get_cnv(runid())
    res <- format_cnv(res,dict, table='genes')
    if(input$goi_only){
      res$table <- res$table[res$table$Gene %in% curve_goi$cnv  & res$table$Event!='Neutral',]
    }
    return(res)
  })

  CNV_SEGS <- reactive({
    validate(need(!is.null(runid()), 'gene table needs runid'))
    out <- db_get_cnv(runid())
    out <- format_cnv(out,dict, table='segments')$table
    out$Chr <- factor(out$Chr,levels=paste0('chr',c(1:22,'X','Y')))
    return(out)
  })

  CNV_PLOT <- reactive({
    validate(need(!is.null(runid()), 'cnv needs runid'))
    validate(need(!is.null(CNV_GENES()$table), 'cnv needs genes'))
    out <- CNV_SEGS() %>% mutate(width=End-Start)
    col <- c('#FB6A4A','#6B58EE','black','#78b09c', "#4500AC", 'gray')
    names(col) <- c('Gain','Loss','Neutral','LOH','Hom. Deletion','Unknown')



    #
    # a quick segment plot
    g <- ggplot(out) + 
            geom_segment(aes(x=Start,xend=End,y=C,yend=C, color=Event),size=1.5, lineend='square') + 
#            geom_point(aes(x=Start,y=C,color=Event,fill=Event),shape=21) + 
#            geom_point(aes(x=End,y=C,color=Event,fill=Event),shape=21) + 
            facet_wrap(.~Chr, nrow=1, scales='free_x') + 
            scale_color_manual(values=col) +
            scale_fill_manual(values=col) +
            theme_minimal() + 
            ylim(c(0,max(out$C)+1)) + 
            theme(
              axis.text.x = element_blank(),
              panel.grid.minor.x=element_blank(),
              panel.grid.major.x=element_blank(),
              text=element_text(size=16),
              legend.position='bottom'
              ) + 
            xlab('') + ylab('Copies')

    #see if a gene is selected
    w <- input$gene_table_rows_selected
    if(!is.null(w)){
      seg <- CNV_GENES()$table[w,]
      height <- max(0.75 * layer_scales(g)$y$get_limits()[2]-seg$C, 3) # put it at the top of the plot at least C=3
      g <- g + geom_point(
                  aes(x=(End+Start)/2, y=C),
                  color='goldenrod',shape=1, size=3,
                  data=out %>% filter(`Seg ID` %in% !!seg$`Seg ID`),
               ) + 
               geom_text_repel(
                   aes(x=(End+Start)/2, y=C, label=seg$Gene, color=Event), size=3,
                   data=out %>% filter(`Seg ID` %in% !!seg$`Seg ID`),
                   min.segment.length=0,
                   nudge_y=height,
                   arrow=arrow(ends='last', length=unit(5,'pt'), type='closed'),
                   show.legend=FALSE
                ) 
    }
  g <- g + guides(fill='none', label='none',color=guide_legend(override.aes=list(size=5,shape=NA, label=NA))) 

  #
  return(g)
  })

  output$cnv_plot <- renderPlot({
    CNV_PLOT()
  })

  output$seg_table <- DT::renderDataTable({
    out <- CNV_SEGS()
    if(input$notable_only){
      out <- out[out$Event!='Neutral',]
    }
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
    ) %>%  DT::formatStyle('Event',target='row',
                            color=DT::styleEqual(
                              c('Gain','Loss','Neutral','LOH','Hom. Deletion','Unknown'),
                              c('#FB6A4A','#6B58EE','black','#78b09c', "#4500AC", 'gray')
                            )
                          )
  })

  output$gene_table <- DT::renderDataTable({
    out <- CNV_GENES()
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
    ) %>%  DT::formatStyle('Event',target='row',
                            color=DT::styleEqual(
                              c('Gain','Loss','Neutral','LOH','Hom. Deletion','Unknown'),
                              c('#FB6A4A','#6B58EE','black','#78b09c', "#4500AC", 'gray')
                            )
                          )
  })

  output$download_cnv_img <- downloadHandler(
    filename=function(){
        paste0('cnv_landscape_image.png')
      },
    content=function(f){ggsave(CNV_PLOT(), file=f, height=6,width=18)}
  )

  output$download_cnv_seg <- downloadHandler(
    filename=function(){
        paste0('cnv_segments.tsv')
      },
    content=function(f){write.table(CNV_SEGS()$table,file=f,sep='\t',quote=F,row.names=F)}
  )

  output$download_cnv_genes <- downloadHandler(
    filename=function(){
        paste0('cnv_genes.csv')
      },
    content=function(f){write.csv(CNV_GENES()$table,file=f,quote=F,row.names=F)}
  )
}
