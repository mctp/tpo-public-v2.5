#' @export
gxpUI <- function(id){
  ns <- NS(id)
  tabsetPanel(type='pill',
    tabPanel('Gene Expression',
        fluidRow(
          column(1),
          column(8,
            hr(),
            actionButton(ns('view_info_modal'),'View'),
            downloadButton(ns('download_gxp'),'Download'),
            withSpinner(DT::dataTableOutput(ns('gxp_table'),width='75%')),
            bsModal(id=ns("info_modal"),title="Gene Info",trigger=ns("view_info_modal"),size='large',
              uiOutput(ns('cancer_cohort_picker')),
              withSpinner(plotOutput(ns('gene_info'))),
              withSpinner(DT::dataTableOutput(ns('gene_table'), width='95%')),
              downloadButton(ns('download_dist'),'Download'),
              downloadButton(ns('download_dist_img'),'Save'),
              hr(),
              radioButtons(ns('plot_type'),NULL, 
                choices=list(Hist='Hist',Box='Box', Point='Point'), 
                selected='Box',
                inline=T
              )
            )
          ),
          column(3)
        )
      )
    #  tabPanel('Signatures',
    #    tagList(
    #      hr(),
    #      uiOutput(ns('dist_picker')),
    #      fluidRow(
    #        column(3,withSpinner(DT::dataTableOutput(ns('sig_gene_table')))),
    #      column(9,align='right',
    #        uiOutput(ns('gene_sig_ht'))
    #      )
    #    )
    #  )
    #)
  )
}

#' @export
gxp <- function(input, output, session, runid){

get_rna <- reactive({
  validate(need(runid()!='NA',paste("No expression data found\n")))
  data <- db_get_gxp(runid())
  validate(need(nrow(data)>0,paste("No expression data found\n")))
  data <- match_ntile(data,gxp_ntile) %>% arrange(desc(rpkm))
  return(data)
  })

output$gxp_table <- DT::renderDataTable({
  out <- format_gxp(get_rna(),dict)

  DT::datatable(out$table, 
    rownames=F,
    selection='single',
    filter='bottom',
    options=list(
      pageLength=22,
      autoWidth=T,
      dom='itlp'
    )
  ) 
})

output$download_gxp <- downloadHandler(
  filename=function(){paste0(runid(),'_gxp.csv')},
  content=function(f){write.csv(get_rna(), f, row.names=F)}
)

#Modal contents
dist <- reactive({
  w <- input$gxp_table_rows_selected
  gene_id <- get_rna()[w,]$gene_id
  dist <- dbGetQuery(
    db_pool,
    sprintf("select * from get_gene('%s')",gene_id)
  ) 
  #
  return(dist)
})
#
#
  output$gene_info <- renderPlot({
    w <- input$gxp_table_rows_selected
    gene <- get_rna()[w,]
    gxp_gene_plot(runid(),gene$gene_name,gene$gene_id,dist(),input$plot_type)
  })

  output$gene_table <- DT::renderDataTable({
      out <- dist()
      out <- format_gxp(out,dict,detail=TRUE)$table
      DT::datatable(out, rownames=F, selection='none',filter='bottom',
        options=list(pageLength=10,dom='itp')) 

  })



  output$download_dist <- downloadHandler(
    filename=function(){
      w <- input$gxp_table_rows_selected
      gene_id <- get_rna()[w,]$gene_id
      paste0(gene_id,'_distribution.csv')
    },
    content=function(f){write.csv(dist(), f, row.names=F)}
  )
#
  output$download_dist_img <- downloadHandler(
    filename=function(){
        w <- input$gxp_table_rows_selected
        gene_id <- get_rna()[w,]$gene_id
        sprintf('%s_%s_expression.pdf',runid(),gene_id)
      },
    content=function(f){
      w <- input$gxp_table_rows_selected
      gene <- get_rna()[w,]
      ggsave(
        gxp_gene_plot(runid(),gene$gene,gene$gene_id,dist(),input$plot_type),
        file=f, 
        height=7,
        width=8
      )
    }
  )
}

#' @export
gxp_gene_plot <- function(runid,gene,id,dist,type){
  mblue <- '#00274c'
  maize <- '#ffcb05'
  dist <- arrange(dist,rpkm)
  dist$x <- seq(nrow(dist))
  sel <- dist %>% filter(runid %in% !!runid) %>% arrange(gene)
  if(type=='Hist'){
    g <- ggplot(dist, aes(x=rpkm)) + 
      geom_density(color=mblue, adjust=1) + 
      scale_x_sqrt() + 
      ylab("Density") + xlab("RPKM") + 
      theme_minimal(base_size=14) + 
      ggtitle(sprintf('Expression distribution for %s (%d samples)', gene, nrow(dist))) 
    #get the smoothed density data and find the yvalue
    d <- ggplot_build(g)$data[[1]]
    val <- sqrt(sel$rpkm)
    w <- which.min(abs(d$x-val))
    yval <- d$y[w]
      
    g <- g +  geom_label_repel(data=sel, aes(x=rpkm, y=yval, label=alignid),
                 size=3, color=maize, fill=mblue, direction='both',
                 force=100,point.padding=unit(2,'pt'),min.segment.length=unit(0,'pt')) + 
              geom_point(data=sel,aes(x=rpkm,y=yval), size=2,color=mblue)

  }else if(type=='Point'){
    g <- ggplot(dist, aes(x=x,y=rpkm)) + 
      geom_point(color=mblue) + 
      ylab('RPKM') + xlab('') + 
      theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) + 
      theme_minimal(base_size=14) + 
      ggtitle(sprintf('Expression distribution for %s (%d samples)', gene, nrow(dist)))  + 
      geom_point(data=sel,aes(x=x,y=rpkm), size=3,shape=23,color=mblue, fill='red') + 
      geom_label_repel(data=sel, aes(x=x, y=rpkm, label=alignid),
                 size=3, color=maize, fill=mblue, direction='both',
                 force=100,point.padding=unit(2,'pt'),min.segment.length=unit(0,'pt'))

  }else if(type=='Box'){
    this_cohort <- sel$cohort[1]
    canc <- dist %>% group_by(cohort) %>% 
      summarize(n=n(),mn=median(rpkm)) %>% 
      filter(n>10 | cohort %in% this_cohort) %>% 
      filter(!is.na(cohort) | cohort %in% this_cohort) %>%
      arrange(desc(mn)) %>% .$cohort

    to_plot <- dist %>% filter(cohort %in% canc) %>%
                          mutate(cohort=factor(cohort,levels=canc))
    #
    g <- ggplot(to_plot, aes(x=cohort,y=rpkm)) + 
      geom_boxplot(color=mblue,outlier.shape=NA) + 
      geom_jitter(
        color=maize,
        alpha=0.6,
        height=0,
        width=1/(2*length(canc))
        ) + 
      ylab('RPKM') + xlab('') + 
      theme_minimal(base_size=14)  + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      geom_point(data=sel, fill='red', size=3, shape=23) + 
      scale_x_discrete(labels=function(x){str_wrap(x,20)}) + 
      ggtitle(sprintf('%s Expression', gene))  + 
      scale_y_sqrt()
    
  }

  return(g)
}

#' Annotate a gxp dataset with the various ntile cutoffs
#' @export
match_ntile <- function(data,gxp_ntile){
  m <- merge(data, gxp_ntile, by='gene_id', all.x=T)
  m$rpkm_ntile <- NA
  m$rpkm_ntile[m$rpkm<m$pct_20] <- '<20%'
  m$rpkm_ntile[m$rpkm>=m$pct_20] <- '>20%'
  m$rpkm_ntile[m$rpkm>=m$pct_30] <- '>30%'
  m$rpkm_ntile[m$rpkm>=m$pct_40] <- '>40%'
  m$rpkm_ntile[m$rpkm>=m$pct_50] <- '>50%'
  m$rpkm_ntile[m$rpkm>=m$pct_75] <- '>75%'
  m$rpkm_ntile[m$rpkm>=m$pct_90] <- '>90%'
  m$rpkm_ntile[m$rpkm>=m$pct_95] <- '>95%'
  m$rpkm_ntile[m$rpkm>=m$pct_99] <- '>99%'
  m$rpkm_ntile[m$rpkm>=m$pct_999] <- '>99.9%'
  return(m %>% select(gene_name=gene_name.x,gene_id,rpkm,rpkm_ntile,c) %>% unique())
}
