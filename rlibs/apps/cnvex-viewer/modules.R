#'@export
gene_tableUI <- function(id){
    ns <- NS(id)
    tabsetPanel(type='pills',
      tabPanel("Gene Table",
        uiOutput(ns("goi_check")),
        withSpinner(DT::dataTableOutput(ns("cnv-gene-table")))
      ),
      tabPanel("Segment Table",
        withSpinner(DT::dataTableOutput(ns("cnv-seg-table")))
      )
    )

}

#'@export
gene_table <- function(input,output,session,dig,pd,goi=NULL,v4=NULL){
  output$`cnv-gene-table` <- DT::renderDataTable({
      digest <- dig()
      plot.data <- pd()
      gene.dt <- unique(plot.data$cov[,.(
          gene_name=unlist(gene_names),
          seg=rep(seg, lengths(gene_names)),
          C=rep(C, lengths(gene_names)),
          K=rep(K, lengths(gene_names)),
          chr=rep(chr, lengths(gene_names))
      )
      ])

      gene.dt$event <- 'Neutral'
      gene.dt[gene.dt$C>2,'event'] <- 'Gain'
      gene.dt[gene.dt$C<2,'event'] <- 'Loss'
      gene.dt[gene.dt$C>=2 & gene.dt$K==0,'event'] <- 'LOH'
      colnames(gene.dt) <- tools::toTitleCase(gsub('_',' ',colnames(gene.dt)))
      colnames(gene.dt) <- gsub('Seg','Segment',colnames(gene.dt))
      if(is.null(goi)){
        gene.dt$notable <- TRUE
      }else{
        req(!is.null(input$goi_only))
        gene.dt$notable <- (gene.dt$`Gene Name` %in% goi) & gene.dt$Event!='Neutral'
        if(input$goi_only){gene.dt <- gene.dt %>% filter(notable)}
      }
      hide <- which(colnames(gene.dt)=='notable')
      DT::datatable(gene.dt %>% arrange(Event),
        rownames = FALSE,
        filter = "top",
        selection='none',
        options = list(
          dom='itlp',
          columnDefs=list(list(visible=FALSE, targets=c(hide-1)))
        )
      ) %>%
        DT::formatStyle('Event',target='row',
          color=DT::styleEqual(c('Gain','Loss','Neutral','LOH'),c('red','blue','black','green'))
        ) %>%
        DT::formatStyle('notable',target='row',fontWeight=DT::styleEqual(TRUE,'bold'))
  })
  
  output$goi_check <- renderUI({
    ns <- session$ns
    checkboxInput(ns("goi_only"), "Notable Events Only", TRUE)
    })
  
  output$`cnv-seg-table` <- DT::renderDT({
      digest <- dig()
      seg <- digest$plot.data$seg[,.(seg, chr, start, end, len,  naf, nlr, C, sC, K, gene_names)] 
      seg$event <- 'Neutral'
      seg[C>2,'event'] <- 'Gain'
      seg[C<2,'event'] <- 'Loss'
      seg[C>=2 & K==0,'event'] <- 'LOH'
      seg <- seg_goi_format(seg,goi,v4) %>% 
        select(-c(sC)) %>% arrange(event) %>%
        mutate(len=round(len/10e6,2)) %>%
        dplyr::rename(Size=len, SNPs=naf, Exons=nlr) %>%
        mutate(Region=sprintf("%s:%d-%d", chr,start,end))
      colnames(seg) <- tools::toTitleCase(gsub('_',' ',colnames(seg)))
      seg <- select(seg, c("Region", "Notable Genes", "C", "K", "Event", "Gene Names", "Exons", "SNPs","Size"))
      if(!is.null(v4)){
        colnames(seg) <- gsub('Gene Names','v4 Genes',colnames(seg))
      }
      DT::datatable(seg, rownames = FALSE, filter = "top", selection='none') %>%
        DT::formatStyle('Event',target='row',
          color=DT::styleEqual(c('Gain','Loss','Neutral','LOH'),c('red','blue','black','green'))
        )
  })

}

#'@export
zoomUI <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(3, selectInput(
        ns("cnv-zoom-sel-chr"),
        "Select Chromosome",
        choices = c(1:22, "X", "Y")
      ))
    ),
    fluidRow(
      style = "margin-top: 15px;",
      withSpinner(plotlyOutput(ns("cnv-zoom-absolute-plot")))
    )
  )
}

#'@export
zoom <- function(input,output,session,dig,pd){

}
