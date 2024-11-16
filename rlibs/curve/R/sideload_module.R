#'@export
sideloadDataUI <- function(id){
  ns <- NS(id)
  tagList(
  uiOutput(ns("sld_choice")),
  actionButton(ns('view_info_modal'),'View'),
  withSpinner(DT::dataTableOutput(ns("sld_display"),width='65%')),
  uiOutput(ns('checkbox')),
  bsModal(id=ns("info_modal"),title="Info",trigger=ns("view_info_modal"),size='large',
    plotOutput(ns('sld_sel_plot')),
    DT::dataTableOutput(ns('sld_table'))
  )
  )
}


#'@export
sideloadData <- function(input, output, session, runid){

  output$checkbox <- renderUI({
    ns <- session$ns
    #
      checkboxGroupInput(ns("include"), "Include:",
        c(
          "NA"="incNA"
          ),
        inline=T,selected=c()
      )
  })

  output$sld_choice <- renderUI({
    ns <- session$ns
    opts <- db_get_sld_list(runid())
    selectInput(ns('sld_choice'),'Dataset',
      choices=opts
    )
  })

  sld <- reactive({
    validate(need(!input$sld_choice=="",""))
    db_get_sld(runid(),input$sld_choice)
  })


  output$sld_display <- DT::renderDataTable({
    ns <- session$ns
    validate(need(!is.null(runid()), ""))
    validate(need(!is.null(input$sld_choice),""))
    ds <- sld()
    if(!'incNA' %in% input$include){
      w <- which(colnames(ds)==input$sld_choice)
      ds <- ds[!is.na(ds[,w]),] 
    }
    DT::datatable(ds, 
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

  output$sld_table <- DT::renderDataTable({
    w <- input$sld_display_rows_selected
    id <- sld()$id[w]
    dist <- db_get_sld_dist(id,input$sld_choice)
    out <- dist %>% arrange(desc(id==runid()),value) %>% select(-id)
    DT::datatable(out, 
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

  output$sld_sel_plot <- renderPlot({
    w <- input$sld_display_rows_selected
    id <- sld()$id[w]
    dist <- db_get_sld_dist(id,input$sld_choice)
    sideload_feature_plot(runid(),id,dist,input$sld_choice)
  })


}



#' @export
sideload_feature_plot <- function(runid,id,dist,nm){
  mblue <- '#00274c'
  maize <- '#ffcb05'
  dist <- arrange(dist,desc(value))
  dist$x <- seq(nrow(dist))
  sel <- dist %>% filter(id %in% !!runid) %>% arrange(desc(value))

  this_cohort <- sel$cohort[1]
  canc <- dist %>% group_by(cohort) %>% 
    summarize(n=n(),mn=median(value)) %>% 
    filter(n>10 | cohort %in% this_cohort) %>% 
    filter(!is.na(cohort) | cohort %in% this_cohort) %>%
    arrange(desc(mn)) %>% .$cohort

  to_plot <- dist %>% filter(cohort %in% canc) %>%
                        mutate(cohort=factor(cohort,levels=canc))
  #
  g <- ggplot(to_plot, aes(x=cohort,y=value)) + 
    geom_boxplot(color=mblue,outlier.shape=NA) + 
    geom_jitter(
      color=maize,
      alpha=0.6,
      height=0,
      width=1/(2*length(canc))
      ) + 
    ylab(nm) + xlab('') + 
    theme_minimal(base_size=14)  + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    geom_point(data=sel, fill='red', size=3, shape=23) + 
    scale_x_discrete(labels=function(x){str_wrap(x,20)}) + 
    ggtitle(sprintf('%s %s', id,nm))  + 
    scale_y_sqrt()
    
  return(g)
  }

