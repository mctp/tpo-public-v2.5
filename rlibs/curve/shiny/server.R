shinyServer(function(input, output, session) {

  #----------- A variants display (from DB) -------------
  callModule(variants,"somatic",
             runid=reactive({input$runid}),
             tp="somatic"
           )

  #----------- A variants display (from DB) -------------
  callModule(variants,"germline",
             runid=reactive({input$runid}),
             tp="germline"
           )

  #----------- An structural variant display (from DB) -------------
  callModule(sv,"struc",
              runid=reactive({input$runid})
            )

  #----------- Expression Display (from DB) -------------
  callModule(gxp,"gxp",
              runid=reactive({input$runid})
            )

  #----------- Fusion Display (from DB) -------------
  callModule(fusion,"fusion",
              runid=reactive({input$runid})
            )

  #----------- CNV Display (from DB) -------------
  callModule(cnv,"som_cnv",
              runid=reactive({input$runid})
            )

  #----------- Sideloaded Data Display (from DB) -------------
  callModule(sideloadData,"sld",
              runid=reactive({input$runid})
            )

  output$pt_info <- renderTable(colnames=FALSE,{
    req(input$runid)
    db_get_meta(input$runid)
  })

  #----------- Debug Table -------------
  output$dev_table <- renderTable({
    out <- rbind(
      c('curve version',as.character(packageVersion('curve'))),
      c('DB Location',db['host']),
      c('DB Name',db['db']),
      c('DB User',db['user']),
      c('cnvex version',as.character(packageVersion('cnvex'))),
      c('codac version',as.character(packageVersion('codac'))),
      c("cohort",input$cohort),
      c("patient",input$patient),
      c("runid",input$dna),
      c("host",system('hostname',intern=T))
      )
      out <- as.data.frame(out)
      colnames(out) <- c('Attribute','Value')
      out
  })


  #----------- Input Boxes -------------
  # cohort
  output$choose_cohort <-renderUI({
    validate(need(nrow(groupings>0),"need groupings"))
    selectInput('cohort', 'Cohort:',
      choices=sort(unique(groupings$cohort)),
      selected='MO'
    )
  })

  # patient
  output$choose_patient <- renderUI({
    validate(need(!is.null(input$cohort),"need cohort"))
    selectInput('patient', 'Patient:',
      choices=groupings %>% 
        filter(cohort %in% input$cohort) %>%
        .$patient %>%
        unique %>%
        sort(decreasing=TRUE)
    )
  })

  # runid
  output$choose_analysis <- renderUI({
    validate(need(!is.null(input$patient),"need patient"))
    choices <- groupings %>%
      filter(patient %in% !!input$patient) %>%
      select(id,alignid_t,alignid_n,alignid_tr,alignid_nr) %>%
      unique()  %>%
      mutate(show=sprintf("%s-%s (DNA) %s-%s (RNA)", alignid_t, alignid_n,alignid_tr,alignid_nr)) %>%
      arrange(show,decreasing=TRUE)
    # remove the patient id for display
    disp <- choices$id
    names(disp) <- choices$show
    radioButtons("runid", "Analysis:", disp)
  })
})
