shinyUI(fluidPage(
  tags$head(HTML("<title>CURVE</title>")),
  sidebarLayout(
    sidebarPanel(width=3,
      uiOutput("choose_cohort"),
      uiOutput("choose_patient"),
      uiOutput("choose_analysis"),
      uiOutput("choose_tr"),
      uiOutput("choose_nr"),
      tags$br(),
      tableOutput('pt_info'),
      HTML(sprintf('<style type="text/css"> .well { background-color: %s; color: %s } </style>',color1,color2))
    ),
    mainPanel(width=9,
      tabsetPanel(id='tabs',
        tabPanel('Somatic Variants',value='som',
          h1('Somatic Variants'),
          variantsUI("somatic")
        ),
        tabPanel('Germline Variants',value='germline',
          h1('Germline Variants'),
          variantsUI("germline")
        ),
        tabPanel('Structural Variants',value='struc',
          h1('Structural Variants'),
          svUI("struc")
        ),
        tabPanel('Expression',value='gxp',
          h1('Gene Expression'),
          gxpUI("gxp")
        ),
        tabPanel('Fusion',value='fus',
          h1('Gene Fusions'),
          fusionUI("fusion")
        ),
        tabPanel('Copy Number',value='cnv',
          tabsetPanel(id='cnv_type',type='pills',
            tabPanel('Somatic',cnvUI('som_cnv'))
            #tabPanel('Germline',cnvUI('ger_cnv'))
          )
        ),
        tabPanel('Additional Data',value='sideload',
          h1('Additional Data'),
          sideloadDataUI('sld')
        ),
        tabPanel('Dev',value='dev',
          h1('Some Relevant Values'),
          fluidRow(column(12,align='center',tableOutput("dev_table"),plotOutput('var_dens',width='40%')))
        )
      )
    )
  )
))
