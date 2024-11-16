source("funcs.R")
source("modules.R")

shinyUI(fluidPage(
  titlePanel("CNVEX browser"),
  tabsetPanel(
    tabPanel(
      "Load Data",
      style = "margin-top: 15px;",
      textInput(
          "cnvex-path",
          "CNVEX Results Path",
          value = CONFIG$default.dir,
          width = "100%"
      ),
      fluidRow(
        column(
          6,
          selectInput("cnvex-run", "CNVEX Run", choices = NULL, width = "100%"),
          checkboxInput("hide-done","Hide Complete",FALSE),
        ),
        column(
          2,
          selectInput("cnvex-model", "Model", choices = c("somatic", "germline"), width = "100%")
        ),
        column(
          2,
          style = "margin-top: 25px;",
          actionButton("cnvex-load-next", "Next", width = "100%", style = "color: black; border-color: #2e6da4")
        ),
        column(
          2,
          style = "margin-top: 25px;",
          actionButton("cnvex-load", "Load", width = "100%", style = "color: black; background-color: #ff8080; border-color: #2e6da4")
        )
      ),
      verbatimTextOutput("cnvex-files", placeholder = TRUE),
      textInput(
        "cnvex-gtf",
        "GTF File",
        value = CONFIG$default.gtf,
        width = "100%"
        
      ),
      p(strong("Stats:")),
      verbatimTextOutput("cnvex-stat")
    ),
    
    tabPanel(
      "QC",
      fluidRow(
        style = "margin-top: 15px; ",
        column(4, withSpinner(plotOutput("qc-gc"))),
        column(4, withSpinner(plotOutput("qc-depth"))),
        column(4, withSpinner(plotOutput("qc-lr-grid")))
      ),
      fluidRow(
          column(2, numericInput("qc-maxn", label = "Number of Points", min=25000, max=1e6, step=25000, value = CONFIG$default.qc.maxn, width="100%"))
      )      
    ),
    
    tabPanel(
      "Select Model",
      fluidRow(
        style = "margin-top: 15px;",
        checkboxGroupInput("model-warn", "Model warnings", choices = 
                               c("over-segmented", "under-segmented", "noisy", 
                                 "BMT", "no-SNPs", "no-CNVs", "no-tumor",
                                 "unstable", "subclonal", "high-ploidy", 
                                 "fit-fail", "GC-bias", "coverage", "indeterminate"), inline = TRUE)
      ),
      fluidRow(
        style = "margin-top: 15px;",
        column(
          6,
          wellPanel(
            style = "height:220px",
            fluidRow(
              column(
                4,
                selectInput(
                  "model", "Model",
                  choices = NULL, width = "100%"
                ),
                checkboxInput("model-adjust", "Custom values", value = TRUE)
              ),
              column(
                4,
                sliderInput(
                  "ploidy",
                  "Tumor Ploidy",
                  2.0,
                  min = 1,
                  max = 9,
                  step = 0.001,
                  width = "100%"
                ),
                sliderInput(
                  "purity",
                  "Tumor Purity",
                  0.5,
                  min = 0,
                  max = 1,
                  step = 0.001,
                  width = "100%"
                )
              ),
              column(
                4,
                sliderInput(
                  "model-mc",
                  "M/C",
                  c(0, 1),
                  min = 0,
                  max = 10,
                  step = 1,
                  width = "100%"
                ),
                sliderInput(
                  "abs-range",
                  "Copy Number Range",
                  10,
                  min = 4,
                  max = 20,
                  step = 1
                )
              )
            )
          )
        ),
        column(6, wellPanel(
          fluidRow(
            column(
              6,
              selectInput("abs-sel-data", "Select Data", choices = c("segment", "tile")),
              selectInput("abs-sel-col", "Select Plot Color", choices = c("segment", "C"))
            ),
            
            column(
              6,
              selectInput(
                "abs-sel-chr",
                "Select Chromosome", choices = NULL
              ),
              sliderInput(
                "abs-hex-bins",
                "Hex Bins",
                100,
                min = 50,
                max = 200,
                step = 10
              )
            )
          )
        ))
      ),
      fluidRow(
        column(
          9,
          withSpinner(
            plotOutput("cnv-absolute-plot")
          )
        ),
        column(
          3,
          withSpinner(
            plotOutput("abs-baf-grid") 
          )
        )
      ),
      fluidRow(column(6, wellPanel(
        dataTableOutput("model-table")
      )))
    ),

    tabPanel(
        "Landscape",
        fluidRow(
            style = "margin-top: 15px;",
            withSpinner(plotOutput("landscape-plot"))
        ),
        fluidRow(
            column(2, selectInput(
                "landscape-sel-col",
                "Select Plot Color",
                choices = c(
                    "segment",
                    "CK"
                )
            )),
            column(2, selectInput(
                "landscape-sel-chr",
                "Select Chromosome", choices = NULL
            )),
            column(
                2,
                sliderInput(
                    "landscape-lr-range-slider",
                    "LogR Range",
                    min = -10,
                    max = 10,
                    value = c(-4, 4)
                ),
                actionButton("landscape-lr-model", label = "from model")
            ),
            column(
                2,
                sliderInput(
                    "landscape-baf-range-slider",
                    "BAF Range",
                    min = 0,
                    max = 1,
                    value = c(0, 1)
                ),
                actionButton("landscape-baf-model", label = "from model")
            ),
            downloadButton('landscape-save-plot', 'Save')
        )
    ),
    tabPanel("Zoom",
             fluidRow(
                 column(2, selectInput(
                     "zoom-sel-chr",
                     "Select Chromosome", choices = NULL
                     ))
             ),
             fluidRow(
                 withSpinner(plotlyOutput("zoom-plot"))
             )
    ),
    tabPanel(
      "Checkout",
      fluidRow(
        column(
          4,
          p(strong("Digest:")),
          verbatimTextOutput("checkout-digest", placeholder = TRUE)
        ),
        column(
          4,
          textAreaInput(
            "checkout-notes",
            label = "Notes:",
            width = "100%",
            height = "100px"
          ),
          actionButton("digest-checkout", "Checkout")
        )
      )
    )
  )
))
