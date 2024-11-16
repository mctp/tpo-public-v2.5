#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(data.table)
library(stringr)
library(DT)

# Define UI for application that draws a histogram
ui <- fluidPage(

    titlePanel("Structural Report Viewer"),

    sidebarLayout(
        sidebarPanel(
            textInput("carat_report", "Report File:", value = "reports13.rds"),
            selectInput("patient", "Patient:", choices = NULL),
            selectInput("tnid", "Sample:", choices = NULL),
            checkboxInput("qual_filter", "Threshold Filter", value=FALSE),
            checkboxInput("manta_filter", "Manta Filter", value=TRUE),
            checkboxInput("splice_filter", "Splice Filter", value=TRUE),
            width = 3
        ),

        mainPanel(
            width = 9,
            tabsetPanel(
                tabPanel("Gene Variants",
                         dataTableOutput("loss_tbl")
                         ),
                tabPanel("Gene Fusions",
                         dataTableOutput("fusion_tbl")
                         ),
                tabPanel("All Variants",
                         dataTableOutput("full_tbl")
                )
            )
        )
    )
)

server <- function(input, output, session) {
    
    REPORTS <- reactive({
        reps <- readRDS(input$carat_report)
        ex1 <- as.integer(str_match(reps$EXON1, "([0-9]*)/")[,2])
        in1 <- as.integer(str_match(reps$INTRON1, "([0-9]*)/")[,2])
        ex2 <- as.integer(str_match(reps$EXON2, "([0-9]*)/")[,2])
        in2 <- as.integer(str_match(reps$INTRON2, "([0-9]*)/")[,2])
        b1 <- ifelse(!is.na(ex1), ex1, in1+1)
        b1[is.na(b1)] <- ""
        b2 <- ifelse(!is.na(ex2), ex2, in2)
        b2[is.na(b2)] <- ""
        reps$exon.range <- ifelse(b1==b2, b1,
                           ifelse(b1=="", paste0("--", b2), 
                           ifelse(b2=="", paste0(b1, "--"),                                  
                           ifelse(as.integer(b1)>as.integer(b2), NA_integer_,
                           paste(b1, b2, sep="--")))))
        setkey(reps, tnid)
        return(reps)
    })
    
    REPORTS.FILT <- reactive({
        reps <- REPORTS()
        if (input$qual_filter) {
            reps <- reps[ADT>6 & AFT>0.025]
        }
        if (input$manta_filter) {
            reps <- reps[MANTA.FILTER=="PASS"]
        }
        return(reps)
    })
    
    TBL <- reactive({
        reps <- REPORTS.FILT()
        reps[J(input$tnid),nomatch=0]        
    })
    
    output$full_tbl <- renderDataTable({
        tbl <- TBL()
        tbl <- tbl[,.(
             `Gene 5'`=SYMBOL1, `Gene 3'`=SYMBOL2, 
             `Tx 5'`=TRANSCRIPT1, `Tx 3'`=TRANSCRIPT2,
             Type=topo,
             Breakends=paste(paste(chr1, pos1, sep=":"), paste(chr2, pos2, sep=":"), sep="<br>"),
             `Exons Affected`=exon.range,
             Reads=ADT,
             Depth=DPT,
             VAF=AFT,
             QUAL,
             END1, END2,
             SCORE=MANTA.SCORE,
             EXON1, INTRON1, STRAND1,
             EXON2, INTRON2, STRAND2,
             Insert=insert,
             Contig=MANTA.CONTIG
        )]
        tbl <- datatable(rbind(tbl), rownames=FALSE, escape=FALSE) %>% formatRound(10, 2)
        
    })
    
    output$fusion_tbl <- renderDataTable({
        tbl <- TBL()
        tbl1 <- tbl[(SYMBOL1!=SYMBOL2) & !is.na(SYMBOL1) & !is.na(SYMBOL2)]
        tbl1 <- tbl1[,.(`Gene 5'`=SYMBOL1, `Gene 3'`=SYMBOL2, 
                        `Tx 5'`=TRANSCRIPT1, `Tx 3'`=TRANSCRIPT2,
                        Type=topo,
                        Breakends=paste(paste(chr1, pos1, sep=":"), paste(chr2, pos2, sep=":"), sep="<br>"),
                        Reads=ADT,
                        Depths=DPT,
                        VAF=AFT,
                        QUAL,
                        SCORE=MANTA.SCORE
        )]
        tbla <- datatable(rbind(tbl1), rownames=FALSE, escape=FALSE) %>% formatRound(9, 2)
        return(tbla)        
    })
    
    output$loss_tbl <- renderDataTable({
        tbl <- TBL()
        tbl1 <- tbl[((SYMBOL1==SYMBOL2) | is.na(SYMBOL1) | is.na(SYMBOL2)) & !is.na(exon.range)]
        tbl1[,SYMBOL:=SYMBOL1]
        tbl1[is.na(SYMBOL1), SYMBOL:=SYMBOL2]
        tbl1[,TRANSCRIPT:=TRANSCRIPT1]
        tbl1[is.na(TRANSCRIPT1), TRANSCRIPT:=TRANSCRIPT2]
        tbl1 <- tbl1[,.(Gene=SYMBOL, Transcript=TRANSCRIPT, Type=topo,
                        Breakends=paste(paste(chr1, pos1, sep=":"), paste(chr2, pos2, sep=":"), sep="<br>"),
                        `Exons Affected`=exon.range,
                        Reads=ADT,
                        Depth=DPT,
                        VAF=AFT,
                        QUAL,
                        SCORE=MANTA.SCORE
        )]
        
        tbla <- datatable(rbind(tbl1), rownames=FALSE, escape=FALSE) %>% formatRound(8, 2)
        return(tbla)

    })

    observe({
        reps <- REPORTS()
        patients <- unique(reps$patient)
        updateSelectInput(session, "patient", choices = patients)
    })
        
    observe({
        reps <- REPORTS()
        tnids <- unique(reps[patient==input$patient]$tnid)
        updateSelectInput(session, "tnid", choices = tnids)
    })
    
    
    


}

shinyApp(ui = ui, server = server)
