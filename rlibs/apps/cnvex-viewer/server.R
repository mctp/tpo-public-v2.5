shinyServer(function(input, output, session) {
  print(sprintf('cnvex version:%s',packageVersion('cnvex')))

  #### Input slows
  input.landscape <- reactive({
        list(
            "landscape-sel-chr"=input$`landscape-sel-chr`,
            "landscape-sel-col"=input$`landscape-sel-col`,
            "landscape-lr-range-slider"=input$`landscape-lr-range-slider`,
            "landscape-baf-range-slider"=input$`landscape-baf-range-slider`
        )
  }) %>% debounce(1000)

  input.model <- reactive({
      list(
          purity = input$purity,
          ploidy = input$ploidy,
          mc = input$`model-mc`,
          abs.range = input$`abs-range`
      )
  }) %>% debounce(1000)

  # Data Objects

  DIR <- reactive({
    dir.names <- list.dirs(input$`cnvex-path`, recursive = FALSE, full.names = FALSE)
    dir.lists <- foreach(dir.name = dir.names) %do% {
      cnv <- paste0(file.path(input$`cnvex-path`, dir.name, basename(dir.name)), ".rds")
      if (!file.exists(cnv)) {
        cnv <- NA_character_
      }
      cnv.somatic.model <- paste0(file.path(input$`cnvex-path`, dir.name, basename(dir.name)), "-somatic-model.rds")
      if (!file.exists(cnv.somatic.model)) {
        cnv.somatic.model <- NA_character_
      }
      cnv.germline.model <- paste0(file.path(input$`cnvex-path`, dir.name, basename(dir.name)), "-germline-model.rds")
      if (!file.exists(cnv.germline.model)) {
          cnv.germline.model <- NA_character_
      }
      cnv.somatic.segment <- paste0(file.path(input$`cnvex-path`, dir.name, basename(dir.name)), "-somatic-segment.rds")
      if (!file.exists(cnv.somatic.segment)) {
        cnv.somatic.segment <- NA_character_
      }
      cnv.germline.segment <- paste0(file.path(input$`cnvex-path`, dir.name, basename(dir.name)), "-germline-segment.rds")
      if (!file.exists(cnv.germline.segment)) {
        cnv.germline.segment <- NA_character_
      }
      cnv.opts <- paste0(file.path(input$`cnvex-path`, dir.name, basename(dir.name)), "-opts.rds")
      if (!file.exists(cnv.opts)) {
        cnv.opts <- NA_character_
      }
      cnv.somatic.digest <- paste0(file.path(input$`cnvex-path`, dir.name, basename(dir.name)), "-somatic-digest.rds")
      if (!file.exists(cnv.somatic.digest)) {
          cnv.somatic.digest <- NA_character_
      }
      cnv.germline.digest <- paste0(file.path(input$`cnvex-path`, dir.name, basename(dir.name)), "-germline-digest.rds")
      if (!file.exists(cnv.germline.digest)) {
        cnv.germline.digest <- NA_character_
      }
      list(cnv = cnv, cnv.opts = cnv.opts,
           cnv.somatic.model = cnv.somatic.model, cnv.germline.model=cnv.germline.model,
           cnv.somatic.segment = cnv.somatic.segment, cnv.germline.segment=cnv.germline.segment,
           cnv.somatic.digest = cnv.somatic.digest, cnv.germline.digest = cnv.germline.digest
           )
    }
    names(dir.lists) <- basename(dir.names)
    return(dir.lists)
  })

  GOBJ <- reactive({
    if (!is.null(CNV()$tile)) {
        gobj <- cnvex::getGobj(unique(genome(CNV()$tile)), NULL, FALSE) 
        return(gobj)
    }
  })

  GENES <- reactive({
    cnvex::robustImport(input$`cnvex-gtf`, GOBJ()$seqi, feature.type = "gene")
  })
  
  CHR <- reactive({
    seqnames(GOBJ()$seqi)
  })

  CNV <- eventReactive(input$`cnvex-load`, {
    cnvex.run <- DIR()[[input$`cnvex-run`]]
    cnv <- readRDS(file.path(cnvex.run$cnv))
    cnv$file_name <- basename(cnvex.run$cnv)
    return(cnv)
  })

  OPTS <- eventReactive(input$`cnvex-load`, {
    cnvex.run <- DIR()[[input$`cnvex-run`]]
    opts <- readRDS(file.path(cnvex.run$cnv.opts))
    return(opts)
  })

  MODEL <- eventReactive(input$`cnvex-load`, {
    cnvex.run <- DIR()[[input$`cnvex-run`]]
    if (input$`cnvex-model`=="somatic") {
        opt.pth <- cnvex.run$cnv.somatic.model
    } else {
        opt.pth <- cnvex.run$cnv.germline.model
    }
    if (!is.na(opt.pth)) {
        opt <- readRDS(opt.pth)
        opt$file_name <- basename(opt.pth)
    } else {
        opt <- list(file_name=NULL)
    }
    return(opt)
  })

  SEG <- eventReactive(input$`cnvex-load`, {
    cnvex.run <- DIR()[[input$`cnvex-run`]]
    if (input$`cnvex-model`=="somatic") {
      seg.pth <- cnvex.run$cnv.somatic.segment
    } else {
      seg.pth <- cnvex.run$cnv.germline.segment
    }
    seg <- readRDS(seg.pth)
    return(seg)
  })

  CAND <- reactive({
    if (!is.null(MODEL()$eval)) {
        cand <- MODEL()$eval
    } else if (!is.null(MODEL()$fine)) {
        fine <- MODEL()$fine
        cand <- fine[order(cand, -iter), .SD[1], cand][order(-aL)]
    }  else {
        cand <- data.table(cand=1, p=0.5, P=2)
    }
    return(cand)
  })

  FIT <- reactive({
      purity <- input.model()$purity
      ploidy <- input.model()$ploidy
      fit <- cnvex::modelGetOne(MCNV(), SEG(), purity = purity, ploidy = ploidy, OPTS())
      return(fit)
  })

  PD.QC <- reactive({
      cnv <- CNV()
      seg <- SEG()
      pd <- plotData(cnv$tile, cnv$var, seg, NULL, GENES(), GOBJ(), OPTS())
      return(pd)
  })

  MCNV <- reactive({
      opts <- OPTS()
      rcnv <- CNV()
      sample <- ifelse(input$`cnvex-model`=="somatic", "tumor", "normal")
      mcnv <- cnvex::modelCnv(sample, rcnv, NULL, GOBJ(), opts)
  })

  DIGEST <- reactive({
    if (input$`model-adjust`) {
      model <- "custom"
      purity <- input.model()$purity
      ploidy <- input.model()$ploidy
    } else {
      model <- input$model
      cand <- CAND()
      purity <- cand[cand == as.integer(input$model), p]
      ploidy <- cand[cand == as.integer(input$model), P]
    }
    digest <- modelDigest(MCNV(), SEG(), FIT(), purity, ploidy, eval=list(p=purity, P=ploidy), opts=OPTS(), log=list(user=session$user, notes=input$`checkout-notes`))
    return(digest)

  })

  PD <- reactive({
    dig <- DIGEST()
    tile <- dig$tile
    var <- dig$var
    seg <- dig$seg
    fit <- dig$fit
    pd <- plotData(tile, var, seg, fit, GENES(), GOBJ(), OPTS())
    return(pd)
  })

  STATS <- reactive({
    abs.range <- input.model()$abs.range
    p <- input.model()$purity
    P <- input.model()$ploidy
    D <- (P * p) + 2 * (1 - p)
    C <- seq(0, abs.range, by = 1)
    Clr <- data.table(C = factor(C), lr = log2((p * C + (1 - p) * 2) / D))
    ymin <- log2((p * min(C) + (1 - p) * 2) / D) - 0.1
    ymax <- log2((p * max(C) + (1 - p) * 2) / D) + 0.1
    beta.grid <- as.data.table(expand.grid(M = 0:abs.range, C = 0:abs.range))[M <= C]
    beta.grid[, lr := log2((p * C + (1 - p) * 2) / D)]
    beta.grid[, AF := (p * M + 1 * (1 - p)) / (p * C + 2 * (1 - p))]
    stats <- list(p = p, D = D, C = C, Clr = Clr, ymin = ymin, ymax = ymax, beta.grid = beta.grid)
    return(stats)
  })

  ######## RENDER

  #### UI updates

  observe({
    if (!is.null(CNV())) {
        updateSelectInput(session, "landscape-sel-chr", choices = c("all",CHR()))
        updateSelectInput(session, "zoom-sel-chr", choices = CHR())
        updateSelectInput(session, "abs-sel-chr", choices = c("all", CHR()))        
    }
  })

  observe({
    l <- DIR()
    w <- sapply(l, function(x) {
      is.na(x$cnv.somatic.digest) & is.na(x$cnv.germline.digest)
    })
    names(w) <- names(l)
    if(input$`hide-done`) {
      l <- l[w]
    }
    updateSelectInput(session, "cnvex-run", choices = names(l))
  })

  observeEvent(input$`cnvex-load`, {
      ## restore default settings
      updateTextInput(session, "checkout-notes", value = "")
  })

  observeEvent(input$`cnvex-load-next`, {
    dir.lists <- DIR()
    print("cnvex-load-next")
    run.names <- names(dir.lists)
    print("cnvex-load-next")
    idx <- which(run.names == input$`cnvex-run`)
    print("cnvex-load-next")
    new.run <- run.names[idx + 1]
    print("cnvex-load-next")
    updateSelectInput(session, "cnvex-run", selected = new.run)
    print("cnvex-load-next")
  })

  #### Data Load

  output$`cnvex-files` <- renderText({
    cnvex.run <- DIR()[[input$`cnvex-run`]]
    if (!is.null(cnvex.run)) {
      sprintf("CNV:\t%s\nOPTS:\t%s\nSOMATIC  MODEL:\t%s\nGERMLINE MODEL:\t%s\nSOMATIC  SEGMENT:\t%s\nGERMLINE SEGMENT:\t%s\nSOMATIC  DIGEST:\t%s\nGERMLINE DIGEST:\t%s",
              cnvex.run$cnv, cnvex.run$cnv.opts,
              cnvex.run$cnv.somatic.model, cnvex.run$cnv.germline.model,
              cnvex.run$cnv.somatic.segment, cnvex.run$cnv.germline.segment,
              cnvex.run$cnv.somatic.digest, cnvex.run$cnv.germline.digest)
    } else {
      sprintf("No CNVEX Run Found.")
    }
  })

  output$`cnvex-stat` <- renderText({
    cnv <- CNV()
    model <- MODEL()
    ##
    tile <- cnv$tile
    n.tgt <- length(tile[tile$target])
    n.off <- length(tile[!tile$target])
    ##
    var <- cnv$var
    n.var <- length(var)
    ##
    sex <- cnv$sex
    sprintf("CNVEX File: %s\nModel File:%s\nOn-target tiles: %s\nOff-target tile: %s\nNumber of SNPs: %s\nSex: %s",
            cnv$file_name, model$file_name, n.tgt, n.off, n.var, sex)
  })

  #### QC

  output$`qc-gc` <- renderPlot({
    plt <- plotGC(PD.QC(), input$`qc-maxn`)
    print(plt)
  })

  output$`qc-depth` <- renderPlot({
    covx <- PD.QC()$cov[, .(target, t.cov.raw, n.cov.raw)]
    covx <- covx[1:min(nrow(covx), input$`qc-maxn`)]
    covm <- melt(covx, id.vars = "target")
    ymax <- quantile(covm$value, 0.995, na.rm=TRUE)
    plt <- ggplot(covm) +
      aes(x = variable, y = value, fill = variable) +
      geom_boxplot() +
      scale_fill_brewer("Dark1", guide = "none") +
      scale_x_discrete(labels = c("tumor", "normal")) +
      coord_cartesian(ylim = c(0, ymax)) +
      theme_pubr() +
      xlab(NULL) +
      ylab("coverage")
    return(plt)
  })

  output$`qc-lr-grid` <- renderPlot({
    fine <- MODEL()$fine
    cand <- fine[order(cand), .SD[1], cand]
    plotGrid(MODEL()$grid, cand=cand, opts = OPTS())
  })

  #### Landscape

  LANDSCAPE.PLOT <- reactive({
      if (input.landscape()$`landscape-sel-chr` == "all") {
          sel.chr <- NULL
      } else {
          sel.chr <- input.landscape()$`landscape-sel-chr`
      }
    pd <- PD()
    
      cov.plt <- plotLogRatioScatterMore(PD(),
                              lr.range = input.landscape()$`landscape-lr-range-slider`,
                              sel.col = input.landscape()$`landscape-sel-col`,
                              sel.chr = sel.chr,
                              sel.data = "tile"
      )
      baf.plt <- plotBaf(PD(),
                         baf.range = input.landscape()$`landscape-baf-range-slider`,
                         sel.col = input.landscape()$`landscape-sel-col`,
                         sel.chr = sel.chr
      )
      plt <- egg::ggarrange(cov.plt, baf.plt, ncol = 1)
  })

  output$`landscape-plot` <- renderPlot({
    plt <- LANDSCAPE.PLOT()
    print(plt)
  })

  output$`landscape-save-plot` <- downloadHandler(
      ## TODO: add sample identifier
      filename = function() {
          sprintf("%s-landscape.pdf", "PATIENT")
          },
      content = function(f) {
          ggsave(LANDSCAPE.PLOT(), file=f, height=6, width=18)
          }
  )

  #### Zoom

  output$`zoom-plot` <- renderPlotly({
      sel.chr <- input$`zoom-sel-chr`
      cov.plt <- plotLogRatio(PD(), lr.range = c(-4,4), sel.col = "CK", sel.data = "tile", max.point = Inf, sel.chr=sel.chr)
      cov.pltly <- ggplotly(hide_legend(cov.plt), tooltip="text")
      baf.plt <- plotBaf(PD(), baf.range = c(0, 1), sel.col = "CK", max.point = Inf, sel.chr = sel.chr)
      baf.pltly <- ggplotly(hide_legend(baf.plt), tooltip="text")
      # This almost works, but has the following problems all plotly bugs:
      # 1. gray rectangles are not displayed (not relevant if used in zoom-chromosome mode)
      # 2. tooltip displays all aesthetics: x, y, color, text although only text is selected
      # this has something to do with only geom_point using it. Probably best to re-write these
      # plotting functions in plot_ly not ggplotly.
      # Warning: Ignoring unknown aesthetics: text
      subplot(cov.pltly, baf.pltly, shareX = TRUE, nrows = 2, titleY = TRUE)
  })

  #### Tune Ploidy

  observe({
    cand <- CAND()
    updateSelectInput(session, "model", choices = cand$cand)
  })

  observe({
    cand <- CAND()
    purity <- cand[cand == input$model, p]
    updateSliderInput(session, "purity", value = purity)
    ploidy <- cand[cand == input$model, P]
    updateSliderInput(session, "ploidy", value = ploidy)
  })

  output$`model-table` <- renderDataTable({
    DT::datatable(CAND()[, .(cand, p, P, aD, pdel, aL, hL, tL, bL)], rownames = FALSE) %>%
      DT::formatRound(columns = c("p", "P", "aD", "pdel", "aL", "hL", "tL", "bL"), digits = 3)
  })

  output$`cnv-absolute-plot` <- renderPlot({
    if (input$`abs-sel-chr` == "all") {
      sel.chr <- NULL
    } else {
      sel.chr <- input$`abs-sel-chr`
    }
    purity <- input.model()$purity
    ploidy <- input.model()$ploidy
    abs.range <- input.model()$abs.range

    if(input$`abs-sel-col` == "segment" & input$`abs-sel-data` == "segment" & input$`cnvex-model`=="somatic") {
      ## just the speed up for this setting
      tile <- CNV()$tile
      seg  <- SEG()
      opts <- OPTS()
      sample <- ifelse(input$`cnvex-model`=="somatic", "tumor", "normal")
      cov.plt <- plotLogRatioQuick(tile, seg, sample,
                                   purity = purity, ploidy = ploidy, C.range = c(0, abs.range),
                                   sel.chr = sel.chr, opts = opts)
    } else {
      if (input$`abs-sel-chr` == "all") {
        sel.chr <- NULL
      } else {
        sel.chr <- input$`abs-sel-chr`
      }
      pd <- PD()
      cov.plt <- plotLogRatioScatterMore(pd,
                              purity = purity, ploidy = ploidy, C.range = c(0, abs.range),
                              sel.col = input$`abs-sel-col`, sel.data = input$`abs-sel-data`, sel.chr = sel.chr
      )
    }
    print(cov.plt)
  })

  output$`abs-baf-grid` <- renderPlot({
    mcnv <- MCNV()
    seg <- SEG()
    opts <- OPTS()
    data <- cnvex:::.opt.data(mcnv = mcnv, seg = seg, opts = opts)
    lr <- data$lr
    lr.seg <- lr[, .(lr = mean(lr, na.rm = TRUE)), seg]
    setkey(lr.seg, seg)
    af <- data$af
    var <- MCNV()$var[af$idx]
    af$lr <- MCNV()$tile[findOverlaps(var, MCNV()$tile + OPTS()$tile.shoulder, select = "first")]$lr
    af$chr <- as.character(seqnames(var))
    af$seg.lr <- lr.seg[af$seg]$lr
    stats <- STATS()
    if (input$`abs-sel-chr` != "all") {
        sel.chr <- input$`abs-sel-chr`
        af <- af[chr %in% sel.chr]
    }
    plt <- with(stats, {
      ggplot(af) +
        aes(x = AF, y = lr) +
        geom_scattermore(aes(y = seg.lr), color = "blue", size = 0.5, alpha = 0.25, na.rm = TRUE) +
        stat_binhex(aes(fill = log(after_stat(count))), bins = input$`abs-hex-bins`, alpha = 0.75) +
        coord_cartesian(ylim = c(ymin, ymax)) +
        ylab("Absolute Copy Number") +
        xlab("Allele Frequency") +
        scale_y_continuous(breaks = Clr$lr, labels = Clr$C) +
        scale_fill_gradient(low = "white", high = "black", guide = "none") +
        geom_point(data = beta.grid, size = 2, color = "red", alpha = 0.75, na.rm = TRUE) +
        theme_pubr()
    })
    return(plt)
  })

  ##### Checkout

  output$`checkout-digest` <- renderText({
    paste(capture.output(str(DIGEST(), max.level = 1)), collapse = "\n")
  })

  observeEvent(input$`digest-checkout`, {
    cnvex.run <- DIR()[[input$`cnvex-run`]]
    out <- str_replace(cnvex.run$cnv, ".rds$", sprintf("-%s-digest.rds", input$`cnvex-model`))
    saveRDS(DIGEST(), out)
    print(out)
  })

})
