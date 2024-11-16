ui <- navbarPage(
  title = "CanIQ",
  selected = "Somatic Mutations",
  collapsible = TRUE,
  theme = bslib::bs_theme(),
  tabPanel(
    title = "Somatic Mutations",
    grid_container(
      layout = c(
        "somatic_filters",
        "somatic_table"
      ),
      row_sizes = c(
        "auto",
        "1fr"
      ),
      col_sizes = c(
        "1fr"
      ),
      gap_size = "10px",
      grid_card(
        area = "somatic_filters",
        card_header("Somatic Variant Selection"),
        card_body(
          grid_container(
            gap_size = "5px",
            layout = c(
              "somatic_filter somatic_select somatic_limit"
            ),
            row_sizes = c(
              "auto"
            ),
            col_sizes = c(
              "1fr",
              "5fr",
              "1fr"
            ),
            grid_card(
              area = "somatic_filter",
              card_header("Filters"),
              card_body_fill(
                checkboxInput(inputId = "somatic_filter_pass", "Passing quality", value=TRUE),
                checkboxInput(inputId = "somatic_filter_impact", "Protein altering", value=TRUE),
              )
            ),
            grid_card(
              area = "somatic_select",
              card_header("Select"),
              card_body_fill(
                textInput(
                  inputId = "somatic_select_symbol",
                  label = "Gene Symbols",
                  value = "",
                  width="100%"
                ),
                textInput(
                  inputId = "somatic_select_sample",
                  label = "Samples",
                  value = "",
                  width="100%"
                )
              )
            ),
            grid_card(
              area = "somatic_limit",
              card_header("Limit"),
              full_screen = FALSE,
              card_body(
                numericInput(
                  inputId = "somatic_limit_maxrow",
                  label = "Maximum Rows",
                  value = 1000,
                  width = "100%"
                )
              )
            )
          )
        )
      ),
      grid_card(
        area = "somatic_table",
        card_body(
          DTOutput(outputId = "somatic_table", width = "100%")
        )
      )
    )
  )
)
