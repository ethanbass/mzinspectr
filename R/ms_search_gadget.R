#' Launch MS search gadget for interactive viewing of spectral matches.
#' @name ms_search_gadget
#' @importFrom graphics rasterImage par plot.new
#' @param data An \code{ms_alignment} object.
#' @export

ms_search_gadget <- function(data){
  check_for_pkg("rcdk")
  if (is.null(data$matches)){
    stop("Matches not found! Please run `search_spectra` before calling the mass search gadget.")
  }

  check_for_pkg("shiny")
  check_for_pkg("DT")
  check_for_pkg("plotly")
  check_for_pkg("rcdk")

  # Define the UI
  ui <- shiny::fluidPage(
    # Sidebar with a list of dataframes to select
    shiny::fluidRow(
      shiny::selectInput("df_select", "Select a dataframe:",
                  choices = names(data$matches))
    ),
    shiny::fluidRow(
      DT::DTOutput("matches",height = "350px", width="50%")
    ),
    shiny::fluidRow(
      shiny::column(6,
             plotly::plotlyOutput("spectrum")
      ),

      shiny::column(6,
             shiny::plotOutput("structure")
      )
    )
  )
  # Define the server
  server <- function(input, output) {
    # Render the selected dataframe as a table
    # escape <- c(escape, include[which(!(include %in% colnames(data$matches[[1]])))])
    include <- c("Name", "Formula", "RI", "Comment", "spectral_match", "ri_match", "total_score")
    output$matches <- DT::renderDT({
      escape.idx <- which(!(colnames(data$matches[[input$df_select]]) %in% include))
      escape = colnames(data$matches[[input$df_select]])[escape.idx]
      DT::datatable(data$matches[[input$df_select]],
                    selection=list(mode="single", selected = 1),
                options=list(columnDefs = list(list(visible=FALSE,
                                                    targets=escape)),
                             dom='tip',
                             paging = TRUE,
                             pageLength = 5)) |>
        DT::formatStyle(c(1:ncol(data$matches[[input$df_select]])), fontSize = '10px')
    })

    output$spectrum <- plotly::renderPlotly({

      # Check if a row has been selected in the table
      shiny::req(input$df_select)
      shiny::req(input$matches_rows_selected)

      # spectrum <- plot_spectrum(data, col = input$df_select)
      ms_mirror_plot(x = ms_get_spectrum(data, input$df_select),
                  y = tidy_eispectrum(data$matches[[input$df_select]][input$matches_rows_selected,"Spectra"]),
                  type="plotly")
    })

    # Generate the chemical structure
    output$structure <- shiny::renderPlot({

      # Check if a row has been selected in the table
      shiny::req(input$matches_rows_selected)

      # Extract the SMILES string
      smiles <- data$matches[[input$df_select]][input$matches_rows_selected, "Smiles"]

      if (smiles == "null"){
        plot.new()
      } else {
        # Generate the chemical structure from the SMILES string
        mol = rcdk::parse.smiles(smiles)[[1]]

        # Generate the 2D depiction of the molecule
        # rcdkplot(mol)
        oldpar <- par(no.readonly = TRUE)
        par(mar = c(3,0,0,5))
        png <- rcdk::view.image.2d(molecule=mol,
                                   depictor = rcdk::get.depictor(width=500, height=500, zoom=10))
        plot(NA, NA, xlim=c(1,10), ylim=c(1,10),
             xaxt='n', yaxt='n', xlab='', ylab='', bty="n")
        rasterImage(png, 1, 1, 10, 10)
        par(oldpar)
      }
    })
  }

  # Run the app
  shiny::shinyApp(ui, server)
}
