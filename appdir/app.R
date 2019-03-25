

solve <- function(input) {
  response <- httr::GET("http://128.174.125.122:8080/cutfree",
                        query=list(unsafe_starting_oligo=input$seq,
                                   unsafe_sites=input$sites,
                                   unsafe_min_blocks=input$min_blocks,
                                   unsafe_re_randomize=input$re_randomize,
                                   raw=TRUE))
  paste(httr::content(response, "text"))
}

server <- function(input, output) {
  result <- shiny::eventReactive(input$goButton, {solve(input)})
  output$result <- shiny::renderText({result()})

  output$codes <- renderTable(tibble::frame_data(~code, ~base,
                                                 "N", "A,C,G,T",
                                                 "D", "A,G,T",
                                                 "V", "A,C,G",
                                                 "B", "C,G,T",
                                                 "H", "A,C,T",
                                                 "R", "A,G",
                                                 "Y", "C,T",
                                                 "M", "A,C",
                                                 "K", "G,T",
                                                 "S", "C,G",
                                                 "W", "A,T"))

  shiny::observeEvent(input$help, {
    shiny::showModal(modalDialog(
      title = "More Information",
      shiny::HTML("<p><b>Starting Sequence:</b> Starting DNA sequence that should be blocked from containing ",
                  "restriction sites. To generate a set of barcodes with the highest diversity, start with a ",
                  "string of N's the length of your oligo.",
                  "In addition to IUB codes, mixed bases can be specified as comma-separated lists ",
                  "surrounded by braces. For example, 'TCCRAC' can be input as 'TCC{A,G}AC'.",
                  "When braces are used for input, the output oligo is also displayed with braces.</p>",
                  "<p>Sequences with up to 40 bases are accepted in the online tool.</p>",
                  "<p><b>Sites to Block:</b> Sequences to block from the oligo pools. Separate multiple sequences by commas.</p>",
                  "<p><b>Minimum # of Blocks:</b> Minimum number of blocks at each site, i.e. the minimum number of changes that need to be made for a cut site to appear anywhere in the oligo.</p>",
                  "<p><b>Increase solution diversity:</b> If this box is checked, re-run the optimization to randomize the codes while maintaining the same number of oligos.</p>",
                  "<hr><p>A new set oligo is generated each time the 'Make Randomer' button is pressed.</p>"),
      easyClose = TRUE,
      footer = shiny::modalButton("Close")
    ))
  })

}

ui <- shiny::fluidPage(
  shiny::titlePanel("CutFree: Randomized DNA barcodes without restriction sites"),
  shiny::fluidRow(
    shiny::column(4, shiny::wellPanel(
      shiny::textInput("seq", "Starting Sequence:", "NNNNNNNNNNNNNNNNNNNN"),
      shiny::textInput("sites", "Sites to Block (comma separated):", "GGTCTC, GGCCGG"),
      shiny::numericInput("min_blocks", "Minimum # of Blocks:", min=1, max=3, value=1, step=1, width='50%'),
      shiny::checkboxInput("re_randomize", "Increase solution diversity", TRUE),
      shiny::actionButton("goButton", "Make Randomer"),
      shiny::actionLink("help", "More Information"),
      shiny::HTML("<hr><center>"),
      shiny::tableOutput("codes"),
      shiny::HTML("</center>")
    )),

    shiny::column(6,
      shiny::verbatimTextOutput(outputId = "result")
    )
  )

)

shiny::shinyApp(ui=ui, server=server)
