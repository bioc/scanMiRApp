#' scanMiRui
#'
#' UI for the scanMiR app.
#'
#' @import shiny shinydashboard
#' @importFrom plotly plotlyOutput ggplotly
#' @importFrom shinycssloaders withSpinner
#' @importFrom waiter use_waiter waiter_show_on_load spin_1
#' @return A shiny ui
#' @export
#' @examples
#' ui <- scanMiRui()
scanMiRui <- function(){
  scanMiRlogo <- paste0("https://raw.git",
                        "hubusercontent.com/ETHZ-INS/scanMiR/",
                        "master/inst/docs/sticker.svg")
  ui <- dashboardPage( skin="black",

    dashboardHeader(title = "scanMiRApp", titleWidth = "300px"),

    ## Sidebar content
    dashboardSidebar( width = "200px",
      sidebarMenu(id="main_tabs",
        menuItem("Species Collection", tabName="tab_collection"),
        menuItem("Search in\ngene/sequence",
          menuSubItem("subject", "tab_subject"),
          menuSubItem("miRNAs", "tab_mirnas"),
          menuSubItem("options", "tab_options"),
          menuSubItem("hits", "tab_hits")
        ),
        menuItem("miRNA-based", tabName="tab_mirna"),
        menuItem("About", tabName="tab_about")
      ),
      tags$a(
        href=paste0("https://git","hub.com/ETHZ-INS/scanMiR"), target="_blank",
        tags$img(src=scanMiRlogo),
        style="display: block; position: fixed; bottom: 5px; left: 20px;")
    ),
    ## Body Content
    dashboardBody(
      use_waiter(spinners = 3),
      waiter_show_on_load(html = tagList(
        tags$img(src=scanMiRlogo), tags$br(), tags$br(),
        tags$h3("Please wait while the application is initialized..."),
        spin_1()
      )),
      tabItems(
        tabItem( tabName = "tab_collection",
          tags$h3("Select a miRNA collection"), tags$br(),
          tabBox(id="collection_type", width=12,
                 tabPanel( title="Pre-built", value="prebuilt",
                   selectInput("mirlist", "miRNA collection", choices=c()))#,
#                 tabPanel(title="Upload", value="upload",
#                          tags$p("Not yet implemented."))
          ),
          box(width=12, withSpinner(verbatimTextOutput("collection_summary")))
        ),
        tabItem(tabName = "tab_subject",
          tabBox(
            title="Search bindings in", side="right",
            id="subjet_type", width=12, selected = "transcript",
            tabPanel(
              title="Custom sequence", value="custom",
              textAreaInput(
                inputId="customseq", label="Sequence (RNA or DNA)",
                placeholder=paste0("Paste in here a sequence in which you want",
                                   " to search for binding sites"),
                width = "100%", height = "250px"),
              checkboxInput(inputId="circular", label="circularize"),
              verbatimTextOutput("custom_info"),
              actionButton("rndseq", "Generate random sequence")
            ),
            tabPanel(
              title="Transcript", value="transcript",
              tags$div(style="float: right; padding: 20px; width: 30%; min-width: 300px;",
                       "At the moment, only the sequences from protein-coding
                       transcripts can be queried in that way. For non-coding
                       transcripts, you'll have to enter the sequence yourself
                       (see 'custom sequence' tab above).", tags$br(),
                       "The next scanMiRApp release in the coming weeks will
                       include also non-coding transcripts!"),
              selectizeInput("annotation", "Genome & Annotation", choices=c()),
              tags$div(selectizeInput("gene", "Gene", choices=c()),
                       style="float: left; padding-right: 10px;"),
              tags$br(), tags$p(style="text-align: left;", "Type in the first ",
                            "few letters of a gene name to populate the list"),
              htmlOutput("gene_link"),
              tags$div(style="clear: left;"),
              selectizeInput("transcript", "Transcript", choices=c()),
              checkboxInput("utr_only", "UTR only", value = TRUE),
              withSpinner(tableOutput("tx_overview"))
            )
          )
        ),
        tabItem(tabName="tab_mirnas",
          box( width=12,
            column(12, selectizeInput("mirnas", choices=c(), multiple=TRUE,
                          label="Selected miRNAs: (type in to filter)")),
            column(4,
              checkboxInput("mirnas_all", "Search for all miRNAs"),
              actionButton("mirnas_confident",
                           "Select confidently annotated miRNAs"),
              actionButton("mirnas_mammals",
                           "Select miRNAs conserved across mammals"),
              actionButton("mirnas_vert",
                           "Select miRNAs conserved across vertebrates"),
              actionButton("mirnas_clear", "Clear all selected miRNAs")
            )
          )
        ),
        tabItem(tabName="tab_options",
          box( width=6,
            sliderInput("minDist",
              label=paste("Minimum distance in NT between two binding sites",
                           "of the same miRNA"),
              min = 0, max = 20, value = 7, step = 1),
            sliderInput(inputId="shadow",
              label=paste("Ribosomal shadow at the beginning of the 3'UTR,",
                          "recommended is 15 NT"),
              min = 0, max = 20, value = 15, step = 1),
            sliderInput(inputId="maxLogKd", label = "Maximum log_kd to report",
                        min = -4, max = 0, value = -1, step = 0.25)
          ),
          box(width=6, checkboxInput("keepmatchseq", "Keep matched sequence"),
              checkboxInput("scanNonCanonical",
                            "Search also for non-canonical sites", value=TRUE))
        ),
        tabItem(tabName="tab_hits",
          column(2, uiOutput("scanBtn")),
          column(10, tags$h5(textOutput("scan_target"))),
          box( width=12, collapsible=TRUE, collapsed=TRUE,
               title="Plot along transcript",
               tags$p("Hover on points to view details, and click to ",
                      "visualize the alignment on the target sequence. You may
                      also select miRNAs to show/hide by clicking on the legend."),
               withSpinner(plotlyOutput("manhattan")),
            column(6, numericInput("manhattan_n", "Max number of miRNAs",
                                   value=10, min=1, max=50)),
            column(6, checkboxInput("manhattan_ordinal", "Ordinal position",
                                    value=FALSE))
          ),
          box(width=12, title="Table", collapsible=TRUE,
            withSpinner(DTOutput("hits_table")),
            downloadLink('dl_hits', label = "Download all"),
            tags$p("Double click on a row to visualize the alignment on the ",
                   "target sequence."),
          ),
          box(width=12, title="Cached hits", collapsible=TRUE, collapsed=TRUE,
            textOutput("cache.info"),
            uiOutput("cached.results"),
            actionButton("loadCache", "Load"),
            actionButton("deleteCache", "Delete")
          )
        ),

        # miRNA-based
        tabItem(tabName="tab_mirna",
          column(5, selectizeInput("mirna", "miRNA", choices=c())),
          column(4, tags$strong("Status"), textOutput("modconservation")),
          column(3, htmlOutput("mirbase_link")),
          box(width=12, title="Affinity plot", collapsible=TRUE, collapsed=TRUE,
            withSpinner(plotOutput("modplot")),
            numericInput("modplot_height", "Plot height (px)", value=400,
                         min=200, max=1000, step=50)
          ),
          box(width=12, title="Targets", collapsible=TRUE,
              uiOutput("targets_ui"))
        ),
        tabItem(tabName = "tab_about",
## TAB ABOUT
    box(width=12, title="About",
        tags$p("The scanMiRApp is an interface to the ",
        tags$a( href=paste0("https://git","hub.com/ETHZ-INS/scanMiR"),
                target="_blank", "scanMiR"),
        "package. The shiny app was developed by Pierre-Luc Germain and ",
        "Michael Soutschek in the context of broader research in the ",
        tags$a( href="http://schrattlab.ethz.ch", "Schratt lab",
                target="_blank"), ".",
        tags$br(), "Bugs reports and feature requests are welcome ",
        tags$a( href=paste0("https://git","hub.com/ETHZ-INS/scanMiRApp/issues"),
                target="_blank", "here"),".", tags$br(),
        style="font-size: 110%;")
    ),
    box(width=12, title="Getting started",
        tags$div( style="font-size: 110%;",
           tags$p("There are two main ways to use scanMiRApp:"),
           tags$br(), tags$h4("Transcript-centered:"),
           tags$p("In the 'Search in gene/sequence' menu, you'll be able to ",
                  "scan the sequence of a given transcript for binding sites ",
                  "of (sets of) miRNA(s). To do so:"),
           tags$ol(
             tags$li("Click on 'Search in gene/sequence' to toggle the ",
                     "visibility of sub-menu items"),
             tags$li("In the 'Subject' tab, first select the sequence you want",
                     " to scan. This can either be a custom sequence (using ",
                     "the 'custom sequence' button on the top-right of the ",
                     "'Subject' tab), or selected from ensembl transcripts."),
             tags$li("In the 'miRNAs' tab, select the miRNAs for which you ",
                     "want to find binding sites."),
             tags$li("When ready, go to the 'hits' tab, and click the 'Scan'",
                     "button at the top to launch the search!")
           ),
           tags$br(), tags$h4("miRNA-centered:"),
           tags$p("In the 'miRNA-based' tab on the left, you'll be able to ",
                  "visualize information relative to a selected miRNA, ",
                  "including for instance it's general binding profile and its",
                  "top targets.")
        )
    )
## END TAB ABOUT
        )
      ),
      tags$head(tags$style(HTML('.content-wrapper { overflow: auto; }')))
    )
  )
}
