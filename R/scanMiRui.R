#' scanMiRui
#' 
#' UI for the scanMiR app.
#' 
#' @import shiny shinydashboard
#' @importFrom plotly plotlyOutput ggplotly
#' @importFrom shinycssloaders withSpinner
#' @export
#' @examples
#' ui <- scanMiRui()
scanMiRui <- function(){

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
        href="https://github.com/ETHZ-INS/scanMiR", target="_blank",
        tags$img(src="https://raw.githubusercontent.com/ETHZ-INS/scanMiR/master/inst/docs/sticker.svg"),
        style="display: block; position: fixed; bottom: 5px; left: 20px;")
    ),
    ## Body Content
    dashboardBody(
      tabItems(
        tabItem( tabName = "tab_collection",
          tags$h3("Select or upload a miRNA collection"), tags$br(),
          tabBox(id="collection_type", width=12,
                 tabPanel( title="Pre-built", value="prebuilt",
                   selectInput("mirlist", "miRNA collection", choices=c())),
                 tabPanel(title="Upload", value="upload", 
                          tags$p("Not yet implemented."))
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
            column(12, selectizeInput("mirnas", label="Selected miRNAs:", 
                                      choices=c(), multiple=TRUE) ),
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
                        min = -4, max = 0, value = -0.5, step = 0.25)
          ),
          box(width=6, checkboxInput("keepmatchseq", "Keep matched sequence"),
              checkboxInput("scanNonCanonical", 
                            "Search also for non-canonical sites", value=TRUE))
        ),
        tabItem(tabName="tab_hits",
          column(2, actionButton("scan", "Scan!", icon = icon("search"))),
          column(10, tags$h5(textOutput("scan_target"))),
          box( width=12, title="Manhattan Plot", collapsible=TRUE, 
               collapsed=TRUE, withSpinner(plotlyOutput("manhattan")),
            column(6, numericInput("manhattan_n", "Max number of miRNAs", 
                                   value=10, min=1, max=50)),
            column(6, checkboxInput("manhattan_ordinal", "Ordinal position", 
                                    value=FALSE))
          ),
          box(width=12, title="Table", collapsible=TRUE,
            withSpinner(DTOutput("hits_table")),
            downloadLink('dl_hits', label = "Download all")
          ),
          box(width=12, id="box_match", title="Match alignment", 
            collapsible=TRUE, collapsed=TRUE,
            textOutput("alignment_header"),
            withSpinner(verbatimTextOutput("alignment"))
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
        tabItem(tabName = "tab_about")
      ),
      tags$head(tags$style(HTML('.content-wrapper { overflow: auto; }')))
    )
  )
}
