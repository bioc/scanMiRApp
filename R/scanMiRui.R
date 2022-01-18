#' scanMiRui
#'
#' UI for the scanMiR app.
#'
#' @import shiny shinydashboard
#' @importFrom plotly plotlyOutput ggplotly
#' @importFrom shinycssloaders withSpinner
#' @importFrom rintrojs introjsUI introBox
#' @importFrom shinyjqui jqui_resizable
#' @importFrom waiter use_waiter waiter_show_on_load spin_1
#' @return A shiny ui
#' @export
#' @examples
#' ui <- scanMiRui()
scanMiRui <- function(){
  scanMiRlogo <- paste0("https://raw.git",
                        "hubusercontent.com/ETHZ-INS/scanMiR/",
                        "master/inst/docs/sticker.svg")
  mer12scheme <- paste0("https://raw.git",
                        "hubusercontent.com/ETHZ-INS/scanMiR/",
                        "master/inst/docs/12mer.png")
  dashboardPage( skin="black",

    dashboardHeader(title = "scanMiRApp", titleWidth = "300px",
                    tags$li(class="dropdown",
                            actionLink("helpBtn", label="Quick start",
                                       icon=icon("question")))),

    ## Sidebar content
    dashboardSidebar( width = "200px",
      sidebarMenu(id="main_tabs",
        menuItem("About", tabName="tab_about"),
        menuItemOutput("menuCollection"),
        menuItem("Search in\ngene/sequence", id="menuSearch",
          menuSubItem("subject", "tab_subject"),
          menuItemOutput("menuMirnas"),
          menuSubItem("options", "tab_options"),
          menuSubItem("hits", "tab_hits"),
          startExpanded=TRUE
        ),
        menuItem(tags$span(id="menuMiRNA", "miRNA-based"), tabName="tab_mirna")
      ),
      tags$a(
        href=paste0("https://git","hub.com/ETHZ-INS/scanMiR"), target="_blank",
        tags$img(src=scanMiRlogo),
        style="display: block; position: fixed; bottom: 5px; left: 20px;")
    ),
    ## Body Content
    dashboardBody(
      introjsUI(),
      use_waiter(),
      waiter_show_on_load(html = tagList(
        tags$img(src=scanMiRlogo), tags$br(), tags$br(),
        tags$h3("Please wait while the application is initialized..."),
        spin_1()
      )),
      tabItems(
        tabItem( tabName = "tab_collection",
          tags$h3("Select a miRNA collection"), tags$br(),
          column(6,
            tabBox(id="collection_type", width=12,
              tabPanel( title="Pre-built", value="prebuilt",
                        actionButton(inputId="help_collections", label="",
                          style="float:right;", icon=icon("question-circle")),
                     selectInput("mirlist", "miRNA collection", choices=c()))
            )),
          column(6, valueBoxOutput("selected_collection", width=12)),
          box(width=12, title="Extra annotation information", 
              collapsible=TRUE, collapsed=TRUE,
              withSpinner(verbatimTextOutput("collection_summary")))
        ),
        tabItem(tabName="tab_subject",
          tags$div(id="subjectsbox", style="min-height: 200px;",
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
              tags$div(id="txbox",
              tags$div(selectizeInput("gene", "Gene", choices=c()),
                       style="float: left; padding-right: 10px;"),
              tags$br(), tags$p(style="text-align: left;", "Type in the first ",
                            "few letters of a gene name to populate the list"),
              htmlOutput("gene_link"),
              tags$div(style="clear: left;"),
              selectizeInput("transcript", "Transcript", choices=c()),
              introBox(selectInput("seqFeature", "Part to scan", 
                                   choices=c("3' UTR only", "CDS+UTR", "whole transcript")),
                       data.hint="You can decide here whether you want to scan 
                       just the 3' UTR or also include the coding sequence,
                       or scan the whole transcript."),
              withSpinner(tableOutput("tx_overview")))
            )
          )
        )),
        tabItem(tabName="tab_mirnas",
          box(id="mirnasbox", width=12,
            column(12, selectizeInput("mirnas", choices=c(), multiple=TRUE,
                                label="Selected miRNAs: (type in to filter)")),
            column(8, actionButton("mirnas_vert",
                          "Select miRNAs conserved across vertebrates"), 
             tags$br(), actionButton("mirnas_mammals",
                          "Select miRNAs conserved across mammals"),
             tags$br(), actionButton("mirnas_confident",
                                     "Select confidently annotated miRNAs")
            ),
            column(4, "The conservation status of miRNAs is obtained from",
                   tags$a(href="http://www.targetscan.org/", target="_blank",
                          "TargetScan")),
            column(8,
              checkboxInput("mirnas_all", "Search for all miRNAs")
            ),
            column(4, actionButton("mirnas_clear",
                                   "Clear all selected miRNAs"))
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
          fluidRow(
            column(2, tags$div(id="scanBtnBox", 
              style="display: inline-block; min-height: 30px; min-width: 100px;",
              uiOutput("scanBtn"), HTML("&nbsp;"))),
            column(10, tags$h5(textOutput("scan_target")))
          ),
          tags$div(id="hitsbox", style="min-height: 300px;",
          box(width=12, collapsible=TRUE, collapsed=TRUE,
              title="Display sites along transcript",
              actionButton("manhattanHelp", label="", style="float: right;", 
                           icon=icon("question-circle")),
              tags$p("Hover on points to view details, and click to ",
                      "visualize the alignment on the target sequence. You may
                      also select miRNAs to show/hide by clicking on the legend."),
              withSpinner(plotlyOutput("manhattan")),
            column(6, numericInput("manhattan_n", "Number of top miRNAs displayed",
                                   value=10, min=1, max=50)),
            column(6, checkboxInput("manhattan_ordinal", "Ordinal position",
                                    value=FALSE))
          ),
          box(width=12, title="Table", collapsible=TRUE,
            fluidRow(
              column(6,
                tags$p("Double click on a row to visualize the alignment",
                        " on the target sequence.")),
              column(6, style="text-align: right;", 
                   actionButton("colHelp","What are those columns?",
                                icon=icon("question-circle")),
                   actionButton("stypeHelp","Site types",
                                icon=icon("question-circle"))
              )
            ),
            withSpinner(DTOutput("hits_table")),
            tags$p(downloadLink('dl_hits', label = "Download all"))
            
          )),
          box(width=12, title="Cached hits", collapsible=TRUE, collapsed=TRUE,
            textOutput("cache.info"),
            uiOutput("cached.results"),
            actionButton("loadCache", "Load"),
            actionButton("deleteCache", "Delete")
          )
        ),

        # miRNA-based
        tabItem(tabName="tab_mirna",
          fluidRow(
            column(5, tags$div(id="mirnainput",
                        selectizeInput("mirna", "miRNA", choices=c()))),
            column(4, tags$strong("TargetScan conservation status"), 
                  textOutput("modconservation")),
            column(3, htmlOutput("mirbase_link"))
          ),
          tags$div(id="affinitybox", style="min-height: 40px;",
            box(width=12, title="Affinity plot", collapsible=TRUE, 
                collapsed=TRUE,
                withSpinner(jqui_resizable(plotOutput("modplot")))
            )
          ),
          tags$div(id="targetsbox", style="min-height: 100px;",
            box(width=12, title="Targets", collapsible=TRUE, 
                uiOutput("targets_ui"))
          )
        ),
        tabItem(tabName = "tab_about",
## TAB ABOUT
    box(width=12, title="miRNA binding prediction",
      tags$div(style="font-size: 110%;",
        tags$img(src=mer12scheme, style="float: right; margin: 20px; width: 25%;"), 
        tags$p("micro-RNAs (miRNAs) regulate transcripts stability and 
          translation by targeting the Argonaute complex to RNAs exhibiting a 
          partial complementarity with, in particular, the miRNA's seed 
          sequence. For a general review on miRNA biology and targeting, see
               ", tags$a(href="https://doi.org/10.1016/j.cell.2018.03.006", 
                         "Bartel (2018)",.noWS = "outside"),".", 
               tags$a(href="https://dx.doi.org/10.1126/science.aav1741", 
            "McGeary, Lin et al. (2019)"), " have additionally shown the 
      relevance of flanking nucleotides by empirically measuring the affinity 
      (specifically the dissociation rate constant, KD) of a set of miRNAs to 
      random 12-mer sequences, and then computationally extrapolating to other 
      miRNAs). The ",
      tags$a( href=paste0("https://git","hub.com/ETHZ-INS/scanMiR"),
              target="_blank", "scanMiR package",.noWS = "outside"), ", to which this app is an 
      interface, builds on this work to offer a powerful and flexible binding 
      and repression framework."))
    ),
    box(width=12, title="Getting started",
        tags$div( style="font-size: 110%;",
           tags$p("For a quick tour of the app, ", actionLink("helpLink", "click here",.noWS = "outside"),"."),
           tags$p("There are two main ways to use scanMiRApp:"),
           tags$br(), tags$h4("Transcript-centered:"),
           tags$p("In the 'Search in gene/sequence' menu, you'll be able to ",
                  "scan the sequence of a given transcript for binding sites ",
                  "of (sets of) miRNA(s)."),
           tags$br(), tags$h4("miRNA-centered:"),
           tags$p("In the 'miRNA-based' tab on the left, you'll be able to ",
                  "visualize information relative to a selected miRNA, ",
                  "including for instance its general binding profile and its",
                  "top targets.")
        )
    ),
    tags$div(style="text-align: right; font-weight: normal;",
      tags$p("The shiny app was developed by Pierre-Luc Germain and ",
           "Michael Soutschek in the context of broader research in the ",
           tags$a( href="http://schrattlab.ethz.ch", "Schratt lab",
                   target="_blank",.noWS = "outside"), ".",
           tags$br(), "Bugs reports and feature requests are welcome on the ",
           tags$a( href=paste0("https://git","hub.com/ETHZ-INS/scanMiRApp/issues"),
                   target="_blank", "github repository",.noWS = "outside"), "."),
           textOutput("pkgVersions")
      )
## END TAB ABOUT
        )
      ),
      tags$head(tags$style(HTML('.content-wrapper { overflow: auto; }')))
    )
  )
}
