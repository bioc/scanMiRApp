#' scanMiRui
#'
#' UI for the scanMiR app.
#'
#' @import shiny shinydashboard
#' @importFrom plotly plotlyOutput ggplotly
#' @importFrom shinycssloaders withSpinner
#' @importFrom rintrojs introjsUI introBox
#' @importFrom waiter use_waiter waiter_show_on_load spin_1
#' @return A shiny ui
#' @export
#' @examples
#' ui <- scanMiRui()
scanMiRui <- function(){
  scanMiRlogo <- paste0("https://raw.git",
                        "hubusercontent.com/ETHZ-INS/scanMiR/",
                        "master/inst/docs/sticker.svg")
  dashboardPage( skin="black",

    dashboardHeader(title = "scanMiRApp", titleWidth = "300px",
                    tags$li(class="dropdown",actionLink("helpBtn", label="Help",
                                                      icon=icon("question")))),

    ## Sidebar content
    dashboardSidebar( width = "200px",
      sidebarMenu(id="main_tabs",
        menuItem(introBox("Species Collection", data.step = 1,
                 data.intro = "First, use this tab to select the miRNA collection
                                 relative to your species of interest"),
                 tabName="tab_collection"),
        menuItem(introBox("Search in\ngene/sequence", data.step=2,
                          data.intro="In this section, you can scan and 
visualize specific transcripts or sequence for miRNA binding sites"),
          menuSubItem("subject", "tab_subject"),
          menuSubItem("miRNAs", "tab_mirnas"),
          menuSubItem("options", "tab_options"),
          menuSubItem("hits", "tab_hits"),
          startExpanded=TRUE
        ),
        menuItem(introBox("miRNA-based", data.step=3, data.intro="
In this section, you can instead start from a specific miRNA, visualize its 
binding profile, and get its top predicted targets.<br/><br/>
We'll start by looking at that part."), tabName="tab_mirna"),
        menuItem("About", tabName="tab_about")
      ),
      tags$a(
        href=paste0("https://git","hub.com/ETHZ-INS/scanMiR"), target="_blank",
        tags$img(src=scanMiRlogo),
        style="display: block; position: fixed; bottom: 5px; left: 20px;")
    ),
    ## Body Content
    dashboardBody(
      introjsUI(),
      use_waiter(spinners = 3),
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
                     selectInput("mirlist", "miRNA collection", choices=c()))#,
  #                 tabPanel(title="Upload", value="upload",
  #                          tags$p("Not yet implemented."))
            )),
          column(6, valueBoxOutput("selected_collection", width=12)),
          box(width=12, title="Details", collapsible=TRUE, collapsed=TRUE,
              withSpinner(verbatimTextOutput("collection_summary")))
        ),
        tabItem(tabName = "tab_subject",
          introBox(data.step=7, data.intro="
Moving on to the <b>transcript-centered</b> section!<br/><br/>
Opening the '<i>Search in gene/sequence</i>' menu provides access to sub-sections,
the first of which is about selecting the sequence we want to scan.<br/><br/>
There are two options here, available through the tabs on the top right:<br>
you can either scan an annotated transcript, or a custom sequence that you'll copy-paste.
We'll look at the former.",
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
              introBox(data.step=8, data.intro="
In this sub-tab, you can first select the gene and then
select the transcript.<br/>Again, you can delete the content
of the box and type the first few words to get the matching options.<br/>
An overview of the selected sequence is then shown at the bottom.",
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
          introBox(data.step=9, data.intro="
The next step is to select the miRNA(s) for which you want to scan binding sites.
<br/><br/>You can click in the box again and type the beginning of the miRNA(s) you wish to add
and the matching miRNAs will show up. Alternatively, you can select a set of miRNAs using the buttons.",
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
          column(2, introBox(data.step=10, uiOutput("scanBtn"), data.intro="
Finally, you can click on the scan button to launch it!<br/>
If it's disabled, it's most likely because you didn't select any miRNA and/or sequence.")),
          column(10, tags$h5(textOutput("scan_target"))),
          introBox(data.step=11, data.intro="
The putative binding sites are shown in the <i>'hits'</i> tab.<br/>
You've got two ways to browse sites: you can either visualize
them on a plot along the length of the sequence, or through the table.<br/>
In both cases, you can <i>double click</i> on a given match to view the supplementary
3' alignment on the target sequence.<br/><br/>
Well, that's about it for the basic functions!",
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
            column(6, tags$p("Double click on a row to visualize the alignment",
                             " on the target sequence."),
                   downloadLink('dl_hits', label = "Download all")),
            column(6, style="text-align: right;", 
                   actionButton("colHelp","What are those columns?"))
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
          column(5, introBox(selectizeInput("mirna", "miRNA", choices=c()),
                             data.step=4, data.intro="
Use this dropdown to select the miRNA you're interested in.<br/><br/>Note that you 
don't have to scroll down until you find it - you can simply erase what's written 
in the box, type the beginning of the miRNA name, and see the matching options.")),
          column(4, tags$strong("Status"), textOutput("modconservation")),
          column(3, htmlOutput("mirbase_link")),
          introBox(data.step=5, data.intro="
This box contains a plot summarizing the binding profile of the miRNA, plotting
the dissociation rate (i.e. affinity) of the top 7mers sequences (with or without
the 'A' at position 1).<br/><br/>If it's not showing, it's because the box is 
collapsed - you can open it by clicking the 'plus' button on the right.",
            box(width=12, title="Affinity plot", collapsible=TRUE, collapsed=TRUE,
              withSpinner(plotOutput("modplot")),
              numericInput("modplot_height", "Plot height (px)", value=400,
                           min=200, max=1000, step=50)
            )
          ),
          introBox(data.step=6, data.intro="
This box contains the predicted repression and number of binding sites for all 
transcripts. You can reorder it anyway you like, and filter or remove any column.<br/><br/>
You can also double-click on one of the row to get the details and visualize 
the individual binding sites - this would automatically move us to the 
transcript-centered section of the app.",
                   box(width=12, title="Targets", collapsible=TRUE,
              uiOutput("targets_ui")))
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
        style="font-size: 110%;"),
        tags$p(textOutput("pkgVersions"))
    ),
    box(width=12, title="Getting started",
        tags$div( style="font-size: 110%;",
           tags$p("For a quick tour of the app, ", actionLink("helpLink", "click here")),
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
    )
## END TAB ABOUT
        )
      ),
      tags$head(tags$style(HTML('.content-wrapper { overflow: auto; }')))
    )
  )
}
