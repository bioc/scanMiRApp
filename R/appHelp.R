.getAppIntro <- function(){
  data.frame(
    element=c(NA, "#menuCollection", "#menuSearch", "#menuMiRNA", "#mirnainput",
              "#affinitybox", "#targetsbox", NA, "#subjectsbox", "#txbox",
              "#mirnasbox", "#scanBtn", "#hitsbox", NA),
    intro=c(
      "This introduction will walk you through the usage of the scanMiR app.",
      "First, use this tab to select the miRNA collection relative to your 
      species of interest",
      "In this section, you can scan and visualize specific transcripts or 
      sequence for miRNA binding sites",
      "In this section, you can instead start from a specific miRNA, visualize 
      its binding profile, and get its top predicted targets.<br/><br/>
      We'll start by looking at that part.",
      "Use this dropdown to select the miRNA you're interested in.<br/><br/>
      Note that you don't have to scroll down until you find it - you can 
      simply erase what's written in the box, type the beginning of the miRNA 
      name, and see the matching options.",
      "This box contains a plot summarizing the binding profile of the miRNA, plotting
the dissociation rate (i.e. affinity) of the top 7mers sequences (with or without
the 'A' at position 1).<br/><br/>If it's not showing, it's because the box is 
collapsed - you can open it by clicking the 'plus' button on the right.",
      "This box contains the predicted repression and number of binding sites for all 
transcripts. You can reorder it anyway you like, and filter or remove any column.<br/><br/>
You can also double-click on one of the row to get the details and visualize 
the individual binding sites - this would automatically move us to the 
transcript-centered section of the app.",
      "Moving on to the <b>transcript-centered</b> section!<br/><br/>
The '<i>Search in gene/sequence</i>' menu on the left provides access to 
sub-sections, the first of which ('subject') is about selecting the sequence we 
want to scan.",
      "There are two ways to select a sequence to scan, available through the 
      tabs on the top right:<br>
      you can either scan an annotated transcript, or a custom sequence that 
      you'll copy-paste. We'll look at the former.",
      "In this sub-tab, you can first select the gene and then select the 
      transcript.<br/>Again, you can delete the content of the box and type the 
      first few letters to get the matching options.<br/>
      An overview of the selected sequence is then shown at the bottom.",
      "The next step is to select the miRNA(s) for which you want to scan 
      binding sites.<br/><br/>
      You can click in the box again and type the beginning of the miRNA(s) you
      wish to add and the matching miRNAs will show up. Alternatively, you can 
      select a set of miRNAs using the buttons.",
      "Finally, you can click on the scan button to launch it!<br/>
      If it's disabled, it's most likely because you didn't select any miRNA 
      and/or sequence.",
      "The putative binding sites are shown in the <i>'hits'</i> tab.<br/>
      You've got two ways to browse sites: you can either visualize them on a 
      plot along the length of the sequence, or through the table.<br/>
      In both cases, you can <i>double click</i> on a given match to view the 
      supplementary 3' alignment on the target sequence.",
      "Well, that's about it for the basic functions! Please let us know if you
      encounter problems!"
    )
  )
}

.getHelpModal <- function(topic){
  switch(topic,
    hitsCol=showModal(modalDialog(title="Columns of the hits table:",
      tags$p("The start and end columns represent the coordinates of the seed
binding site (i.e. corresponding to positions 1-8 of the miRNA) on the
target sequence. The coordinates are based on the beginning of the sequence
scanned (i.e. if you scanned only the UTR region, 1 represents the first
nucleotide of the UTR)."), tags$p("
The 'type' column indicates the type of match in the seed region."), tags$p("
The 'log_kd' column indicates the log of the estimated dissociation rate. A
lower log_kd value is indicative of a stronger affinity."), tags$p("
The 'p3.score' represent the strength of the 3' supplementary alignment. In
most cases, it roughly corresponds to the number of consecutive complementary
bases (exluding the seed region). However, because of tolerance for G-U
bindings or mismatches in large regions of complementarity, it can depart from
this rule of thumb. To visualize this alignment, you can simply double click
on one of the rows of the table."), tags$p("
The 'note' column contains eventual special features of the binding site (e.g.
prediction of TDMD sites)."), easyClose=TRUE, footer=NULL)),
    stypes=showModal(modalDialog(easyClose=TRUE, footer=NULL,
      title="Description of the site types",
      tags$p("MicroRNAs can bind to mRNAs in different kinds, generally defined by the 
             'binding site types'. See the figure below from"), 
      tags$a(href="https://doi.org/10.1016/j.cell.2009.01.002", "Bartel (2009)"),
      tags$p("for a first overview and the corresponding review for detailed 
             information."),
      tags$img(src="Bartel_2009.png"),
      tags$p("In brief, there are canonical and non-canonical binding sites. The canonical
             binding sites are primarily classified by the number of consecutive nucleotides of the 
             extended seed region of a miRNA that are bound to their complements. '7mer-m8'
             indicates that the 8th nucleotide (nt) of the miRNA is bound at that binding site,
             '7mer-A1' that there is an 'A' opposite the first nt of the miRNA, which is beneficial
             at that specific position. Non-canoncial binding sites most commonly include wobble bindings
             or bulged out nucleotides, though often do not cause effective silencing of transcripts."),
      tags$a(href="https://dx.doi.org/10.1126/science.aav1741", 
             "McGeary, Lin et al. (2019)"), tags$p("show that in order to determine the efficacy of a binding
             site, it's important to not only consider the type of the site but also take the flanking 
             di-nucleotides into account. By neglecting this immediate surrounding information one can
             though generally state that 8mer binding sites are the most effective ones, followed by 
             '7mer-m8' and '7mer-A1' binding sites.")
    )),
    manhattan=showModal(modalDialog(easyClose=TRUE, footer=NULL,
      title="Plotting sites along transcript",
      tags$p("This plot shows the detected miRNA binding sites at their positions within the
             selected sequence. The '-log_kd' (= -log(KD)) information on the y-axis is a proxy for the 
             affinity of a specific miRNA to its complementary sequence at that position. The higher 
             a '-log_kd' value is, the more time the Ago-miRNA complex will spend at that site and the 
             more likely it gets that this particular binding site leads to post-transcriptional 
             repression of the transcript. HOWEVER: Due to the nature of the biochemical model,
             this is not true for comparisons of binding sites bound by different miRNAs!! In case the 
             scan has been performed with a single miRNA, a green area highlighting the typical 8mer affinity 
             range (meaning the -log_kd values which on average correspond to the most effective binding sites
             of the chosen miRNA) can give further information on the strenght of individual binding sites."),
      tags$p("If the scan has been performed with several miRNAs, those with the highest affinity values
             (see also the order in the hits-table) will be displayed by default in the plot. The number
             of shown miRNAs can be altered in the 'max number of miRNA' input"),
      tags$p("A 'log_kd' of -1 (indicated by the dashed red line) indicates background binding.")
    )),
    modalDialog(title=topic, "No help currently available for this topic.")
  )
}
