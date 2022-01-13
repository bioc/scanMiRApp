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
    modalDialog(title=topic, "No help currently available for this topic.")
  )
}
