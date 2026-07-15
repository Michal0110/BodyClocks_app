ui <- fluidPage(
  tags$head(
    includeCSS("www/styles.css"),
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
    tags$link(rel = "shortcut icon", href = "favicon.ico"),
    tags$script(src = "sine-wave.js"), # Sine wave animation
    tags$script(src = "custom-tooltip.js"), # Custom tooltip handling for visNetwork

    # Script to handle completely clearing the checkbox group when new dataset is selected
    tags$script(HTML("
      Shiny.addCustomMessageHandler('clearCheckboxGroup', function(message) {
        // Get the container with all checkboxes and completely clear it
        $('#' + message).find('.shiny-options-group').empty();
      });
    ")),

    # Script to make enrichment checkboxes behave like radio buttons
    tags$script(HTML("
      var enrichmentTermSelector = '#selected_term .shiny-options-group input[type=\"checkbox\"], [id$=\"-selected_term\"] .shiny-options-group input[type=\"checkbox\"]';

      function syncSingleEnrichmentTerm($checkbox) {
        var $group = $checkbox.closest('#selected_term, [id$=\"-selected_term\"]');

        if ($checkbox.prop('checked')) {
          $group.find('.shiny-options-group input[type=\"checkbox\"]').not($checkbox).prop('checked', false);
        }

        return $group;
      }

      $(document).on('click', enrichmentTermSelector, function(e) {
        syncSingleEnrichmentTerm($(this));
      });

      $(document).on('change', enrichmentTermSelector, function(e) {
        syncSingleEnrichmentTerm($(this));
      });
    ")),

),
  
  
  div(class = "bc-app-shell",
  # Wrap the title in a div with class for easier styling
  div(class = "title-container",
    tags$div(id = "sine-wave-container", class = "sine-wave-bg"),
    div(class = "title-text-container",
      h1("BodyClocks", class = "main-title"),
      h3("Circadian Gene Expression and Network Visualization v1.0", class = "subtitle")
    )
  ),
  
  tabsetPanel(
    tabPanel(
      "Circadian Explorer",
      explorerUI(
        id                         = "circ",
        title                      = "Circadian",
        species_choices            = CIRC_SPECIES_CHOICES,
        tissue_choices_by_species  = CIRC_TISSUE_CHOICES_BY_SPECIES,
        dataset_choices_by_species = CIRC_DATASET_CHOICES_BY_SPECIES,
        legend_html                = "www/legend.html",
        show_plot_controls         = TRUE
      )
    ),

    tabPanel(
      "Dataset Comparison",
      comparisonUI(
        id                         = "comp",
        species_choices            = CIRC_SPECIES_CHOICES,
        tissue_choices_by_species  = CIRC_TISSUE_CHOICES_BY_SPECIES,
        dataset_choices_by_species = CIRC_DATASET_CHOICES_BY_SPECIES
      )
    )
  ),
  ),
  
  # Duplicate tooltips for each module instance (namespaced ids)
  bsTooltip(id = "circ-tissue_help", title = "Select tissue type first.", placement = "right", trigger = "click"),
  bsTooltip(id = "circ-dataset_help", title = "Select dataset from the dropdown list.", placement = "right", trigger = "click"),

  bsTooltip(id = "circ-gene_help", title = "Enter the gene name or a list of gene names separated by commas to plot expression. Make sure you are using the official gene name.", placement = "right", trigger = "click"),
  bsTooltip(id = "circ-enrichment_help",
            title = paste(
              "Choose whether enrichment is calculated against genes tested",
              "for rhythmicity in the selected dataset or against the",
              "species-wide STRING proteome. This changes the term list and enrichment",
              "shown in node tooltips and which nodes are highlighted, but not",
              "the STRING interaction edges.",
              "Then select a category and term to highlight its genes; click",
              "'Plot selected genes' to view their expression patterns."
            ),
            placement = "right", trigger = "click"),
  bsTooltip(id = "circ-data_table_help", title = "This table displays the whole dataset. Use the filters at the top to refine the results. Click on rows to display gene expression plots. Click the button at the button to clear the selection and plots.", placement = "right", trigger = "click"),
  bsTooltip(id = "circ-network_help", title = "This is the STRING protein interaction network graph for the selected dataset. Hover over nodes to see details.Click the node to display gene expression plot.", placement = "left", trigger = "click"),
  bsTooltip(id = "circ-graph_help", title = "Use these controls to toggle the sine fitting and connecting line on the gene expression plots. BH.Q value is derived from RAIN rhythmicity test. Fitted sine is for visualisation only.", placement = "middle", trigger = "click"),


)
