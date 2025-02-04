library(shiny)
library(bslib)
library(reticulate)
library(DT)  # For interactive tables
library(Biostrings)  # For parsing FASTA files
library(jsonlite)
library(ggplot2)
library(colourpicker)
library(shinyWidgets)
library(shinytitle)
library(shinyjs)
library(httr)
library(stringr)
library(readr)
library(jsonlite)
library(r3dmol)
library(xml2)

# Configure reticulate options and virtual environment
options(reticulate.output_handler = function(x) cat(x, "\n"))

# Ensure virtual environment exists
envs <- reticulate::virtualenv_list()
if (!'venv_seqnovelty_shiny' %in% envs) {
  reticulate::virtualenv_create(
    envname = 'venv_seqnovelty_shiny'
    # python = '/opt/homebrew/bin/python3'
  )
  reticulate::virtualenv_install(
    'venv_seqnovelty_shiny',
    packages = c("biopython")
  )
}
reticulate::use_virtualenv('venv_seqnovelty_shiny', required = TRUE)

shinyUI(page_navbar(
  tags$script(HTML("
// Listen for clicks on table cells
$(document).on('click', '#alignment_table table.dataTable tbody td', function() {
    var table = $('#alignment_table table.dataTable').DataTable();
    var cell = table.cell(this);
    var colIndex = cell.index().column;
    var rowIndex = cell.index().row;
    
    // Remove previous column highlights
    $('#alignment_table table.dataTable tbody td').removeClass('selected-column');

    // Highlight all cells in the selected column
    $('#alignment_table table.dataTable tbody tr').each(function() {
        $(this).find('td').eq(colIndex).addClass('selected-column');
    });

    // Remove previous row highlights
    $('#alignment_table table.dataTable tbody tr td').removeClass('selected-row');

    // Highlight the clicked row
    var rowElement = $(this).closest('tr');
    rowElement.find('td').addClass('selected-row');

    // Remove previous cell border
    $('#alignment_table table.dataTable tbody td').removeClass('selected-cell');

    // Add border to the selected cell
    $(this).addClass('selected-cell');

    // Send column index to Shiny
    Shiny.setInputValue('col_selected', colIndex , {priority: 'event'});
    Shiny.setInputValue('row_selected', rowIndex , {priority: 'event'});
  });

    
  Shiny.addCustomMessageHandler('scrollTableToColumn', function(columnIndex) {
    console.log('scrollTableToColumn message received with columnIndex:', columnIndex);

// Find the DataTable's scrollable container
var scrollContainer = document.querySelector('#alignment_table .dataTables_scrollBody');
var table = document.querySelector('#alignment_table table.dataTable');

if (scrollContainer && table) {
    console.log('✅ Found table');

    // Use DataTables API to get rows
    var dt = $('#alignment_table table.dataTable').DataTable();
    var rowNodes = dt.rows().nodes();
    
    // Get columns from the first row instead of the header
    var columns = rowNodes.length > 0 ? rowNodes[0].querySelectorAll('td') : [];
    
    console.log('✅ Column count based on first row:', columns.length);

    if (columnIndex < columns.length) {
        var offset = 0;
        
        for (var i = 0; i < columnIndex; i++) {
            offset += columns[i].offsetWidth;  // Sum the widths of preceding columns
        }
        
        console.log('✅ Calculated offset for column:', offset);
        scrollContainer.scrollLeft = offset;  // Scroll to the calculated position
    } else {
        console.log('❌ Invalid columnIndex:', columnIndex);
    }
} else {
    console.log('❌ Scrollable container or table not found!');
}

  });
"))
  ,
  tags$head(
    tags$style(HTML("
    .selected-column {
      background-color: rgba(0, 123, 255, 0.2) !important; /* Bootstrap blue with transparency */
    }
     /* Override DT default row selection color */
    .selected-row {
      background-color: rgba(0, 222, 188, 0.2) !important;
    }
    
    /* Selected Cell Border (Black) */
.selected-cell {
    box-shadow: #de00c6 0px 0px 1px 1px inset !important;
    background-color: #9900FF69 !important;
}

/* 🔹 Add a dashed border around the hovered cell */
#alignment_table table.dataTable tbody tr:not(:first-child) td:not(:nth-child(1)):hover {
    box-shadow: #de00c6 0px 0px 1px 1px inset !important;
    cursor: crosshair;
}
    
    .error_out {
        height: 100%;
        border: 1px solid transparent;
        box-shadow: 2px 2px 8px rgba(50, 180, 200, 0.3);
        border-image: repeating-linear-gradient(to right,
  #c4e17f 7%, #c4e17f 14%,
  #c9ff56 21%, #f2fa71 27%,
  #fad071 33%, #fad071 40%,
  #f0766b 47%, #f0766b 54%,
  #db9dbe 60%, #db9dbe 66%,
  #c49cdf 71%, #c49cdf 78%,
  #6599e2 84%, #6599e2 90%,
  #61c2e4 95%, #61c2e4 100%);
        border-image-slice: 1;
        background: linear-gradient(to right, rgba(0, 100, 50, 0.01), rgba(0, 200, 150, 0.05)); 
        border-radius: 0.5rem;
        padding: 0.5rem;
    }
    
.flex-container {
  display: flex;
  align-items: center; /* Ensures items align on the same row */
  gap: 20px; /* Spacing between items */
}

.flex-item {
  display: flex; /* Allows nested elements to stretch correctly */
  align-items: center; /* Ensures inner elements align properly */
}
.colourpicker-input {
padding-top:4px;
padding-bottom:4px;
}
.action-button {
padding-top:6px;
padding-bottom:6px;
}
.btn-default {
    --bs-btn-bg: #cceef8 !important; /* White border */
    margin-right: 3px;
    margin-left: 3px;
}
[data-bs-theme='dark'] .btn-default {
    --bs-btn-bg: #0e495a !important; /* White border */
}
    
    table.dataTable td {
      padding-left: 0.55px !important;
      padding-right: 0.55px !important;
      padding-top: 0 !important;
      padding-bottom: 0 !important;
      margin-top: 0.25rem;
      margin-bottom: 0.25rem;
      max-width: 9px !important;
      font-family: 'Consolas', 'Menlo', 'Monaco', 'Roboto Mono', 'Courier New', Courier, monospace !important;
      font-size: 15px;
      text-align: center; /* Center align other cells */
      vertical-align: middle;
      overflow: hidden; /* Ensure text does not overflow out of the cell */
    }
    
/* Make the rowname column sticky */
table.dataTable .rowname-col {
  text-align: left !important;
  width: 152px !important;
  max-width: 152px !important;
  overflow: hidden !important;
  margin-right: 10px;
  padding-right: 10px;
  align-items: center !important; /* Vertically center content */

  /* Sticky column settings */
  position: sticky !important;
  left: 0 !important;
  z-index: 2 !important; /* Ensure it's above other columns */
}

/* Ensure the span fully stretches */
table.dataTable .rowname-col div {
  display: flex !important;
  align-items: center !important;
  width: 100% !important;
  line-height: 100%;
  overflow-y: hidden; /* Disable vertical scrolling */
  scrollbar-width: thin; /* Firefox */
  -ms-overflow-style: none; /* Hide scrollbar for Edge */
  overflow-x: auto !important; /* Enable scrolling */
  white-space: nowrap !important;
}

/* Hide scrollbar for WebKit browsers (Chrome, Safari) */
table.dataTable .rowname-col div::-webkit-scrollbar {
  display: none;
  height: 6px; /* Adjust scrollbar thickness */
}

/* Ensure the span inside the column fully expands */
table.dataTable .rowname-col span {
  display: block !important;
  width: 100% !important;
}
.spacer-col2 {
  border-right-width: 3px;
  border-right-style: solid;
  border-right-color: #8A8A8A1D;
}

    table.dataTable th {
      font-weight: normal !important; /* Override <strong> styling */
      font-family: 'Consolas', 'Menlo', 'Monaco', 'Roboto Mono', 'Courier New', Courier, monospace !important; /* Modern fonts with fallback */      
      font-size: 9px;
      padding-bottom: 0;
      border-bottom-width: 0.5px;
      border-bottom-style: solid;
      border-bottom-color: grey;
      letter-spacing: -1px; /* Decrease character spacing */
      padding-left: 0.5px !important;
      padding-right: 0.5px !important;
      text-align: center !important;

    }
    
      #seq_input {
        padding: 6px !important;
        font-size: 13px !important; /* Adjust font size */
      }

  ")),
    tags$link(
      href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0-beta3/css/all.min.css", 
      rel = "stylesheet"
    )
  ),
  use_shiny_title(),
  useShinyjs(),
  title = HTML("<strong><i>de novo</i> Sequence Novelty</strong>"),
  nav_spacer(), # push nav items to the right
  nav_item(actionBttn('info', HTML('&nbsp help/info &nbsp'), color = 'royal', style = 'stretch', size = 'xs')),
  nav_item(HTML('<span style="white-space: nowrap; font-size: 11px;">made by <a href ="https://bsky.app/profile/maxfus.bsky.social" target="_blank">Max Fürst</a></span>'), ),
  nav_item(),
  nav_item(
    input_dark_mode(id = "dark_mode")
  ),
  sidebarLayout(
    sidebarPanel(
      radioGroupButtons(
        'mode',
        NULL,
        choices = c(
          `<i class='fa fa-search'></i><br>Search` = "search",
          `<i class='fa fa-align-left'></i><br>align` = "align",
          `<i class='fa fa-align-justify'></i><br>MSA` = "msa"),
        size = 'sm',
      justified = TRUE
    ),
    conditionalPanel(
      "input.mode=='search'",
      textAreaInput("seq_input", NULL, height= "15vh", placeholder = "Paste single sequence (FASTA or seq only)"),
      radioGroupButtons(
        'searchmode',
        NULL,
        choices = c(
          'BLAST',
          "mmseqs"
          ),
        size = 's',
        justified = TRUE
      ),
      conditionalPanel(
        "input.searchmode=='BLAST'",
          fluidRow(column(6,
                            pickerInput(
                              inputId = "blast_db",
                              label = "database:", 
                              choices = c("nr", "swissprot", "pdb", "pataa", "refseq_select", "refseq_protein", "landmark", "env_nr", "tsa_nr"),
                              selected = 'nr',
                              width = "100%")),
                   column(6,
                            sliderTextInput(
                              inputId = "blast_seq_n",
                              label = "max num seq:", 
                              choices = c(50, 100, 250, 500, 1000),
                              selected = 250,
                              grid = TRUE
                            ))),
            actionButton("run_blast", "BLAST")
      ),
      conditionalPanel(
        "input.searchmode=='mmseqs'",
        pickerInput(
          inputId = "mmseqs_db",
          label = HTML("database:&nbsp"), inline = TRUE,
          choices = c("afdb50","afdb-swissprot","BFVD","afdb-proteome","bfmd","cath50","mgnify_esm30","pdb100","gmgcl_id"),
          selected = 'afdb50',
          width = "fit"),
        actionButton("run_mmseqs", "Run mmseqs")
      )
    ),
    conditionalPanel(
      "input.mode=='align'", fileInput("fasta_file", "Upload multi sequence .fasta file", accept = c(".fasta", ".fa"))
    ),
    conditionalPanel(
      "input.mode=='msa'", fileInput("fasta_file_msa", "Upload MSA .fasta file", accept = c(".fasta", ".fa"))
    ),
      uiOutput("seq_id_dropdown"),  # Dynamically generated dropdown for sequence IDs
    br(),
    tags$details(
      tags$summary(HTML("Visualization Settings<br><br>")),
      div(
        HTML('<small>'),
        HTML('<strong>Ref sequence</strong>'),
          HTML('<hr style="margin: 0px 0 10px; border-color: #656565; "/>'),
        div(class = "flex-container",
            div(class = "flex-item", HTML("novel:&nbsp<br>&nbsp")),
            div(class = "flex-item", prettySwitch('boldn', 'bold', value = TRUE, status = 'primary')),
            div(class = "flex-item", colourpicker::colourInput(
              "novel_col", NULL, "black",
              allowTransparent = TRUE,
              closeOnClick = TRUE))),
        div(class = "flex-container",
            div(class = "flex-item", HTML("matched:&nbsp&nbsp<br>&nbsp"),
                   prettySwitch('boldm', 'bold', value = FALSE, status = 'primary')),
            div(class = "flex-item", colourpicker::colourInput(
              "matched_col", NULL, "black",
              allowTransparent = TRUE,
              closeOnClick = TRUE))),
          prettySwitch('ref_nums', 'show ref seq numbering', value = FALSE, status = 'primary'),
          HTML('<strong>Aligned sequences</strong>'),
          HTML('<hr style="margin: 0px 0 10px; border-color: #656565; "/>'),
        div(class = "flex-container",
            div(class = "flex-item", "consumed match"),
            div(class = "flex-item", colourpicker::colourInput(
              "match_con_col", NULL, "black",
              allowTransparent = TRUE,
              closeOnClick = TRUE))),
        div(class = "flex-container",
            div(class = "flex-item", "unconsumed match"),
            div(class = "flex-item", 
                colourpicker::colourInput(
              "match_noncon_col", NULL, "#94949480",
              allowTransparent = TRUE,
              closeOnClick = TRUE))),
        div(class = "flex-container",
            div(class = "flex-item", "no match"),
            div(class = "flex-item", colourpicker::colourInput(
             "nomatch_col", NULL, "#9494942C",
             allowTransparent = TRUE,
             closeOnClick = TRUE))),
        HTML('</small')
      )
    ),
    actionButton("run_script", "Run Analysis"),
    actionButton("run_testm", "Demo Run"),
    # actionButton("run_testm", "Test Run M"),
    # prettySwitch('rem_gaps', 'rem_gaps', value = TRUE, status = 'primary'),
    width=3
    ),
    
    mainPanel(
      div(id = 'hide_scroll', style="text-align:center;",
          actionLink("xxx", NULL)
      ),
      DT::dataTableOutput("alignment_table"),  # Render the output table here
      HTML('<hr style="margin: 0px 0 5px; border-color: #656565; "/>'),
      uiOutput('novelty_text'), br(),
      uiOutput('novelty_media'),
      br(),
      width = 9
    )
  )
))
