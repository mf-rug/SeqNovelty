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
library(httr2)
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
    envname = 'venv_seqnovelty_shiny',
    python = '/opt/homebrew/bin/python3'
  )
  reticulate::virtualenv_install(
    'venv_seqnovelty_shiny',
    packages = c("biopython")
  )
}
reticulate::use_virtualenv('venv_seqnovelty_shiny', required = TRUE)

shinyUI(page_navbar(
  tags$head(
    tags$link(href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0-beta3/css/all.min.css", rel = "stylesheet")
  ), 
  includeCSS("www/theme.css"),
  includeScript("www/custom.js"),
  useShinyjs(),
  use_shiny_title(),
  title = HTML("<strong><i>de novo</i> Sequence Novelty</strong>"),
  nav_spacer(), # push nav items to the right
  nav_item(actionBttn('help', HTML('help'), color = 'royal', style = 'stretch', size = 'xs')),
  nav_item(actionBttn('about', HTML('about'), color = 'royal', style = 'stretch', size = 'xs')),
  nav_item(HTML('<span style="white-space: nowrap; font-size: 11px;">&nbsp made by <a href ="https://bsky.app/profile/maxfus.bsky.social" target="_blank">Max Fürst</a></span>'), ),
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
          `<i class='fa fa-search'></i><br>search` = "search",
          `<i class='fa fa-align-left'></i><br>align` = "align",
          `<i class='fa fa-align-justify'></i><br>MSA` = "msa",
          `<i class='fa fa-folder-open'></i><br>restore` = "restore"),
      size = 'sm',
      justified = TRUE
    ),
    conditionalPanel(
      "input.mode=='search'",
      textAreaInput(
        "seq_input",
        NULL,
        height = "15vh",
#         value = '>esmGFP
# MSKVEELIKPDMKMKLEMEGEVNGHKFSIEAEGEGKPYEGKQTIKAWSTTGKLPFAWDILSTSLTYGNRAFTKYPEGLEQHDFFKQSFPEGYSWERTITYEDGATVKVTADISLEDGVLINKVKFKGENFPSDGPVMQKKTTGWEASTELITPDPATGGLKGEVKMRLKLEGGGHLLADFKTTYRSKKKEKLPLPGVHYVDHRIVNEKATHPEGKEYMIQYEHAVARL',
        placeholder = "Paste single sequence (FASTA or seq only)"
      ), radioGroupButtons(
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
                              choices = c(100, 250, 500, 1000),
                              selected = 500,
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
    conditionalPanel(
      "input.mode=='restore'", 
      textAreaInput("restore_id", 
                    "Paste the restore code", 
                    placeholder = 'e.g. 0123456789aAbBcCdDeEfFgGhHiIjJkKlLmMnNoOpPqQrRsStTuUvVwWxXyYzZ',
                    height = "15vh"),
      actionButton('restore_btn', 'Restore')
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
    conditionalPanel(
      "input.mode!='restore'",
        HTML('<hr style="margin: 0px 0 10px; border-color: #656565; "/>'),
        radioGroupButtons(
          'align_mode',
          NULL,
          choices = c(
            'pairwise',
            "MSA"
          ),
          size = 's',
          justified = TRUE
        ),
        conditionalPanel(
          "input.align_mode=='MSA'",
          pickerInput(
            inputId = "msa_tool",
            label = HTML("MSA algorithm:&nbsp"), inline = TRUE,
            choices = c("clustalo","kalign","mafft","muscle","muscle5","tcoffee"),
            selected = 'mafft',
            width = "fit")
        ),
        conditionalPanel(
          "input.align_mode=='pairwise'",
          sliderInput('ident_cutoff', '% identity cutoff', 0, 100, 33)
        ),
        HTML('<hr style="margin: 0px 0 15px; border-color: #656565; "/>'),
        actionButton("run_script", "Run Analysis"),
        actionButton("run_testm", "Demo Run"),
        # actionButton("run_testm", "Test Run M"),
        # prettySwitch('rem_gaps', 'rem_gaps', value = TRUE, status = 'primary'),
      ),
    width=3
    ),
    
    mainPanel(
      div(
        style = "display: flex; justify-content: space-between; width: 100%;",
        
        div(id = 'hide_scroll',
            actionButton("scroll", NULL, class='btn-jump')
        ),
        div(id = 'save_hide',
            actionButton("save", 'save result', class='btn-jump')
        )
      ),
      DT::dataTableOutput("alignment_table"),  # Render the output table here
      HTML('<hr style="margin: 0px 0 5px; border-color: #656565; "/>'),
      fluidRow(
        column(6, uiOutput('novelty_text')),
        column(6, uiOutput('res_analysis'))
      ), br(),
      uiOutput('novelty_media'),
      br(),
      width = 9
    )
  )
))
