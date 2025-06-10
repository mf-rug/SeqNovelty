sstop <- function(session, html, title, text, type) {
  sendSweetAlert(
    session = session,
    html = html,
    title = title,
    text = text,
    type = type
  )
}


# blast_json <- fromJSON("blast_out.blast")
# hits <- blast_json$BlastOutput2$report$results$search$hits
# hit_df <- map_dfr(hits, function(hit) {
#   descriptions <- hit$description
#   hsps <- hit$hsps
#   
#   # If there are multiple HSPs, take the first one (most significant)
#   evalue <- hsps[[1]]$evalue
#   
#   # For each description, extract id and evalue
#   map_dfr(descriptions, function(desc) {
#     tibble(
#       id = desc$id,
#       accession = desc$accession,
#       title = desc$title,
#       evalue = evalue
#     )
#   })
# })

shinyServer(function(input, output, session) {
  
  change_window_title(session, 'SeqNovelty')
  # Reactive values to store results and sequence IDs
  results <- reactiveValues(log = NULL, seq_ids = NULL, alignment_table = NULL, seq_info = NULL, coverage_df = NULL)
  ref_name <- reactiveVal(NULL)
  test <- reactiveVal(NULL)
  ref_seq_start_col <- reactiveVal(NULL)
  structure <- reactiveVal('')
  saveID <- reactiveVal(NULL)
  hit_df <- reactiveVal(NULL)
  
  hl <- function(...) {
    paste0('<span style="color:', ifelse(input$dark_mode == 'dark', 'aqua', 'teal'), 
           ';"><big>', ...,'</span></big>')
  }
  
  parse_fasta_in <- function(input_string) {
    if (str_detect(input_string, '^>') & str_detect(input_string, '\n')) {
      input_id <- str_extract(input_string, '(?<=>).*(?=\n)') 
      input_seq <- str_remove(input_string, fixed(input_id)) %>% str_remove_all(., '^>|\n| ')
      input_id <- input_id %>% str_replace_all(., ' ', '_')
    } else {
      input_id <- 'UnknownInput'
      input_seq <- str_remove_all(input_string, '\n| ')
    }
    if (is.na(input_seq) || is.null(input_seq) || str_detect(input_seq, '[^ACDEFGHIKLMNPQRSTVWYZ]')) {
      invalid_str <- paste(str_extract_all(input_seq, '[^ACDEFGHIKLMNPQRSTVWYZ]')[[1]], collapse = ',')
      if (!is.na(invalid_str)) {
        invalid_info <- paste0(': ', invalid_str)
      } else {
        invalid_info <- '.'
      }
      cat('Invalid format / invalid characters in input seq:', input_seq, '\n')
      sendSweetAlert(
        session = session,
        html = TRUE,
        title = "Input error",
        text = tags$span(
          paste0("Your input sequence contains invalid characters", invalid_info)
        ), 
        type = "error"
      )
      return(NULL)
    }
    
    ref_name(input_id)
    cat(paste0('input id is: ', input_id, '\ninput seq is: ', input_seq, ' . Starting Blast\n'))
    return(c(input_id, input_seq))
  }
  
  convert_col_to_residue <- function(sequence, col_indices) {
    # Find positions of actual residues (non-gap)
    non_gap_positions <- cumsum(sequence != "-") * (sequence != "-")
    
    # Extract corresponding residue numbers
    residue_numbers <- non_gap_positions[col_indices]
    
    return(residue_numbers)
  }
  

  # Function to submit and retrieve MAFFT alignment
  run_msa <- function(fasta_file, msa_tool) {
    cutoff <- ifelse(input$msa_tool == 'clustalo',
                     3999,
                     ifelse(input$msa_tool == 'kalign', 1999, 499))
    
    fasta_content <- readLines(fasta_file)
    # browser()
    
    nseqs <- sum(str_count(fasta_content, '^>'))
    if (nseqs > cutoff) {
      showNotification(
        paste0(
          'This tool is limited to ',
          cutoff,
          ' sequences, removing ',
          nseqs - cutoff,
          ' sequences from analysis.'
        )
      )
      Sys.sleep(4)
      fasta_content <- fasta_content[1:which(str_detect(fasta_content, '^>'))[cutoff +1]-1]
    } 
    # Read the FASTA file
    
    # fasta_id_lines <- which(str_detect(fasta_content, '^>'))
    # max_seqs <- min(nseqs, cutoff + 1)
    # cutoff_line <- which(str_detect(fasta_content, '^>'))[max_seqs]
    # fasta_content <- fasta_content[1:cutoff_line-1]
    

    fasta_content <- paste(fasta_content, collapse = '\n')
    
    # Step 1: Submit the job
    response <- POST(
      url = paste0("https://www.ebi.ac.uk/Tools/services/rest/",msa_tool ,"/run"),
      body = list(
        email = "your@email.com",
        title = "R_seqnovelty_Job",
        stype = "protein",
        sequence = fasta_content,
        gapopen = 1.6,
        gapext = 0.15
      ),
      encode = "form"
    )
    
    if (str_detect(content(response, as = "text"), 'error')) {
      showNotification(paste('error', content(response, as = "text")))
      cat('error', content(response, as = "text"))
      return(NULL)
    } else {
      job_id <- content(response, as = "text")
    }
    # Extract Job ID
    live_status <- paste0("https://www.ebi.ac.uk/Tools/services/rest/",msa_tool ,"/status/", job_id)
    cat("Job submitted, Job ID:", job_id, "\n")
    
    # Step 2: Check job status
    repeat {
      Sys.sleep(10)  # Wait 5 seconds before checking status
      status_response <- GET(live_status)
      status <- content(status_response, as = "text")
      status_out <- paste0("Job status: ", status, "\n", live_status, '\n')
      cat(status_out)
      incProgress(amount = 0.001, detail = status_out)
      
      if (status == "FINISHED") {
        break
      }
    }
    
    # Step 3: Retrieve alignment in FASTA format
    result_response <- GET(paste0('https://www.ebi.ac.uk/Tools/services/rest/', msa_tool, '/result/', job_id, "/fa"))
    aligned_sequences <- content(result_response, as = "text")
    
    return(aligned_sequences)
  }
  
  
  # Function to retrieve FASTA sequences in batches
  get_fasta_from_uniprot <- function(ids, batch_size = 500) {
    base_url <- "https://rest.uniprot.org/uniprotkb/stream"
    all_fasta <- ""
    
    # Process IDs in batches
    for (batch in split(ids, ceiling(seq_along(ids) / batch_size))) {
      query <- paste0("query=(", paste(batch, collapse = "%20OR%20"), ")")
      url <- paste0(base_url, "?", query, "&format=fasta")
      
      response <- request(url) |>
        req_perform()
      
      fasta_data <- resp_body_string(response)  # Read as plain text
      all_fasta <- paste0(all_fasta, fasta_data, "\n")  # Append results
    }
    cat('Retrieved ', str_count(all_fasta, '^>'), ' fasta sequences.\n')
    return(all_fasta)
  }
  
  run_mmseqs <- function(file) {
    # Define API URLs
    api_base <- "https://search.foldseek.com/api"
    query_endpoint <- paste0(api_base, "/ticket")
    status_endpoint <- paste0(api_base, "/ticket/")
    result_endpoint <- paste0(api_base, "/result/download/")
    
    # Read original sequence from file
    fasta_lines <- readLines(file)
    if (length(fasta_lines) < 2) stop("Invalid FASTA file format.")
    sequence_name <- sub("^>", "", fasta_lines[1])  # Extract sequence ID
    sequence <- paste(fasta_lines[-1], collapse = "")  # Combine sequence lines
    
    # Step 1: Fetch the 3DI sequence
    three_di_url <- paste0("https://3di.foldseek.com/predict/", sequence)
    three_di_response <- GET(three_di_url)
    if (status_code(three_di_response) != 200) stop("Failed to fetch 3DI sequence.")
    
    three_di_sequence <- content(three_di_response, "text", encoding = "UTF-8") %>%
      str_replace_all('^"|"$', '') %>%  # Remove leading/trailing quotes
      str_trim()  # Trim any extra newlines or spaces
    
    # Step 2: Create a new FASTA file with both sequences
    modified_fasta <- "modified_ex.fasta"
    writeLines(c(
      paste0(">", sequence_name),
      sequence,
      ">3DI",
      three_di_sequence
    ), modified_fasta)
    
    cat("Modified FASTA file created:", modified_fasta, "\n")
    
    # Step 3: Submit the request (upload modified_ex.fasta)
    query_response <- POST(
      query_endpoint,
      body = list(
        q = upload_file(modified_fasta),  # Upload new FASTA
        mode = "3diaa",  # Specify mode
        iterative = 1,  # Search AFDB50 database 
        `database[]` = "afdb50"  # Search AFDB50 database 
      ),
      encode = "multipart"
    )
    
    # Parse response correctly
    query_result <- fromJSON(content(query_response, "text", encoding = "UTF-8"))
    
    if (!"id" %in% names(query_result)) {
      stop("Error: Failed to submit query.")
    }
    job_id <- query_result$id
    cat("Job submitted to foldseek server. ID:", job_id, "\n")
    
    # Step 4: Poll for job completion
    repeat {
      Sys.sleep(5)  # Wait 5 seconds before checking status
      status_response <- GET(paste0(status_endpoint, job_id))
      status_result <- fromJSON(content(status_response, "text", encoding = "UTF-8"))
      
      cat("Current Status:", status_result$status, "\n")
      
      if (status_result$status == "COMPLETE") break
      if (status_result$status == "ERROR") stop("Error: Job failed.")
    }
    
    # Step 5: Download the results
    output_file <- "msa_results.tar.gz"
    download_url <- paste0(result_endpoint, job_id)
    download.file(download_url, destfile = output_file, mode = "wb")
    cat("Results downloaded as", output_file, "\n")
    
    # Step 6: Extract the tar.gz file
    untar(output_file, exdir = "msa_results") 
    file.remove(output_file)
    m8_file <- list.files("msa_results", pattern = "\\.m8$", full.names = TRUE)
    if (length(m8_file) == 0) stop("No .m8 file found in extracted results.")
    
    # Step 7: Process `.m8` file and extract hit IDs
    m8_data <- read_tsv(m8_file, col_names = FALSE, show_col_types = FALSE)
    hit_ids <- str_extract(m8_data$X2, "(?<=^AF-)[^-]+")
    return(hit_ids)
  }
  
  
  
  
  
  run_blast <- function(sequence, program = "blastp", database = "nr") {
    # NCBI BLAST API URL
    base_url <- "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    
    # Step 1: Submit the BLAST job
    response <- POST(
      url = base_url,
      body = list(
        CMD = "Put",               # Submit job
        PROGRAM = program,         # Specify BLAST program (e.g., blastn, blastp)
        DATABASE = input$blast_db,       # Specify database (e.g., nr, pdb, swissprot)
        QUERY = sequence,          # Input sequence
        FORMAT_TYPE = "JSON2_S",   # Output format
        MAX_NUM_SEQ = input$blast_seq_n     # Maximum number of sequences
      ), encode = "form"
    )
    
    if (response$status_code != 200) {
      stop("Failed to submit BLAST job.")
    }
    # Extract Request ID (RID) from response
    rid <- sub(".*RID = ([^\\n]+).*", "\\1", content(response, as = "text"))
    rid <- trimws(rid)
    rid <- str_extract(rid, '.*(?=\n)')
    
    cat("BLAST job submitted. RID:", rid, "\nhttps://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=", rid)
    
    # Step 2: Poll for results
    for (i in 1:30) {
      incProgress(0.01, message = 'BLAST running') #, detail = paste0('https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=', rid))
      url <- paste0('https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=', rid)
      session$sendCustomMessage(
        "updateProgressDetail",
        paste0('<a href="', url, '" target="_blank">', url, ' &nbsp <i class="fas fa-external-link-alt"></i></a>')
      )
      Sys.sleep(30) # Wait 30 seconds before checking
      status_response <- GET(base_url, query = list(CMD = "Get", RID = rid, FORMAT_TYPE = "JSON2_S"))
      status_text <- content(status_response, as = "text")
      
      if (grepl("Status=(?:WAITING|UNKNOWN)", status_text)) {
        cat("Still waiting for results of", rid, ":", str_extract(status_text, 'Status=.*'), ",\n")
      } else if (grepl("Status=FAILED", status_text)) {
        stop("BLAST job failed.")
      } else if (!grepl('Status=', status_text)) {
        cat("BLAST results ready.\n")
        
        # Retrieve results in text format
        results_response <- GET(base_url, query = list(
          CMD = "Get",
          RID = rid,
          FORMAT_TYPE = "JSON2_S",          # Options: JSON2, XML, Text
          FORMAT_OBJECT = "Alignment"   # Retrieve alignments
        ))
        return(list(rid, content(results_response, as = "text")))
      } else {
        print('something odd')
        # browser()
      }
    }
    sstop(session, TRUE, "Error in Blast", 
          text = tags$span(
            "BLAST job did not complete in time. You might want to run",  
            tags$a(href = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins', target='_blank', 'BLAST'),
            'manually and paste the resulting link in here instead.'),
          'error')
  }
  
  extract_blast_df <- function(json_file) {

    # Read JSON file
    blast_data <- fromJSON(json_file)
    
    # Get all hits
    hits <- blast_data$BlastOutput2$report$results$search$hits
    
    df <- map_dfr(hits, function(hit) {
      descriptions <- hit$description
      hsps <- hit$hsps
      
      # If there are multiple HSPs, take the first one (most significant)
      evalue <- hsps[[1]]$evalue
      
      # For each description, extract id and evalue
      map_dfr(descriptions, function(desc) {
        tibble(
          id = desc$id,
          evalue = evalue
        )
      })
    })
    
    return(df)
  }
  
  extract_blast_ids <- function(json_file, single_id=TRUE) {
    # Read the JSON file
    blast_data <- fromJSON(json_file)
    
    # Navigate to the hits section
    hits <- blast_data$BlastOutput2$report$results$search$hits
    
    # Extract IDs of each hit
    hit_ids <- sapply(hits[[1]]$description, function(desc) {
      desc$id  # Extract the 'id' field
    })
    if (single_id) {
      sapply(hit_ids, function(x) x[1]) 
    } else{
      return(hit_ids)
    }
  }
  
  fetch_sequences <- function(ids, output_file = "search_sequences.fasta", chunk_size=100) {
    # Flatten all IDs into a single vector
    all_ids <- unlist(ids)  # Flatten if there are multiple IDs per hit
    
    # Split IDs into chunks of 200
    id_chunks <- split(all_ids, ceiling(seq_along(all_ids) / chunk_size))
    
    # Base URL for efetch
    base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    
    # Initialize a variable to store the merged FASTA data
    all_fasta <- ""
    
    # Loop through each chunk and fetch sequences
    for (chunk in id_chunks) {
      print(paste('fetching chunk'))
      # Create a comma-separated string of IDs for the current chunk
      # here used to be a [[1]] list, why?
      # seqs <- paste(chunk[[1]], collapse = ",")
      seqs <- paste(chunk, collapse = ",")
      
      # Construct the URL for the current chunk
      query_url <- modify_url(base_url, query = list(
        db = "protein",
        id = paste0(seqs, collapse = ''),
        rettype = "fasta",
        retmode = "text"
      ))

      # Perform the GET request
      response <- GET(query_url)
      
      # Check for successful response
      if (status_code(response) != 200) {
        cat("Failed to fetch sequences for chunk. HTTP status code: ", status_code(response), '\n')
        
        # try once more
        Sys.sleep(5)
        response <- GET(query_url)
        if (status_code(response) != 200) {
          cat("Failed to fetch sequences for chunk. HTTP status code: ", status_code(response), '\n')
          next
        }
      }
      # Append the fetched FASTA data to the merged result
      all_fasta <- paste0(all_fasta, content(response, as = "text", encoding = "UTF-8"))
      Sys.sleep(2)
    }
    
    return(all_fasta)
  }
  
  
  observe({
    if (input$seq_input == '') {
      shinyjs::disable('run_blast') 
      shinyjs::disable('load_blast') 
    } else {
      shinyjs::enable('run_blast')
      if (input$blast_url != '') {
        shinyjs::enable('load_blast')
      } else {
        shinyjs::disable('load_blast') 
      }
    }
  })
  
  
  observeEvent(input$run_mmseqs, {
    input_parse <- parse_fasta_in(input$seq_input)
    req(input_parse)
    input_id <- input_parse[1]
    input_seq <- input_parse[2]
    
    if (str_count(input_seq) > 400) {
      sendSweetAlert(
        session = session,
        html = TRUE,
        title = "mmseqs only supports up to 400 amino acids.",
        text = tags$span(
          "Please reduce length or run BLAST instead."
        ), 
        type = "error"
      )
      return()
    }
    
    write_file(input$seq_input, 'input_fasta.fasta')
    
    withProgress(message = "Running mmseqs", value = 0.15, {
      ids <- run_mmseqs('input_fasta.fasta')
    
      incProgress(message = "Getting sequences", amount = 0.6)
      tryCatch({
        fasta_data <- get_fasta_from_uniprot(ids)
      }, error = function(e) {
        showNotification(e$message)
        cat(e$message, '\n')
        return()
      })
    })
    # browser()
    write_file(paste0('>', input_id, '\n', input_seq, 
                      '\n', fasta_data),
               file = 'search_sequences.fasta')
    sendSweetAlert(
      session = session,
      html = TRUE,
      title = "mmseqs results ready.",
      text = tags$span(
        "mmseqs has finished successfully, ", length(ids), " sequences were found.", tags$br(),
        "To analyse the sequences press the 'Run Analysis' button."
      ), 
      type = "info"
    )
  })
  
  
  observe({
    if (!is.null(ref_name())) {
      print('not null')
    }
  })
  
  
  observeEvent(input$load_blast, {
    ref_name(NULL)
    results$alignment_table <- NULL
    results$coverage_df <- NULL
    results$seq_ids <- NULL
    results$seq_ids <- NULL
    print(paste('ref_name() and results now', ref_name()))
    input_parse <- parse_fasta_in(input$seq_input)
    req(input_parse)
    input_id <- input_parse[1]
    input_seq <- input_parse[2]
    ref_name(input_id)

    if (startsWith(input$blast_url, 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=')) {
      rid <- str_remove(input$blast_url, fixed('https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=')) %>% str_remove(., fixed('%20'))
    } else if (str_detect(input$blast_url,'^[A-Z0-9]{10,12}$')) {
      rid <- input$blast_url
    } else {
      sstop(session, TRUE, 'Error', 'You have to provide a valid URL or RID', 'error')
      return()
    }
    cat('rid is: ', rid, '\n')
    base_url <- "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    response <- GET(
      url = base_url,
      query = list(
        CMD = "Get",
        RID = rid,
        FORMAT_TYPE = "JSON2_S",          # Options: JSON2, XML, Text
        FORMAT_OBJECT = "Alignment"   # Retrieve alignments
      )
    )
    
    if (status_code(response) != 200 ) {
      sstop(session, TRUE, 'Error', paste("Failed to fetch BLAST result. Status: ", status_code(response)), 'error')
      return()
    } else {
      if (str_detect(content(response, 'text'), 'Error: Problems retrieving results for RID')) {
        sstop(session, TRUE, 'Error', paste("Failed to fetch BLAST result. Problems retrieving results for RID"), 'error')
        return()
      }
      else {
        blast_out <- content(response, "text")
      }
    }
    

    withProgress(message = "Processing BLAST", value = 0.6, {
      write_file(blast_out, 'blast_out.blast')
      ids <- extract_blast_ids('blast_out.blast', single_id=TRUE)
      hit_df <- extract_blast_df('blast_out.blast')
      hit_df(hit_df)
    })

    withProgress(message = "Getting sequences", value = 0.7, {
      fasta_data <- fetch_sequences(ids)
      # if (success) {
      writeLines(fasta_data, con = 'search_sequences.fasta')
      write_file(paste0('>',input_id,'\n',input_seq,'\n'), 
                 'search_sequences.fasta', 
                 append = TRUE)
      # } else { stop() }
    })
    
    blast_link <- paste0("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=", rid)

    sendSweetAlert(
      session = session,
      html = TRUE,
      title = "BLAST ready.",
      text = tags$span(
        "Blast finished successfully, ", length(ids), " sequences found. You can look at the results under ",
        tags$a(href = blast_link, target='_blank', blast_link), tags$br(), tags$br(),
        "To analyse the sequences press the 'Run Analysis' button."
      ), 
      type = "info"
    )
  })
  
  observeEvent(input$run_blast, {
    ref_name(NULL)
    results$alignment_table <- NULL
    results$coverage_df <- NULL
    results$seq_ids <- NULL
    results$seq_ids <- NULL
    print(paste('ref_name() and results now', ref_name()))
    input_parse <- parse_fasta_in(input$seq_input)
    req(input_parse)
    input_id <- input_parse[1]
    input_seq <- input_parse[2]
    
    withProgress(message = "Starting BLAST", value = 0.15, {
      blast_results <- run_blast(input_seq, program = 'blastp', database = 'nr')
      rid <- blast_results[[1]]
      blast_out <- blast_results[[2]]
    })
    # if (success) {
    withProgress(message = "Processing BLAST", value = 0.6, {
      write_file(blast_out, 'blast_out.blast')
      ids <- extract_blast_ids('blast_out.blast', single_id=TRUE)
    })
    # } else { stop() }
    
    withProgress(message = "Getting sequences", value = 0.7, {
      fasta_data <- fetch_sequences(ids)
      # if (success) {
      writeLines(fasta_data, con = 'search_sequences.fasta')
      write_file(paste0('>',input_id,'\n',input_seq,'\n'), 
                 'search_sequences.fasta', 
                 append = TRUE)
      # } else { stop() }
    })
    
    blast_link <- paste0("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=", rid)
    sendSweetAlert(
      session = session,
      html = TRUE,
      title = "BLAST ready.",
      text = tags$span(
        "Blast finished successfully, ", length(ids), " sequences found. You can look at the results under ",
        tags$a(href = blast_link, target='_blank', blast_link), tags$br(), tags$br(),
        "To analyse the sequences press the 'Run Analysis' button."
      ), 
      type = "info"
    )
})
  
  
  
  # Function to run the Python script and process results
  run_analysis <- function(input_file, ref_seq_id) {
    # if (is.null(ref_seq_id)) {
    #   ref_seq_id <- 'esmGFP'
    #   browser()
    # }
    updateActionButton(session, 'save', label = 'save result')
    saveID(NULL)
    
    tryCatch({
      # Import Python modules
      sys <- import("sys")
      if (is.null(test())) {
        # print('no test, performing msa')
        
        if (input$mode != 'msa' && input$align_mode == 'MSA') {
          
          # run local
          # # mafft alignment script
          # source_python("run_msa.py")  
          # 
          # # Set up `sys.argv` for the script
          # sys_args <- c(
          #   "run_msa.py",
          #   input_file
          # )
          # sys$argv <- sys_args
          # cat('Running this command:\npython3', paste(sys_args, collapse = " "),'\n')
          # # Run the Python script
          # withProgress(message = "Running sequence alignment", value = 0.15, {
          #   capture.output({
          #     py$main()
          #   })
          # })
          
          # run MSA online
          withProgress(message = paste0('Running ', input$msa_tool, ' alignment on EBI server'), value = 0.15, {
            print('EBi server')
            msa_out <- run_msa(input_file, input$msa_tool)
          })
          
          if (is.null(msa_out)) {
            cat('msa failed\n')
            return(NULL)
          }
          writeLines(msa_out, 'initial_alignment.fasta')
          cat('Success')

        } else {
          print(paste('its an msa / pairwise requested, so skipping ', input$msa_tool, ' and using', input_file, 'as initial_alignment.fasta'))
          if (input$align_mode == 'MSA') {
            file.copy(input_file, 'initial_alignment.fasta')
          } else {
            # perform all pairwise alignments
            tryCatch({
              # pairwise alignment script
              source_python("run_pairwise_align.py")  
              # Set up `sys.argv` for the script
              sys_args <- c(
                'run_pairwise_align.py',
                input_file,
                ref_seq_id,
                input$ident_cutoff
              )
              sys$argv <- sys_args
              cat('Running this command:\npython3', paste(sys_args, collapse = " "),'\n')
              # Run the Python script
              withProgress(message = "Running pairwise alignments", value = 0.22, {
                capture.output({
                  py$main()
                })
              })
            }, error = function(e) {
              results$log <- paste("Error:", e$message)
              print(paste("Error running pairwise alignment script. ", e$message), type = "error")
              showNotification(paste("Error running pairwise alignment script. ", e$message), type = "error")
            })
            
            tryCatch({
              # pairwise alignment merge script
              source_python("merge_pairwise.py")  
              # Set up `sys.argv` for the script
              sys_args <- c(
                'merge_pairwise2.py',
                'pairwise_alignments/',
                'initial_alignment.fasta'
              )
              sys$argv <- sys_args
              cat('Running this command:\npython3', paste(sys_args, collapse = " "),'\n')
              # Run the Python script
              withProgress(message = "Merging pairwise alignments", value = 0.3, {
                capture.output({
                  py$main()
                })
              })
            }, error = function(e) {
              results$log <- paste("Error:", e$message)
              print(paste("Error merging pairwise alignments ", e$message), type = "error")
              showNotification(paste("Error merging pairwise alignments ", e$message), type = "error")
            })
            
          }
        }
      } else {
        print('its a demo')
      }
    }, error = function(e) {
      results$log <- paste("Error:", e$message)
      print(paste('Error running ', input$msa_tool, ': ', e$message), type = "error")
      showNotification(paste('Error running ', input$msa_tool, ': ', e$message), type = "error")
      browser()
    })
    if (is.null(test()) && !file.exists('initial_alignment.fasta')) {
      sendSweetAlert(
        session = session,
        html = TRUE,
        title = "No Matching Sequences Found!",
        text = tags$span("Either your sequence is completely novel (congrats!) or you have to broaden the search criteria",
                         tags$br(), "Consider using BLAST to search, and when using pairwise alignments, set the identity cutoff to a lower %."),
        type = "info"
      )
      return()
      }
    
    tryCatch({
      if (is.null(test())) {
        # print('no test seqnovelty')
        print(ref_seq_id)
        # browser()
        # aln_df <- data.frame(id = str_remove(aln[seq(1,length(aln), by=2)], '^>'),
        #                      seq = aln[seq(2,length(aln), by=2)])
        # 
        # 
        # matched_df <- map_dfr(aln_df$id, function(qid) {
        #   match <- hit_df()[grepl(qid, hit_df()$id), ]
        #   
        #   if (nrow(match) > 0) {
        #     match[1, ]  # Take the first matching row
        #   } else {
        #     # Return a row with NAs and the same structure as hit_df
        #     tibble(id = NA_character_, evalue = NA_real_)
        #   }
        # })
        # aln_df <- bind_cols(aln_df, matched_df$evalue)  
        # 
        # # CONTINUE HERE< NOT DONE
        
        
        
        # sequence analysis script
        source_python("SeqNoveltyMin.py")  
        # Set up `sys.argv` for the script
        sys_args <- c(
          "SeqNoveltyMin.py",
          'initial_alignment.fasta',
          ref_seq_id
        )
        sys$argv <- sys_args
        cat('Running this command:\npython3', paste(sys_args, collapse = " "),'\n')
        # Run the Python script
        withProgress(message = "Running sequence analysis", value = 0.3, {
          capture.output({
            py$main()
          })
        })
      } else {
        print('its a test')
      }
    }, error = function(e) {
      results$log <- paste("Error:", e$message)
      print(paste("Error running Python seqnovelty script. ", e$message), type = "error")
      showNotification(paste("Error running Python seqnovelty script. ", e$message), type = "error")
    })
    
    tryCatch({
      # Read the alignment output (e.g., FASTA format with matched residues)
      print('still the same')
      if (!is.null(test()) && test() == 'many') {
        print('still many')
        aligned_sequences <- readAAStringSet("final_alignment_demo.fasta")  # Adjust file name
        # print('done')
      } else if (!is.null(test()) && test() == 'many') {
        print('still few')
        aligned_sequences <- readAAStringSet("final_alignment_few.fasta")  # Adjust file name
      } else {
        aligned_sequences <- readAAStringSet("final_alignment.fasta")  # Adjust file name
      }
      # Build the alignment table for visualization
      ref_name(names(aligned_sequences)[1])
      sequences <- as.character(aligned_sequences)
      seq_names <- names(aligned_sequences)
      
      # Create a data frame where each residue is in its own cell
      alignment_table <- data.frame(
        Name = seq_names,
        do.call(rbind, strsplit(sequences, ""))
      )
      
      # Store the table for rendering
      print('filling results as part of run_analysis')
      results$alignment_table <- alignment_table
      
      # Read JSON for sequence information
      print('still same 2')
      if (!is.null(test()) && test() == 'many') {
        print('still many 2')
        # browser()
        json_file <- "alignment_summary_demo.json"
        # json_file <- file.path(output_dir, "alignment_summary.json")
      } else if (!is.null(test()) && test() == 'few') {
        print('still few 2')
        json_file <- "alignment_summary_few.json"
      } else {
        json_file <- "alignment_summary.json"
        
      }
      if (file.exists(json_file)) {
        results$seq_info <- jsonlite::fromJSON(json_file)
      } else {
        print("JSON file not found. Ensure the Python script generated the summary correctly.", type = "error")
        showNotification("JSON file not found. Ensure the Python script generated the summary correctly.", type = "error")
      }
      # Compute cumulative coverage for each sequence
      ref_length <- sum(results$alignment_table[1, -1] != "-")  # Non-gap residues in reference
      seq_coverage <- data.frame(Sequence = character(), CumulativeCoverage = numeric(), stringsAsFactors = FALSE)
      cumulative_matched <- 0
      
      for (seq_name in seq_names[-1]) {  # Skip reference sequence
        if (!is.null(results$seq_info[[seq_name]]) && !is.null(results$seq_info[[seq_name]]$matched_indices)) {
          matched_indices <- results$seq_info[[seq_name]]$matched_indices +1
          cumulative_matched <- cumulative_matched + length(matched_indices)
          coverage <- (cumulative_matched / ref_length) * 100
          seq_coverage <- rbind(seq_coverage, data.frame(Sequence = seq_name, CumulativeCoverage = coverage))
        }
      }
      
      # Add an initial value for 0% coverage at the start
      # Add an initial value for 0% coverage at the start
      seq_coverage <- rbind(data.frame(Sequence = " ", CumulativeCoverage = 0), seq_coverage)
      
      # Set levels for the Sequence column to ensure the correct order
      seq_coverage$Sequence <- factor(seq_coverage$Sequence, levels = c(" ", seq_names[-1]))
      
      results$coverage_df <- seq_coverage
      
    }, error = function(e) {
      results$log <- paste("Error:", e$message)
      print(paste("Error parsing alignment ", e$message), type = "error")
      showNotification(paste("Error parsing alignment. ", e$message), type = "error")
    })
  }
  
  # Extract sequence IDs from uploaded FASTA file
  observeEvent(input$fasta_file, {
    req(input$fasta_file)
    fasta_file <- input$fasta_file$datapath
    
    # Read the FASTA file and extract sequence IDs
    tryCatch({
      fasta <- readAAStringSet(fasta_file)
      results$seq_ids <- names(fasta)  # Extract sequence IDs
    }, error = function(e) {
      results$seq_ids <- NULL
      print("Error reading FASTA file. Ensure it is correctly formatted.", type = "error")
      showNotification("Error reading FASTA file. Ensure it is correctly formatted.", type = "error")
    })
  })
  observeEvent(input$fasta_file_msa, {
    req(input$fasta_file_msa)
    fasta_file <- input$fasta_file_msa$datapath
    
    # Read the FASTA file and extract sequence IDs
    tryCatch({
      fasta <- readAAStringSet(fasta_file)
      results$seq_ids <- names(fasta)  # Extract sequence IDs
    }, error = function(e) {
      results$seq_ids <- NULL
      print("Error reading FASTA file. Ensure it is correctly formatted.", type = "error")
      showNotification("Error reading FASTA file. Ensure it is correctly formatted.", type = "error")
    })
  })
  
  # Dynamically populate dropdown for sequence IDs
  output$seq_id_dropdown <- renderUI({
    req(results$seq_ids)
    selectInput(
      "seq_id",
      paste(
        length(results$seq_ids),
        "sequences found. Select Reference:"),
      choices = results$seq_ids,
      selected = NULL
    )
  })
  
  # Run the Python script and generate alignment visualization
  observeEvent(input$run_script, {
    test(NULL)
    if (input$mode != 'search') {
      req(input$seq_id)  # Ensure inputs are provided
    }
    # browser()
    if (is.null(input$fasta_file)  && is.null(input$fasta_file_msa)) {
      if (is.null(ref_name())) {
        showNotification("In search mode, you have to run a search (blast / mmseqs) first.", type = "warning")
        return()
      }
      print('looks like blast')
      file = 'search_sequences.fasta'
      # browser()
      print(paste('ref_name() is', ref_name()))
      analysis_return <- run_analysis(file, ref_name())
      if (is.null(analysis_return)) {
        return()
      }
    } else {
      print('file provided')
      if (input$mode == 'msa') {
        file = input$fasta_file_msa$datapath
      } else {
        file = input$fasta_file$datapath
      }
      run_analysis(file, input$seq_id)
    }
  })
  
  # Test Run logic
  observeEvent(input$run_test, {
    print('few')
    # test('few')
    run_analysis("test.fasta", "esmGFP")
  })
  
  observeEvent(input$run_testm, {
    print('many')
    test('many')
    updateTextAreaInput(session, 'seq_input', value = '>esmGFP
MSKVEELIKPDMKMKLEMEGEVNGHKFSIEAEGEGKPYEGKQTIKAWSTTGKLPFAWDILSTSLTYGNRAFTKYPEGLEQHDFFKQSFPEGYSWERTITYEDGATVKVTADISLEDGVLINKVKFKGENFPSDGPVMQKKTTGWEASTELITPDPATGGLKGEVKMRLKLEGGGHLLADFKTTYRSKKKEKLPLPGVHYVDHRIVNEKATHPEGKEYMIQYEHAVARL')
    run_analysis("initial_alignment_demo.fasta", "esmGFP")
  })
  
  # Display log output in the UI
  output$output_log <- renderText({
    req(results$log)
    paste(results$log, collapse = "\n")
  })
  
  observe({
    if (input$dark_mode == 'light') {
      colourpicker::updateColourInput(session = session, inputId = 'novel_col', value = '#610000')
      colourpicker::updateColourInput(session = session, inputId = 'matched_col', value = '#004558')
      colourpicker::updateColourInput(session = session, inputId = 'match_con_col', value = '#000000')
    } else {
      colourpicker::updateColourInput(session = session, inputId = 'novel_col', value = '#ffc3c3')
      colourpicker::updateColourInput(session = session, inputId = 'matched_col', value = '#c3f2ff')
      colourpicker::updateColourInput(session = session, inputId = 'match_con_col', value = '#FFFFFF')
    }
  })
  
  # Render the alignment table
  output$alignment_table <- DT::renderDataTable({
    req(results$alignment_table, results$seq_info, results$seq_info[[ref_name()]])
    withProgress(message = "parsing df", value = 0.7, {
        
      # Make a new df object for modification
      df <- results$alignment_table[, -1]
      
      # Create a logical vector to track non-all-gap columns
      non_gap_columns <- !apply(df, 2, function(x) all(str_detect(x, '-')))
      
      # Keep track of the original indices of the non-gap columns
      original_indices <- which(non_gap_columns)
      
      # if (input$rem_gaps) {
        # Remove the all-gap columns
        df <- df[, non_gap_columns]
      # }
      # Reference sequence name dynamically
      if (!is.null(isolate(input$seq_id))) {
        print(paste('manually setting ref_name to', isolate(input$seq_id)))
        ref_name(isolate(input$seq_id))
      } 
      
      # Apply coloring to non-reference sequences
      seq_color <- input$col
      for (seq_name in rownames(df)[rownames(df) != ref_name()]) {  
        if (!is.null(results$seq_info[[seq_name]]) && !is.null(results$seq_info[[seq_name]]$matched_indices)) {
          
          matched_indices <- results$seq_info[[seq_name]]$matched_indices +1
          
          # if (input$rem_gaps) {
            # Map original indices to new indices
            matched_indices <- which(original_indices %in% matched_indices)       
          # }
          # ref_seq <- results$alignment_table[ref_name(), -1]
          ref_seq <- df[1,]
          # Apply color to matched residues
          # browser()
          df[seq_name, ] <- sapply(seq_len(ncol(df)), function(x) {
            if (!x %in% (matched_indices)) {
              if (df[seq_name, x] == ref_seq[x]) {
                paste0('<span style="color:', input$match_noncon_col, '">', df[seq_name, x], '</span>')
              } else {
                paste0('<span style="color:', input$nomatch_col, '">', df[seq_name, x], '</span>')
              }
            } else {
              paste0('<span style="color:', input$match_con_col, '">', df[seq_name, x], '</span>')
            }
          })
        }
      }
    })
    withProgress(message = "parsing df...", value = 0.8, {
      # Check if seq_info contains the reference name and matched_indices
      if (!is.null(results$seq_info[[ref_name()]]) && !is.null(results$seq_info[[ref_name()]]$matched_indices)) {
        matched_indices <- results$seq_info[[ref_name()]]$matched_indices +1
        
        # if (input$rem_gaps) {
          # Map original indices to new indices
          matched_indices <- which(original_indices %in% matched_indices)      
        # }
        # save the column where the ref seq starts (no gap) for prescrolling
        ref_seq_start_col(min(which(df[1,] != '-')))

        # Apply modifications: Surround matched residues with <b></b>
        
        if (input$boldn) {
          tags_bef <- paste0('<strong><span style="color:', input$novel_col, '">')
          tags_aft <- '</span></strong>'
        } else {
          tags_bef <- paste0('<span style="color:', input$novel_col, '">')
          tags_aft <- '</span>'
        }
        
        ref_bold <- sapply(seq_len(ncol(df)), function(x) {
          if (x %in% (matched_indices)) {  
            # paste0('<strong>', df[1, x], '</strong>')
              paste0('<span style="color:', input$matched_col, '">', df[1, x], '</span')
          } else {
            if (df[1, x] == '-') {
              df[1, x]
            } else {
              paste0(tags_bef, df[1, x], tags_aft)
            }
          }
        })

        
        if (input$ref_nums) {
          ref_seq_wg <- df[1, ]  # input seq
          counter <- 0
          ref_nums <- sapply(ref_seq_wg, function(x) {
            if (x != "-") {  # Check if it's not a "-"
              counter <<- counter + 1  # Increment the counter globally
              if (counter %% 5 == 0) {
                return(paste0('<span style="font-size: 9px;">', counter, '</span>'))  # Return the counter for every 5th position
              } else {
                return(paste0('<span style="font-size: 9px;"> </span>'))  
              }
            } else {
              return(" ")  # Keep "-" for gaps
            }
          })
          
          df[1,] <- paste0(ref_nums,'<br>', ref_bold)
        } else {
          df[1,] <- ref_bold
        }
      } else {
        # browser()
        # If seq_info or matched_indices is missing, log a warning
        print("No matched indices found for the reference sequence.", type = "warning")
        showNotification("No matched indices found for the reference sequence.", type = "warning")
      }
    })
    colnames(df) <- sapply(seq_along(colnames(df)), function(i) {
      if (i %% 5 == 0) {
        as.character(i)  # Show the header text for every 5th column
      } else {
        " "  # Leave other headers blank
      }
    })  
    # browser()
    # df$Sequence <- paste0('<a href="www.test.com" target="_blank">', df$Sequence, '</a>')
    rownames(df)[-1] <- paste0('<a href="https://www.ncbi.nlm.nih.gov/protein/', rownames(df)[-1] ,'" target="_blank">', rownames(df)[-1], '</a>')
    
    df <- cbind('&nbsp', '&nbsp', df)
    df <- rbind(colnames(df), df)
    df[1,1:2] <- ''
    rownames(df)[1] <- ''
    
    # Render the DataTable
    datatable(df,
              escape = FALSE,  # Allows HTML rendering for bold residues
              extensions = c("FixedColumns", "Buttons"),
              selection = 'none',
              # selection = list(mode = 'single', target = "row+column"),
              options = list(
                scrollX = TRUE,  # Enable horizontal scrolling
                paging = FALSE,# Disable paging for a continuous view
                deferRender = TRUE,
                ordering = FALSE, 
                buttons = list('colvis'),  # Adds column visibility buttons
                info = FALSE,
                fixedColumns = list(leftColumns = 3),  # Fix the first column
                searching = FALSE,  # Disable search for simplicity
                columnDefs = list(
                  # dom = 't',
                  list(targets = 3:ncol(df), width = "22px"),  # Set residue columns to 30px
                  list(targets = 2, className = "spacer-col"),  # Set residue columns to 30px
                  list(targets = 1, className = "spacer-col2"),  # Set residue columns to 30px
                  list(
                    targets = 0,  
                    className = "rowname-col",
                    render = JS("function(data, type, row) { return '<div><span>' + data + '</span></div>'; }")
                  )
                ),
                headerCallback = JS(
                  "function(thead, data, start, end, display){",
                  "  $(thead).remove();",
                  "}"),
                rowCallback = JS(
                  "function(row, data, index) {",
                  "  if (index === 0) {",  # Target the first row (index starts at 0)
                  "    $(row).find('td').css({",
                  paste0("      'border-bottom': '0.5px solid ", ifelse(input$dark_mode == 'dark', 'white', 'black'),"',"),
                  "      'border-top': '0px solid transparent !important',",
                  "      'letter-spacing': '-1px',",
                  "      'overflow': 'visible',",
                  "      'font-family': 'Consolas',",
                  "      'font-size': '9px'",
                  "    });",
                  "  }",
                  "  if (index === 1) {",  # Target the second row (index starts at 0)
                  "    $(row).find('td').css({",
                  paste0("      'border-bottom': '0.5px solid ", ifelse(input$dark_mode == 'dark', 'white', 'black'),"',"),
                  "    });",
                  "  }",
                  "}"
                )
              )
    ) 
  })
  
  observe({
    session$sendCustomMessage("scrollTableToColumn", ref_seq_start_col() -20)  # Adjust column index
  })
  
  observe({
    if (!is.null(results$alignment_table) && is.data.frame(results$alignment_table) && !is.null(results$seq_info[[ref_name()]])) {
      # print('clearly, results are not null:')
      # print(dim(results$alignment_table))
      shinyjs::show('hide_scroll')
      shinyjs::show('save_hide')
      updateActionButton(session, 'scroll', label = HTML(paste0("<small><i>jump to <strong>", ref_name() ,"</strong> start →</i></small>")))
    } else {
      # print('clearly, results are NULL:')
      # print(dim(results$alignment_table))
      shinyjs::hide('save_hide')
      shinyjs::hide('hide_scroll')
    }
  })
  observeEvent(input$scroll,{
    print(ref_seq_start_col() -20)
    session$sendCustomMessage("scrollTableToColumn", ref_seq_start_col() -20)  # Adjust column index
  })
  
  save_btn <- function(id) {
    tags$div(tags$span(
      'To retrieve the result later, use this code in the \'restore\' mode: ',
      
      tags$div(class = "code-box",
               tags$button(id = "copy-btn", class = "copy-btn", onclick="copyToClipboard()", "Copy"),
               tags$code(id="code-text", id)
      ), 
      
      
      # JavaScript function to copy text
      tags$script(HTML('
        function copyToClipboard() {
            var codeText = document.getElementById("code-text").innerText;
            var copyButton = document.getElementById("copy-btn");
            
            navigator.clipboard.writeText(codeText).then(() => {
                copyButton.innerText = "Copied!";
                copyButton.style.backgroundColor = "#28a745"; // Green color
                setTimeout(() => {
                    copyButton.innerText = "Copy";
                    copyButton.style.backgroundColor = "#007bff"; // Reset color
                }, 2000);
            }).catch(err => {
                console.error("Failed to copy: ", err);
            });
        }
    '))
    ))
  }
  
  observe({
    if (input$restore_id != '') {
      shinyjs::enable('restore_btn')
    } else {
      shinyjs::disable('restore_btn')
    }
  })
  
  observeEvent(input$restore_btn, {
    # Load it back
    if (file.exists(paste0(input$restore_id, ".rds"))) {
      import_data <- readRDS(paste0(input$restore_id, ".rds"))
      import_results <- import_data[[1]]
      ref_name(import_data[[2]])
      print('filling results with restore')
      results$alignment_table <- import_results$alignment_table
      results$coverage_df <- import_results$coverage_df
      results$seq_info <- import_results$seq_info
      results$seq_ids <- import_results$seq_ids
      updateActionButton(session, 'save', label = 'result saved')
      saveID(input$restore_id)
    } else {
      sendSweetAlert(session, title = 'Not found!', text = 'This code is not recognized.')
    }
    })
    
  observeEvent(input$save, {
    if (is.null(saveID())) {
      id <- paste0(sample(c(LETTERS, letters, 0:9,0:9)), collapse = '')
      saveID(id)
      sendSweetAlert(
        title ="Results saved!",
        html = TRUE,
        text = save_btn(id)
      )
      updateActionButton(session, 'save', label = 'result saved')
      
      #save logic
      # Save an object
      saveRDS(list(results, ref_name()), file = paste0(id, ".rds"))

    } else {
      sendSweetAlert(
        title = NULL,
        html = TRUE,
        text = save_btn(saveID())
      )
    }
  })
  
  observeEvent(input$help, {
    runjs(
      '$("#mode").tooltip({title: "First choose the input mode. If you want the app to look for similar sequences for you choose <strong>search</strong> mode. If you have your own set of similar sequences, unaligned, select <strong>align</strong> mode. If you have your own MSA, select <strong>MSA</strong> mode.", 
      placement: "right", html: true}).tooltip("show");'
    )
    if (input$mode == 'search') {
      runjs(
        '$("#seq_input").tooltip({title: "Enter a protein sequence to search. Either Fasta format with a >header on the first line and amino acid residues on subsequent lines, or sequence only.", 
      placement: "right", html: true}).tooltip("show");'
      )
      runjs(
        '$("#searchmode").tooltip({title: "Choose a search tool and adapt search parameters. For the most extensive search, choose BLAST and the nr database. For faster, less comprehensive search in more restricted database, you cab also choose mmseqs. Note that mmseqs may return sequences quite different from the input, occassionally resulting in unrealistic alignments.", 
      placement: "right", html: true}).tooltip("show");'
      )
    }
    runjs(
      '$("#align_mode").tooltip({title: "Choose an alignment mode. In pairwise mode, each hit sequence is aligned separately to the target, in MSA mode, one multiple sequence alignment with target and all hits is generated. MSAs run externally on <a href=\'https://www.ebi.ac.uk/jdispatcher/msa\' target=\'_blank\'>EBI servers which imposes particular limits</a> but can sometimes give better alignments.", 
      placement: "right", html: true}).tooltip("show");'
    )
  })
  observeEvent(input$about, {
    sendSweetAlert(
      session = session,
      html = TRUE,
      title = "de novo Sequence Novelty",
      text = tags$span(style="text-align:left", "This app determines the novelty of a given protein sequence by aligning it to homologs and reporting identical residues.",
                       tags$br(),tags$br(),"In Search mode, the app will find sequences from the selected database for you.",
                       tags$br(),tags$br(),"In align mode, the app will align and only use the sequences provided in the input.",
                       tags$br(),tags$br(),"In MSA mode, the app will use a provided MSA and only perform the analysis.",
                       tags$br(),tags$br(), "Bugs? Problems? Requests?", tags$a(href = 'https://github.com/mf-rug/SeqNovelty/issues', target='_blank', 'Raise an issue on Github!')
                       ),
      type = "info"
    )
  })
  
  # Plot cumulative coverage
  output$cum_coverage <- renderPlot({
    req(results$alignment_table, results$seq_info, results$seq_info[[ref_name()]])
    df <- results$coverage_df
    dark <- input$dark_mode == 'dark'

    ref_seq <- results$alignment_table[ref_name(), -1]
    ref_len <- length(ref_seq[!str_detect(ref_seq, '-')])
    match_len <- length(results$seq_info[[ref_name()]]$matched_indices)
    
    orig_levels <- levels(df$Sequence)
    
    # Convert to character and shorten strings
    df$Sequence <- paste0(substr(as.character(df$Sequence), 1, 20), 
                          ifelse(nchar(as.character(df$Sequence)) > 20, '...', ''))
    
    shortened_levels <- unique(paste0(substr(orig_levels, 1, 20), 
                                      ifelse(nchar(orig_levels) > 20, '...', '')))
    
    # Convert back to a factor, maintaining original level order
    df$Sequence <- factor(df$Sequence, levels = shortened_levels)
    
    
    p1 <- ggplot(df, aes(x = Sequence, y = CumulativeCoverage)) +
      geom_hline(yintercept = max(df$CumulativeCoverage), linetype = 2, color='#58aebf', linewidth= 1) +
      geom_line(group = 1, color = ifelse(dark, '#d9d1ff', "#0f0059"), linewidth=1.5) +
      geom_point(data = df[-1,], color = "#6f51ff", size=10.5) +
      geom_text(data = df[-1,], aes(y = CumulativeCoverage, label = round(CumulativeCoverage,1)), color = "white", size=4) +
      # geom_text(x= df$Sequence[length(df$Sequence)], 
      #           y=max(df$CumulativeCoverage), 
      #           color= '#58aebf',
      #           fontface = 'bold',
      #           label = paste0('novelty:\n',  
      #                          ref_len - match_len, ' / ', ref_len, '\n',
      #                          '(', round(100-max(df$CumulativeCoverage), 1), '%)'), 
      #           size=5.5,
      #           vjust = 0.5, hjust=-0.1) +
      coord_cartesian(clip = "off") +  # Prevent clipping of text outside the plot area
      scale_x_discrete(expand = c(0,0)) +
      scale_y_continuous(limits = c(0,100), breaks = c(seq(0, 100, 20)), expand = c(0,0)) +
      labs(x = "Sequences", y = "Cumulative Coverage (%)") +
      theme_bw() +
      theme(text = element_text(size = 16),
            plot.margin = margin(15, 15, 15, 15),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
      )
    if (dark) {
      p1 <- p1 + theme(text = element_text(color = 'white'),
                       axis.text = element_text(color = 'white'),
                       panel.background = element_rect(fill = 'black'),
                       plot.background = element_rect(fill = '#1d1f21', size = 0),
                       panel.grid.major = element_line(color='#94949480'),
                       panel.grid.minor = element_line(color='#94949480'))
    }
    p1
  })
  
  
  output$novelty_media <- renderUI({
    req(results$alignment_table, results$seq_info, results$seq_info[[ref_name()]])
    ref_seq <- results$alignment_table[ref_name(), -1]
    sequence <- paste0(ref_seq, collapse = '') %>% str_remove_all(., '-')
    if (str_count(sequence) <= 400) {
      structureUI <- actionButton("structure", HTML("&nbspPredict structure & highlight"), icon = icon("robot"))
    } else {
      structureUI <- fileInput(
        "structure_upload",
        HTML('&nbspToo large to esmFold, upload own model <i class="fa-solid fa-file"></i>'),
        placeholder = 'select a .pdf or .cif file',
        accept = c('.pdb', '.cif'),
        width = '100%'
      )
    }
    
    renderUI(
      fluidRow(column(6,
                      div(class = 'error_out', plotOutput('cum_coverage', height = '60vh'))),
               column(6,
                      div(class = 'error_out', style = "text-align:center;", 
                          structureUI,
                          br(),br(),
                          conditionalPanel(
                            condition = "input.structure > 0",
                            style = "display: none;",
                            uiOutput("structure_view", height = "100%")
                          )
                         )
               )
      )
    )
  })
  
  output$novelty_text <- renderUI({
    req(results$alignment_table, results$seq_info, results$seq_info[[ref_name()]])
    ref_seq <- results$alignment_table[ref_name(), -1]
    ref_len <- length(ref_seq[!str_detect(ref_seq, '-')])
    match_len <- length(results$seq_info[[ref_name()]]$matched_indices)

    label = paste0('Sequence novelty:<br>', hl(ref_len - match_len, ' / ', ref_len),
                   ' amino acids',
                   hl(' (',round(100 - max(results$coverage_df$CumulativeCoverage), 1), '%)'))
    renderUI(HTML(
      paste0(
        '<div class="error_out" style="height:100%;">',
        label,
        '</div>'
      )
    ))
  })  
  
  
  output$res_analysis <- renderUI({
    req(results$alignment_table, results$seq_info, results$seq_info[[ref_name()]])
    
    selrow <- input$row_selected
    selcol <- input$col_selected -2 # correct indexing and first col is name

    df <- results$alignment_table
    non_gap_columns <- !apply(df, 2, function(x) all(str_detect(x, '-')))
    
    # Remove the all-gap columns
    df <- df[, non_gap_columns]
    
    # browser()
    if (is.null(selrow) || length(df$Name[selrow]) == 0) {
      info <- '<i>Select a residue for more info.</i>'
    } else {
      if (selcol < 1) {
        info <- paste0(
          "Selected: ",
          ifelse(df$Name[selrow] == ref_name(), "analysed sequence ", "hit sequence "),
          tags$strong(df$Name[selrow])
        )
      } else {
        selseq <- df[selrow, -1]
        selres <- selseq[selcol]
        other_seqs_res <- df[-selrow, selcol + 1]
        num_identical <- length(which(other_seqs_res %in% selres))
        perc_identical <- round(num_identical / length(other_seqs_res) * 100,1)
        
        # browser()
        if (selres != '-') {
          selected_col_name <- names(selres)
          match_index <- which(names(selseq) == selected_col_name)
          non_gap_indices <- which(selseq != "-")
          selected_non_gap <- which(non_gap_indices == match_index)
          selection_str <- paste0(', residue ', tags$strong(selected_non_gap))
        } else {
          selection_str <- ''
        }
        info <- paste0(
          "Selected: ",
          ifelse(df$Name[selrow] == ref_name(), "analysed sequence ", "hit sequence "),
          tags$strong(df$Name[selrow]),
          selection_str,
          '<br>Selected ',
          ifelse(selres == '-', 'gap', paste0('amino acid ', hl(selres))),
          ', found in ',
          hl(
            num_identical,
            ' / ',
            length(other_seqs_res),
            ' (',
            perc_identical,
            '%)'
          ),
          ' of aligned sequences',
          ifelse(df$Name[selrow] == ref_name(), 
                 ifelse(
                   perc_identical == 0,
                   paste0('<br>→ <span style="color:green;"> novel',
                          ifelse(all(other_seqs_res == '-'), 
                                 ' insertion', 
                                 ifelse(selres == '-' && all(other_seqs_res != '-'), 
                                        ' deletion', 
                                        ' mutation'))),
                   '<br>→ <span style="color:red;">not novel'
                 ), '')
        )
        }
      }
    
    renderUI(HTML(
      paste0('<div class="error_out" style="height:100%;">', info, '</div>')
    ))
    
  })
  
  
  observeEvent(input$structure_upload, {
    print('struc upload')
    ref_seq <- results$alignment_table[1,-1]
    sequence <- paste0(ref_seq, collapse = '') %>% str_remove_all(., '-')
    
    matched_col <-  results$seq_info[[ref_name()]]$matched_indices +1
    matched_resi <- convert_col_to_residue(ref_seq, matched_col)
    
    output$structure_view <- renderUI(r3dmolOutput('molView'))
    
    output$molView <- renderR3dmol({
      pdb_result <- input$structure_upload$datapath
      print(pdb_result)
      r3dmol() %>%
        
        m_add_model(pdb_result) %>%
        m_set_style(style = m_style_cartoon()) %>%
        m_zoom_to() %>% 
        m_add_style(
          style = m_style_cartoon(color = input$matched_col),
          sel = m_sel(resi = matched_resi)
        ) %>% 
        m_add_style(
          style = m_style_cartoon(color = input$novel_col),
          sel = m_sel(resi = setdiff(1:str_count(sequence), matched_resi))
        ) %>% 
        m_set_background_color(hex = ifelse(input$dark_mode == 'dark', '#000000', '#ffffff'))
    })
  })
  
  observeEvent(input$structure, {
    print('predicting structure')
    shinyjs::disable('structure')
    ref_seq <- results$alignment_table[1,-1]
    sequence <- paste0(ref_seq, collapse = '') %>% str_remove_all(., '-')
    
    matched_col <-  results$seq_info[[ref_name()]]$matched_indices +1
    matched_resi <- convert_col_to_residue(ref_seq, matched_col)
    
    # browser()
    
    if (str_count(sequence) <= 400) {
      # Call esmfold API
      response <- POST("https://api.esmatlas.com/foldSequence/v1/pdb/",
                       body = sequence, encode = "raw", timeout(60), progress())
      if (status_code(response) == 200) {
        print('success! got a response from esm server')
        pdb_result <- content(response, "text")
        print('downloaded pdb')
      
        # Save PDB structure temporarily
        pdb_file <- tempfile(fileext = ".pdb")
        writeLines(pdb_result, pdb_file)
        
        # Render the PDB in 3Dmol.js
        output$structure_view <- renderUI(r3dmolOutput('molView'))
        
        structure(pdb_file)

      } else {
        cat('oh no something went wrong with esmfold!:', status_code(response), '\n')
        output$molView <- renderText("Error fetching structure.")
        sstop(session,TRUE, 'Error modeling structure', 'The esmFold server may be down. Try again later.', 'error')
        shinyjs::hide('structure_view')
        shinyjs::hide('molView')
      }
    } else {
      print('nahh')
      output$structure_view <- renderUI(uiOutput('molView'))
      output$molView <- renderUI(
        renderUI(
          HTML(
            "<br><strong>Error:</strong><br>Structure prediction in this app is limited to sequences below 400 amino acids."
          )
        )
      )
    }
  })


  output$molView <- renderR3dmol({
    if (structure() != '') {
      pdb_result <- structure()
      ref_seq <- results$alignment_table[1,-1]
      sequence <- paste0(ref_seq, collapse = '') %>% str_remove_all(., '-')
      matched_col <-  results$seq_info[[ref_name()]]$matched_indices +1
      matched_resi <- convert_col_to_residue(ref_seq, matched_col)
      
      selrow <- input$row_selected
      selcol <- input$col_selected -2 # correct indexing and first col is name
      
      df <- results$alignment_table
      non_gap_columns <- !apply(df, 2, function(x) all(str_detect(x, '-')))
      
      # Remove the all-gap columns
      df <- df[, non_gap_columns]
      
      
      ref_seq_wg <- df[1, ]  # input seq
      counter <- 0
      ref_nums <- sapply(ref_seq_wg, function(x) {
        if (x != "-") {  # Check if it's not a "-"
          counter <<- counter + 1  # Increment the counter globally
            return(counter)  # Return the counter for every 5th position
        } else {
          return(" ")  # Keep "-" for gaps
        }
      })
      
      # browser()
      out <<- r3dmol() %>%
        m_add_model(pdb_result, format = "pdb") %>%
        m_set_style(style = m_style_cartoon()) %>%
        m_zoom_to() %>% 
        m_add_style(
          style = m_style_cartoon(color = input$matched_col),
          sel = m_sel(resi = matched_resi)
        ) %>% 
        m_add_style(
          style = m_style_cartoon(color = input$novel_col),
          sel = m_sel(resi = setdiff(1:str_count(sequence), matched_resi))
        ) %>% 
        m_set_background_color(hex = ifelse(input$dark_mode == 'dark', '#000000', '#ffffff'))
      
      if (length(selrow) > 0) {
        out <- out %>% 
          m_add_style(
            style = m_style_stick(color = '#de00c6'),
            sel = m_sel(resi = as.numeric(ref_nums[selcol]))
          )
      }
      out
      }
  })
})