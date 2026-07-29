library(shiny)
library(shinydashboard)
library(DT)
library(dplyr)
library(tidyr)
library(stringr)

# Allows file upload up to (até 30 MB)
options(shiny.maxRequestSize = 30 * 1024^2)

# ==============================================================================
# UI (User Interface)
# ==============================================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Epitope Sniper"),
  dashboardSidebar(
    sidebarMenu(
      id = "sidebarMenu",
      menuItem("Binder Analysis", tabName = "analysis", icon = icon("dna")),
      menuItem("How to Cite", tabName = "cite", icon = icon("quote-right"))
    ),
    hr(),
    conditionalPanel(
      condition = "input.sidebarMenu == 'analysis'",
      fileInput("file_xml", "Upload NetMHCpan File (.xls / .csv / .txt):",
                accept = c(".xml", ".xls", ".txt", ".csv", ".out")),
      hr(),
      h4(" Rank Thresholds (%)", style = "margin-left: 15px; color: #fff;"),
      numericInput("sb_cutoff", "Strong Binder (Rank ≤ %):", value = 0.500, min = 0, max = 100, step = 0.1),
      numericInput("wb_cutoff", "Weak Binder (Rank ≤ %):", value = 2.000, min = 0, max = 100, step = 0.1),
      hr(),
      uiOutput("allele_filter_ui"),
      uiOutput("filter_ui"),
      hr(),
      h4(" Export Options", style = "margin-left: 15px; color: #fff;"),
      downloadButton("download_csv", "Download CSV", class = "btn-success", 
                     style = "margin-left: 15px; margin-bottom: 5px; width: 85%;"),
      downloadButton("download_fasta", "Download FASTA", class = "btn-info", 
                     style = "margin-left: 15px; margin-bottom: 15px; width: 85%;")
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "analysis",
              fluidRow(
                valueBoxOutput("total_box", width = 3),
                valueBoxOutput("sb_box", width = 3),
                valueBoxOutput("wb_box", width = 3),
                valueBoxOutput("promiscuous_box", width = 3)
              ),
              fluidRow(
                box(
                  title = "Epitope Table & Binding Affinity", 
                  status = "primary", 
                  solidHeader = TRUE, 
                  width = 12,
                  DTOutput("netmhc_table")
                )
              )
      ),
      tabItem(tabName = "cite",
              fluidRow(
                box(
                  title = "How to Cite NetMHCPan Epitope Sniper",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  h3("Citation Information"),
                  p("NetMHCPan Epitope Sniper: An Interactive Shiny Application for High-Throughput Epitope Screening."),
                  wellPanel(tags$code("Queiroz Jr., et al. (2026). NetMHCPan Epitope Sniper."))
                )
              )
      )
    )
  )
)

# ==============================================================================
# SERVER (Application Logic)
# ==============================================================================
server <- function(input, output, session) {
  
  # 1. Parsing of file NetMHCpan
  parsed_long_data <- reactive({
    req(input$file_xml)
    path <- input$file_xml$datapath
    
    tryCatch({
      lines <- readLines(path, warn = FALSE)
      
      # Search lines on header with col (Pos, Peptide, ID, core...)
      header_line_idx <- grep("^Pos\tPeptide|^Pos ", lines)[1]
      
      if (is.na(header_line_idx)) {
        stop("Formato não reconhecido. Certifique-se de carregar um arquivo de saída válido do NetMHCpan.")
      }
      
      # Allele extraction by anterior header line 
      allele_line_idx <- header_line_idx - 1
      extracted_alleles <- c()
      
      if (allele_line_idx >= 1) {
        raw_tokens <- unlist(strsplit(lines[allele_line_idx], "\t"))
        clean_tokens <- trimws(raw_tokens)
        extracted_alleles <- clean_tokens[clean_tokens != ""]
      }
      
      # Fallback if direct extracion of line 2 fails
      if (length(extracted_alleles) == 0) {
        header_text <- paste(lines[1:header_line_idx], collapse = " ")
        extracted_alleles <- unlist(str_extract_all(
          header_text, 
          "(?:HLA|SLA|BoLA|Eqca|Patr|Mamu|Gogo|H-2)[-A-Za-z0-9:*]+"
        ))
        extracted_alleles <- unique(extracted_alleles)
      }
      
      # Table reading by tabulation
      df <- read.delim(
        text = paste(lines[header_line_idx:length(lines)], collapse = "\n"),
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
      colnames(df) <- make.unique(colnames(df))
      
      # Keep only numeric valid pos
      df <- df[suppressWarnings(!is.na(as.numeric(as.character(df$Pos)))), ]
      
      # Peptides Col Identification
      pep_col_idx <- which(grepl("^Peptide$", colnames(df), ignore.case = TRUE))[1]
      if (is.na(pep_col_idx)) pep_col_idx <- 2
      
      # Rank Col Identification
      rank_col_indices <- which(grepl("^EL_rank", colnames(df), ignore.case = TRUE))
      if (length(rank_col_indices) == 0) {
        rank_col_indices <- which(grepl("^BA_rank|^Rank", colnames(df), ignore.case = TRUE))
      }
      
      # Rank col Mapping and Allele Mapping
      long_list <- list()
      for (i in seq_along(rank_col_indices)) {
        col_idx <- rank_col_indices[i]
        
        current_allele <- if (i <= length(extracted_alleles)) {
          extracted_alleles[i]
        } else {
          paste0("Allele_", i)
        }
        
        sub_df <- data.frame(
          Pos = suppressWarnings(as.numeric(as.character(df$Pos))),
          Peptide = df[[pep_col_idx]],
          Allele = current_allele,
          Rank = suppressWarnings(as.numeric(as.character(df[[col_idx]]))),
          stringsAsFactors = FALSE
        )
        long_list[[i]] <- sub_df
      }
      
      df_out <- bind_rows(long_list) %>%
        filter(!is.na(Rank) & !is.na(Peptide) & Peptide != "")
      
      return(df_out)
      
    }, error = function(e) {
      showNotification(paste("Erro ao ler o arquivo:", e$message), type = "error")
      return(NULL)
    })
  })
  
  # 2. Ligands Class and Promiscous Evaluation
  processed_data <- reactive({
    req(parsed_long_data(), input$sb_cutoff, input$wb_cutoff)
    
    df <- parsed_long_data()
    sb_thresh <- as.numeric(input$sb_cutoff)
    wb_thresh <- as.numeric(input$wb_cutoff)
    
    df_classified <- df %>%
      mutate(
        BindType = case_when(
          Rank <= sb_thresh ~ "Strong Binder",
          Rank > sb_thresh & Rank <= wb_thresh ~ "Weak Binder",
          TRUE ~ "Non-Binder"
        )
      )
    
    res_df <- df_classified %>%
      group_by(Peptide) %>%
      mutate(
        Bound_Alleles_Count = n_distinct(Allele[BindType %in% c("Strong Binder", "Weak Binder")]),
        Promiscuous = ifelse(Bound_Alleles_Count >= 2, "Yes", "No")
      ) %>%
      ungroup()
    
    return(res_df)
  })
  
  # 3. Dynamic Components (UI)
  output$allele_filter_ui <- renderUI({
    req(processed_data())
    df <- processed_data()
    unique_alleles <- sort(unique(df$Allele))
    
    selectizeInput(
      "selected_alleles",
      "Filter by Allele(s):",
      choices = c("All Alleles", unique_alleles),
      selected = "All Alleles",
      multiple = TRUE
    )
  })
  
  output$filter_ui <- renderUI({
    req(processed_data())
    radioButtons(
      "bind_filter",
      "Display Binders:",
      choices = c("Only Binders (SB + WB)" = "binders",
                  "Strong Binders Only (SB)" = "SB",
                  "Weak Binders Only (WB)" = "WB",
                  "Promiscuous Only (≥2 Alleles)" = "promiscuous",
                  "All Peptides (SB + WB + NB)" = "all"),
      selected = "binders"
    )
  })
  
  # 4. Data Filtering 
  filtered_data <- reactive({
    req(processed_data(), input$bind_filter)
    df <- processed_data()
    
    if (!is.null(input$selected_alleles) && !"All Alleles" %in% input$selected_alleles) {
      df <- df %>% filter(Allele %in% input$selected_alleles)
    }
    
    if (input$bind_filter == "SB") {
      df <- df %>% filter(BindType == "Strong Binder")
    } else if (input$bind_filter == "WB") {
      df <- df %>% filter(BindType == "Weak Binder")
    } else if (input$bind_filter == "binders") {
      df <- df %>% filter(BindType %in% c("Strong Binder", "Weak Binder"))
    } else if (input$bind_filter == "promiscuous") {
      df <- df %>% filter(Promiscuous == "Yes" & BindType %in% c("Strong Binder", "Weak Binder"))
    }
    
    return(df)
  })
  
  # 5. Main Metrics (ValueBoxes)
  output$total_box <- renderValueBox({
    total <- if (!is.null(filtered_data())) nrow(filtered_data()) else 0
    valueBox(total, "Total Entries", icon = icon("list"), color = "black")
  })
  
  output$sb_box <- renderValueBox({
    sb_count <- if (!is.null(filtered_data())) sum(filtered_data()$BindType == "Strong Binder", na.rm = TRUE) else 0
    valueBox(sb_count, "Strong Binders (SB)", icon = icon("fire"), color = "red")
  })
  
  output$wb_box <- renderValueBox({
    wb_count <- if (!is.null(filtered_data())) sum(filtered_data()$BindType == "Weak Binder", na.rm = TRUE) else 0
    valueBox(wb_count, "Weak Binders (WB)", icon = icon("leaf"), color = "green")
  })
  
  output$promiscuous_box <- renderValueBox({
    prom_count <- if (!is.null(filtered_data())) n_distinct(filtered_data()$Peptide[filtered_data()$Promiscuous == "Yes"]) else 0
    valueBox(prom_count, "Promiscuous Epitopes", icon = icon("star"), color = "purple")
  })
  
  # 6. Interactive DT Table
  output$netmhc_table <- renderDT({
    req(filtered_data())
    
    datatable(
      filtered_data(),
      extensions = 'Buttons',
      options = list(
        pageLength = 25, 
        scrollX = TRUE, 
        autoWidth = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      rownames = FALSE
    ) %>%
      formatStyle(
        'BindType',
        target = 'row',
        backgroundColor = styleEqual(c("Strong Binder", "Weak Binder"), c('#ffcccc', '#d4edda')),
        fontWeight = styleEqual(c("Strong Binder", "Weak Binder"), c('bold', 'normal'))
      )
  })
  
  # 7. Exporting Results (CSV AND FASTA)
  output$download_csv <- downloadHandler(
    filename = function() { paste0("Epitope_Sniper_Results_", Sys.Date(), ".csv") },
    content = function(file) { write.csv(filtered_data(), file, row.names = FALSE) }
  )
  
  output$download_fasta <- downloadHandler(
    filename = function() { paste0("Epitope_Sniper_Sequences_", Sys.Date(), ".fasta") },
    content = function(file) {
      df <- filtered_data()
      req(nrow(df) > 0)
      
      headers <- sprintf(">Epitope_%d|Pos_%d|Peptide_%s|Allele_%s|Type_%s|Rank_%.3f",
                         seq_len(nrow(df)), df$Pos, df$Peptide, df$Allele, df$BindType, df$Rank)
      
      fasta_lines <- as.vector(rbind(headers, df$Peptide))
      writeLines(fasta_lines, file)
    }
  )
}

# ==============================================================================
# Launch App
# ==============================================================================
shinyApp(ui = ui, server = server)
