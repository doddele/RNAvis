library(shiny)
library(shinydashboard)
library(data.table)
library(shinyFiles)
library(tidyverse)
library(dplyr)
library(tidyverse)
library(plotly)





ui <- dashboardPage(
  dashboardHeader(title = "RNA visualization",
                  
                  # Permet d'afficher des notifications réactives.
                  dropdownMenuOutput("notificationMenu")
                  
                  ), #dashboardHeader
  
  dashboardSidebar(
    sidebarMenu( ##Le texte ne dépasse pas de la sidebar.
      
      
      br(),
      
      h4(HTML("&nbsp; 1 | Select a 4_Tail folder")),
      
      #HTML("&nbsp; Select the 4_Tail folder where are<br>&nbsp; located your files of interest."),
      
      
      shinyDirButton("dir", 'Select a folder', 'Please select a 4_Tail folder  (use the arrows)',
                     icon=icon(name="", class="fa-regular fa-folder-open")),
    
      br(),
      
      h4(HTML("&nbsp; 1.5 | Create a tabix index")),
      h5(HTML("&nbsp;&nbsp; Only run if the file is not yet indexed<br>&nbsp;&nbsp; by tabix"), style = "color:grey"),
      
      actionButton("tabix", class="btn-warning",
                   "Create Tabix index"),
    
      br(),
      
      h4(HTML("&nbsp; 2 | Enter an AGI to search for")),
      textInput("AGI",
                label=NULL
      ), #textInput
      
      helpText("AT4G31830.1"),
      
      br(),
      
      h4(HTML("&nbsp; 3 | Select genotypes")),
      h5(HTML("&nbsp;&nbsp; If genotypes do not appear, you may<br>&nbsp;&nbsp;  need to create a tabix index"), style = "color:grey"),
      
      checkboxGroupInput("check_barcode",
                       "Genotypes:",
                       choices = NULL
      ), #checkboxGroupInput
      
     br(),
     
     h4(HTML("&nbsp; 4 | Update")),
     
     actionButton("update_df",
                  "Update",
                  icon=icon(name="", class="fa-solid fa-arrows-rotate")),
     
     br(),
     
     h4(HTML("&nbsp; 5 | Upload the DataTable")),
     
     actionButton("save_df",
                  "Upload DataTable",
                  icon = icon(name="", class="fa-solid fa-download"))
    ) #sidebarMenu
  ), #dashboardSidebar
  
  dashboardBody(
    fluidRow(
      tabsetPanel(
                  tabPanel("mRNA statistics",
      tabsetPanel(
        title = NULL,
        # The id lets us use input$tabset1 on the server to find the current tab
        id = "tabset1",
        
        tabPanel("DataTable",
                 dataTableOutput("df_file")
                 ),
        
        tabPanel("Poly(A) distribution",
                 plotOutput("plot_6"),
                 br(),
                 plotOutput("plot_polyA"),
                 br(),
                 sliderInput("slider", "x-axis", min=0, max=100,
                             value=c(0,100))
        ),
        
        tabPanel("Poly(U) distribution",
                 plotOutput("plot_10"),
                 br(),
                 plotOutput("plot_11"),
                 br(),
                 fluidRow(
                   column(3,
                          sliderInput("slider_polyU", "x-axis", min=0, max=100,
                                      value=c(0,100))
                          ),
                   column(3,
                          numericInput("cutoff_U_distrib", "U cut off", 70)
                   )
                 )
        ),
        
        tabPanel("Percentages of nucleic acids in additional tails",
                 plotOutput("plot_3")
        ),
        
        tabPanel("Genome coordinates",
                 plotOutput("plot_1")
                 ),
        
        tabPanel("5' position",
                 plotOutput("plot_gen_start")
        ),
        
        tabPanel("3' position",
                 plotOutput("plot_1A")
                 ),

        tabPanel("polyA/polyU",
                 fluidRow(
                   column(6,
                          plotOutput("plot_A"),
                          #numericInput("cutoff_A", "Cut off A", 70)

                   ),
                   column(6,
                          plotOutput("plot_U"),
                          #numericInput("cutoff_U", "Cut off U", 70)
                          
                   )
                 ),
                 br(),
                 fluidRow(
                   column(6,
                          plotOutput("plot_C"),
                          #numericInput("cutoff_C", "Cut off C", 70)
                          ),
                   column(6,
                          plotOutput("plot_G"),
                          #numericInput("cutoff_G", "Cut off G", 70),
                          
                   )
                   ),
                 box(column(6,
                            numericInput("cutoff_A", "Cut off A", 70),
                            numericInput("cutoff_C", "Cut off C", 70)
                            ),
                     column(6,
                            numericInput("cutoff_U", "Cut off U", 70),
                            numericInput("cutoff_G", "Cut off G", 70)
                            ),
                     sliderInput("AU_slider", "x-axis", min=0, max=100,
                             value=c(0,100)))
                 ) #tabPanel
      ) #tabsetPanel
      ) #tabPanel
      ) #tabsetPanel
    ),
    fluidRow(infoBoxOutput("tabset1Selected"))
  )
) #ui





server <- function(input, output, session) {
  

  # Fait débuter le browser de fichiers dans le home.
  shinyDirChoose(
    input,
    'dir',
    roots = c(home = '~')
  ) #shinyDirChoose
  
  
  global <- reactiveValues(datapath = getwd())
  

  dir <- reactive(input$dir)
  
  
  output$dir <- renderText({ 
    global$datapath 
    }) #output$dir
  
  
  # Met a jour le path du dossier choisi.
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {input$dir},
               handlerExpr = {
                 if (!"path" %in% names(dir())) return()
                 home <- normalizePath("~")
                 global$datapath <- file.path(home,
                                              paste(unlist(dir()$path[-1]),
                                                    collapse = .Platform$file.sep)
                                              )
               }) #observeEvent

  
  # PAS ENCORE UTILISE
  barcode_correspondance <- reactive({
    # Path du dossier contenant le fichier de correspondance des barcodes.
    path_minus_4_Tail <- gsub("/4_Tail", "", global$datapath);
    
    # Stocke le fichier de correspondance des barcodes.
    check_if_barcode_correspondance <- list.files(path=path_minus_4_Tail,
                                                  pattern = ".tsv$");
    
    # Lis les génotypes associés aux barcodes.
    barcode_correspondance <- fread(paste0(path_minus_4_Tail,
                                           "/",
                                           check_if_barcode_correspondance),
                                    header=FALSE,
                                    col.names=c("Barcodes", "Genotypes"));
  }) #barcode_correspondance
  
  
  # Actualise les checkboxes avec les noms des barcodes présents dans le dossier.
  observeEvent(
    input$dir, {
      
      # Stocke uniquement les fichiers ".gz" dans un dataframe.
      barcode_file <- list.files(path=global$datapath, pattern = ".gz$");
      barcode_file_df <- data.frame(barcode_file);
      
      # Affiche les barcodes sans le reste du nom du fichier.
      barcode_choices <-  strsplit(barcode_file_df[[1]], split = "$");
      barcode_choices <- gsub("\\..*", "", barcode_choices);
      
      # Met à jour la checkbox.
      updateCheckboxGroupInput(session,
                               "check_barcode",
                               "Genotypes:",
                               choiceNames = barcode_choices,
                               choiceValues = barcode_file_df[[1]],
                               selected = barcode_file_df[[1]])
    }
  )
  
  
  # Demande la confirmation de l'indexation du dossier en TABIX.
  observeEvent(input$tabix, {
    showModal(modalDialog(
      tagList(
        p("Do you want to index files in the current folder ?")
      ), 
      title="Tabix index",
      footer = tagList(actionButton("confirm_tabix", "Confirm"),
                       modalButton("Cancel")
      )
    ))
  })
  
  
  # Indexe le dossier en TABIX.
  observeEvent(input$confirm_tabix, { 
    
    # Close the confirmation window.
    removeModal()
    
    # Create 0-row data frame which will be used to store data
    dat <- data.frame(x = numeric(0), y = numeric(0))
    
    withProgress(message = 'Creating tabix index', value = 0, {
      # Number of times we'll go through the loop
      n <- 1
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Doing part", i))
        
        # Commande bash executant le script d'index tabix sur le dossier choisi.
        system(paste0("cd ", global$datapath, " && bash ~/scripts/tabix_script.sh"))
      }
    })
  })
  
  
  reactive_icon <- reactive({
    req(input$dir)
    icon <- paste0("/biotools/htslib/1.9/bin/tabix ", global$datapath,
                   "/", input$check_barcode, " ", input$AGI)
  }) #reactiveIcon
  
  
  read_file <- reactive({
    # Liste contenant les listes des freads des barcodes.
    read_file <- lapply(reactive_icon(), fread) ; 
    
    # Liste contenant les noms complets des fichiers barcodes du dossier 
    # ainsi que l'AGI.
    list_short <- lapply(reactive_icon(), basename);
    
    # Liste contenant les numéros des barcodes (ex: barcode09), et les associe 
    # a read_file.
    liste_file_split_short <- gsub("\\..*", "", list_short) ;
    names(read_file) <- liste_file_split_short;
    
    read_file <- read_file[lengths(read_file) != 0];
  })
  

  # Renvoie une liste de dataframe des barcodes utilisés dans la recherche.
  reactive_read_file <- reactive ({
    
    req(length(read_file())>0)
    
    # Liste contenant les listes des freads des barcodes.
    #read_file <- lapply(reactive_icon(), fread) ; 
    
    # Liste contenant les noms complets des fichiers barcodes du dossier 
    # ainsi que l'AGI.
    #list_short <- lapply(reactive_icon(), basename);
    
    # Liste contenant les numéros des barcodes (ex: barcode09), et les associe 
    # a read_file.
    #liste_file_split_short <- gsub("\\..*", "", list_short) ;
    #names(read_file) <- liste_file_split_short;
  
    #read_file <- read_file[lengths(read_file) != 0];
    
    #print(length(read_file));
    #print(lengths(read_file));
    
    read_file <- read_file()
    
    csv_barcode <- gsub("_sorted.csv.gz", "", input$check_barcode[1]);
    colnames <- fread(paste0('head -n 1 ', global$datapath, "/",csv_barcode));
    

    # Ajoute le header aux dataframes de la liste.
    for(i in 1:length(read_file())) {
      colnames(read_file[[i]]) <- names(colnames)
    };
      return(read_file)

  }) #reactiveReadFile
  
  

  
  
  # Renvoie le dataframe affiché en mergeant tous les dataframes recherchés, et
  # ajoute la colonne barcode.
  reactive_df <- eventReactive(
    input$update_df,
    {
      
      # Create 0-row data frame which will be used to store data
      dat <- data.frame(x = numeric(0), y = numeric(0))
      
      withProgress(message = 'Making DataTable', value = 0, {
        # Number of times we'll go through the loop
        n <- 10
        
        for (i in 1:n) {
          # Each time through the loop, add another row of data. This is
          # a stand-in for a long-running computation.
          dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
          
          # Increment the progress bar, and update the detail text.
          incProgress(1/n, detail = paste("Doing part", i))
          
          # Rajoute la colonne Barcode aux dataframes et garde le nom des lignes.
          df <- rbindlist(reactive_read_file(), idcol="Barcode", unname(FALSE));
        }
      })
      
      return(df)
    }
  ) #reactiveDf
  
  
  # Affiche la datatable.
  output$df_file <- renderDataTable({
    reactive_df()
  },
  options = list(scrollX = TRUE)
  ) #output$df-file
  
  
  # Fait débuter le browser de fichiers de sauvegarde de df dans le home.
  shinyDirChoose(
    input,
    'save_dir',
    roots = c(home = '~')
  ) #shinyDirChoose
  
  
  save_global <- reactiveValues(datapath = getwd())
  
  
  save_dir <- reactive(input$save_dir)
  
  
  output$save_dir <- renderText({ 
    save_global$datapath 
  }) #output$dir
  
  
  # Met a jour le path du dossier choisi.
  observeEvent(ignoreNULL = TRUE, eventExpr = {input$save_dir},
               handlerExpr = {
                 if (!"path" %in% names(save_dir())) return()
                 home <- normalizePath("~")
                 save_global$datapath <- file.path(home,
                                              paste(unlist(save_dir()$path[-1]),
                                                    collapse = .Platform$file.sep))
               }) #observeEvent
  
  
  
  # Fenêtre de confirmation de sauvegarde du dataframe.
  observeEvent(input$save_df, {
    showModal(modalDialog(
      tagList(
        textInput("filename", label = "Filename", placeholder = "my_file.csv"),
        
        shinyDirButton("save_dir", 'Select where to save the dataframe', 'Please select a folder', FALSE),
        verbatimTextOutput("save_dir", placeholder = TRUE)
      ), 
      title="Save the dataframe as .csv",
      footer = tagList(actionButton("confirmCreate", "Create"),
                       modalButton("Cancel")
      )
    ))
  })
  
  
  # Enregistre le dataframe en .csv.
  observeEvent(input$confirmCreate, { 
    
    # Requiert un nom de fichier pour créer la sauvegarde.
    req(input$filename)
    
    # Close the confirmation window.
    removeModal()
    
    # Create 0-row data frame which will be used to store data
    dat <- data.frame(x = numeric(0), y = numeric(0))
    
    withProgress(message = 'Saving file', value = 0, {
      # Number of times we'll go through the loop
      n <- 1
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Doing part", i))
        
        # Sauvegarde le dataframe en .csv.
        write.csv(reactive_df(), paste0(save_global$datapath, "/", input$filename), row.names=FALSE, quote=FALSE)
      }
    })
  })

  

  # Extrait les coords et renvoie le dataframe utilisé pour créer les plots.
  reactive_selected_df <- reactive({
      
      # Sélectionne les données qu'on va utiliser dans le plot.
      df_plot <- select(reactive_df(), 
                        Barcode,
                        mRNA,      
                        read_core_id,      
                        coords_in_read,    
                        sense,
                        add_tail_pct_T,
                        add_tail_pct_A,
                        add_tail_pct_G,
                        add_tail_pct_C);     
      
      # Ajoute un index.
      df_plot <- df_plot %>% mutate(read=seq(1:nrow(df_plot)));
      
      # Sépare la colonne read_core_id.
      df_plot <- separate(df_plot,
                          col=read_core_id,
                          into=c("readname", "Chr",
                                 "genome_start_coord", "genome_end_coord"), 
                          sep=',')
      
      # Sépare la colonne coords_in_read.
      df_plot <- separate(df_plot, 
                          col=coords_in_read,
                          into=c("gen_start", "gen_end", 
                                 "polya_start", "polya_end", 
                                 "add_start", "add_end", 
                                 "adapt_start", "adapt_end"), 
                          sep=':')
      
      # Transforme les colonnes en numeric.
      df_plot <- df_plot %>%
        transform(genome_start_coord=as.numeric(genome_start_coord)) %>%
        transform(genome_end_coord=as.numeric(genome_end_coord)) %>%
        transform(gen_start=as.numeric(gen_start)) %>%
        transform(gen_end=as.numeric(gen_end)) %>%
        transform(polya_start=as.numeric(polya_start)) %>%
        transform(polya_end=as.numeric(polya_end)) %>%
        transform(add_start=as.numeric(add_start)) %>%
        transform(add_end=as.numeric(add_end)) %>%
        transform(adapt_start=as.numeric(adapt_start)) %>%
        transform(adapt_end=as.numeric(adapt_end));
      
      df_plot$polyA_length <- df_plot$polya_end-df_plot$polya_start;
      
      df_plot$polyU_length <- df_plot$add_end-df_plot$add_start;
      
      # Trie le df en fonction de la colonne gen_end et du sens.
      df_plot <- df_plot[order(df_plot$gen_end),];
      
      # Index.
      df_plot <- df_plot %>% mutate(read=seq(1:nrow(df_plot)))
      
      
    }
  )
  
  # Crée le dataframe requis pour créer le BoxPlot du pourcentage des bases
  # nucléiques.
  reactive_df_add_tail <- reactive({
      
      df_add_tail <- select(reactive_df(),
                            mRNA,      
                            readname,      
                            add_tail_pct_A,
                            add_tail_pct_T,
                            add_tail_pct_G,
                            add_tail_pct_C,
                            sense);
      
      # Supprime les NA.
      df_add_tail <- df_add_tail[complete.cases(df_add_tail)] ##### MARCHE PAS ? A VERIFIER
      
      # Index -> count affiché sur le boxplot.
      df_add_tail <- df_add_tail %>% mutate(read=seq(1:nrow(df_add_tail)))
    }
  )
  
  
  # Change le format du df pour le plot 1.
  reactive_df_plot <- reactive({

      # Fait passer le dataframe en long.
      df_plot <- reactive_selected_df() %>%
        pivot_longer(ends_with(c("start", "end")),
                     names_to = "type_of_sequence",
                     values_to = "coord");
      
      # Sépare la colonne "type_of_sequence" en deux colonnes "type" et
      # "start_end".
      df_plot <- separate(df_plot,
                          col=type_of_sequence,
                          into=c("type", "start_end"),
                          sep="_");
      
      # Fait passer le dataframe en wide selon la colonne start_end et les 
      # valeurs de coordonnées.
      df_plot <- df_plot %>%
        pivot_wider(names_from=start_end, values_from=coord);
    }
  )
  
  
  # Update le slider Poly(A) distribution.
  observeEvent(input$update_df, {
    
    # Ordonne le dataframe selon les tailles de polyA décroissantes.
    df_plotted <- reactive_selected_df() %>%
      arrange(-polyA_length)
    
    
    # Stocke la plus grande valeur de polyA dans max_value
    max_value <- df_plotted$polyA_length[1]
    
    updateSliderInput(session, "slider", 
                      min=0, max = max_value,
                      value=c(0,max_value))
  })
  
  
  # Update le slider Poly(U) distribution.
  observeEvent(input$update_df, {
    
    # Ordonne le dataframe selon les tailles de polyU décroissantes.
    df_plotted <- reactive_selected_df() %>%
      arrange(-polyU_length)
    
    
    # Stocke la plus grande valeur de polyU dans max_value
    max_value <- df_plotted$polyU_length[1]
    
    updateSliderInput(session, "slider_polyU", 
                      min=0, max = max_value,
                      value=c(0,max_value))
  })
  
  
  # Update le slider PolyA/polyU
  observeEvent(input$update_df, {
    df_cutoff <- reactive_selected_df()
    
    # Ordonne le dataframe selon les tailles de polyA décroissantes.
    df_plotted <- df_cutoff %>%
      arrange(-polyU_length)
    
    # Stocke la plus grande valeur de polyA dans max_value
    max_value <- df_plotted$polyU_length[1]
    
    updateSliderInput(session, "AU_slider", 
                      min=0, max = max_value,
                      value=c(0,max_value))
    
  })
  

  
  
  # PLOT POLYA
  output$plot_6 <- renderPlot({
    req(input$update_df)
    
    ggplot(reactive_selected_df(), aes(polyA_length))+
      geom_bar(aes(fill=Barcode))+
      facet_wrap(~ Barcode) +
      scale_x_continuous(limits=input$slider) +
      xlab("Poly(A) tail length of genes") +
      ylab("Number of genes") +
      labs(title=paste("Number of poly(A) tail length per genotype of AGI", input$AGI)) +
      guides(fill="none")
  })
  
  output$plot_polyA <- renderPlot({
    req(input$update_df)
    
    ggplot(reactive_selected_df(), aes(polyA_length))+
      geom_density(aes(fill=Barcode))+
      facet_wrap(~ Barcode) +
      scale_x_continuous(limits=input$slider) +
      xlab("Poly(A) tail length of genes") +
      ylab("Ratio of genes") +
      labs(title=paste("Distribution of poly(A) tail length per genotype of AGI", input$AGI)) +
      guides(fill="none")
  })
  
  # PLOTS POLYU
  output$plot_10 <- renderPlot({
    req(input$update_df)
    
    df_cutoff <- reactive_selected_df()
    # Supprime les NA.
    df_cutoff[is.na(df_cutoff)] <- 0
    
    # Sélectionne les valeurs >= au cut off.
    df_cutoff <- df_cutoff[df_cutoff$add_tail_pct_T >= input$cutoff_U_distrib, ]
    
    ggplot(df_cutoff, aes(polyU_length))+
      geom_bar(aes(fill=Barcode))+
      facet_wrap(~ Barcode) +
      scale_x_continuous(limits=input$slider_polyU) +
      xlab("Poly(U) tail length of genes") +
      ylab("Number of genes") +
      labs(title=paste("Number of poly(U) tail length per genotype of AGI", input$AGI, ", tail composition >=", input$cutoff_U_distrib, "% uracil")) +
      guides(fill="none")
  })
  
  output$plot_11 <- renderPlot({
    req(input$update_df)
    
    df_cutoff <- reactive_selected_df()
    # Supprime les NA.
    df_cutoff[is.na(df_cutoff)] <- 0
    
    # Sélectionne les valeurs >= au cut off.
    df_cutoff <- df_cutoff[df_cutoff$add_tail_pct_T >= input$cutoff_U_distrib, ]
    
    ggplot(df_cutoff, aes(polyU_length))+
      geom_density(aes(fill=Barcode))+
      facet_wrap(~ Barcode) +
      scale_x_continuous(limits=input$slider_polyU) +
      xlab("Poly(U) tail length of genes") +
      ylab("Ratio of genes") +
      labs(title=paste("Distribution of poly(U) tail length per genotype of AGI", input$AGI, ", tail composition >=", input$cutoff_U_distrib, "% uracil")) +
      guides(fill="none")
  })
  
  # Affiche le plot position extremité 3'.
  output$plot_1 <- renderPlot({
    req(input$update_df)
    
    ggplot(reactive_df_plot(),
           aes(xmin= start+genome_start_coord, xmax= end+genome_end_coord, 
               ymin= read, ymax= read+0.8,
               text=Barcode)) +
      geom_rect(aes(fill=type)) +
      labs(fill="Type") +
      scale_fill_discrete(labels=c("Adapter",
                                   "Additional tail", 
                                   "Genome read", 
                                   "PolyA tail"))
  }) #output$plot
  
  
  output$plot_gen_start <- renderPlot({
    req(input$update_df)
    
    ggplot(reactive_selected_df(), aes(genome_start_coord, y=..count../sum(..count..), fill=Barcode)) +
      geom_bar() +
      facet_wrap(~ Barcode) +
      labs(title=paste("Percentage of reads based on 5' ends of AGI", input$AGI, ", n=", nrow(reactive_selected_df()))) +
      ylab("Percentage of reads") +
      xlab("5' ends coordinates") +
      guides(fill="none")
  })
  
  
  output$plot_1A <- renderPlot({
    req(input$update_df)
    
    ggplot(reactive_selected_df(), aes(genome_start_coord + gen_end, y=..count../sum(..count..), fill=Barcode)) +
      geom_histogram() +
      facet_wrap(~ Barcode) +
      labs(title=paste("Percentage of reads based on 3' ends of AGI", input$AGI, ", n=", nrow(reactive_selected_df()))) +
      ylab("Percentage of reads") +
      xlab("3' ends coordinates") +
      guides(fill="none")
  })
  
  
  
  # Affiche le BoxPlot des bases nucléiques.
  output$plot_3 <- renderPlot({
    req(input$update_df)

    ggplot(reactive_df(), aes("A", y=add_tail_pct_A, color=Barcode)) +
      geom_boxplot() +
      geom_boxplot(data=reactive_df(), aes("T", y=add_tail_pct_T)) +
      geom_boxplot(data=reactive_df(), aes("G", y=add_tail_pct_G)) +
      geom_boxplot(data=reactive_df(), aes("C", y=add_tail_pct_C)) +
      ylab("Percentage")+
      xlab("Base")+
      labs(title=paste("Percentage of C, G, A and U in additional tails of AGI", input$AGI, ", n=", nrow(reactive_df())))
  })
  
  # Affiche le BoxPlot des bases nucléiques.
  output$boxplot_cutoff <- renderPlot({
    req(input$update_df)
    
    df_cutoff <- reactive_df()
    df_cutoff <- df_cutoff[df_cutoff$add_tail_pct_T >= input$cutoff]
    
    ggplot(df_cutoff, aes("A", y=add_tail_pct_A, color=Barcode)) +
      geom_boxplot() +
      geom_boxplot(data=df_cutoff, aes("T", y=add_tail_pct_T))+
      geom_boxplot(data=df_cutoff, aes("G", y=add_tail_pct_G))+
      geom_boxplot(data=df_cutoff, aes("C", y=add_tail_pct_C))+
      ylab("Percentage")+
      xlab("Base")+
      labs(title=paste("Percentage of C, G, A and U in additional tails of AGI", input$AGI, ", n=", nrow(df_cutoff)))
  })
  
  output$plot_A <- renderPlot({
    req(input$update_df)
    
    df_cutoff <- reactive_selected_df()
    # Supprime les NA.
    df_cutoff[is.na(df_cutoff)] <- 0
    
    # Sélectionne les valeurs >= au cut off.
    df_cutoff <- df_cutoff[df_cutoff$add_tail_pct_A >= input$cutoff_A, ]
    
    df_cutoff <- df_cutoff %>%
      group_by(polyU_length, polyA_length) %>%
      summarise(n=n())
    
    
    ggplot(df_cutoff, aes(x=polyU_length, y=polyA_length))+
      geom_tile(aes(fill= n)) +
      scale_fill_distiller(palette = "Spectral", direction = -1) +
      labs(x="Poly(U) length",
           y="Poly(A) length",
           title= paste("Poly(U) length in relation to poly(A) length of AGI", input$AGI, ", tail composition >=", input$cutoff_U, "% adenine")) +
      scale_x_continuous(limits=input$AU_slider)
  })
  
  output$plot_U <- renderPlot({
    req(input$update_df)
    
    df_cutoff <- reactive_selected_df()
    # Supprime les NA.
    df_cutoff[is.na(df_cutoff)] <- 0
    
    # Sélectionne les valeurs >= au cut off.
    df_cutoff <- df_cutoff[df_cutoff$add_tail_pct_T >= input$cutoff_U, ]
    
    df_cutoff <- df_cutoff %>%
      group_by(polyU_length, polyA_length) %>%
      summarise(n=n())
    
    
    ggplot(df_cutoff, aes(x=polyU_length, y=polyA_length))+
      geom_tile(aes(fill= n)) +
      scale_fill_distiller(palette = "Spectral", direction = -1) +
      labs(x="Poly(U) length",
           y="Poly(A) length",
           title= paste("Poly(U) length in relation to poly(A) length of AGI", input$AGI, ", tail composition >=", input$cutoff_U, "% uracil")) +
      scale_x_continuous(limits=input$AU_slider)
  })
  
  output$plot_C <- renderPlot({
    req(input$update_df)
    
    df_cutoff <- reactive_selected_df()
    # Supprime les NA.
    df_cutoff[is.na(df_cutoff)] <- 0
    
    # Sélectionne les valeurs >= au cut off.
    df_cutoff <- df_cutoff[df_cutoff$add_tail_pct_C >= input$cutoff_C, ]
    
    df_cutoff <- df_cutoff %>%
      group_by(polyU_length, polyA_length) %>%
      summarise(n=n())
    
    
    ggplot(df_cutoff, aes(x=polyU_length, y=polyA_length))+
      geom_tile(aes(fill= n)) +
      scale_fill_distiller(palette = "Spectral", direction = -1) +
      labs(x="Poly(U) length",
           y="Poly(A) length",
           title= paste("Poly(U) length in relation to poly(A) length of AGI", input$AGI, ", tail composition >=", input$cutoff_U, "% cytosine")) +
      scale_x_continuous(limits=input$AU_slider)
  })
  
  output$plot_G <- renderPlot({
    req(input$update_df)
    
    df_cutoff <- reactive_selected_df()
    # Supprime les NA.
    df_cutoff[is.na(df_cutoff)] <- 0
    
    # Sélectionne les valeurs >= au cut off.
    df_cutoff <- df_cutoff[df_cutoff$add_tail_pct_G >= input$cutoff_G, ]
    
    df_cutoff <- df_cutoff %>%
      group_by(polyU_length, polyA_length) %>%
      summarise(n=n())
    
    
    ggplot(df_cutoff, aes(x=polyU_length, y=polyA_length))+
      geom_tile(aes(fill= n)) +
      scale_fill_distiller(palette = "Spectral", direction = -1) +
      labs(x="Poly(U) length",
           y="Poly(A) length",
           title= paste("Poly(U) length in relation to poly(A) length of AGI", input$AGI, ", tail composition >=", input$cutoff_U, "% guanine")) +
      scale_x_continuous(limits=input$AU_slider)
  })
  
  
  output$plot_7 <- renderPlot({
    req(input$update_df)
    
    ggplot(reactive_selected_df(), aes(x=polya_end-polya_start, y=add_end-add_start, fill=Barcode))+
      geom_tile(alpha=0.5) +
      labs(x="Poly(A) length",
           y="Poly(U) length",
           title= paste("Poly(U) length in relation to poly(A) length of AGI", input$AGI)) +
      facet_wrap(~ Barcode)
  })
  
  output$plot_8 <- renderPlot({
    req(input$update_df)

    ggplot(reactive_df(), aes(add_tail_T, fill=Barcode)) +
      geom_bar() +
      labs(title = paste("Distribution of the number of Uridines in additional tails of AGI", input$AGI)) +
      xlab("Number of uridines") +
      facet_wrap(~ Barcode) +
      guides(fill="none")
  })

} #server
  


shinyApp(ui, server)