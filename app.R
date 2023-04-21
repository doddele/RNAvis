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

      shinyDirButton("dir", 'Select a folder', 'Please select a folder',
                     icon=icon(name="", class="fa-regular fa-folder-open")),
      
      verbatimTextOutput("dir", placeholder = TRUE),
      
      br(),
      
      actionButton("tabix",
                   "Create Tabix index", class = "btn-warning"),
    
      textInput("AGI",
              label = h3("AGI input")
      ), #textInput
      
      helpText("AT1G01010.1"),
      helpText("AT1G75830.1"),
      helpText("ATMG01350.1"),
    
      checkboxGroupInput("check_barcode",
                       "Barcodes:",
                       choices = NULL
      ), #checkboxGroupInput
      
     br(),
     
     actionButton("update_df",
                  "Update",
                  icon=icon(name="", class="fa-solid fa-arrows-rotate")),
     actionButton("save_df",
                  "Download DataTable",
                  icon = icon(name="", class="fa-solid fa-download"))
    ) #sidebarMenu
  ), #dashboardSidebar
  
  dashboardBody(
    fluidRow(
      tabsetPanel(tabPanel("Global statistics",
                           tabBox(
                             title = "First tabBox",
                             # The id lets us use input$tabset1 on the server to find the current tab
                             id = "tabset1", height = "250px"
                           )),
                  tabPanel("mRNA statistics",
      tabsetPanel(
        title = NULL,
        # The id lets us use input$tabset1 on the server to find the current tab
        id = "tabset1",
        
        tabPanel("DataTable",
                 dataTableOutput("df_file")
                 ),
        
        tabPanel("Plot",
                 h4("Brush and double click to zoom"),
                 plotlyOutput(outputId="plotly")
                 ),
        
        tabPanel("Histogram",
                 plotOutput("histogram")),
        
        tabPanel("BoxPlot",
                 h4("Percentage of each nucleic acid in additional tails"),
                 plotOutput("plot_3")
                 ),
        
        tabPanel("PolyA distribution",
                 h4("PolyA length distribution"),
                 plotOutput("plot_4")
                 ),
        
        tabPanel("PolyA distribution 2",
                 h4("Mean PolyA length per Barcode"),
                 plotOutput("plot_5")
                 ),
        
        tabPanel("PolyA length",
                 h4("Another mean PolyA length per Barcode"),
                 plotOutput("plot_6"),
                 sliderInput("slider", "x-axis", min=0, max=100,
                             value=c(0,100))),
        tabPanel("polyA/polyU",
                 plotOutput("plot_7"))
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
                               "Barcodes:",
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
    

    # Ajoute le header aux dataframes de la liste. <- TODO: CRASH QUAND UPDATE
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
                        readname,      
                        coords_in_read,    
                        sense);     
      
      # Ajoute un index.
      df_plot <- df_plot %>% mutate(read=seq(1:nrow(df_plot)));
      
      # Sépare la colonne V24 (coords_in_read).
      df_plot <- separate(df_plot, 
                          col=coords_in_read,
                          into=c("gen_start", "gen_end", 
                                 "polya_start", "polya_end", 
                                 "add_start", "add_end", 
                                 "adapt_start", "adapt_end"), 
                          sep=':')
      
      # Transforme les colonnes en numeric.
      df_plot <- df_plot %>%
        transform(gen_start=as.numeric(gen_start)) %>%
        transform(gen_end=as.numeric(gen_end)) %>%
        transform(polya_start=as.numeric(polya_start)) %>%
        transform(polya_end=as.numeric(polya_end)) %>%
        transform(add_start=as.numeric(add_start)) %>%
        transform(add_end=as.numeric(add_end)) %>%
        transform(adapt_start=as.numeric(adapt_start)) %>%
        transform(adapt_end=as.numeric(adapt_end));
      
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
      df_add_tail <- df_add_tail[complete.cases(df_add_tail)]
      
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

  
  # Affiche le plot 1.
  output$plotly <- renderPlotly({
    req(input$update_df)
    
    # Create 0-row data frame which will be used to store data
    dat <- data.frame(x = numeric(0), y = numeric(0));
    
    withProgress(message = 'Making Plot', value = 0, {
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Doing part", i))
        
        # Making the plot.
        p <- ggplot(reactive_df_plot(),
               aes(xmin= as.integer(start), xmax= as.integer(end), 
                   ymin= read, ymax= read+0.8,
                   text=Barcode)) +
          geom_rect(aes(fill=type)) +
          labs(fill="Type") +
          scale_fill_discrete(labels=c("Adapter",
                                       "Additional tail", 
                                       "Genome read", 
                                       "PolyA tail"));
        fig <- ggplotly(p, tooltip = "text")
        
      }
    })
    
    # Print the plot.
    fig
  }) #output$plot
  
  
 # Affiche l'histogramme.
  output$histogram <- renderPlot({
    req(input$update_df)
    
    # Create 0-row data frame which will be used to store data
    dat <- data.frame(x = numeric(0), y = numeric(0))
    
    withProgress(message = 'Making Histogram', value = 0, {
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Doing part", i))
        
        # Making the histogram.
        p <- ggplot(reactive_selected_df(), aes(x=gen_end, fill=Barcode)) +
          geom_histogram() +
          xlab("PolyA start coordinates")
      }
    })
    
    # Print the plot.
    p
  }) #output$plot2
  
  
  # Affiche le BoxPlot des bases nucléiques.
  output$plot_3 <- renderPlot({
    req(input$update_df)
    
    # Create 0-row data frame which will be used to store data
    dat <- data.frame(x = numeric(0), y = numeric(0))
    
    withProgress(message = 'Making BoxPlot', value = 0, {
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Doing part", i))
        
        # Making the plot.
        p <- ggplot(reactive_df(), aes("A", y=add_tail_pct_A, color=Barcode)) +
          geom_boxplot() +
          geom_boxplot(data=reactive_df(), aes("T", y=add_tail_pct_T))+
          geom_boxplot(data=reactive_df(), aes("G", y=add_tail_pct_G))+
          geom_boxplot(data=reactive_df(), aes("C", y=add_tail_pct_C))+
          ylab("Percentage")+
          xlab("Base")+
          labs(title=paste("n=", nrow(reactive_df())))
      }
    })
    p
  })
  
  
  # Affiche le plot PolyA length distribution.
  output$plot_4 <- renderPlot({
    req(input$update_df)
    
    # Create 0-row data frame which will be used to store data
    dat <- data.frame(x = numeric(0), y = numeric(0))
    
    withProgress(message = 'Making BoxPlot', value = 0, {
      # Number of times we'll go through the loop
      n <- 10
      
      for (i in 1:n) {
        # Each time through the loop, add another row of data. This is
        # a stand-in for a long-running computation.
        dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
        
        # Increment the progress bar, and update the detail text.
        incProgress(1/n, detail = paste("Doing part", i))
        
        # Making the plot.
        p <- ggplot(reactive_selected_df(), aes(x=polya_start, fill=Barcode, color=Barcode)) +
          geom_density(alpha=0.5) +
          xlab("PolyA start coordinates")
      }
    })
    p
  })
  
  
  output$plot_5 <- renderPlot({
    req(input$update_df)
    
    df_plot <- aggregate(reactive_df()$polya_length, by=list(Barcode=reactive_df()$Barcode), FUN=mean);
    ggplot(df_plot, aes(x=Barcode, y=x, color=Barcode, fill=Barcode))+
      geom_col()+
      theme(legend.position="none")
    })
  
  
  observeEvent(input$update_df, {
    
    # Ordonne le dataframe selon les tailles de polyA décroissantes.
    df_plotted <- reactive_df() %>%
      arrange(-polya_length)
    
    # Stocke la plus grande valeur de polyA dans max_value
    max_value <- df_plotted$polya_length[1]
    
    updateSliderInput(session, "slider", 
                      min=0, max = max_value,
                      value=c(0,max_value))
  })
  
  
  output$plot_6 <- renderPlot({
    req(input$update_df)
    
    ggplot(reactive_df(), aes(polya_length, color=Barcode))+
      geom_density()+
      coord_cartesian(xlim=input$slider)
  })
  
  
  output$plot_7 <- renderPlot({
    req(input$update_df)
    
    ggplot(reactive_selected_df(), aes(x=polya_end-polya_start, y=add_end-add_start, color=Barcode))+
      geom_point()+
      geom_line()
  })


} #server
  


shinyApp(ui, server)