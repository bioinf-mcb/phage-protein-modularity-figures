#
# This is a Shiny web application. You can run the application by clicking
# httr requires curl 5.0.2 which is not available
# install older version of hhtr: install_version("httr", version = "1.4.2")



# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shinylogs)
library(shiny)
library(R.utils)
# adding it at the bottom may brak igraph plotting
library(ggiraph)  
library(stringr)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggnetwork)
library(plotly)
set.seed(1)
options(stringsAsFactors = F)


source("R/helpers.R")
ecod.levels = data.frame(ecod.group = c("H", "T"), stringsAsFactors = FALSE) %>%
  mutate(level = 1:2)

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
    withMathJax(),
    # Application title
    titlePanel("Phage structural part architecture lookup"),
    
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
          tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 0px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #CCFF66;
               z-index: 105;
             }
          ")),
            h4("Parameters"),
            br(),
            selectInput("ecod_group",
                     label = "Choose the ECOD aggregation group:",
                     ecod.levels %>%  pull(ecod.group),
                     selected = "T"),
            uiOutput("category_selection"),
            uiOutput("annotation_selection"),
            downloadButton("downloadData", "Download"),
            #conditionalPanel(condition = "isvalidCategory",
            #                 uiOutput("annotation_selection")),
            conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                             tags$div("Please wait, the server is loading data...",id="loadmessage")),
          width = 3
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(type = "tabs",
                      tabPanel("All Domains per annotation",
                               h4("Domain Architecture of all representative proteins from selected category that have at least one domain.
                                   Hover over the plot to see the details."),
                               girafeOutput("domainplot", height = "500px"),# making it bigger doesn't help, the plot is too crowded
                               plotOutput("domainplotlegend"),
                               girafeOutput("domainplotsimplified", height = "500px"),
                      ),
                      tabPanel("Domains shared with other annotations",
                               uiOutput("second_annotation_selection"),
                               h4("Sequences are filtered to those containing domains common to both annotation"),
                               girafeOutput("domainplotfromotherannotations"),
                               girafeOutput("domainplotsharedwithotherannotations"),
                      )
          ),
          
        )
    )
))


# Define server logic required to draw a histogram
server <- shinyServer(function(input, output, session) {
  #print(sessionInfo())
  #print(version)
  
  
  # LOAD ALL DATASET IF POSSIBLE
  # DON'T KEEP TOO MICH MEMORY IN THE DATA FOLDEE OR THE SERVER WILL BE SUPER SLOW
  
  track_usage(storage_mode = store_json(path = "logs/"))
  # Create a Progress object
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  #tile.data = data.table::fread('data/domain.positions.txt')
  annnotations.and.categories = readRDS("data/annnotations.and.categories.Rds")
  tile.data = readRDS('data/domain.positions.Rds')  
  progress$inc(10, detail = "Loading data |")
  #tile.data.multi = data.table::fread('data/multiple_domains.domain.positions.txt')
  tile.data.multi = readRDS('data/multiple_domains.domain.positions.Rds')
  progress$inc(20, detail = "Loading data II")
  text.size = 12
  theme.no.verical = theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
  

  
  output$category_selection <- renderUI({
    selection = annnotations.and.categories %>% filter(domain.level == input$ecod_group) %>% distinct(category) %>% pull(category)
    #print(selection)
    selectInput(inputId = "Category",
                label = "Select PHROG category:",
                choices = selection,
                selected = "tail") 
  })
  output$annotation_selection <- renderUI({
    if (isvalidCategory() & !(is.null(input$ecod_group))) {
    annot.selection = annnotations.and.categories %>% filter(category == input$Category & domain.level == input$ecod_group) %>% pull(annotation) %>% unique() %>% sort()
    #print(annot.selection)
    selected.annot =  "tail fiber"
      selectInput(inputId = "Annotation",
                  label = "Choose specific annotation :",
                  choices = annot.selection,
                  selected = selected.annot) 
    } else {
      #print("Invalid Category")
    }
  })
  
  output$second_annotation_selection <- renderUI({
    if(isValid_input()) {
    selectInput("second.Annotation",
                label = paste0("Choose an annotation from those sharing a domain with : ",input$Annotation),
                all.tiles.sharing.domains.with.this.annotation() %>% pull(annotation) %>% unique(),
                selected = all.tiles.sharing.domains.with.this.annotation() %>% head(1) %>% pull(annotation) %>% unique()) 
    }
  })
  

  progress$inc(30, detail = "Loading data this annot")
  this.category.tile.data  = reactive({tile.data %>% filter(category == input$Category)})
  this.category.tile.data.multi  = reactive({tile.data.multi %>% filter(category == input$Category)})
  
  #this.category.tile.data  = reactive({readRDS(sprintf('data/domain.positions_%s.Rds', input$Category))})
  #this.category.tile.data.multi  = reactive({readRDS(sprintf('data/multiple_domains.domain.positions_%s.Rds', input$Category))})
  this.tile.data = reactive({this.category.tile.data() %>% 
      filter(annotation == input$Annotation & category == input$Category & domain.level == input$ecod_group) })
  second.annotation.tile.data = reactive({tile.data %>% filter(annotation == input$second.Annotation & domain.level == input$ecod_group)})
  this.tile.data.multi = reactive({this.category.tile.data.multi() %>%
      filter(annotation == input$Annotation & category == input$Category & domain.level == input$ecod_group) })
  progress$inc(100, detail = "Done")
  #print(head(this.tile.data))
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      # Use the selected dataset as the suggested file name
      paste0("domain_positions_",input$Annotation, ".csv")
    },
    content = function(file) {
      # Write the dataset to the `file` that will be downloaded
      write.csv(this.tile.data() %>% select(-row, -row.plus, -family), file)
    }
  )
  
  
  isValid_input <- reactive({
    not.nulls = !is.null(input$Category) &  !is.null(input$Annotation) & !is.null(input$ecod_group)
    return(not.nulls)
    })
  
  isvalidCategory <- reactive({
    return(!is.null(input$Category))
  })
  
  
    output$errortext.any <- renderText({
        if(isValid_input()){ invisible(NULL) 
        }else{ #
            "Wrong input"
        }
    })
    output$errortextmulti <- renderText({
      if(isValid_input()){ invisible(NULL) 
      }else{ #
        "Wrong input"
      }
    })
    

    colm.any.domain = reactive({         
      colormap = this.tile.data() %>% distinct(domain) %>% filter(!(domain %in% c("undetected", "multiple domains")))
      colormap$color = randomcoloR::distinctColorPalette(nrow(colormap))
      colormap = rbind(colormap, data.frame(domain = c("undetected", "multiple domains"), 
                                            color = c( "#FFFFFF",  "#000000"), stringsAsFactors = FALSE))
      colm = colormap$color
      names(colm) = colormap$domain
      return(colm)
    })
    

    
    
    output$domainplot <- renderGirafe({ 
      if(isValid_input()){
        gg_rect = Show.Domain.Position.Within.Data(
          this.tile.data(),
          colm.any.domain(), 
          theme.no.verical) +
          theme(strip.text.x = element_text(size = 1.5*text.size), 
                legend.position = "none") +
          theme(
            axis.title.x = element_text(size = text.size),
            axis.title.y = element_text(size = text.size),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank())
        
        return(girafe(ggobj = gg_rect))
      }else{ #
        return(invisible(NULL))
      }
    })
    

    
    output$domainplotsimplified <- renderGirafe({ 
      if(isValid_input()){ 
        #print(this.tile.data.multi())
        gg_rect = Show.Domain.Position.Within.Data(
          this.tile.data.multi(),
          colm.any.domain(), 
          theme.no.verical) +
          theme(strip.text.x = element_text(size = 1.5*text.size), 
                legend.position = "none") +
          theme(
            axis.title.x = element_text(size = text.size),
            axis.title.y = element_text(size = text.size),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank())
        
        return(girafe(ggobj = gg_rect))
      }else{ #
        return(invisible(NULL))
      }
    })
    
    
    output$domainplotlegend <- renderPlot({
      if(isValid_input()){ 
        gg_rect = Show.Domain.Position.Within.Data(
          this.tile.data(),
          colm.any.domain(), 
          theme.no.verical) +
          theme(strip.text.x = element_text(size = 1.5*text.size)) +
          theme(
            axis.title.x = element_text(size = text.size),
            axis.title.y = element_text(size = text.size),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank())
        
        leg = ggpubr::get_legend(gg_rect)
        l=ggpubr::as_ggplot(leg)
        return(l)
      }else{ #
        return(invisible(NULL))
      }
    })
    


    all.tiles.sharing.domains.with.this.annotation = reactive({
      domains.in.this.annotation = this.tile.data() %>% filter(domain != "undetected" & domain != "multiple domains") %>% distinct(domain) 
      qnames.with.any.domain.from.this.annotation = tile.data %>% 
        filter(domain.level == input$ecod_group) %>% 
        inner_join(domains.in.this.annotation) %>%
        distinct(qname)
      all.tiles.sharing.domains.with.this.annotation = tile.data %>% 
        filter(domain.level == input$ecod_group) %>%
        inner_join(qnames.with.any.domain.from.this.annotation) %>%
        filter(annotation != input$Annotation)
      return(all.tiles.sharing.domains.with.this.annotation)
    })
    
    domains.shared.between.both.annotations = reactive({
      domains.in.this.annotation = this.tile.data() %>% 
        filter(domain != "undetected" & domain != "multiple domains") %>% 
        distinct(domain)       
      domains.in.the.second.annotation = second.annotation.tile.data() %>% 
        filter(domain != "undetected" & domain != "multiple domains") %>% 
        distinct(domain) 
      domains.shared.between.both.annotations = domains.in.this.annotation %>%
        inner_join(domains.in.the.second.annotation)
      return(domains.shared.between.both.annotations)
    })
    
    first.annotation.qname.data.where.domain.shared.with.second.annotation= reactive({
      qnames.with.any.domain.shared.between.annotations = this.tile.data() %>% 
        inner_join(domains.shared.between.both.annotations()) %>%
        distinct(qname)
      tiles =  this.tile.data() %>%
        inner_join(qnames.with.any.domain.shared.between.annotations)
      return(tiles)
    })
    
    second.annotation.qname.data.where.domain.shared.with.first.annotation= reactive({
      qnames.with.any.domain.shared.between.annotations = second.annotation.tile.data() %>% 
        inner_join(domains.shared.between.both.annotations()) %>%
        distinct(qname)
      tiles =  second.annotation.tile.data() %>%
        inner_join(qnames.with.any.domain.shared.between.annotations)
      return(tiles)
    })
    
    
    colormapbothannotations = reactive({
      colormap = rbind(
        second.annotation.qname.data.where.domain.shared.with.first.annotation() %>%
          distinct(domain) %>% 
          filter(!(domain %in% c("undetected", "multiple domains"))),
        first.annotation.qname.data.where.domain.shared.with.second.annotation() %>% 
          distinct(domain) %>% 
          filter(!(domain %in% c("undetected", "multiple domains")))) %>% 
        distinct()
      
      
      colormap$color = randomcoloR::distinctColorPalette(nrow(colormap))
      colormap = rbind(colormap, data.frame(domain = c("undetected", "multiple domains"), 
                                            color = c( "#FFFFFF",  "#000000"), stringsAsFactors = FALSE))
      colm = colormap$color
      names(colm) = colormap$domain
      return(colm)
      
    })
    
    
    
    output$domainplotfromotherannotations <- renderGirafe({ 
      if(isValid_input() & !(is.null(input$second.Annotation))){ 

        gg_rect = Show.Domain.Position.Within.Data(
          second.annotation.qname.data.where.domain.shared.with.first.annotation(),
          colormapbothannotations(), 
          theme.no.verical) +
          theme(strip.text.x = element_text(size = 1.5*text.size), 
                legend.position = "none") +
          theme(
            axis.title.x = element_text(size = text.size),
            axis.title.y = element_text(size = text.size),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank())
        
        return(girafe(ggobj = gg_rect))
      }else{ #
        return(invisible(NULL))
      }
    })



    output$domainplotsharedwithotherannotations <- renderGirafe({ 
      if(isValid_input() & !is.null(input$second.Annotation)){ 

        
        gg_rect = Show.Domain.Position.Within.Data(
          first.annotation.qname.data.where.domain.shared.with.second.annotation(),
          colormapbothannotations(), 
          theme.no.verical) +
          theme(strip.text.x = element_text(size = 1.5*text.size), 
                legend.position = "none") +
          theme(
            axis.title.x = element_text(size = text.size),
            axis.title.y = element_text(size = text.size),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank())
        
        return(girafe(ggobj = gg_rect))
      }else{ #
        return(invisible(NULL))
      }
    })
    
    
    session$onSessionEnded(function() {
      stopApp()
    })

})

# Run the application 
shinyApp(ui = ui, server = server)
