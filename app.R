library(shiny)
library(reactable)
library(reactablefmtr)
library(readxl)
library(dplyr)
library(htmltools)
library(tidyverse)
library(rsconnect)



library(readxl)
T1II_summstats_jan21 <- read_excel("data/T1II_summstats_jan21.xlsx", 
                                   col_types = c("text", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "text"))

T4I_summstats_may21 <- read_excel("data/T4I_summstats_may21.xlsx", 
                                   col_types = c("text", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "text"))

T4II_summstats_may21 <- read_excel("data/T4II_summstats_may21.xlsx", 
                                   col_types = c("text", "numeric", "numeric", 
                                                 "numeric", "numeric", "numeric", 
                                                 "numeric", "numeric", "text"))

feb_summstats_jan21 <- read_excel("data/feb_summstats_jan21.xlsx", 
                                  col_types = c("text", "numeric", "numeric", 
                                                "numeric", "numeric", "numeric", 
                                                "numeric", "numeric", "text"))

ALL <- read_excel("data/ALL.xlsx", 
                  col_types = c("text", "numeric", "numeric", 
                                "numeric", "numeric", "numeric", 
                                "numeric", "numeric", "text"))

samp_summ <- read_excel("data/sample_summ.xlsx", col_types = c("text", "text","text", "numeric","numeric")) 

colorscale <- function(x) rgb(colorRamp(c("#d5e5e7", "#86cbff", "#50aef5"))(x), maxColorValue = 255)

# Define UI ----
ui <- fluidPage(
  fluidRow(
    column(12, style = "background-color:#FFFFFF;",
           h2(strong("Metagenomic analysis of bacterial communities from an oil refinery wastewater treatment plant"), align = "center", style = "color:#2971FF"), 
           h4("Henry Say | Gloor Lab", align = "center"),
           p("Department of Biochemistry, Western University, London, Canada", align = "center", style = "color:#B0B0B0"),
           fluidRow(
             #Left column
             column(4,style = "background-color:#FFFFFF;",
                    h3(strong("Introduction"), style = "background-color:#2971FF;color:#FFFFFF"),
                    p("The process of refining crude oil produces large amounts of wastewater containing toxic compounds including heavy metals, various hydrocarbons and aliphatic carboxylic compounds such as naphthenic acids. 
                      Given the grand scale of the oil production industry, wastewater is being accumulated at increasing rates and must be treated in order to limit its environmental impact. Purification of this wastewater is 
                      done at various treatment facilities, and typically involves multiple steps that involve physical separation of solids, and even bioremediation by resident bacterial communities found in all stages of purification downstream of primary treatment."),
                    #WASTEWATER FACILITY DIAGRAM
                    img(src = "suncor_sarnia_schematicalt.png", width = "100%", align = "center"),
                    p("Bacterial communities have been observed to influence the success of the purification process and are correlated to times of facility failure, which occur unpredictably several times per year. However, they are almost completely uncharacterized. Characterizing the metagenome of bacterial communities throughout the wastewater treatment facility could provide insight into
                       the relationship between bacterial constituents and capacity to remove certain toxic compounds at various stages. Furthermore, time series sampling could reveal how changes in the facility's microbiome are related to the functionality of the purification plant at each stage. Finally, meta-transcriptome sequencing could reveal genes that are activated at various treatment stages.", style = "font-size:14px;"),
                    h3(strong("Methods"), style = "background-color:#2971FF;color:#FFFFFF"),
                    #p(""),
                    #ASSEMBLY LONG READS DIAGRAM
                    img(src = "methods.png", width = "90%", align = "center")
             ),
             #Center column
             column(4,style = "background-color:#FFFFFF;",
                    h3(strong("Results"), style = "background-color:#2971FF;color:#FFFFFF"),
                    #SUMMARY STATS TABLE
                    #INSERT RUN STATS 
                    #tags$div(tags$ul(
                    #  tags$li(tags$span(h4())),
                    #)),
                    
                    reactableOutput("table2"),
                    p("\n"),
                    p(strong(h4("16/27 complete metagenome assembled genomes (>90% completeness, <10% redundancy)"), align = "center")),
                    selectInput("contig",
                                label = "",
                                choices = list('Accumulibacter sp005524045 (Feb)',
                                               'Accumulibacter sp005524045 (T4)',
                                               'Bacteriovoracales (T4)',
                                               'CG2-30-37-29 (T4)',
                                               'DRLM01 (T1)',
                                               'DRLM01 (T4)',
                                               'DRLM01(Feb)',
                                               'Fen-1342 (T1)',
                                               'Fen-1342 (T4)',
                                               'Ga0077554 (T1)',
                                               'Ga0077554(Feb)',
                                               'Ga0077554(Feb)',
                                               'Micavibrionales_A (T1)',
                                               'Micavibrionales_A (T4)',
                                               'OLB17 (T1)',
                                               'OLB17(Feb)',
                                               'RAS1 (T1)',
                                               'RBC074 (T1)',
                                               'SXVE01 (T1)',
                                               'SZUA-149 (T1)',
                                               'UBA2112 (T4)',
                                               'UBA2386 (Feb)',
                                               'UBA7966 (Feb)',
                                               'UBA9160 (Feb)',
                                               'UBA9628(Feb)'
                                ),
                                selected = 'Accumulibacter sp005524045 (T4)'
                    ),
                    uiOutput("circos")
             ), 
             #Right column
             column(4,style = "background-color:#FFFFFF;",
                    h3("x ", style = "background-color:#2971FF;color:#2971FF"),

                    selectInput("samples",
                                label = "",
                                choices = list("Trailer 1 GAC (10-Aug-21)",
                                               "Trailer 4 GAC - sample 1 (14-May-21)",
                                               "Trailer 4 GAC - sample 2 (14-May-21)",
                                               #"Trailer 2  GAC (03-May-21)",
                                               #"Trailer 3 GAC (07-May-21)",
                                               "GAC (XX-Feb-20)",
                                               #"GAC (01-Jan-20)")
                                               "ALL"
                                )
                    ),
                    reactableOutput("table"),
                    
                    h3(strong("Discussion"), style = "background-color:#2971FF;color:#FFFFFF"),
                    tags$div(tags$ul(
                      tags$li(tags$span("Relatively low coverage for most assembled genomes, many are incomplete. Even for complete ones, coverage could be improved - per base accuracy is probably low")),
                      tags$li(tags$span("Nanopore with their new basecalling/assembly model can reach Q50 at 100x coverage.")),
                      tags$li(tags$span("Even with 20-30 gb of data collected, coverage is not ideal. Indicate that maybe even twice as much data is needed to be collected from one gac sample to achieve that (Q50)"))
                    )),
                    h4(strong("Next steps:")),
                    tags$div(tags$ul(
                      tags$li(tags$span("Sequence the metatranscriptome: use complete annotated genomes as a reference to detect upregulated/downregulated genes (i.e. hydrocarbon degradation genes")),
                      tags$li(tags$span("Observe differences in gene expression between stages of purification, and potentially relate gene expression to the wastewater constituents before and after each stage to infer functional characteristics of bacterial communities.")),
                      tags$li(tags$span(" Observe similarities/differences in bacterial constituents/abundance in samples collected at different times - infer how it may play a role in times of failure.")),
                      tags$li(tags$span(" i.e. one GAC sample collected from Feb 2020 had almost 500 fold coverage from a single run - indicates that this bacteria was extremely abundant vs other species at this time.")),
                      tags$li(tags$span(" Use 16s rRNA to potentially get a more specific taxonomic classification.")),
                      tags$li(tags$span(" Current taxonomic predications aren't very specifc - anticipate that most bacterial species identified will be novel.")),
                      tags$li(tags$span(" Potentially combine reads from runs for a pangenome analysis ")),
                      tags$li(tags$span(" Potentially filter and combine reads for species found in multiple samples to improve coverage/assembly coverage")),
                      tags$li(tags$span(" Identify subgenomic sequences.")),
                      style = "font-size:14px;"
                      )),
                    h3(strong("Conclusion"), style = "background-color:#2971FF;color:#FFFFFF"),
                    p("Nanopore technology allowed for assembly of complete genomes from a metagenome. Ultimately, this allows us to identify bacteria and their functional characteristics in this community"),
                    h3(strong("Acknowledgements"), style = "background-color:#2971FF;color:#FFFFFF"),
                    p("Thank you to Dr. Gloor, Dr. Edgell, Dr. Yeung, Daniel Giguere, and Ben Joris.")
             ))
)
))

# Define server logic ----
server <- function(input, output) {
  output$table <- renderReactable({
    #REACTABLE SAMPLE SELECT
    data <- switch(input$samples, 
                   "Trailer 1 GAC (10-Aug-21)" = T1II_summstats_jan21, 
                   "Trailer 4 GAC - sample 1 (14-May-21)" = T4I_summstats_may21,
                   "Trailer 4 GAC - sample 2 (14-May-21)" = T4II_summstats_may21,
                   "GAC (XX-Feb-20)" = feb_summstats_jan21,
                   "ALL" = ALL
    )
    
    reactable(data[,1:8],
              height = 350,
              compact = TRUE,
              outlined = TRUE,
              highlight = TRUE,
              pagination = TRUE,
              style = list(fontFamily = "Work Sans, sans-serif", fontSize = "14px"),
              defaultColDef = colDef(align = "center", format = colFormat(digits = 2)),
              columns = list(
                Size = colDef(name = "Size (bp)",
                              format = colFormat(digits = 0)),
                Classification = colDef(align = "left", minWidth = 120, cell = function(value, index) {
                  species <- data$Class2[index]
                  div(
                    div(style = list(fontWeight = 600), value),
                    div(style = list(fontSize = 10), species))}
                ),
                Contamination = colDef(minWidth = 100,
                                       name = "R(%)",
                                       format = colFormat(digits = 0),
                                       cell = data_bars(
                                         data, round_edges = TRUE,
                                         force_outside = c(0,3),
                                         max_value = 10,
                                         fill_color = "#f1c40f",
                                         text_position = "outside-base",
                                         number_fmt = scales::number
                                         ),
                                       
                                       
                ),
                Coverage = colDef(style = function(value) {
                  normal <- (value - min(data$Coverage)) / (max(data$Coverage) - min(data$Coverage))
                  color <- colorscale(normal)
                  list(background = color)
                }
                ),
                GC = colDef(name = "GC (%)"),
                Completeness = colDef(minWidth = 100,
                                      name = "C(%)",
                                      cell = data_bars(data,
                                                       round_edges = TRUE,
                                                       max_value = 100,
                                                       fill_color = c("#1abc9c"),
                                                       text_position = "outside-base",
                                                       number_fmt = scales::number
                                                       ),
                              
                                      ),
                rRNA_genes = colDef(name = "rRNA genes", minWidth = 80, format = colFormat(digits = 0)),
                tRNA_genes = colDef(name = "tRNA genes", minWidth = 80, format = colFormat(digits = 0)),
                species = colDef(show = FALSE)
              )
    )
  })
  output$table2 <- renderReactable({
    reactable(samp_summ,
              compact = TRUE,
              outlined = TRUE,
              highlight = TRUE,
              style = list(fontFamily = "Work Sans, sans-serif", fontSize = "14px"),
              defaultColDef = colDef(align = "center", format = colFormat(digits = 2))
              )
  })
  
  output$circos <- renderUI({
    if(input$contig == "Accumulibacter sp005524045 (T4)"){            
      img(width = "100%", src = 'Accumulibacter sp005524045 (T4).png', position = "center")}
    else if (input$contig == "Fen-1342 (T4)"){            
      img(width = "100%", src = 'Fen-1342 (T4).png', position = "center")}
    else if (input$contig == "DRLM01 (T4)"){            
      img(width = "100%", src = 'DRLM01 (T4).png', position = "center")}
    else if (input$contig == "Bacteriovoricales (T4)"){            
      img(width = "100%", src = 'Bacteriovoricales (T4).png', position = "center")}    
    else if (input$contig == "RBC074 (T1)"){            
      img(width = "100%", src = 'RBC074 (T1).png', position = "center")} 
    else if (input$contig == "SXVE01 (T1)"){            
      img(width = "100%", src = 'SXVE01 (T1).png', position = "center")} 
    else if (input$contig == "Fen-1342 (T1)"){            
      img(width = "100%", src = 'Fen-1342 (T1).png', position = "center")} 
    else if (input$contig == "RAS1 (T1)"){            
      img(width = "100%", src = 'RAS1 (T1).png', position = "center")} 
    else if (input$contig == "DRLM01 (T1)"){            
      img(width = "100%", src = 'DRLM01 (T1).png', position = "center")} 
    else if (input$contig == "Ga0077554 (T1)"){            
      img(width = "100%", src = 'Ga0077554 (T1).png', position = "center")} 
    else if (input$contig == "SZUA-149 (T1)"){            
      img(width = "100%", src = 'SZUA-149 (T1).png', position = "center")}
    else if (input$contig == "Micavibrionales_A (T1)"){            
      img(width = "100%", src = 'Micavibrionales_A (T1).png', position = "center")}
    else if (input$contig == "OLB17 (T1)"){            
      img(width = "100%", src = 'OLB17 (T1).png', position = "center")}
  })
}
  
# Run the app ----
shinyApp(ui = ui, server = server)