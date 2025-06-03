#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#load packages
library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(ComplexHeatmap)
library(Hmisc)
library(corrplot)
library(factoextra)
library(ggvenn)
library(RColorBrewer)
# library(palettetown)
library(forcats)
library(shinybrowser)
library(colorRamp2)
library(shinyWidgets)
library(shinyBS)
library(ggcorrplot)
library(reshape2)
library(viridis)


#read in metadata files
metadata <- read.csv("data/pankbase_human_donor_report_2025_5_27_17h_8m.csv")
metadata$Description.of.diabetes.status <- ifelse(metadata$Description.of.diabetes.status == "non-diabetic", "No diabetes", ifelse(metadata$Description.of.diabetes.status == "type 1 diabetes", "Type 1 diabetes", ifelse(metadata$Description.of.diabetes.status == "type 2 diabetes", "Type 2 diabetes", metadata$Description.of.diabetes.status)))
metadata$Ethnicities <- ifelse(metadata$Ethnicities == "Caucasian", "White", metadata$Ethnicities)

#filter out "test" sample
metadata <- metadata[metadata$Collections != "",]
#change names of diabetes status
metadata <- metadata %>% rename("Program" = Collections,
                                "Age (years)" = Age..years.,
                                "C. Peptide (ng/ml)" = C.Peptide..ng.ml.,
                                "AAB-GADA value (unit/ml)" = AAB.GADA.value..unit.ml.,
                                "AAB-IA2 value (unit/ml)" = AAB.IA2.value..unit.ml.,
                                "AAB-IAA value (unit/ml)" = AAB.IAA.value..unit.ml.,
                                "AAB-ZNT8 value (unit/ml)" = AAB.ZNT8.value..unit.ml.,
                                "HbA1C percentage" = HbA1C..percentage.,
                                "Hospital stay (hours)" = Hospital.Stay..hours.,
                                "Description of diabetes status" = Description.of.diabetes.status,
                                "Cause of death" = Cause.of.Death,
                                "Ethnicity" = Ethnicities)

metadata$Program <- ifelse(metadata$Program == "IIDP,Prodo", "IIDP", metadata$Program)

#rename values that have different capitalization so they get grouped together
metadata$`Cause of death` <- ifelse(metadata$`Cause of death` == "Cerebrovascular/stroke", "Cerebrovascular/Stroke", metadata$`Cause of death`)
metadata$`Cause of death` <- ifelse(metadata$`Cause of death` == "Head trauma", "Head Trauma", metadata$`Cause of death`)

#change capitalization of some fields
metadata$Sex <- ifelse(metadata$Sex == "female", "Female",
                       ifelse(metadata$Sex == "male", "Male",
                              metadata$Sex))

metadata$`Description of diabetes status` <- ifelse(metadata$`Description of diabetes status` == "cystic fibrosis diabetes", "Cystic fibrosis diabetes", metadata$`Description of diabetes status`)
metadata$`Description of diabetes status` <- ifelse(metadata$`Description of diabetes status` == "gestational diabetes", "Gestational diabetes", metadata$`Description of diabetes status`)
metadata$`Description of diabetes status` <- ifelse(metadata$`Description of diabetes status` == "monogenic diabetes", "Monogenic diabetes", metadata$`Description of diabetes status`)
metadata$`Description of diabetes status` <- ifelse(metadata$`Description of diabetes status` == "diabetes unspecified", "Diabetes unspecified", metadata$`Description of diabetes status`)
metadata$`Description of diabetes status` <- ifelse(metadata$`Description of diabetes status` == "steroid-induced diabetes", "Steroid-induced diabetes", metadata$`Description of diabetes status`)

#rename some NAs to "unknown"
metadata$Ethnicity <- ifelse(metadata$Ethnicity == "", "Unknown", metadata$Ethnicity)
metadata$`Cause of death` <- ifelse(metadata$`Cause of death` == "", "Unknown", metadata$`Cause of death`)
metadata <- metadata %>% mutate(across(starts_with("AAB."), ~ifelse( is.na(.x), "Unknown", .x)))
metadata <- metadata %>% rename("AAB-GADA Positive" = AAB.GADA.POSITIVE,
                                "AAB-IA2 Positive" = AAB.IA2.POSITIVE,
                                "AAB-IAA Positive" = AAB.IAA.POSITIVE,
                                "AAB-ZNT8 Positive" = AAB.ZNT8.POSITIVE)


#get all categorical variables. we actually don't want to plot all of these though
categorical_vars <- c("Program", "Description of diabetes status", "Cause of death", "Sex", "Ethnicity",
                      "AAB-GADA Positive", "AAB-IA2 Positive", "AAB-IAA Positive", "AAB-ZNT8 Positive")
all_continuous_vars <- metadata %>% select_if(is.numeric) %>% colnames()
continuous_vars <- c("Age (years)", "BMI", "C. Peptide (ng/ml)", "HbA1C percentage", "Hospital stay (hours)",
                     "AAB-GADA value (unit/ml)", "AAB-IA2 value (unit/ml)", 
                     "AAB-IAA value (unit/ml)",
                     "AAB-ZNT8 value (unit/ml)")

metadata_mini <- metadata[,c("Program", "Description of diabetes status", "Cause of death", "Sex", "Ethnicity",
                             "AAB-GADA Positive", "AAB-IA2 Positive", "AAB-IAA Positive", "AAB-ZNT8 Positive")]

## Wrangle donor-by-assay data here
data_avail <- metadata %>% select(c("ID", "Data.Available")) %>% filter(Data.Available != "")
data_avail$Data.Available2 <- gsub("\'","\"", data_avail$Data.Available, fixed=TRUE)
data_avail$Data.Available3 <- gsub("[[:space:]]", "", data_avail$Data.Available2, fixed=TRUE)

#get out df that has donor by assay by tissue info
all_donor_df <- data.frame()
for (d in unique(data_avail$ID)) {
  mini_df <- data_avail %>% filter(ID == d)
  
  donor_mini_df <- do.call(plyr::rbind.fill, lapply(paste0("[",mini_df$Data.Available3,"]"), function(x) jsonlite::fromJSON(x)))
  donor_mini_df$ID <- d
  
  all_donor_df <- rbind(all_donor_df, donor_mini_df)
  
}

all_donor_df$dataset <- ifelse(all_donor_df$dataset == "HLA_typing", "HLA typing", 
                           ifelse(all_donor_df$dataset == "RNAseq", "RNA-seq",
                                  ifelse(all_donor_df$dataset == "scRNAseq", "scRNA-seq",
                                         ifelse(all_donor_df$dataset == "snATACseq", "snATAC-seq", all_donor_df$dataset))))

islet_df <- all_donor_df %>% filter(dataset_tissue %in% c("Islet",  "-"))


#set up universal color palette
all_palette <- colorRampPalette(c("#FFBE0B", "#FB5607", "#FF006E", "#8338EC", "#3A86FF"))

#set default font size for ggplot
ggplot2::theme_set(theme_classic(base_size=18))




# Define UI for application
ui <- fluidPage(    
  shinybrowser::detect(),
  tags$head(
            tags$style(HTML("@import url('https://fonts.googleapis.com/css2?family=Open+Sans:ital,wght@0,300..800;1,300..800&display=swap');
                            .tabbable > .nav > li > a                  {color:black; border-top: 5px solid transparent; border-left: 0px solid transparent; border-right: 0px solid transparent; border-bottom: 5px solid transparent}
                            .tabbable > .nav > li[class=active]    > a {color:rgba(33, 145, 151, 1); border-top: 5px solid transparent; border-left: 0px solid transparent; border-right: 0px solid transparent; border-bottom: 5px solid rgba(33, 145, 151, 1)}
                            .tabbable > .nav > li > a:hover {color:rgba(33, 145, 151, 1); background-color:transparent}
                            .well                {background-color: rgba(33, 145, 151, 0.8); color:black}
                            input[type='checkbox'] {accent-color:rgba(148, 201, 94, 1); color:black}
                            .selectize-dropdown-content {background-color: rgba(148, 201, 94, 1); color:black}
                            .selectize-dropdown-content .active {background-color:rgba(148, 201, 94, 1) !important; color:white !important}
                            .selectize-dropdown-content .selected {background-color:rgba(148, 201, 94, 1); color:rgba(33, 145, 151, 0.8)}
                            * {font-family:'google', 'Open Sans', sans-serif; font-weight: 600}
                            a {color:rgba(33, 145, 151, 1)}
}
  "))),
  
  #make tabset on main panel 
  tabsetPanel(

    tabPanel("Donor Metadata Summary",
             tabsetPanel(
               tabPanel("Bar Plots",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("DataFilter", "Program Filter:",
                                        choices = c("All", unique(metadata$Program))),
                            bsTooltip(id = "DataFilter",
                                      title = "Select a program to subset the metadata"),
                            selectInput("Variable", "Grouping:", 
                                        choices=colnames(metadata_mini)),
                            bsTooltip(id = "Variable",
                                      title = "Select a grouping by which to separate donors"),
                            downloadButton("DownloadBar",
                                           "Download Barplot"),
                            selectInput("RemoveNA", "Remove 'unknowns' from stacked bar?",
                                        choices=c("No", "Yes")),
                            bsTooltip(id = "RemoveNA",
                                      title = "Choose whether to remove donors whose metadata for that variable is not currently in PanKbase"),
                            downloadButton("DownloadStacked", 
                                           "Download Stacked Barplot")
                            
                          ),
                          mainPanel(
                            p(),
                            p("Select inputs on left to show the number of donors that fall into different categories. By default, donors are shown from all programs listed on the x-axis.", style = "font-weight:400"),
                            plotOutput("barPlot"),
                            plotOutput("stackedBar"))
                        )),
               tabPanel("Ridgeline Plots",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("DataFilter2", "Program Filter:",
                                        choices = c("All", unique(metadata$Program))),
                            bsTooltip(id = "DataFilter2",
                                      title = "Select a program to subset the metadata"),
                            selectInput("Variable2", "Grouping:",
                                        choices=colnames(metadata_mini)),
                            bsTooltip(id = "Variable2",
                                      title = "Select a grouping by which to separate donors"),
                            selectInput("MetrictoPlot2", "Metric:",
                                        choices=continuous_vars),
                            bsTooltip(id = "MetrictoPlot2",
                                      title = "Select a metric to plot the distribution of"),
                            downloadButton("DownloadRidges",
                                           "Download")
                          ),
                          mainPanel(
                            p(),
                            p("Select inputs on left to show ridgeline plots of a continuous metric, separated by grouping.", style = "font-weight:400"),
                            plotOutput("ridgelines"))
                        )),
               tabPanel("Clustering Heatmap",
                        sidebarLayout(
                          sidebarPanel(
                            checkboxGroupInput("checkbox", "Metrics to include:",
                                               choices=continuous_vars,
                                               selected=continuous_vars)
                          ),
                          mainPanel(
                            p(),
                            p("Select checkboxes on the left to show a clustered heatmap of donors based on different continuous metrics (e.g. age, BMI). 
                              Metrics are provided for all donor programs for which they are available. 
                              For more information about each program, see ",  tags$a(href = "https://pankbase.org/programs.html", "here", .noWS = "outside"), ". For more information about each donor, explore the ", tags$a(href = "https://data.pankbase.org", "Data Library", .noWS = "outside"), ".", style = "font-weight:400"),
                            plotOutput("heatmap"))
                        )),
               tabPanel("Correlation Matrix",
                        sidebarLayout(
                          sidebarPanel( 
                            selectInput("corr_type", "Correlation method:",
                                        choices=setNames(c("pearson", "spearman"), c("Pearson", "Spearman"))),
                            bsTooltip(id = "corr_type",
                                      title = "Select a correlation method"),
                            checkboxGroupInput("checkbox3", "Metrics to include:",
                                               choices=continuous_vars,
                                               selected=continuous_vars)),
                          mainPanel(
                            p(),
                            p("Select checkboxes and correlation type on the left to show the correlation between different metrics for PanKbase donors. The correlation coefficient (R for Pearson's correlation and Rho for Spearman's correlation) is denoted by the color and number inside the box, and the p-value is indicated by the asterisks. Significance levels are as follows: * p < 0.05, ** p < 0.01, *** p < 0.001. P-values are re-calculated for multiple comparisons based on the number of comparisons selected using the checkboxes.", style = "font-weight:400"),
                            plotOutput("corr_plot"))
                        )),
               tabPanel("Scatter Plot",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("Var1", "x-axis metric",
                                        choices=continuous_vars),
                            bsTooltip(id = "Var1",
                                      title = "Select a metric to plot on the x-axis"),
                            selectInput("Var2", "y-axis metric",
                                        choices=continuous_vars),
                            bsTooltip(id = "Var2",
                                      title = "Select a metric to plot on the y-axis"),
                            selectInput("Color2", "Color points by:",
                                        choices=categorical_vars),
                            bsTooltip(id = "Color2",
                                      title = "Select a metric to color points by"),
                            selectInput("corr_type2", "Correlation method:",
                                        choices=setNames(c("pearson", "spearman"), c("Pearson", "Spearman"))),
                            bsTooltip(id = "corr_type2",
                                      title = "Select a correlation method"),
                            selectInput("plot_lm", "Plot line of best fit with 95% CI?",
                                        choices=c("Yes", "No")),
                            bsTooltip(id = "plot_lm",
                                      title = "Add fit line with confidence interval?"),
                            downloadButton("DownloadScatter",
                                           "Download")),
                          mainPanel(
                            p(),
                            p("Select inputs on the left to create a scatter plot showing the correlation between two metrics.", style = "font-weight:400"),
                            plotOutput("scatterPlot"))
                        )
               ),
               tabPanel("PCA",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("Color", "Color by:",
                                        choices=categorical_vars,
                                        selected = "Description of diabetes status"),
                            bsTooltip(id = "Color",
                                      title = "Select a metric to color points by"),
                            checkboxGroupInput("checkbox2", "Metrics to include:",
                                               choices=continuous_vars,
                                               selected=continuous_vars),
                            selectInput("Ellipses", "Add 95% CI per-group ellipses?",
                                        choices = c("No", "Yes")),
                            bsTooltip(id = "Ellipses",
                                      title = "Add ellipses showing 95% confidence interval per group?")),
                          mainPanel(
                            p(),
                            p("Select metrics on the left to calculate and plot a principal component analysis (PCA) of the donor metadata. Each point in the plot represents one donor. For more information about each donor, explore the ", tags$a(href = "https://data.pankbase.org", "Data Library", .noWS = "outside"), ".", style = "font-weight:400"),
                            plotOutput("pca_plot"))
                          
                        ),
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("PC", "Principal Component",
                                        choices= setNames(1:9,paste0("PC", 1:9))),
                            bsTooltip(id = "PC",
                                      title = "Choose which PC to view")),
                          mainPanel(
                            p(),
                            p("Select a principal component (PC) to show the contribution of different variables to that PC (i.e. how much does that metric influence the PC?) based on the plot above.", style = "font-weight:400"),
                            plotOutput("pca_contribs"))
                        )
               )
               
               
             )
      
    ),
    
    tabPanel("Assays-by-Donors",
             tabsetPanel(
               tabPanel("Matrix",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("Tissue", "Tissue",
                                        choices = c("Islet"),
                                        selected = "Islet")
                            # downloadButton("DownloadDonorMat", "Download")
                          ),
                          mainPanel(
                            p(),
                            p("Browse which assays are available for which donors. For more information about each donor, explore the ", tags$a(href = "https://data.pankbase.org", "Data Library", .noWS = "outside"), ".", style = "font-weight:400"),
                            plotOutput("donor_matrix"))
                        )),
               tabPanel("UpSet Plot",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("Tissue3", "Tissue",
                                        choices = c("Islet"),
                                        selected = "Islet"),
                            checkboxGroupInput("checkbox_upset", "Assays to include:",
                                               choices=unique(islet_df$dataset),
                                               selected=c("Genotyping", "RNA-seq", "ATAC-seq", "scRNA-seq", "snATAC-seq"))
                            # downloadButton("DownloadUpset", "Download")
                          ),
                          mainPanel(
                            p(),
                            p("Select checkboxes on the left to see the number of donors with data available for specific combinations of assays. ", tags$b("Important!"), " These intersections are ", tags$i("inclusive", style = "font-weight:400"), "(e.g. if a donor has Genotyping, scRNAseq, and RNAseq data available, they are included in the count for Genotyping+scRNAseq as well as the count for Genotyping+scRNAseq+RNAseq. This means that an individual donor may be represented in multiple intersection groups. For more information about each donor, explore the ", tags$a(href = "https://data.pankbase.org", "Data Library", .noWS = "outside"), ".", style = "font-weight:400"),
                            plotOutput("upset"))
                        )
                        
               ),
               tabPanel("Stacked Bar Plot",
                        sidebarLayout(
                          sidebarPanel(
                            selectInput("Tissue2", "Tissue",
                                        choices = c("Islet"), #unique(all_donor_df$dataset_tissue)
                                        selected = "Islet"),
                            selectInput("Grouping2", "Grouping",
                                        choices = categorical_vars,
                                        selected = "Program"
                            ),
                            downloadButton("DownloadDonorStacked", "Download")
                          ),
                          mainPanel(
                            p(),
                            p("Select inputs on the left to see the total number of donors who have data for each assay type by the grouping variable. For more information about each donor, explore the ", tags$a(href = "https://data.pankbase.org", "Data Library", .noWS = "outside"), ".", style = "font-weight:400"),
                            plotOutput("donor_stacked_barplot")))
               ),
               tabPanel("Venn Diagram of AAB status",
                        mainPanel(
                          p(),
                          p("Venn Diagram of auto-antibody (AAB) positivity for donors positive for at least one AAB. Currently displayed donors are from the ", 
                            tags$a(href = "https://pankbase.org/hpap-program.html", "Human Pancreas Analysis Program (HPAP)"), 
                            ". For more information about each donor, explore the ", tags$a(href = "https://data.pankbase.org", "Data Library", .noWS = "outside"), ".", style = "font-weight:400"),
                          p(),
                          plotOutput("venndiag")))
               
             )
             ),
    tabPanel("Assay Descriptions",
             h3("Assays included in the Donor Summary Tool are defined as follows:"),
             h1(),
             h5(strong("scRNA-seq:"), "single cell measurement of gene transcript abundance"),
             h5(strong("snATAC-seq:"), "single cell measurement of accessible/open chromatin"),
             h5(strong("HLA typing:"), "serological or genetic measurement of human leukocyte antigen types"),
             h5(strong("Genotyping:"), "measurement of genetic information either my microarray or sequencing"),
             h5(strong("Function:"), "measurement of pancreatic islet function, for example perifusion"),
             h5(strong("RNA-seq:"), "bulk measurement of gene transcript abundance, either from tissue or sorted cells"),
             h5(strong("ATAC-seq:"), "bulk measurement of accessible/open chromatin, either from tissue or sorted cells"),
             h5(strong("WGBS:"), "measurement of methylated DNA based on bisulfate sequencing"),
             h5(strong("Morphology:"), "measurement describing the size, shape, and structure of a tissue and/or cells within a tissue"),
             h5(strong("CyTOF:"), "single cell measurement of protein markers based on mass cytometry"),
             h5(strong("Imaging:"), "measurement describing the location and abundance of specific molecules within a tissue by imaging techniques"),
             h5(strong("Histology:"), "measurement describing the anatomy of a tissue based on microscopy"),
             h5(strong("CODEX:"), "measurement of protein markers in tissue samples using spatial multiplexed imaging"),
             h5(strong("Imaging mass spec:"), "measurement of molecules in tissue samples using mass spectrometry"),
             h5(strong("Single cell multiome:"), "paired single cell measurement of gene transcript abundance and accessible/open chromatin"),
             h5(strong("BCR-seq:"), "paired single cell measurement of B cell receptor sequence and gene transcript abundance"),
             h5(strong("TCR-seq:"), "paired single cell measurement of T cell receptor sequence and gene transcript abundance"),
             h5(strong("CITE-seq:"), "paired single cell measurement of surface protein markers and gene transcript abundance"),
             h5(strong("Flow cytometry:"), "measurement of physical and/or chemical characteristics of cells in a sample, for example to quantify cell type abundance")))
    
    
)


             


# Define server logic required to draw all plots
server <- function(input, output, session) {

  #create plot by itself first
  barPlot_fxn <- function() {
    #filter if we want
    if (input$DataFilter != "All") {
      metadata_mini <- metadata_mini[metadata_mini[,"Program"] == input$DataFilter,]
    }
    
    #subset data to be the data we want
    metadata_summ <- metadata_mini %>% summarise(counts = n(), .by = input$Variable)
    metadata_summ <- metadata_summ[order(metadata_summ[,1]),]
    
    #set up color palettes
    collections <- unique(metadata[,input$Variable]) %>% sort()
    collection_pal <- all_palette(length(collections))
    names(collection_pal) <- collections

    #make a plot
    ggplot(metadata_summ, aes(fill=metadata_summ[,input$Variable], y=counts, x = metadata_summ[,input$Variable])) + 
      geom_bar(stat="identity", position = "dodge2") +
      geom_text(aes(label=metadata_summ[,"counts"]), position=position_dodge(width=0.9), vjust=-0.25, size = 5) +
      scale_fill_manual(values = collection_pal) +
      xlab(input$Variable) + ylab("Number of Donors") + 
      labs(fill = input$Variable) + 
      coord_cartesian(clip = "off") +
      theme(axis.text.x = element_text(angle = 45, hjust=1)) +
      theme(plot.margin = margin(1.5,1,1.5,1, "cm"))
  }
  
  # Fill in the spot we created for a plot
  output$barPlot <- renderPlot({
    barPlot_fxn()
  })
  
  stackedBar_fxn <- function() {
    #set up color palettes
    collections <- unique(metadata[,input$Variable]) %>% sort()
    collection_pal <- all_palette(length(collections))
    names(collection_pal) <- collections

    #filter if we want
    if (input$DataFilter != "All") {
      metadata_mini <- metadata_mini[metadata_mini[,"Program"] == input$DataFilter,]
    }
    
    #subset data to be the data we want
    metadata_summ <- metadata_mini %>% summarise(counts = n(), .by = input$Variable)
    metadata_summ$perc <- (metadata_summ$counts/sum(metadata_summ$counts))*100
    metadata_summ <- metadata_summ[order(metadata_summ[,1]),]
    
    #remove unknown values if you want to
    if (input$RemoveNA == "Yes") {
      metadata_summ <- metadata_summ[metadata_summ[,1] != "unknown" & metadata_summ[,1] != "Unknown" & metadata_summ[,1] != "NA" & metadata_summ[,1] != "",]
    }
    
    #make a plot
    ggplot(metadata_summ, aes(fill=metadata_summ[,input$Variable], y=perc, x = input$DataFilter)) + 
      geom_bar(stat="identity", color = "white") +
      scale_fill_manual(values = collection_pal) +
      ylab("Percent of Donors") + 
      labs(fill = input$Variable) + 
      # theme_classic() +
      theme(plot.margin = margin(1.5,1,1.5,1, "cm"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  output$stackedBar <- renderPlot({
    stackedBar_fxn()
  })
  
  ridgelines_fxn <- function() {
    #set up colors for plot
    collections <- unique(metadata[,input$Variable2]) %>% sort()
    collection_pal <- all_palette(length(collections))
    names(collection_pal) <- collections
    
    #filter if we want
    if (input$DataFilter2 != "All") {
      metadata_filt <- metadata[metadata[,"Program"] == input$DataFilter2,]
    }
    else {
      metadata_filt <- metadata
    }
    
    #remove values that are NA
    metadata_filt <- metadata_filt[!is.na(metadata_filt[,input$Variable2]) & metadata_filt[,input$Variable2] != "unknown" & metadata_filt[,input$Variable2] != "",]
    metadata_filt <- metadata_filt[order(metadata_filt[,input$Variable2]),]
    metadata_filt[,input$Variable2] <- metadata_filt[,input$Variable2] %>% fct_rev()
    
    #need to filter for values where they only occur once so they don't get plotted - otherwise it plots them as empty lines
    metadata_filt <- metadata_filt[metadata_filt[,input$Variable2] %in% unique(metadata_filt[,input$Variable2][duplicated(metadata_filt[,input$Variable2])]),]
    
    #make a plot
    ggplot(metadata_filt, aes(x = metadata_filt[,input$MetrictoPlot2], y = metadata_filt[,input$Variable2], fill = fct_rev(metadata_filt[,input$Variable2]))) + 
      geom_density_ridges(scale = 0.9, 
                          # alpha = 0.8,
                          # rel_min_height = 0.01,
                          quantile_lines = TRUE, quantiles = 2) +
      scale_fill_manual(values = collection_pal) +
      xlim(c(0, max(metadata_filt[,input$MetrictoPlot2]))) + #limit minimum of ridge to zero
      xlab(input$MetrictoPlot2) +
      ylab(input$Variable2) +
      guides(fill=guide_legend(title=input$Variable2)) +
      coord_cartesian(clip = "off") +
      # theme_classic() +
      theme(plot.margin = unit(c(1.5,1,1,1), "cm"))
  }
  
  output$ridgelines <- renderPlot({
    ridgelines_fxn()
  })
  
  output$venndiag <- renderPlot({

    #get lists of positive donors to intersect
    AAB_GADA_donors <- metadata %>% filter(`AAB-GADA Positive` == TRUE) %>% .$ID
    AAB_IA2_donors <- metadata %>% filter(`AAB-IA2 Positive` == TRUE) %>% .$ID
    AAB_IAA_donors <- metadata %>% filter(`AAB-IAA Positive` == TRUE) %>% .$ID
    AAB_ZNT8_donors <- metadata %>% filter(`AAB-ZNT8 Positive` == TRUE) %>% .$ID

    venn_list <- list("GADA positive" = AAB_GADA_donors,
                      "IA2 positive" = AAB_IA2_donors,
                      "IAA positive" = AAB_IAA_donors,
                      "ZNT8 positive" = AAB_ZNT8_donors)

    ggvenn(venn_list,
           fill_color = all_palette(length(venn_list))
          ) +
      coord_cartesian(clip = "off")
  })
  
  output$heatmap <- renderPlot({
    #set up colors for plot
    collections <- unique(metadata$Program) %>% sort()
    collection_pal <- all_palette(length(collections))
    names(collection_pal) <- collections
    
    collections2 <- unique(metadata$`Description of diabetes status`) %>% sort()
    collection_pal2 <- all_palette(length(collections2))
    names(collection_pal2) <- collections2

    set.seed(123)
    metadata_vars <- metadata[,c("Program", "Description of diabetes status", input$checkbox)]
    metadata_vars <- metadata_vars %>% na.omit()
    metadata_mat <- metadata_vars[,input$checkbox] %>% as.matrix() %>% scale()
    
    #set up color palette
    col_fun <- colorRamp2(c(-3, 0, 3), hcl_palette = "Viridis")
    
    row_ha <- rowAnnotation(Program = metadata_vars$Program, `Diabetes status` = metadata_vars$`Description of diabetes status`, col = list(Program = collection_pal, `Diabetes status` = collection_pal2))
    draw(Heatmap(metadata_mat, right_annotation = row_ha, name = "Scaled value", 
                 show_row_names = FALSE, row_title = "Donors", 
                 row_title_gp = gpar(fontsize = 16), row_title_side = "left",
                 col = col_fun), 
         padding = unit(c(5, 5, 5, 5), "mm")) #, width = unit(4, "cm"), height = unit(6, "cm")
  },
  height = 600)
  
  output$corr_plot <- renderPlot({
    #get values you want
    metadata_vars <- metadata[,c("Program", input$checkbox3)]

    if (input$corr_type == "pearson") {
      plot_title <- "Pearson Correlation Matrix"
      legend_title <- "Pearson R"
    }
    else if (input$corr_type == "spearman") {
      plot_title <- "Spearman Correlation Matrix"
      legend_title <- "Spearman Rho"
    }

    corr_all <- rcorr(as.matrix(metadata_vars[,-1]),type=input$corr_type)
    corr <- corr_all$r
    p_mat <- corr_all$P
    n_tests <- (length(p_mat) - length(diag(p_mat)))/2

    # Get p-value matrix
    p.df = as.data.frame(p_mat)
    p.df <-p.df*n_tests
    # Function to get asteriks
    labs.function = function(x){
      case_when(x >= 0.05 ~ "",
                x < 0.05 & x >= 0.01 ~ "*",
                x < 0.01 & x >= 0.001 ~ "**",
                x < 0.001 ~ "***")
    }
    
    # Get asterisks matrix based on p-values
    p.labs = p.df  %>%                      
      mutate_all(labs.function)
    
    # Reshaping asteriks matrix to match ggcorrplot data output
    p.labs$Var1 = as.factor(rownames(p.labs))
    p.labs = melt(p.labs, id.vars = "Var1", variable.name = "Var2", value.name = "lab")

    # Initial ggcorrplot
    cor.plot = ggcorrplot(corr, hc.order = FALSE, type = "lower",
                          lab = TRUE,
                          title = plot_title) +
      scale_fill_gradient2(legend_title,
                           low = "#3A86FF", 
                           mid = "white", 
                           high = "#FD2244",
                           midpoint = 0,
                           limits = c(-1, 1))
    
    # Subsetting asteriks matrix to only those rows within ggcorrplot data
    p.labs$in.df = ifelse(is.na(match(paste0(p.labs$Var1, p.labs$Var2),
                                      paste0(cor.plot[["data"]]$Var1, cor.plot[["data"]]$Var2))),
                          "No", "Yes")

    p.labs = select(filter(p.labs, in.df == "Yes"), -in.df)

    # Add asteriks to ggcorrplot
    cor.plot.labs = cor.plot +
      geom_text(aes(x = p.labs$Var1,
                    y = p.labs$Var2),
                label = p.labs$lab,
                nudge_y = 0.25,
                size = 5) 
    cor.plot.labs

    
  })
  
  output$pca_plot <- renderPlot({
    #set up color palettes
    collections <- unique(metadata[,input$Color]) %>% sort()
    collection_pal <- all_palette(length(collections)) #maybe try #70A7FF for blue or #85B4FF
    names(collection_pal) <- factor(collections)

    df_pca <- metadata[complete.cases(metadata[ , input$checkbox2]),]
    res.pca <- prcomp(df_pca[,input$checkbox2], scale = TRUE)

    #what do you want to color by
    groups <- factor(df_pca[,input$Color])
    unq_pal <- collection_pal[names(collection_pal) %in% unique(groups)]

    #plot pca with color
    if (input$Ellipses == "No") {
      fviz_pca_ind(res.pca,
                   col.ind = groups,
                   palette = collection_pal,
                   label = FALSE,
                   legend.title = input$Color)  +
      xlab("PC1") +
      ylab("PC2") +
      theme(axis.text = element_text(size=12),
            axis.title = element_text(size=14),
            plot.title = element_text(size=16))
    }
    else if (input$Ellipses == "Yes") {
      
      if (length(unique(df_pca[,input$Color])) < 2) {
        nam_vec <- df_pca[,input$Color]
        col_vec <- rep(unq_pal[[1]], nrow(df_pca)) #%>% as.factor()
        names(col_vec) <- nam_vec
        names(unq_pal) <- as.character(names(unq_pal))
        unq_col <- unname(unq_pal)
        fviz_pca_ind(res.pca,
                     col.var = df_pca[,input$Color],
                     col.ind = unq_pal,
                     label = FALSE,
                     addEllipses=TRUE,
                     ellipse.level=0.95,
                     legend.title = input$Color,
                     legend.position = "right",
                     add.legend = TRUE) +
          xlab("PC1") +
          ylab("PC2") +
          theme(axis.text = element_text(size=12),
                axis.title = element_text(size=14),
                plot.title = element_text(size=16))
      }
      
      else {
        fviz_pca_ind(res.pca,
                     col.ind = groups,
                     palette = collection_pal,
                     label = FALSE,
                     addEllipses=TRUE,
                     ellipse.level=0.95,
                     legend.title = input$Color) +
          xlab("PC1") +
          ylab("PC2") +
          theme(axis.text = element_text(size=12),
                axis.title = element_text(size=14),
                plot.title = element_text(size=16))
      }
      
      
    }

  })
  
  output$pca_contribs <- renderPlot({
    df_pca <- metadata[complete.cases(metadata[ , continuous_vars]),]
    res.pca <- prcomp(df_pca[,input$checkbox2], scale = TRUE)
    #plot
    fviz_contrib(res.pca, choice = "var", axes = as.numeric(input$PC), top = 10, fill = "#219197", color = "#219197", ggtheme = theme_classic()) +
      ggtitle(paste0("Contributions of variables to PC", input$PC)) +
      theme(axis.text = element_text(size=12),
            axis.title = element_text(size=14),
            plot.title = element_text(size=16))
  })
  
  scatterPlot_fxn <- function() {
    #set up colors for plot
    collections <- unique(metadata[,input$Color2]) %>% sort()
    collection_pal <- all_palette(length(collections))
    names(collection_pal) <- collections
    
    #filtered df that removes NAs
    metadata_filt <- metadata[!is.na(metadata[,input$Var1]) & 
                                metadata[,input$Var1] != "" & 
                                !is.na(metadata[,input$Var2]) &
                                metadata[,input$Var2] != "",]
    
    #calculate correlation
    c <- cor.test(metadata_filt[,input$Var1], metadata_filt[,input$Var2], method = input$corr_type2) #changed from rcorr
    
    #generate a label for the correlation value shown on the graph
    corr_lab <- ifelse(input$corr_type2 == "pearson", paste0("R = ", signif(c$estimate, 4)), paste0("Rho = ", signif(c$estimate, 4)))
    
    #filter values for NA in NEITHER column to get good dimensions for plot - probably don't need these anymore
    max_filt_val1 <- subset(metadata[,input$Var1], !is.na(metadata[,input$Var2]) & !is.na(metadata[,input$Var1])) %>% max()
    max_filt_val2 <- subset(metadata[,input$Var2], !is.na(metadata[,input$Var1]) & !is.na(metadata[,input$Var2])) %>% max()
    
    min_filt_val1 <- subset(metadata[,input$Var1], !is.na(metadata[,input$Var2]) & !is.na(metadata[,input$Var1])) %>% min()
    min_filt_val2 <- subset(metadata[,input$Var2], !is.na(metadata[,input$Var1]) & !is.na(metadata[,input$Var2])) %>% min()
    
    #plot the data
    p <- ggplot(metadata_filt, aes(x = metadata_filt[,input$Var1], y = metadata_filt[,input$Var2], color = metadata_filt[,input$Color2])) + 
      geom_point() + 
      scale_color_manual(values = collection_pal) +
      annotate("text", 
               x = max_filt_val1*0.9, 
               y=max_filt_val2, 
               label = corr_lab,
               size=8) +
      annotate("text", 
               x = max_filt_val1*0.9, 
               y=max_filt_val2*0.9, 
               label = paste0("p = ", signif(c$p.value, digits = 4)),
               size=8) +
      xlab(input$Var1) +
      ylab(input$Var2) +
      guides(color=guide_legend(title=input$Color2))
    
    #add trend line or not
    if (input$plot_lm == "Yes") {
      p <- p + geom_smooth(method = "lm", color = "black")
    }
    
    print(p)
    

  }
  
  output$scatterPlot <- renderPlot ({
    scatterPlot_fxn()
  })
  
  donor_matrix_fxn <- function() {
    tissue_df <- all_donor_df %>% filter(dataset_tissue %in% c(input$Tissue,  "-"))

    lt <- list()
    for (i in unique(tissue_df$dataset)) {
      dons <- tissue_df$ID[tissue_df$dataset == i]
      lt[[i]] <- dons
    }
    
    df_wide <- tissue_df %>% pivot_wider(names_from = dataset, id_cols = ID, values_from = dataset_tissue)
    cols <- colnames(df_wide)[-1]
    df_wide[cols] <- +(!is.na(df_wide[cols]))
    df_wide <- df_wide %>% merge(metadata, on = "ID")
    df_wide$sum <- rowSums(df_wide[cols])
    df_wide <- df_wide[order(df_wide$`Description of diabetes status`, -df_wide$sum),]

    to_plot <- as.matrix(df_wide[,cols]) %>% t()

    col_fun <- colorRamp2(c(0, 1), c("#fafafa", "#219197"))
    
    collections <- unique(metadata[,"Description of diabetes status"]) %>% sort()
    collection_pal <- all_palette(length(collections))
    names(collection_pal) <- collections
    
    color_this_plot <- collection_pal[c("No diabetes", "Type 1 diabetes", "Type 2 diabetes")]

    ha <- HeatmapAnnotation(`Diabetes status` = df_wide$`Description of diabetes status`,
                             col = list(`Diabetes status` = color_this_plot),
                             show_annotation_name = FALSE,
                             show_legend = FALSE)
    ra <- rowAnnotation(`Total # Donors` = anno_barplot(rowSums(to_plot),
                                                  gp = gpar(fill = "#94c95e", lty="blank"),
                                                  border = FALSE,
                                                  axis_param = list(labels_rot = 0,
                                                                    side = "bottom"),
                                                  width = unit(4, "cm")),
                        `# Donors2` = anno_text(rowSums(to_plot)))
    
    ht = Heatmap(to_plot, name = "Donor Assay Availability", top_annotation = ha, right_annotation = ra,
            cluster_rows = FALSE, cluster_columns = FALSE, 
            row_order = c("scRNA-seq", "RNA-seq", "snATAC-seq", "ATAC-seq", "Genotyping", "HLA typing", "Function", "Morphology", "CyTof", "WGBS", "Imaging"),
            row_names_side = "left",
            col = col_fun, 
            rect_gp = gpar(col = "white", lwd = 1),
            column_split = df_wide$`Description of diabetes status`,
            show_heatmap_legend = FALSE,
            show_column_names = FALSE,
            row_names_gp = grid::gpar(fontsize = 14))
    draw(ht, column_title = "Donors", column_title_gp = gpar(fontsize = 16), column_title_side = "bottom",
         padding = unit(c(1, 1, 1, 1), "cm"))

  }
  
  output$donor_matrix <- renderPlot ({
    donor_matrix_fxn()
  })
  
  donor_stacked_bar_fxn <- function() {
    tissue_df <- all_donor_df %>% filter(dataset_tissue %in% c(input$Tissue2,  "-"))
    tissue_df <- tissue_df %>% merge(metadata, on = "ID")
    
    lt <- list()
    for (i in unique(tissue_df$dataset)) {
      dons <- tissue_df$ID[tissue_df$dataset == i]
      lt[[i]] <- dons
    }
    
    calc_df <- tissue_df %>% merge(metadata, on = "ID") %>% summarise(n_counts = n(), .by = c(input$Grouping2, "dataset"))

    #to get order of stacked bar
    summ_df <- calc_df %>% summarise(total_samps = sum(n_counts), .by = "dataset")
    summ_df <- summ_df[order(summ_df$total_samps, decreasing = TRUE),]
    data_factors <- summ_df$dataset
      
    collections <- unique(metadata[,input$Grouping2]) %>% sort()
    collection_pal <- all_palette(length(collections))
    names(collection_pal) <- collections

    calc_df$dataset <- factor(calc_df$dataset, levels = data_factors)
    ggplot(calc_df, aes(y = n_counts, x = dataset)) + 
      geom_bar(position = "stack", stat="identity", aes(fill = calc_df[,input$Grouping2]), color = "white") + 
      geom_text(data = summ_df, aes(x = dataset, y = total_samps+6, label = total_samps)) +
      scale_fill_manual(values = collection_pal) +
      ylab("Number of Donors") + xlab("Assay") +
      ggtitle(input$Tissue2) +
      guides(fill=guide_legend(title=input$Grouping2)) +
      theme(axis.text.x = element_text(angle = 45, hjust=1), plot.title = element_text(hjust = 0.5))

  }
  
  output$donor_stacked_barplot <- renderPlot ({
    donor_stacked_bar_fxn()
  })
  
  upset_fxn <- function() {
    tissue_df <- all_donor_df %>% filter(dataset_tissue %in% c(input$Tissue3,  "-"))
    tissue_df <- tissue_df %>% filter(dataset %in% input$checkbox_upset)
    
    lt <- list()
    for (i in unique(tissue_df$dataset)) {
      dons <- tissue_df$ID[tissue_df$dataset == i]
      lt[[i]] <- dons
    }
    
    m1 <- make_comb_mat(lt, mode = "intersect")
    ss = set_size(m1)
    cs = comb_size(m1)
    
    ht = UpSet(m1,
               right_annotation = NULL,
               comb_order = order(comb_degree(m1), -cs),
               top_annotation = HeatmapAnnotation(
                 "Number of Donors" = anno_barplot(cs,
                                             ylim = c(0, max(cs)*1.1),
                                             border = FALSE,
                                             height = unit(4, "cm"),
                                             gp = gpar(fill = "#219197", lty="blank"),
                                             axis_param = list(gp = gpar(fontsize=16))),
                                          annotation_name_gp = gpar(fontsize=16),
                 
                 annotation_name_side = "left"
               ),
               row_names_gp = gpar(fontsize=16),
               comb_col = "#94c95e",
               pt_size = unit(8, "mm"))
    ht = draw(ht, padding = unit(c(1, 1, 1, 1), "cm"))
    od = column_order(ht)
    decorate_annotation("Number of Donors", {
      grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(5, "pt"),
                default.units = "native", just = c("center", "bottom"),
                gp = gpar(fontsize = 12, col = "black"), rot = 0)
    })
  }
  
  output$upset <- renderPlot ({
    upset_fxn()
  })
  
  output$scRNAavail <- function() {
    
  }
  
  stacked_assays_fxn <- function() {
    collections <- unique(metadata$Program) %>% sort()
    collection_pal <- all_palette(length(collections))
    names(collection_pal) <- collections
    
    all_donor_program <- all_donor_df %>% merge(metadata[,c("ID", "Program")], on = "ID") %>% filter(dataset_tissue %in% c(input$Tissue4, "-"))
    
    #subset data to be the data we want
    program_summ <- all_donor_program %>% summarise(counts = n(), .by = c("dataset", "Program"))

    dataset_summ <- all_donor_program %>% summarise(total_dataset_counts = n(), .by = "dataset")
    order_df <- dataset_summ[order(dataset_summ[,"total_dataset_counts"]),]
    order <- order_df$dataset
    
    ggplot(program_summ, aes(y = counts, x = factor(dataset, levels = rev(order)))) + 
      geom_bar(position = "stack", stat="identity", aes(fill = Program)) + 
      geom_text(data = dataset_summ, aes(x = dataset, y = total_dataset_counts+6, label = total_dataset_counts)) +
      scale_fill_manual(values = collection_pal) +
      ylab("Number of Donors") + xlab("Assay") +
      guides(fill=guide_legend(title="Program")) +
      theme(axis.text.x = element_text(angle = 45, hjust=1),
            plot.margin = margin(1.5,1,1.5,1, "cm"))
    
  }
  
  output$stacked_assays <- renderPlot ({
    stacked_assays_fxn()
  })

  #make plot downloaders
  output$DownloadBar <- downloadHandler(
    filename = function() { paste0("barplot_", Sys.time(), ".png") },
    content = function(file) {
      # device <- function(..., width, height) grDevices::png(..., width = 8, height = 7, res = 300, units = "in")
      ggsave(file, plot = barPlot_fxn(), width = shinybrowser::get_width()*4, height = shinybrowser::get_height()*4, units="px") #plot=last_plot()
    }
  )
  
  output$DownloadStacked <- downloadHandler(
    filename = function() {paste0("stacked_barplot_", Sys.time(), ".png")},
    content = function(file) {
      ggsave(file, plot = stackedBar_fxn(), width = shinybrowser::get_width()*3, height = shinybrowser::get_height()*4, units = "px") #plot=last_plot()
      # ggsave(file, plot = stackedBar_fxn(), width = 16, height = 12, units = "in") #plot=last_plot()
      
    }
  )
  
  output$DownloadRidges <- downloadHandler(
    filename = function() {paste0("ridgeline_plot_", Sys.time(), ".png")},
    content = function(file) {
      ggsave(file, plot = ridgelines_fxn(), width = shinybrowser::get_width()*3, height = shinybrowser::get_height()*4, units = "px")
      
    }
  )
  
  output$DownloadScatter <- downloadHandler(
    filename = function() {paste0("scatter_plot_", Sys.time(), ".png")},
    content = function(file) {
      ggsave(file, plot = scatterPlot_fxn(), width = shinybrowser::get_width()*3, height = shinybrowser::get_height()*4, units = "px")
    }
  )
  
  output$DownloadDonorMat <- downloadHandler(
    filename = function() {paste0("donor_assay_matrix_", Sys.time(), ".png")},
    content = function(file) {
      png(file, width = 12, height = 4, units = "in", res = 300)
      donor_matrix_fxn()
      dev.off()
    }
  )
  
  output$DownloadUpset <- downloadHandler(
    filename = function() {paste0("upset_plot_", Sys.time(), ".png")},
    content = function(file) {
      png(file, width = 12, height = 6, units = "in", res = 300)
      upset_fxn()
      dev.off()
    }
  )
  
  output$DownloadDonorStacked <- downloadHandler(
    filename = function() {paste0("donor_stacked_barplot_", Sys.time(), ".png")},
    content = function(file) {
      ggsave(file, plot = donor_stacked_bar_fxn(), width = 12, height = 8, units = "in")
    }
  )
  
}


# Run the application 
shinyApp(ui = ui, server = server)
