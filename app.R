#
# This is the Shiny web application that runs the donor summary browser for PanKbase
# See it live at https://dev.pankbase.org/donor-metadata.html
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
library(forcats)
library(shinybrowser)
library(colorRamp2)
library(shinyWidgets)
library(shinyBS)
library(ggcorrplot)
library(reshape2)
library(viridis)

set.seed(123)

#read in metadata files
metadata <- read.csv("data/pankbase_human_donor_report_2025_5_27_17h_8m.csv")

#filter out "test" sample
metadata <- metadata[metadata$ID != "<NA>",]

#change capitalization of some metadata fields
metadata$Description.of.diabetes.status <- recode(metadata$Description.of.diabetes.status,
                                                  "non-diabetic" = "No diabetes",
                                                  "type 1 diabetes" = "Type 1 diabetes",
                                                  "type 2 diabetes" = "Type 2 diabetes",
                                                  "cystic fibrosis diabetes" = "Cystic fibrosis diabetes",
                                                  "gestational diabetes" = "Gestational diabetes",
                                                  "diabetes unspecified" = "Diabetes unspecified",
                                                  "steroid-induced diabetes" = "Steroid-induced diabetes",
                                                  "monogenic diabetes" = "Monogenic diabetes")

metadata$Ethnicities <- recode(na_if(metadata$Ethnicities, ""),
                               "Caucasian" = "White",
                               .missing = "Unknown")
metadata$Sex <- recode(metadata$Sex,
                       female = "Female",
                       male = "Male")
metadata$Collections <- recode(metadata$Collections,
                               "IIDP,Prodo" = "IIDP")

#rename values that have different capitalization so they get grouped together and change blanks to unknown
metadata$Cause.of.Death <- recode(na_if(metadata$Cause.of.Death, ""),
                                  "Cerebrovascular/stroke" = "Cerebrovascular/Stroke",
                                  "Head Trauma" = "Head trauma",
                                  "ICH/stroke" = "ICH/Stroke",
                                  "Cerebral Edema (DKA)" = "Cerebral edema (DKA)",
                                  "IHC" = "ICH",
                                  .missing = "Unknown")



#rename columns to friendly names
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

#change NA AAB status to unknown
metadata <- metadata %>% mutate(across(starts_with("AAB."), ~ifelse( is.na(.x), "Unknown", .x)))

metadata <- metadata %>% rename("AAB-GADA Positive" = AAB.GADA.POSITIVE,
                                "AAB-IA2 Positive" = AAB.IA2.POSITIVE,
                                "AAB-IAA Positive" = AAB.IAA.POSITIVE,
                                "AAB-ZNT8 Positive" = AAB.ZNT8.POSITIVE)

#get all categorical variables we want to plot
categorical_vars <- c("Program", "Description of diabetes status", "Cause of death", "Sex", "Ethnicity",
                      "AAB-GADA Positive", "AAB-IA2 Positive", "AAB-IAA Positive", "AAB-ZNT8 Positive")
#get all continuous variables we want to plot
continuous_vars <- c("Age (years)", "BMI", "C. Peptide (ng/ml)", "HbA1C percentage", "Hospital stay (hours)",
                     "AAB-GADA value (unit/ml)", "AAB-IA2 value (unit/ml)", "AAB-IAA value (unit/ml)", "AAB-ZNT8 value (unit/ml)")

#subset df to most important cols
metadata_mini <- metadata[,categorical_vars]

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

#make dataset column look pretty
all_donor_df$dataset <- recode(all_donor_df$dataset,
                               "HLA_typing" = "HLA typing",
                               "RNAseq" = "RNA-seq",
                               "scRNAseq" = "scRNA-seq",
                               "snATACseq" = "snATAC-seq")

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
                                           "Download Bar Plot"),
                            selectInput("RemoveNA", "Remove 'unknowns' from stacked bar?",
                                        choices=c("No", "Yes")),
                            bsTooltip(id = "RemoveNA",
                                      title = "Choose whether to remove donors whose metadata for that variable is not currently in PanKbase"),
                            downloadButton("DownloadStacked", 
                                           "Download Stacked Bar Plot")
                            
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
                                               selected=continuous_vars),
                            downloadButton("DownloadHeatmap", 
                                           "Download")
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
                                               selected=continuous_vars),
                            downloadButton("DownloadCorrMat", 
                                           "Download")),
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
                                      title = "Add ellipses showing 95% confidence interval per group?"),
                            downloadButton("DownloadPCA", "Download")),
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
                                      title = "Choose which PC to view"),
                            downloadButton("DownloadContribs", "Download")),
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
                                        selected = "Islet"),
                            downloadButton("DownloadDonorMat", "Download")
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
                                               selected=c("Genotyping", "RNA-seq", "ATAC-seq", "scRNA-seq", "snATAC-seq")),
                            downloadButton("DownloadUpset", "Download")
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
                        sidebarLayout(
                          sidebarPanel(
                            downloadButton("DownloadVenn", "Download")
                          ),
                        mainPanel(
                          p(),
                          p("Venn Diagram of auto-antibody (AAB) positivity for donors positive for at least one AAB. Currently displayed donors are from the ", 
                            tags$a(href = "https://pankbase.org/hpap-program.html", "Human Pancreas Analysis Program (HPAP)"), 
                            ". For more information about each donor, explore the ", tags$a(href = "https://data.pankbase.org", "Data Library", .noWS = "outside"), ".", style = "font-weight:400"),
                          p(),
                          plotOutput("venndiag")))
               
             )
             )),
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

  ##for bar plots tab
  
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
      labs(x = input$Variable,
           y = "Number of Donors",
           fill = input$Variable) + 
      coord_cartesian(clip = "off") +
      theme(axis.text.x = element_text(angle = 45, hjust=1), plot.margin = margin(1.5,1,1.5,1, "cm"))
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
    
    if (nrow(metadata_summ) > 0) {
      #make a plot if there is data to plot
      p <- ggplot(metadata_summ, aes(fill=metadata_summ[,input$Variable], y=perc, x = input$DataFilter)) + 
            geom_bar(stat="identity", color = "white") +
            scale_fill_manual(values = collection_pal) +
            labs(y = "Percent of Donors",
                 fill = input$Variable) + 
            theme(plot.margin = margin(1.5,1,1.5,1, "cm"),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank())
    }
    
    else {
      #if no data, print error message
      par(mar = c(0,0,0,0))
      p <- plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      p <- p + text(x = 0.5, y = 0.5, paste("No data fit these criteria.\n",
                                            "Please select other parameters.\n"), 
                    cex = 1.6, col = "black")
    }
    print(p)
    
  }
  
  output$stackedBar <- renderPlot({
    stackedBar_fxn()
  })

  
  ##for ridgelines tab
  
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
    metadata_filt <- metadata_filt[!is.na(metadata_filt[,input$Variable2]) & !metadata_filt[,input$Variable2] %in% c("unknown", "NA", ""),]
    metadata_filt <- metadata_filt[!is.na(metadata_filt[,input$MetrictoPlot2]),]
    metadata_filt <- metadata_filt[order(metadata_filt[,input$Variable2]),]
    metadata_filt[,input$Variable2] <- metadata_filt[,input$Variable2] %>% fct_rev()
    
    #filter because must have at least 3 data points to plot a ridge
    metadata_filt <- metadata_filt %>% filter(n() > 2, .by = input$Variable2)
    
    
    if (nrow(metadata_filt) > 1) {
      #make a plot if there is enough data to make a plot
      p <- ggplot(metadata_filt, aes(x = metadata_filt[,input$MetrictoPlot2], y = metadata_filt[,input$Variable2], fill = fct_rev(metadata_filt[,input$Variable2]))) + 
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
          theme(plot.margin = unit(c(1.5,1,1,1), "cm"))
    }
    
    else {
      #if no data, print error message
      par(mar = c(0,0,0,0))
      p <- plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      p <- p + text(x = 0.5, y = 0.5, paste("No data fit these criteria.\n",
                                            "Please select other parameters.\n"), 
                    cex = 1.6, col = "black")

    }
    print(p)
    
  }
  
  output$ridgelines <- renderPlot({
    ridgelines_fxn()
  })
  
  
  ##for clustering heatmap tab
  
  heatmap_fxn <- function() {
    #set up colors for plot
    collections <- unique(metadata$Program) %>% sort()
    collection_pal <- all_palette(length(collections))
    names(collection_pal) <- collections
    
    collections2 <- unique(metadata$`Description of diabetes status`) %>% sort()
    collection_pal2 <- all_palette(length(collections2))
    names(collection_pal2) <- collections2
    
    if (length(input$checkbox) >= 2) {
      #create matrix to plot
      metadata_vars <- metadata[,c("Program", "Description of diabetes status", input$checkbox)]
      metadata_vars <- metadata_vars[complete.cases(metadata[,continuous_vars]),]
      metadata_vars[,input$checkbox] <- lapply(metadata_vars[,input$checkbox], as.numeric)
      metadata_mat <- metadata_vars[,input$checkbox] %>% as.matrix() %>% scale()
      
      #set up color palette
      col_fun <- colorRamp2(c(-3, 0, 3), hcl_palette = "Viridis")
      
      #draw heatmap
      row_ha <- rowAnnotation(Program = metadata_vars$Program, `Diabetes status` = metadata_vars$`Description of diabetes status`, col = list(Program = collection_pal, `Diabetes status` = collection_pal2))
      draw(Heatmap(metadata_mat, right_annotation = row_ha, name = "Scaled value", 
                   show_row_names = FALSE, row_title = "Donors", 
                   row_title_gp = gpar(fontsize = 16), row_title_side = "left",
                   col = col_fun), 
           padding = unit(c(5, 5, 5, 5), "mm")) 
    }
    
    else {
      #if no data, print error message
      par(mar = c(0,0,0,0))
      p <- plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      p <- p + text(x = 0.5, y = 0.5, paste("Please select at least 2 variables.\n"), 
                    cex = 1.6, col = "black")
      print(p)
    }
    
  }
  
  output$heatmap <- renderPlot({
    heatmap_fxn()
  },
  height = 600)
  
  
  ##for correlation matrix tab
  corr_plot_fxn <- function() {
    
    if (length(input$checkbox3) >= 2) {
      #get values you want
      metadata_vars <- metadata[,c("Program", input$checkbox3)]
      metadata_vars[,input$checkbox3] <- lapply(metadata_vars[,input$checkbox3], as.numeric)
      
      if (input$corr_type == "pearson") {
        plot_title <- "Pearson Correlation Matrix"
        legend_title <- "Pearson R"
      }
      else if (input$corr_type == "spearman") {
        plot_title <- "Spearman Correlation Matrix"
        legend_title <- "Spearman Rho"
      }
      
      #get correlations
      corr_all <- rcorr(as.matrix(metadata_vars[,-1]),type=input$corr_type)
      corr <- corr_all$r
      p_mat <- corr_all$P
      n_tests <- (length(p_mat) - length(diag(p_mat)))/2
      
      # Get p-value matrix
      p.df = as.data.frame(p_mat)
      p.df <-p.df*n_tests
      # Function to get asterisks
      labs.function = function(x){
        case_when(x >= 0.05 ~ "",
                  x < 0.05 & x >= 0.01 ~ "*",
                  x < 0.01 & x >= 0.001 ~ "**",
                  x < 0.001 ~ "***")
      }
      
      # Get asterisks matrix based on p-values
      p.labs = p.df  %>%                      
        mutate_all(labs.function)
      
      # Reshaping asterisks matrix to match ggcorrplot data output
      p.labs$Var1 = as.factor(rownames(p.labs))
      p.labs = melt(p.labs, id.vars = "Var1", variable.name = "Var2", value.name = "lab")
      
      # Initial ggcorrplot
      cor.plot = ggcorrplot(corr, hc.order = FALSE, type = "lower",
                            lab = TRUE,
                            lab_size = 6,
                            tl.cex = 16,
                            title = plot_title) +
        scale_fill_gradient2(legend_title,
                             low = "#3A86FF", 
                             mid = "white", 
                             high = "#FD2244",
                             midpoint = 0,
                             limits = c(-1, 1)) +
        theme(plot.title = element_text(size = 16))
      
      # Subsetting asterisks matrix to only those rows within ggcorrplot data
      p.labs$in.df = ifelse(is.na(match(paste0(p.labs$Var1, p.labs$Var2),
                                        paste0(cor.plot[["data"]]$Var1, cor.plot[["data"]]$Var2))),
                            "No", "Yes")
      
      p.labs = select(filter(p.labs, in.df == "Yes"), -in.df)
      
      # Add asterisks to ggcorrplot
      cor.plot.labs = cor.plot +
        geom_text(aes(x = p.labs$Var1,
                      y = p.labs$Var2),
                  label = p.labs$lab,
                  nudge_y = 0.25,
                  size = 8) 
      print(cor.plot.labs)
    }
    
    else {
      #if no data, print error message
      par(mar = c(0,0,0,0))
      p <- plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      p <- p + text(x = 0.5, y = 0.5, paste("Please select at least 2 variables.\n"), 
                    cex = 1.6, col = "black")
      print(p)
    }
    
  }
  
  output$corr_plot <- renderPlot({
    corr_plot_fxn()
  },
  height = 600)
  
  ##for scatter plot tab
  
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
    c <- cor.test(metadata_filt[,input$Var1], metadata_filt[,input$Var2], method = input$corr_type2)
    
    #generate a label for the correlation value shown on the graph
    corr_lab <- ifelse(input$corr_type2 == "pearson", paste0("R = ", signif(c$estimate, 4)), paste0("Rho = ", signif(c$estimate, 4)))
    
    #filter values for NA in NEITHER column to get good location to put stats on plot
    max_filt_val1 <- subset(metadata[,input$Var1], !is.na(metadata[,input$Var2]) & !is.na(metadata[,input$Var1])) %>% max()
    max_filt_val2 <- subset(metadata[,input$Var2], !is.na(metadata[,input$Var1]) & !is.na(metadata[,input$Var2])) %>% max()
    
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
  
  
  ##for PCA tab
  
  pca_fxn <- function() {
    #set up color palettes
    collections <- unique(metadata[,input$Color]) %>% sort()
    collection_pal <- all_palette(length(collections))
    names(collection_pal) <- factor(collections)
    
    
    if (length(input$checkbox2) >= 2) {
      #do PCA
      df_pca <- metadata[complete.cases(metadata[ , continuous_vars]),]
      df_pca[,input$checkbox2] <- lapply(df_pca[,input$checkbox2], as.numeric)
      res.pca <- prcomp(df_pca[,input$checkbox2], scale = TRUE)
      
      #what do you want to color by
      groups <- factor(df_pca[,input$Color])
      unq_pal <- collection_pal[names(collection_pal) %in% unique(groups)]
      
      #plot pca with color
      if (input$Ellipses == "No") {
        p <- fviz_pca_ind(res.pca,
                          col.ind = groups,
                          palette = collection_pal,
                          label = FALSE,
                          legend.title = input$Color)  +
          labs(x = "PC1",
               y = "PC2") +
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
          p <- fviz_pca_ind(res.pca,
                            col.var = df_pca[,input$Color],
                            col.ind = unq_pal,
                            label = FALSE,
                            addEllipses=TRUE,
                            ellipse.level=0.95,
                            legend.title = input$Color,
                            legend.position = "right",
                            add.legend = TRUE) +
            labs(x = "PC1",
                 y = "PC2") +
            theme(axis.text = element_text(size=12),
                  axis.title = element_text(size=14),
                  plot.title = element_text(size=16))
        }
        
        else {
          p <- fviz_pca_ind(res.pca,
                            col.ind = groups,
                            palette = collection_pal,
                            label = FALSE,
                            addEllipses=TRUE,
                            ellipse.level=0.95,
                            legend.title = input$Color) +
            labs(x = "PC1", 
                 y = "PC2") +
            theme(axis.text = element_text(size=12),
                  axis.title = element_text(size=14),
                  plot.title = element_text(size=16))
        }
        
      }
    }
    
    else {
      #if no data, print error message
      par(mar = c(0,0,0,0))
      p <- plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      p <- p + text(x = 0.5, y = 0.5, paste("Please select at least 2 variables.\n"), 
                    cex = 1.6, col = "black")
    }
    
    print(p)
    
  }
  
  output$pca_plot <- renderPlot({
    pca_fxn()
    
  })
  
  pca_contribs_fxn <- function() {
    
    if (length(input$checkbox2) >= 2) {
      #do PCA
      df_pca <- metadata[complete.cases(metadata[ , continuous_vars]),]
      df_pca[,input$checkbox2] <- lapply(df_pca[,input$checkbox2], as.numeric)
      res.pca <- prcomp(df_pca[,input$checkbox2], scale = TRUE)
      
      #plot variable contributions
      p <- fviz_contrib(res.pca, choice = "var", axes = as.numeric(input$PC), top = 10, fill = "#219197", color = "#219197", ggtheme = theme_classic()) +
        ggtitle(paste0("Contributions of variables to PC", input$PC)) +
        theme(axis.text = element_text(size=12),
              axis.title = element_text(size=14),
              plot.title = element_text(size=16))
          }
    
    else {
      #if no data, print error message
      par(mar = c(0,0,0,0))
      p <- plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      p <- p + text(x = 0.5, y = 0.5, paste("Please select at least 2 variables.\n"), 
                    cex = 1.6, col = "black")
    }
    
    print(p)
    
  }
  
  output$pca_contribs <- renderPlot({
    pca_contribs_fxn()
  })
  
  
  ##for donor matrix tab
  
  donor_matrix_fxn <- function() {
    #filter df for tissue of interest
    tissue_df <- all_donor_df %>% filter(dataset_tissue %in% c(input$Tissue,  "-"))
    
    #get number of donors for each assay, then sort descending
    df_wide <- tissue_df %>% pivot_wider(names_from = dataset, id_cols = ID, values_from = dataset_tissue)
    cols <- colnames(df_wide)[-1]
    df_wide[cols] <- +(!is.na(df_wide[cols]))
    df_wide <- df_wide %>% merge(metadata, on = "ID")
    df_wide$sum <- rowSums(df_wide[cols])
    df_wide <- df_wide[order(df_wide$`Description of diabetes status`, -df_wide$sum),]
    
    #get out columns you want to plot in matrix
    to_plot <- as.matrix(df_wide[,cols]) %>% t()
    
    #set up colors
    col_fun <- colorRamp2(c(0, 1), c("#fafafa", "#219197"))
    collections <- unique(metadata[,"Description of diabetes status"]) %>% sort()
    collection_pal <- all_palette(length(collections))
    names(collection_pal) <- collections
    color_this_plot <- collection_pal[c("No diabetes", "Type 1 diabetes", "Type 2 diabetes")]
    
    #get top heatmap annotation
    ha <- HeatmapAnnotation(`Diabetes status` = df_wide$`Description of diabetes status`,
                            col = list(`Diabetes status` = color_this_plot),
                            show_annotation_name = FALSE,
                            show_legend = FALSE)
    #add barplot of sums on the right
    ra <- rowAnnotation(`Total # Donors` = anno_barplot(rowSums(to_plot),
                                                        gp = gpar(fill = "#94c95e", lty="blank"),
                                                        border = FALSE,
                                                        axis_param = list(labels_rot = 0,
                                                                          side = "bottom"),
                                                        width = unit(4, "cm")),
                        `# Donors2` = anno_text(rowSums(to_plot)))
    
    #draw matrix
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
  
  
  ##for upset plot tab
  
  upset_fxn <- function() {
    #make sure there is at least one box checked
    if (length(input$checkbox_upset) >= 1) {
      #filter to tissue of interest
      tissue_df <- all_donor_df %>% filter(dataset_tissue %in% c(input$Tissue3,  "-"))
      tissue_df <- tissue_df %>% filter(dataset %in% input$checkbox_upset)
      
      #get list of donors in each dataset
      lt <- list()
      for (i in unique(tissue_df$dataset)) {
        dons <- tissue_df$ID[tissue_df$dataset == i]
        lt[[i]] <- dons
      }
      
      #make combination matrix
      m1 <- make_comb_mat(lt, mode = "intersect")
      ss = set_size(m1)
      cs = comb_size(m1)
      
      #plot upset plot
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
    
    else {
      #if no data, print error message
      print("no data!")
      par(mar = c(0,0,0,0))
      p <- plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      p <- p + text(x = 0.5, y = 0.5, paste("Please select at least 1 variable.\n"), 
                    cex = 1.6, col = "black")
      print(p)
    }
    
    
    
  }
  
  output$upset <- renderPlot ({
    upset_fxn()
  })
  
  
  ##for stacked bar plot tab
  
  donor_stacked_bar_fxn <- function() {
    #filter df for tissue of interest
    tissue_df <- all_donor_df %>% filter(dataset_tissue %in% c(input$Tissue2,  "-"))
    tissue_df <- tissue_df %>% merge(metadata, on = "ID")
    
    #count number of donors in each group
    calc_df <- tissue_df %>% merge(metadata, on = "ID") %>% summarise(n_counts = n(), .by = c(input$Grouping2, "dataset"))
    
    #to get order of stacked bar
    summ_df <- calc_df %>% summarise(total_samps = sum(n_counts), .by = "dataset")
    summ_df <- summ_df[order(summ_df$total_samps, decreasing = TRUE),]
    data_factors <- summ_df$dataset
    
    #set up color palette
    collections <- unique(metadata[,input$Grouping2]) %>% sort()
    collection_pal <- all_palette(length(collections))
    names(collection_pal) <- collections
    
    #re-factor to get correct order
    calc_df$dataset <- factor(calc_df$dataset, levels = data_factors)
    
    #plot
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
  
  
  
  ##for venn diagram tab
  
  venn_diag_fxn <- function() {
    #get lists of positive donors to intersect
    AAB_GADA_donors <- metadata %>% filter(`AAB-GADA Positive` == TRUE) %>% .$ID
    AAB_IA2_donors <- metadata %>% filter(`AAB-IA2 Positive` == TRUE) %>% .$ID
    AAB_IAA_donors <- metadata %>% filter(`AAB-IAA Positive` == TRUE) %>% .$ID
    AAB_ZNT8_donors <- metadata %>% filter(`AAB-ZNT8 Positive` == TRUE) %>% .$ID
    
    venn_list <- list("GADA positive" = AAB_GADA_donors,
                      "IA2 positive" = AAB_IA2_donors,
                      "IAA positive" = AAB_IAA_donors,
                      "ZNT8 positive" = AAB_ZNT8_donors)
    
    p <- ggvenn(venn_list,
           fill_color = all_palette(length(venn_list))) +
      coord_cartesian(clip = "off")
    print(p)
  }
  
  output$venndiag <- renderPlot({
    venn_diag_fxn()
  })
  



  #make plot downloaders
  output$DownloadBar <- downloadHandler(
    filename = function() { paste0("PanKbase_BarPlot_", Sys.time(), ".png") },
    content = function(file) {
      ggsave(file, plot = barPlot_fxn(), width = shinybrowser::get_width()*4, height = shinybrowser::get_height()*4, units="px", bg = "white")
    }
  )
  
  output$DownloadStacked <- downloadHandler(
    filename = function() {paste0("PanKbase_StackedBarPlot_", Sys.time(), ".png")},
    content = function(file) {
      ggsave(file, plot = stackedBar_fxn(), width = shinybrowser::get_width()*3, height = shinybrowser::get_height()*4, units = "px", bg = "white")
    }
  )
  
  output$DownloadRidges <- downloadHandler(
    filename = function() {paste0("PanKbase_RidgelinePlot_", Sys.time(), ".png")},
    content = function(file) {
      ggsave(file, plot = ridgelines_fxn(), width = shinybrowser::get_width()*3, height = shinybrowser::get_height()*4, units = "px", bg = "white")
    }
  )
  
  output$DownloadHeatmap <- downloadHandler(
    filename = function() {paste0("PanKbase_ClusterHeatmap_", Sys.time(), ".png")},
    content = function(file) {
      png(file, width = 12, height = 8, units = "in", res = 300) 
      heatmap_fxn()
      dev.off()
    }
  )
  
  output$DownloadCorrMat <- downloadHandler(
    filename = function() {paste0("PanKbase_CorrelationMatrix_", Sys.time(), ".png")},
    content = function(file) {
      ggsave(file, plot = corr_plot_fxn(), width = shinybrowser::get_width()*3, height = shinybrowser::get_height()*3, units = "px", bg = "white")
    }
  )
  
  output$DownloadScatter <- downloadHandler(
    filename = function() {paste0("PanKbase_ScatterPlot_", Sys.time(), ".png")},
    content = function(file) {
      ggsave(file, plot = scatterPlot_fxn(), width = shinybrowser::get_width()*3, height = shinybrowser::get_height()*4, units = "px", bg = "white")
    }
  )
  
  output$DownloadPCA <- downloadHandler(
    filename = function() {paste0("PanKbase_PCA_", Sys.time(), ".png")},
    content = function(file) {
      png(file, width = 10, height = 8, units = "in", res = 300)
      pca_fxn()
      dev.off()
    }
  )
  
  output$DownloadContribs <- downloadHandler(
    filename = function() {paste0("PanKbase_PCA_Contributions_", Sys.time(), ".png")},
    content = function(file) {
      png(file, width = 12, height = 8, units = "in", res = 300) 
      pca_contribs_fxn()
      dev.off()
    }
  )
  
  output$DownloadDonorMat <- downloadHandler(
    filename = function() {paste0("PanKbase_DonorAssayMatrix_", Sys.time(), ".png")},
    content = function(file) {
      png(file, width = 12, height = 4, units = "in", res = 300) 
      donor_matrix_fxn()
      dev.off()
    }
  )
  
  output$DownloadUpset <- downloadHandler(
    filename = function() {paste0("PanKbase_UpSetPlot_", Sys.time(), ".png")},
    content = function(file) {
      png(file, width = 12, height = 6, units = "in", res = 300)
      upset_fxn()
      dev.off()
    }
  )
  
  output$DownloadDonorStacked <- downloadHandler(
    filename = function() {paste0("PanKbase_DonorStackedBarPlot_", Sys.time(), ".png")},
    content = function(file) {
      ggsave(file, plot = donor_stacked_bar_fxn(), width = 12, height = 8, units = "in", bg = "white")
    }
  )
  
  output$DownloadVenn <- downloadHandler(
    filename = function() {paste0("PanKbase_VennDiagram_", Sys.time(), ".png")},
    content = function(file) {
      ggsave(file, plot = venn_diag_fxn(), width = 12, height = 8, units = "in", bg = "white")
    }
  )
  
}


# Run the application 
shinyApp(ui = ui, server = server)
