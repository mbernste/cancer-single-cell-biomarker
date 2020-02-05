library(shiny)
library(ggplot2)
library(rhdf5)
library(viridis)
library(heatmaply)

#TUMORS = c('PJ016', 'PJ018', 'PJ025', 'PJ048', 'PJ030', 'PJ035', 'PJ017', 'PJ032')
TUMORS = c('PJ025')
DATA_DIR = '/Users/matthewbernstein/Development/single-cell-hackathon/data'
TMP_DIR = '/Users/matthewbernstein/Development/single-cell-hackathon/tmp'

#print('Loading counts data...')
counts = h5read(paste0(DATA_DIR, "/GSE103224_normalized.h5"),"count")
cells = h5read(paste0(DATA_DIR, "/GSE103224_normalized.h5"), "cell")
tumor_ids = h5read(paste0(DATA_DIR, "/GSE103224_normalized.h5"), "tumor")
genes = h5read(paste0(DATA_DIR, "/GSE103224_normalized.h5"), "gene_name")

counts <- t(counts)
counts <- data.frame(counts)
rownames(counts) <- cells
colnames(counts) <- genes

tumor_dfs <- c()
for (tumor in TUMORS) {
    # Load data matrix
    tumor_indices <- c()
    for (i in 1:length(cells)) {
        if (tumor_ids[[i]] == tumor) {
            tumor_indices <- append(tumor_indices, i)
        }
    }
    tumor_cells <- cells[tumor_indices]
    tumor_counts <- counts[tumor_indices,]
    rownames(tumor_counts) <- tumor_cells
    colnames(tumor_counts) <- genes
    
    # Read dimension-reduction coordinates
    df <- read.csv(
        file=paste0(TMP_DIR, '/', tumor, '_PHATE.tsv'), 
        sep = '\t', 
        header = FALSE
    )
    colnames(df) <- c('PHATE 1', 'PHATE 2')
    rownames(df) <- tumor_cells
    df <- data.frame(df)
    print('Merging...')
    df <- merge(df, tumor_counts,  by=0)
    print('done.')
    tumor_dfs[[tumor]] <- df
}

aligned_umap_df <- read.csv(
    paste0(TMP_DIR, '/all_tumors_aligned_UMAP.tsv'),
    sep = '\t',
    header = TRUE,
    row.names = 1
)

gene_df <- counts['SOX2']
rownames(gene_df) <- cells
curr_df <- merge(aligned_umap_df, gene_df, by=0)
print(curr_df)

server <- function(input, output) {
    output$plot1<-renderPlot({
        ggplot(tumor_dfs[[input$tumor1]], aes(x=PHATE.1,y=PHATE.2, color=get(input$colorby1))) + 
        scale_color_viridis() + 
        geom_point() + 
        theme_classic()
    })
    output$plot2<-renderPlot({
        ggplot(tumor_dfs[[input$tumor2]], aes(x=PHATE.1,y=PHATE.2, color=get(input$colorby2))) + 
        scale_color_viridis() + 
        geom_point() + 
        theme_classic()
    })
    output$plot3<-renderPlot({
        ggplot(tumor_dfs[[input$tumor3]], aes(x=PHATE.1,y=PHATE.2, color=get(input$colorby3))) + 
        scale_color_viridis() + 
        geom_point() + 
        theme_classic()
    })
    output$plot4<-renderPlot({
        ggplot(tumor_dfs[[input$tumor4]], aes(x=PHATE.1,y=PHATE.2, color=get(input$colorby4))) + 
        scale_color_viridis() + 
        geom_point() + 
        theme_classic()
    })
    

    output$aligned_plot1 <- renderPlot({
        colby <- input$aligned_colorby1
        if (colby != 'tumor') {
            gene_df <- counts[toString(input$aligned_colorby1)]
            rownames(gene_df) <- cells
            curr_df <- merge(aligned_umap_df, gene_df, by=0)
            ggplot(curr_df, aes(x=UMAP1, y=UMAP2, color=get(input$aligned_colorby1))) +
            scale_color_viridis() +
            geom_point(size = 0.1) +
            theme_classic()
        }
        else {
            curr_df <- aligned_umap_df
            ggplot(curr_df, aes(x=UMAP1, y=UMAP2, color=get(input$aligned_colorby1))) +
            geom_point(size = 0.1) +
            theme_classic()
        }
    })
    output$aligned_plot2 <- renderPlot({
        colby <- input$aligned_colorby2
        if (colby != 'tumor') {
            gene_df <- counts[toString(input$aligned_colorby2)]
            rownames(gene_df) <- cells
            curr_df <- merge(aligned_umap_df, gene_df, by=0)
            ggplot(curr_df, aes(x=UMAP1, y=UMAP2, color=get(input$aligned_colorby2))) +
            scale_color_viridis() +
            geom_point(size = 0.1) +
            theme_classic()
        }
        else {
            curr_df <- aligned_umap_df
            ggplot(curr_df, aes(x=UMAP1, y=UMAP2, color=get(input$aligned_colorby2))) +
            geom_point(size = 0.1) +
            theme_classic()
        }
    })
}   

ui <- fluidPage(

  # App title ----
  titlePanel(title=div("CACTUS: CharACterizing TUmor Subpopulations", img(src="Cactus.png", height=60, width=40))),

  tabsetPanel(
        tabPanel("Individual tumors", 
            sidebarLayout(
                sidebarPanel(
                    selectInput("tumor1", "Plot 1 tumor:", TUMORS),
                    textInput("colorby1", "Select gene:"), 
                    selectInput("tumor2", "Plot 2 tumor:", TUMORS),
                    textInput("colorby2", "Select gene:"), 
                    selectInput("tumor3", "Plot 3 tumor:", TUMORS),
                    textInput("colorby3", "Select gene:"),
                    selectInput("tumor4", "Plot 4 tumor:", TUMORS),
                    textInput("colorby4", "Select gene:")
                ),               
                mainPanel(
                    "main panel",
                    fluidRow(
                        splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot1"), plotOutput("plot2"))
                    ),
                    fluidRow(
                        splitLayout(cellWidths = c("50%", "50%"), plotOutput("plot3"), plotOutput("plot4"))
                    )
                )
            )
        ),
        tabPanel("Cluster GSEA", 
            heatmaply(mtcars) 
        ),
        tabPanel("Tumors integrated",
           sidebarLayout(
                sidebarPanel(
                    textInput("aligned_colorby1", "Plot 1, color by:"),
                    textInput("aligned_colorby2", "Plot 2, color by:")
                ),
                mainPanel(
                    fluidRow(
                        splitLayout(cellWidths = c("50%", "50%"), plotOutput("aligned_plot1"), plotOutput("aligned_plot2"))
                    )
                )
            )
        )
  )
)

shinyApp(ui = ui, server = server)
