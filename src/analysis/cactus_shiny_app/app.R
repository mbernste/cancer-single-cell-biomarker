library(shiny)
library(ggplot2)
library(rhdf5)
library(viridis)
library(heatmaply)
library(dplyr)
library(shinycssloaders)
library(RColorBrewer)

#TUMORS = c('PJ016', 'PJ018', 'PJ025', 'PJ048', 'PJ030', 'PJ035', 'PJ017', 'PJ032')
TUMORS = c('PJ048')
DATA_DIR = '/Users/matthewbernstein/Development/single-cell-hackathon/data'
TMP_DIR = '/Users/matthewbernstein/Development/single-cell-hackathon/tmp'

print('Loading counts data...')
counts = h5read(paste0(DATA_DIR, "/GSE103224_normalized.h5"),"count")
cells = h5read(paste0(DATA_DIR, "/GSE103224_normalized.h5"), "cell")
tumor_ids = h5read(paste0(DATA_DIR, "/GSE103224_normalized.h5"), "tumor")
genes = h5read(paste0(DATA_DIR, "/GSE103224_normalized.h5"), "gene_name")
counts <- t(counts)
counts <- data.frame(counts)
rownames(counts) <- cells
colnames(counts) <- genes
print('done.')

# Load the GSEA data
gsea_df <- read.csv(
    paste0(TMP_DIR, "/gsea_results.tsv"),
    sep = '\t',
    header = TRUE
)

# Load tumor-cluster 'expression' data
tumor_clust_expr_df <- read.csv(
    paste0(TMP_DIR, "/tumor_cluster_gene_expression.tsv"),
    sep = '\t',
    header = TRUE
)

# Load each tumor's expression data and dimension reduction data
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
    
    # Read dimension reduction
    df <- read.csv(
        file=paste0(TMP_DIR, '/', tumor, '_PHATE_3.tsv'), 
        sep = '\t', 
        header = FALSE
    )
    colnames(df) <- c('PHATE 1', 'PHATE 2', 'PHATE 3')
    rownames(df) <- tumor_cells
    df <- data.frame(df)
    print('Merging...')
    df <- merge(df, tumor_counts,  by=0)
    print('done.')
    tumor_dfs[[tumor]] <- df
}

# Load the UMAP coordinates for the aligned tumors
aligned_umap_df <- read.csv(
    paste0(TMP_DIR, '/all_tumors_aligned_PHATE_3.tsv'),
    sep = '\t',
    header = TRUE,
    row.names = 1
)

server <- function(input, output) {
    
    output$plot1 <- renderPlotly({plot_ly(
        tumor_dfs[[input$tumor1]], 
        x = tumor_dfs[[input$tumor1]]$PHATE.1, 
        y = tumor_dfs[[input$tumor1]]$PHATE.2, 
        z = tumor_dfs[[input$tumor1]]$PHATE.3, 
        height = 700,
        marker = list(color = tumor_dfs[[input$tumor1]][[input$colorby1]], colorscale = 'Viridis', showscale = TRUE, size = 5)
        ) %>% add_markers() %>% layout(
            title = "\nPlot 1",
            titlefont = list(family = "arial", size = 25),
            scene = list(xaxis = list(title = 'PHATE 1'),
            yaxis = list(title = 'PHATE 2'),
            zaxis = list(title = 'PHATE 3')))
    })
    #output$plot1<-renderPlot({
    #    ggplot(tumor_dfs[[input$tumor1]], aes(x=PHATE.1,y=PHATE.2, color=get(input$colorby1))) + 
    #    scale_color_viridis() + 
    #    geom_point() + 
    #    theme_classic()
    #})
    output$plot2 <- renderPlotly({plot_ly(
        tumor_dfs[[input$tumor2]],
        x = tumor_dfs[[input$tumor2]]$PHATE.1,
        y = tumor_dfs[[input$tumor2]]$PHATE.2,
        z = tumor_dfs[[input$tumor2]]$PHATE.3,
        height = 700,
        marker = list(color = tumor_dfs[[input$tumor2]][[input$colorby2]], colorscale = 'Viridis', showscale = TRUE, size = 5)
        ) %>% add_markers() %>% layout(
            title = "Plot 2",
            scene = list(xaxis = list(title = 'PHATE 1'),
                     yaxis = list(title = 'PHATE 2'),
                     zaxis = list(title = 'PHATE 3')))
    })
    #output$plot2<-renderPlot({
    #    ggplot(tumor_dfs[[input$tumor2]], aes(x=PHATE.1,y=PHATE.2, color=get(input$colorby2))) + 
    #    scale_color_viridis() + 
    #    geom_point() + 
    #    theme_classic()
    #})
    output$plot3 <- renderPlotly({plot_ly(
        tumor_dfs[[input$tumor3]],
        x = tumor_dfs[[input$tumor3]]$PHATE.1,
        y = tumor_dfs[[input$tumor3]]$PHATE.2,
        z = tumor_dfs[[input$tumor3]]$PHATE.3,
        height = 700,
        marker = list(color = tumor_dfs[[input$tumor3]][[input$colorby3]], colorscale = 'Viridis', showscale = TRUE, size = 5)
        ) %>% add_markers() %>% layout(scene = list(xaxis = list(title = 'PHATE 1'),
                     yaxis = list(title = 'PHATE 2'),
                     zaxis = list(title = 'PHATE 3')))
    })
    #output$plot3<-renderPlot({
    #    ggplot(tumor_dfs[[input$tumor3]], aes(x=PHATE.1,y=PHATE.2, color=get(input$colorby3))) + 
    #    scale_color_viridis() + 
    #    geom_point() + 
    #    theme_classic()
    #})
    output$plot4 <- renderPlotly({plot_ly(
        tumor_dfs[[input$tumor4]],
        x = tumor_dfs[[input$tumor4]]$PHATE.1,
        y = tumor_dfs[[input$tumor4]]$PHATE.2,
        z = tumor_dfs[[input$tumor4]]$PHATE.3,
        height = 700,
        marker = list(color = tumor_dfs[[input$tumor4]][[input$colorby4]], colorscale = 'Viridis', showscale = TRUE, size = 5)
        ) %>% add_markers() %>% layout(scene = list(xaxis = list(title = 'PHATE 1'),
                     yaxis = list(title = 'PHATE 2'),
                     zaxis = list(title = 'PHATE 3')))
    })
    #output$plot4<-renderPlot({
    #    ggplot(tumor_dfs[[input$tumor4]], aes(x=PHATE.1,y=PHATE.2, color=get(input$colorby4))) + 
    #    scale_color_viridis() + 
    #    geom_point() + 
    #    theme_classic()
    #})
    
    # GSEA
    output$heatmap <- renderPlotly({
        tumor_clusts <- select(filter(tumor_clust_expr_df, tumor_clust_expr_df[input$gsea1] > 0), c('tumor_cluster'))    
        all_go_terms = select(filter(gsea_df, gsea_df$tumor_cluster == tumor_clusts[[1]]), c('GO_term'))
        for (tumor_clust in tumor_clusts[2:length(rownames(tumor_clusts)),]) {
          h <- select(filter(gsea_df, gsea_df$tumor_cluster == tumor_clust), c('GO_term'))
          all_go_terms <- rbind(
            all_go_terms,
            select(
              filter(
                gsea_df, gsea_df$tumor_cluster == tumor_clust
              ), 
              c('GO_term')
            )
          )
        }
        all_go_terms <- unique(all_go_terms["GO_term"])

        matrix = c()
        print(dim(all_go_terms['GO_term']))
        for (go_term in all_go_terms[['GO_term']][1:200]) {
          row <- c()
          for (tumor_clust in tumor_clusts[['tumor_cluster']]) {
            curr_terms <- select(filter(gsea_df, gsea_df$tumor_cluster == tumor_clust), c('GO_term'))['GO_term']
            if(go_term %in% curr_terms[['GO_term']]) {
              row <- append(row, 1.0)
            }
            else{
              row <- append(row, 0.0)
            }
          }
          print(row)
          matrix <- append(matrix, row)
        } 

        heatmap <- matrix(matrix, ncol = length(tumor_clusts[['tumor_cluster']]), nrow = length(all_go_terms[['GO_term']][1:200]))
        rownames(heatmap) <- all_go_terms[['GO_term']][1:200]
        colnames(heatmap) <- tumor_clusts[['tumor_cluster']]
        heatmaply(
            heatmap, 
            showticklabels = c(TRUE, FALSE), 
            colors = gray.colors(100, rev = TRUE) 
        ) %>% layout(height = 2000, width=2000)
    })

    # Tumor integration
    output$aligned_plot1 <- renderPlotly({
        colby <- input$aligned_colorby1
        if (colby != 'tumor') {
            gene_df <- counts[toString(input$aligned_colorby1)]
            rownames(gene_df) <- cells
            curr_df <- merge(aligned_umap_df, gene_df, by=0)
            plot_ly(
                curr_df,
                x = curr_df$PHATE1,
                y = curr_df$PHATE2,
                z = curr_df$PHATE3,
                height = 700,
                marker = list(color = curr_df[[colby]], colorscale = 'Viridis', showscale = TRUE, size = 2)
                ) %>% add_markers() %>% layout(scene = list(xaxis = list(title = 'PHATE 1'),
                             yaxis = list(title = 'PHATE 2'),
                             zaxis = list(title = 'PHATE 3'))
            )
        }
        else {
            curr_df <- aligned_umap_df
            plot_ly(
                curr_df,
                x = curr_df$PHATE1,
                y = curr_df$PHATE2,
                z = curr_df$PHATE3,
                color = curr_df[['tumor']],
                colors = brewer.pal(length(unique(curr_df[['tumor']])), "Set1"),
                height = 700,
                marker = list(size = 2)
                ) %>% add_markers() %>% layout(scene = list(xaxis = list(title = 'PHATE 1'),
                             yaxis = list(title = 'PHATE 2'),
                             zaxis = list(title = 'PHATE 3'))
            )
        }
    })
    output$aligned_plot2 <- renderPlotly({
        colby <- input$aligned_colorby2
        if (colby != 'tumor') {
            gene_df <- counts[toString(input$aligned_colorby2)]
            rownames(gene_df) <- cells
            curr_df <- merge(aligned_umap_df, gene_df, by=0)
            plot_ly(
                curr_df,
                x = curr_df$PHATE1,
                y = curr_df$PHATE2,
                z = curr_df$PHATE3,
                height = 700,
                marker = list(color = curr_df[[colby]], colorscale = 'Viridis', showscale = TRUE, size = 2)
                ) %>% add_markers() %>% layout(scene = list(xaxis = list(title = 'PHATE 1'),
                             yaxis = list(title = 'PHATE 2'),
                             zaxis = list(title = 'PHATE 3'))
            )
        }
        else {
            curr_df <- aligned_umap_df
            plot_ly(
                curr_df,
                x = curr_df$PHATE1,
                y = curr_df$PHATE2,
                z = curr_df$PHATE3,
                color = curr_df[['tumor']],
                colors = brewer.pal(length(unique(curr_df[['tumor']])), "Set1"),
                height = 700,
                marker = list(size = 2) 
                ) %>% add_markers() %>% layout(scene = list(xaxis = list(title = 'PHATE 1'),
                             yaxis = list(title = 'PHATE 2'),
                             zaxis = list(title = 'PHATE 3'))
            )
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
                    selectInput("tumor1", "Plot 1 dataset:", TUMORS),
                    textInput("colorby1", "Select gene:", value = "OLIG1"), 
                    selectInput("tumor2", "Plot 2 dataset:", TUMORS),
                    textInput("colorby2", "Select gene:", value = "STMN2"), 
                    selectInput("tumor3", "Plot 3 dataset:", TUMORS),
                    textInput("colorby3", "Select gene:", value = "MOG"),
                    selectInput("tumor4", "Plot 4 dataset:", TUMORS),
                    textInput("colorby4", "Select gene:", value = "GFAP")
                ),               
                mainPanel(
                    "",
                    fluidRow(
                        splitLayout(
                            cellWidths = c("50%", "50%"), 
                            plotlyOutput("plot1", height = "100%") %>% withSpinner(color="#000000"), 
                            plotlyOutput("plot2", height = "100%") %>% withSpinner(color="#000000"), 
                            style = "height:700px;"
                        )
                    ),
                    fluidRow(
                        splitLayout(
                            cellWidths = c("50%", "50%"),
                            plotlyOutput("plot3", height = "100%") %>% withSpinner(color="#000000"),
                            plotlyOutput("plot4", height = "100%") %>% withSpinner(color="#000000"),
                            style = "height:700px;" 
                        )
                    )
                )
            )
        ),
        tabPanel("Cluster GSEA",
            sidebarLayout(
                sidebarPanel(
                    textInput("gsea1", "Select gene:"),
                    width = 2
                ),
                mainPanel("",
                    plotlyOutput("heatmap", width="100%") %>% withSpinner(color="#000000"),
                    width = 6,
                    style = "height:700px;"
                ),
            ),
            class = 'leftAlign'
            #heatmaply(mtcars, dendrogram = "column", showticklabels = c(TRUE, FALSE)) 
        ),
        tabPanel("Tumors integrated",
           sidebarLayout(
                sidebarPanel(
                    textInput("aligned_colorby1", "Plot 1, color by:", value = "tumor"),
                    textInput("aligned_colorby2", "Plot 2, color by:", value = "GFAP"),
                    width=2
                ),
                mainPanel(
                    fluidRow(
                        splitLayout(
                            cellWidths = c("50%", "50%"), 
                            plotlyOutput("aligned_plot1", height = "100%") %>% withSpinner(color="#000000"), 
                            plotlyOutput("aligned_plot2", height = "100%") %>% withSpinner(color="#000000")
                        )
                    )
                )
            )
        )
  ),
  theme = "bootstrap.css"
)

shinyApp(ui = ui, server = server)
