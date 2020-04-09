

library(shiny)
library(ggplot2)
library(rhdf5)
library(viridis)
library(heatmaply)
library(dplyr)
library(shinycssloaders)
library(RColorBrewer)
library(rjson)
library(DT)
#library(bigmemory)

##########################   Configuration ########################################
#DATA_DIR = '/Users/matthewbernstein/Development/single-cell-hackathon/tmp_test'
###################################################################################

options(bigmemory.allow.dimnames=TRUE)

COORDS = c('UMAP', 'PHATE')
#DB_FILENAME = 'GSE103224_GSE72056_normalized.h5'
DB_FILENAME = 'GSE103224_GSE123904_normalized.h5'

# Load meatadata
meta_df <- fromJSON(file = paste0(DATA_DIR, '/tumor_metadata.json'))
TUMORS <- meta_df[['datasets']]

de_genes_json <- paste0(DATA_DIR, "/cluster/cluster_de_genes.json")
print(de_genes_json)
de_data <- fromJSON(file = de_genes_json)

# Load the GSEA data
gsea_df <- read.csv(
    paste0(DATA_DIR, "/cluster/gsea_results.tsv"),
    sep = '\t',
    header = TRUE
)   
gsea_full_df <- read.csv(
    paste0(DATA_DIR, "/cluster/gsea_binary_matrix.tsv"),
    sep = '\t',
    header = TRUE,
    row.names = 1
)

# Load tumor-cluster 'expression' data
tumor_clust_expr_df <- read.csv(
    paste0(DATA_DIR, "/cluster/tumor_cluster_gene_expression.tsv"),
    sep = '\t',
    header = TRUE
)
    

print('Loading counts data...')
counts = h5read(paste0(DATA_DIR, '/', DB_FILENAME),"count")
cells = h5read(paste0(DATA_DIR, '/', DB_FILENAME), "cell")
tumor_ids = h5read(paste0(DATA_DIR, '/', DB_FILENAME), "dataset")
genes = h5read(paste0(DATA_DIR, '/', DB_FILENAME), "gene_name")
counts <- data.frame(t(counts))
#counts <- as.big.matrix(counts)
#gc()
rownames(counts) <- cells
colnames(counts) <- genes
print('done.')


# Load each tumor's expression data and dimension reduction data
tumor_dfs <- c()
for (tumor in TUMORS) {
    print(paste('Creating dataset for tumor', tumor))
    # Load data matrix
    tumor_indices <- c()
    for (i in 1:length(cells)) {
        if (tumor_ids[[i]] == tumor) {
            tumor_indices <- append(tumor_indices, i)
        }
    }
    tumor_cells <- cells[tumor_indices]
    #tumor_counts <- as.big.matrix(counts[tumor_indices,])
    tumor_counts <- counts[tumor_indices,]
    rownames(tumor_counts) <- tumor_cells
    colnames(tumor_counts) <- genes
    print('done.')    

    # Read dimension reduction
    print(paste('Loading PHATE data for tumor', tumor))
    phate_df <- read.csv(
        file=paste0(DATA_DIR, '/dim_reduc/', tumor, '_PHATE_3.tsv'), 
        sep = '\t', 
        header = TRUE,
        row.names = 1
    )
    colnames(phate_df) <- c('PHATE1', 'PHATE2', 'PHATE3')
    umap_df <- read.csv(
        file=paste0(DATA_DIR, '/dim_reduc/', tumor, '_UMAP_3.tsv'),
        sep = '\t',
        header = TRUE,
        row.names = 1
    )
    colnames(umap_df) <- c('UMAP1', 'UMAP2', 'UMAP3')
    print(length(tumor_cells))
    print(nrow(phate_df))
    rownames(phate_df) <- tumor_cells
    print('done.')

    clust_df <- read.csv(
        file=paste0(
            DATA_DIR, 
            '/cluster/', 
            tumor, 
            '_clusters.res_0_8.tsv'
        ),
        sep = '\t',
        header = TRUE
    )
    colnames(clust_df) <- c('cell', 'cluster')
    rownames(clust_df) <- tumor_cells

    phate_df <- data.frame(phate_df)
    #print('Merging...')
    #df <- merge(df, tumor_counts,  by=0)
    
    tumor_counts$PHATE1 <- phate_df$PHATE1
    tumor_counts$PHATE2 <- phate_df$PHATE2
    tumor_counts$PHATE3 <- phate_df$PHATE3
    tumor_counts$UMAP1 <- umap_df$UMAP1
    tumor_counts$UMAP2 <- umap_df$UMAP2
    tumor_counts$UMAP3 <- umap_df$UMAP3
    tumor_counts$cluster <- as.character(clust_df$cluster)

    #print('done.')
    tumor_dfs[[tumor]] <- tumor_counts
}

##################################################################
# Load the PHATE coordinates for the aligned tumors
##################################################################

aligned_dfs <- list()
aligned_counts_df <- list()
for (integ_set_name in names(meta_df$integration_sets)) {
    print(paste('Loading integrated data for integration set ', integ_set_name)) 
    phate_df <- read.csv(
        paste0(DATA_DIR, '/dim_reduc/', integ_set_name, '_aligned_PHATE_3.tsv'),
        sep = '\t',
        header = TRUE,
        row.names = 1
    )
    umap_df <- read.csv(    
        paste0(DATA_DIR, '/dim_reduc/', integ_set_name, '_aligned_UMAP_3.tsv'),
        sep = '\t',
        header = TRUE,
        row.names = 1
    )
    aligned_dfs[[integ_set_name]] <- transform(merge(phate_df, umap_df, by=0), row.names=Row.names, Row.names=NULL)
    aligned_dfs[[integ_set_name]] <- transform(merge(aligned_dfs[[integ_set_name]], counts, by=0), row.names=Row.names, Row.names=NULL)
    names(aligned_dfs[[integ_set_name]])[names(aligned_dfs[[integ_set_name]]) == "tumor.x"] <- "tumor"
    names(aligned_dfs[[integ_set_name]])[names(aligned_dfs[[integ_set_name]]) == "subtype.x"] <- "subtype" 
    aligned_counts_df[[integ_set_name]] <- counts[rownames(aligned_dfs[[integ_set_name]]),]
}

#aligned_dfs[['All']] <- all_aligned_df
#for (tumor_i in names(tumor_dfs)) {
#    for (tumor_j in names(tumor_dfs)) {
#        if (tumor_i == tumor_j) {
#            break
#        }
#        print(paste(
#            "Loading alignment data for tumor", 
#            tumor_i, 
#            "and", 
#            tumor_j
#        ))
#        a_df <- read.csv(
#            paste0(
#                DATA_DIR, 
#                '/pairwise_integrations_PHATE/', 
#                tumor_i, 
#                '_', 
#                tumor_j, 
#                '_aligned_PHATE_3.tsv'
#            ),
#            sep = '\t',
#            header = TRUE,
#            row.names = 1
#        )
#        aligned_dfs[[paste0(tumor_i,'_', tumor_j)]] <- a_df
#    }
#}

server <- function(input, output) {
    output$plot1 <- renderPlotly({
        if (input$coords1 == 'PHATE') {
            px <- tumor_dfs[[input$tumor1]]$PHATE1
            py <- tumor_dfs[[input$tumor1]]$PHATE2
            pz <- tumor_dfs[[input$tumor1]]$PHATE3
        }
        else if (input$coords1 == 'UMAP') {
            px <- tumor_dfs[[input$tumor1]]$UMAP1
            py <- tumor_dfs[[input$tumor1]]$UMAP2
            pz <- tumor_dfs[[input$tumor1]]$UMAP3            
        }
        if (input$colorby1 != 'cluster') {
            plot_ly(
                tumor_dfs[[input$tumor1]], 
                x = px, 
                y = py, 
                z = pz, 
                height = 700,
                marker = list(
                    color = tumor_dfs[[input$tumor1]][[input$colorby1]], 
                    colorscale = 'Viridis', 
                    showscale = TRUE, 
                    size = input$size_plot1
                )
            ) %>% 
            add_markers() %>% 
            layout(
                title = "\nPlot 1",
                titlefont = list(family = "arial", size = 25),
                scene = list(xaxis = list(title = 'PHATE 1')),
                yaxis = list(title = 'PHATE 2'),
                zaxis = list(title = 'PHATE 3')
            )
        }
        else {
            plot_ly(
                tumor_dfs[[input$tumor1]],
                x = px,
                y = py,
                z = pz,
                color = tumor_dfs[[input$tumor1]][['cluster']],
                colors = brewer.pal(length(unique(tumor_dfs[[input$tumor1]][['cluster']])), "Set1"),
                height = 700,
                marker = list(size = input$size_plot1)
            ) %>% 
            add_markers() %>% 
            layout(
                title = "\nPlot 1",
                titlefont = list(family = "arial", size = 25),
                scene = list(xaxis = list(title = 'PHATE 1'),
                yaxis = list(title = 'PHATE 2'),
                zaxis = list(title = 'PHATE 3'))
            )
        }
    })
    output$plot2 <- renderPlotly({
        units <- input$coords2
        if (units == 'PHATE') {
            px <- tumor_dfs[[input$tumor2]]$PHATE1
            py <- tumor_dfs[[input$tumor2]]$PHATE2
            pz <- tumor_dfs[[input$tumor2]]$PHATE3
        }
        else if (units == 'UMAP') {
            px <- tumor_dfs[[input$tumor2]]$UMAP1
            py <- tumor_dfs[[input$tumor2]]$UMAP2
            pz <- tumor_dfs[[input$tumor2]]$UMAP3
        }
        if (input$colorby2 != 'cluster') {
            plot_ly(
                tumor_dfs[[input$tumor2]],
                x = px,
                y = py,
                z = pz,
                height = 700,
                marker = list(
                    color = tumor_dfs[[input$tumor2]][[input$colorby2]], 
                    colorscale = 'Viridis', 
                    showscale = TRUE, 
                    size = input$size_plot2
                )
            ) %>% 
            add_markers() %>% 
            layout(
                title = "\nPlot 2",
                titlefont = list(family = "arial", size = 25),
                scene = list(xaxis = list(title = paste(units, '1'))),
                yaxis = list(title = paste(units, '2')),
                zaxis = list(title = paste(units, '3'))
            )
        }
        else {
            plot_ly(
                tumor_dfs[[input$tumor2]],
                x = px,
                y = py,
                z = pz,
                color = tumor_dfs[[input$tumor2]][['cluster']],
                colors = brewer.pal(length(unique(tumor_dfs[[input$tumor2]][['cluster']])), "Set1"),
                height = 700,
                marker = list(size = input$size_plot2)
            ) %>% 
            add_markers() %>% 
            layout(
                title = "\nPlot 2",
                titlefont = list(family = "arial", size = 25),
                scene = list(
                    xaxis = list(title = paste(units, '1')),
                    yaxis = list(title = paste(units, '2')),
                    zaxis = list(title = paste(units, '3'))
                )
            )
        }

    })
    output$plot3 <- renderPlotly({
        if (input$coords3 == 'PHATE') {
            px <- tumor_dfs[[input$tumor3]]$PHATE1
            py <- tumor_dfs[[input$tumor3]]$PHATE2
            pz <- tumor_dfs[[input$tumor3]]$PHATE3
        }
        else if (input$coords3 == 'UMAP') {
            px <- tumor_dfs[[input$tumor3]]$UMAP1
            py <- tumor_dfs[[input$tumor3]]$UMAP2
            pz <- tumor_dfs[[input$tumor3]]$UMAP3
        }
        if (input$colorby3 != 'cluster') {
            plot_ly(
                tumor_dfs[[input$tumor3]],
                x = px,
                y = py,
                z = pz,
                height = 700,
                marker = list(
                    color = tumor_dfs[[input$tumor3]][[input$colorby3]], 
                    colorscale = 'Viridis', 
                    showscale = TRUE, 
                    size = input$size_plot3
                )
            ) %>% 
            add_markers() %>% 
            layout(
                title = "\nPlot 3",
                titlefont = list(family = "arial", size = 25),
                scene = list(
                xaxis = list(title = 'PHATE 1')),
                yaxis = list(title = 'PHATE 2'),
                zaxis = list(title = 'PHATE 3')
            )
        }
        else {
            plot_ly(
                tumor_dfs[[input$tumor3]],
                x = px,
                y = py,
                z = pz,
                color = tumor_dfs[[input$tumor3]][['cluster']],
                colors = brewer.pal(length(unique(tumor_dfs[[input$tumor3]][['cluster']])), "Set1"),
                height = 700,
                marker = list(size = input$size_plot3)
            ) %>% 
            add_markers() %>% 
            layout(
                title = "\nPlot 3",
                titlefont = list(family = "arial", size = 25),
                scene = list(
                    xaxis = list(title = 'PHATE 1'),
                    yaxis = list(title = 'PHATE 2'),
                    zaxis = list(title = 'PHATE 3')
                )
            )
        }

    })
    output$plot4 <- renderPlotly({
        units <- input$coords4
        if (units == 'PHATE') {
            px <- tumor_dfs[[input$tumor4]]$PHATE1
            py <- tumor_dfs[[input$tumor4]]$PHATE2
            pz <- tumor_dfs[[input$tumor4]]$PHATE3
        }
        else if (units == 'UMAP') {
            px <- tumor_dfs[[input$tumor4]]$UMAP1
            py <- tumor_dfs[[input$tumor4]]$UMAP2
            pz <- tumor_dfs[[input$tumor4]]$UMAP3
        }
        if (input$colorby4 != 'cluster') {
            plot_ly(
                tumor_dfs[[input$tumor4]],
                x = px,
                y = py,
                z = pz,
                height = 700,
                marker = list(
                    opacity = 1.0,
                    color = tumor_dfs[[input$tumor4]][[input$colorby4]], 
                    colorscale = 'Viridis', 
                    showscale = TRUE, 
                    size = input$size_plot4
                )
            ) %>% 
            add_markers() %>% 
            layout(
                title = "\nPlot 4",
                titlefont = list(family = "arial", size = 25),
                scene = list(
                    xaxis = list(title = 'PHATE 1'),
                    yaxis = list(title = 'PHATE 2'),
                    zaxis = list(title = 'PHATE 3')
                )
            )
        }
        else {
            plot_ly(
                tumor_dfs[[input$tumor4]],
                x = px,
                y = py,
                z = pz,
                color = tumor_dfs[[input$tumor4]][['cluster']],
                colors = brewer.pal(length(unique(tumor_dfs[[input$tumor4]][['cluster']])), "Set1"),
                height = 700,
                marker = list(size = input$size_plot4)
            ) %>% 
            add_markers() %>% 
            layout(
                title = "\nPlot 4",
                titlefont = list(family = "arial", size = 25),
                scene = list(
                    xaxis = list(title = 'PHATE 1'),
                    yaxis = list(title = 'PHATE 2'),
                    zaxis = list(title = 'PHATE 3')
                )
            )
        }
    })
    
    # GSEA
    output$heatmap <- renderPlotly({
        tumor_clusts <- select(
            filter(
                tumor_clust_expr_df, 
                tumor_clust_expr_df[input$gsea1] > 0
            ), 
            c('tumor_cluster')
        )
        select_cols <- c()
        for (tc in colnames(gsea_full_df)) {
            if(tc %in% tumor_clusts[['tumor_cluster']]) {
                select_cols <- c(select_cols, 1.0)
            }
            else {
                select_cols <- c(select_cols, 0.0)
            }
        }

        # String manipulation to get the tumor names from 
        # the tumor-cluster strings
        splits <- strsplit(colnames(gsea_full_df), "_")
        tumor_annot <- c()
        for (split in splits) {
            tumor_annot <- c(tumor_annot, paste(split[1:length(split)-1], collapse = '_'))
        }   
        row_annot <- data.frame(
            tumor = tumor_annot,
            gene_is_DE = select_cols
        )
        heatmaply(
            gsea_full_df,
            col_side_colors = row_annot, 
            showticklabels = c(TRUE, FALSE),
            margins=c(0,100), 
            colors = gray.colors(100, rev = TRUE) 
        ) %>% layout(height = 1000, width=700)
    })

    # Tumor integration
    output$aligned_plot1 <- renderPlotly({
        colby <- input$aligned_colorby1
        aligned_df <- aligned_dfs[[input$aligned1]]
        if (input$coords_aligned1 == 'PHATE') {
            px <- aligned_df$PHATE1
            py <- aligned_df$PHATE2
            pz <- aligned_df$PHATE3
        }
        else if (input$coords_aligned1 == 'UMAP') {
            px <- aligned_df$UMAP1
            py <- aligned_df$UMAP2
            pz <- aligned_df$UMAP3
        }
        if (colby == 'tumor' | colby == 'subtype' | colby == 'cluster') {
            curr_df <- aligned_df
            plot_ly(
                curr_df,
                x = px,
                y = py,
                z = pz,
                color = curr_df[[colby]],
                colors = brewer.pal(length(unique(curr_df[[colby]])), "Set1"),
                #colors = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999"),
                height = 700,
                marker = list(size = input$aligned_size_plot1)
                ) %>% add_markers() %>% layout(
                    scene = list(xaxis = list(title = 'PHATE 1'),
                    yaxis = list(title = 'PHATE 2'),
                    zaxis = list(title = 'PHATE 3')
                )
            )
        }
        else {
            aligned_df <- aligned_dfs[[input$aligned1]]
            #gene_df <- counts[toString(input$aligned_colorby1)]
            #rownames(gene_df) <- cells
            #curr_df <- merge(aligned_df, gene_df, by=0)
            plot_ly(
                aligned_df,
                x = px,
                y = py,
                z = pz,
                height = 700,
                marker = list(
                    #color = curr_df[[colby]],
                    color = aligned_dfs[[input$aligned1]][[colby]], 
                    colorscale = 'Viridis', 
                    showscale = TRUE, 
                    size = input$aligned_size_plot1
                )
            ) %>% add_markers() %>% layout(
                scene = list(
                    xaxis = list(title = 'PHATE 1'),
                    yaxis = list(title = 'PHATE 2'),
                    zaxis = list(title = 'PHATE 3')
                )
            )
        }
    })
    output$aligned_plot2 <- renderPlotly({
        colby <- input$aligned_colorby2
        aligned_df <- aligned_dfs[[input$aligned2]]
        if (input$coords_aligned2 == 'PHATE') {
            px <- aligned_df$PHATE1
            py <- aligned_df$PHATE2
            pz <- aligned_df$PHATE3
        }
        else if (input$coords_aligned2 == 'UMAP') {
            px <- aligned_df$UMAP1
            py <- aligned_df$UMAP2
            pz <- aligned_df$UMAP3
        }
        if (colby == 'tumor' | colby == 'subtype' | colby == 'cluster') {
            curr_df <- aligned_df
            plot_ly(
                curr_df,
                x = curr_df$PHATE1,
                y = curr_df$PHATE2,
                z = curr_df$PHATE3,
                color = curr_df[[colby]],
                colors = brewer.pal(length(unique(curr_df[[colby]])), "Set1"),
                height = 700,
                marker = list(size = input$aligned_size_plot2)
                ) %>% add_markers() %>% layout(
                    scene = list(xaxis = list(title = 'PHATE 1'),
                    yaxis = list(title = 'PHATE 2'),
                    zaxis = list(title = 'PHATE 3')
                )
            )
        }
        else{
            aligned_df <- aligned_dfs[[input$aligned2]]
            #gene_df <- counts[toString(input$aligned_colorby2)]
            #rownames(gene_df) <- cells
            #curr_df <- merge(aligned_df, gene_df, by=0)
            plot_ly(
                aligned_df,
                x = px,
                y = py,
                z = pz,
                height = 700,
                marker = list(
                    #color = curr_df[[colby]],
                    color = aligned_dfs[[input$aligned2]][[colby]],
                    colorscale = 'Viridis', 
                    showscale = TRUE, 
                    size = input$aligned_size_plot2
                )
            ) %>% add_markers() %>% layout(
                scene = list(xaxis = list(title = 'PHATE 1'),
                    yaxis = list(title = 'PHATE 2'),
                    zaxis = list(title = 'PHATE 3')
                )
            )
        }
    })


    # DE tables
    output$de_table1 <- renderDataTable({
        table <- de_data[[input$de_table_tumor1]][[input$de_table_clust1]]
        table <- data.frame(table)
        colnames(table) <- c('DE Genes')
        datatable(table)
    })
    output$de_table2 <- renderDataTable({
        table <- de_data[[input$de_table_tumor2]][[input$de_table_clust2]]
        table <- data.frame(table)
        colnames(table) <- c('DE Genes')
        datatable(table)
    })
}   

ui <- fluidPage(

  # App title ----
  titlePanel(title=div(img(src="charts_logo.png", height=100, width=100), "CHARTS: CHARacterizing Tumor Subpopulations")),

  tabsetPanel(
        tabPanel("Individual tumors", 
            sidebarLayout(
                sidebarPanel(
                    selectInput("tumor1", "Plot 1 dataset:", TUMORS),
                    #textInput("colorby1", "Select gene:", value = "OLIG1"),
                    selectInput("colorby1", "Select gene or variable:", c('cluster', genes)),
                    selectInput("coords1", "Select coordinates:", COORDS), 
                    sliderInput("size_plot1", "Dot size:",
                        min = 1, max = 10,
                        value = 3
                    ),
                    selectInput("tumor2", "Plot 2 dataset:", TUMORS),
                    #textInput("colorby2", "Select gene:", value = "STMN2"),
                    selectInput("colorby2", "Select gene or variable:", c('cluster', genes)),
                    selectInput("coords2", "Select coordinates:", COORDS),
                    sliderInput("size_plot2", "Dot size:",
                        min = 1, max = 10,
                        value = 3
                    ), 
                    selectInput("tumor3", "Plot 3 dataset:", TUMORS),
                    #textInput("colorby3", "Select gene:", value = "MOG"),
                    selectInput("colorby3", "Select gene or variable:", c('cluster', genes)),
                    selectInput("coords3", "Select coordinates:", COORDS),
                    sliderInput("size_plot3", "Dot size:",
                        min = 1, max = 10,
                        value = 3
                    ),
                    selectInput("tumor4", "Plot 4 dataset:", TUMORS),
                    #textInput("colorby4", "Select gene:", value = "GFAP"),
                    selectInput("colorby4", "Select gene or variable:", c('cluster', genes)),
                    selectInput("coords4", "Select coordinates:", COORDS),
                    sliderInput("size_plot4", "Dot size:",
                        min = 1, max = 10,
                        value = 3
                    ),
                    width=2
                ),               
                mainPanel(
                    "",
                    fluidRow(
                        splitLayout(
                            cellWidths = c("50%", "50%"), 
                            plotlyOutput("plot1", height = "100%") %>% withSpinner(color="#000000"), 
                            plotlyOutput("plot2", height = "100%") %>% withSpinner(color="#000000"), 
                            style = "height:600px;"
                        )
                    ),
                    fluidRow(
                        splitLayout(
                            cellWidths = c("50%", "50%"),
                            plotlyOutput("plot3", height = "100%") %>% withSpinner(color="#000000"),
                            plotlyOutput("plot4", height = "100%") %>% withSpinner(color="#000000"),
                            style = "height:600px;" 
                        )
                    )
                )
            )
        ),
        tabPanel("Integrated Tumors",
           sidebarLayout(
                sidebarPanel(
                    selectInput("aligned1", "Plot 1 dataset:", names(aligned_dfs)),
                    selectInput("aligned_colorby1", "Plot 1, color by:", c('tumor', 'subtype', genes)),
                    selectInput("coords_aligned1", "Select coordinates:", COORDS),
                    sliderInput("aligned_size_plot1", "Dot size:",
                        min = 1, max = 10,
                        value = 2
                    ),
                    selectInput("aligned2", "Plot 2 dataset:", names(aligned_dfs)),
                    selectInput("aligned_colorby2", "Plot 2, color by:", c('tumor', 'subtype', genes)),
                    selectInput("coords_aligned2", "Select coordinates:", COORDS),
                    sliderInput("aligned_size_plot2", "Dot size:",
                        min = 1, max = 10,
                        value = 2
                    ),
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
        ),
        tabPanel("Cluster GSEA",
            sidebarLayout(
                sidebarPanel(
                    textInput("gsea1", "Select gene:", value = "SLC16A3"),
                    width = 2
                ),
                mainPanel("",
                    plotlyOutput("heatmap", width="100%") %>% withSpinner(color="#000000"),
                    width = 6,
                    style = "height:700px;"
                ),
            ),
            class = 'leftAlign'
        ),
        tabPanel("Tumor-Cluster DE",
            sidebarLayout(
                sidebarPanel(
                    textInput("de_table_tumor1", "Table 1, select tumor:", value = 'PJ016'),
                    textInput("de_table_clust1", "Table 1, select cluster:", value = '0'),
                    textInput("de_table_tumor2", "Table 2, select tumor:", value = 'PJ017'),
                    textInput("de_table_clust2", "Table 2, select cluster:", value = '0'),
                    width = 2
                ),
                mainPanel("",
                    fluidRow(
                        splitLayout(
                            cellWidths = c("50%", "50%"),
                            h3('Table 1'),
                            h3('Table 2')
                        )
                    ),
                    fluidRow(
                        splitLayout(
                            cellWidths = c("50%", "50%"),
                            dataTableOutput("de_table1") %>% withSpinner(color="#000000"),
                            dataTableOutput("de_table2") %>% withSpinner(color="#000000")
                        )
                    ),
                    width = 6,
                    style = "height:700px;"
                )
            )
        )
  ),
  theme = "bootstrap.css"
)

shinyApp(ui = ui, server = server)
