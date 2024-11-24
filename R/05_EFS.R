barplot_efs <- function(efs_table, filter_by = 0.5, order = TRUE) {
  # Filter and order features
  if (order) {
    efs_table <- efs_table[, colSums(efs_table) > filter_by]
    efs_table <- efs_table[, order(colSums(efs_table))]
  }
  
  # Get the number of features
  num_features <- ncol(efs_table)
  
  # Check the number of features to adjust plot behavior
  if (num_features > 100) {
    # Create a barplot for large datasets
    bar_heights <- colSums(efs_table)
    barplot(
      bar_heights,
      ylim = c(0, 1),
      main = 'Ensemble Feature Selection',
      xlab = "Features",
      ylab = "Importance values",
      axisnames = FALSE
    )
    return(invisible()) # Stop further processing for large datasets
  }
  
  # Adjust legend height based on the number of features
  legend_height <- if (num_features < 35) 10 else num_features / 5
  
  # Define colors
  colors <- c(
    'goldenrod1', 'navy', 'royalblue', 'indianred3',
    'darkolivegreen1', 'darkgreen', 'darkolivegreen3', 'chartreuse4'
  )
  
  # Prepare barplot settings
  par(mar = c(5, 4, 4, 10), xpd = TRUE)
  
  # Create the barplot
  bar_positions <- barplot(
    efs_table,
    xlim = c(0, 1),
    main = 'Ensemble Feature Selection',
    horiz = TRUE,
    las = 2,
    names.arg = abbreviate(colnames(efs_table)),
    col = colors
  )
  
  # Add a legend
  legend(
    "topright", inset = c(-0.2, 0),
    legend = rownames(efs_table),
    col = colors, lty = 1, lwd = 12
  )
  
  # Add text labels for bar values
  text(
    colSums(efs_table) + 0.065, bar_positions,
    format(round(colSums(efs_table), 2), nsmall = 2)
  )
  
  # Add a vertical reference line
  segments(1, 0, 1, 1.25 * num_features, lty = 3, col = "gray40")
}

getVolcanoPlotEFS <- function(DAP_table_factor, FactorLevel.1, FactorLevel.2, mrlimit = 1.2, pValue_Limit= 0.05, legend_title = NULL, additional_data = NULL, column_gene_name = "Gene_primary"){
  if(is.null(additional_data)){
    message("additional_data is mandatory")
    break
  }else{
    differential_data <- data.frame(gene_name = DAP_table_factor$Gene_primary , 
                                    meanRatio = 2^DAP_table_factor[[colnames(DAP_table_factor)[grepl("logFC",colnames(DAP_table_factor))]]], 
                                    pValue = DAP_table_factor[[colnames(DAP_table_factor)[grepl("adj.P.Val",colnames(DAP_table_factor))]]], 
                                    EFS = "No")# forming a new differential data for plotting
    
    # Color
    differential_data$EFS <- "No"
    differential_data[(differential_data$meanRatio> mrlimit) & (differential_data$pValue< pValue_Limit), "EFS"] = "Sig"
    differential_data[(differential_data$meanRatio< 1/mrlimit) & (differential_data$pValue< pValue_Limit), "EFS"] = "Sig"
    differential_data[differential_data$gene_name %in% additional_data[[column_gene_name]], "EFS"] <- "Yes"
    
    differential_data$EFS <- factor(differential_data$EFS, levels = c("Yes", "Sig", "No"))
    regulation_color <- setNames(c("red", "black", "grey"),
                                 levels(differential_data$EFS))
    # Text 
    differential_data$delabel <- ifelse(differential_data$EFS == "Yes",
                                        differential_data$gene_name, "")
    # Plotting
    volcano_plot = ggplot(data = differential_data,
                          mapping = aes(log2(meanRatio),
                                        -log10(pValue),
                                        label = delabel)) +
      # point
      geom_point(
        aes(colour = EFS), 
        # shape= 21,
        size = 2,
        alpha = 5/10
      )  +
      # line
      geom_vline(xintercept = c(-log2(mrlimit),
                                log2(mrlimit)),
                 col = "black",
                 linetype = 'dashed') +
      geom_hline(yintercept = -log10(pValue_Limit),
                 col = "black",
                 linetype = 'dashed')+
      #scale_fill_manual(values = rep("white", 3))   +
      # text: to overcome the text overlap
      geom_text_repel(max.overlaps = Inf, 
                      color =  "red",
                      position = 
                        position_nudge_to(x = 2.3),
                      min.segment.length = 0,
                      segment.color = "black",
                      arrow = arrow(length = unit(0.015, "npc")),
                      direction = "y") +
      labs(#color = legend_title,
        x = substitute(log[2] ~ "Fold Change (" * F/C * ")", list(F = FactorLevel.1, C = FactorLevel.2)) ,
        y = expression("-log"[10]*"(Pvalue)"),
        caption = glue("Pvalue < {pValue_Limit}; Fold Change > {mrlimit}")) +
      scale_color_manual(values = regulation_color) + 
      
      
      theme_minimal() +
      theme(text = element_text(family = "serif"),
            axis.line = element_line(),
            axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
            axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
            plot.title = element_text(hjust = 0.5)
      )
    
    return(volcano_plot)
    
  }
  
}
