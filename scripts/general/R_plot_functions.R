# libraries 

# check available fonts: names(pdfFonts()) 

library(ggplot2)
library(tidyr)
library(dplyr)
library(glue)
library(patchwork)
library(ggridges)
library(plotly)
library(cowplot)
library(gridExtra)
library(stringr)



# Define root dir for plots 
plotRoot <- "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/all_plots"
# plotRoot <- "/home/xil4009/store_tilgnerlab/Exome_Enrich_additional/plots"

# Custome theme (theme_ipsum) plot_title_face = "bold"
theme_jennie <- function(base_family = "ArialMT", base_size = 12,
                        plot_title_family = base_family, plot_title_size = 18,
                        plot_title_face = "plain", plot_title_margin = 10,
                        subtitle_family = base_family, subtitle_size = 12,
                        subtitle_face = "plain", subtitle_margin = 15,
                        strip_text_family = base_family, strip_text_size = 12,
                        strip_text_face = "plain", caption_family = base_family,
                        caption_size = 9, caption_face = "italic", caption_margin = 10,
                        axis_text_size = base_size, axis_title_family = subtitle_family,
                        axis_title_size = base_size, axis_title_face = "plain", axis_title_just = "cc", #rt
                        plot_margin = margin(30, 30, 30, 30), grid_col = "#cccccc",
                        grid = FALSE, axis_col = "#2b2b2b", axis = TRUE, ticks = TRUE,
                        tick_label_gap = 2) {

  # Initialize theme with minimal base
  ret <- ggplot2::theme_minimal(base_size = base_size, base_family = base_family) + 
    theme(legend.background = element_blank(),
          legend.key = element_blank(), 
          legend.text = element_text(size = base_size-1))

  # Grid lines settings
  if (isTRUE(grid) || is.character(grid)) {
    ret <- ret + theme(panel.grid.major.y = element_line(color = grid_col, linewidth = 0.2), # horizontal only 
                       panel.grid.minor.y = element_line(color = grid_col, linewidth = 0.15))
      # panel.grid = element_line(color = grid_col, linewidth = 0.2),
      #                  panel.grid.major = element_line(color = grid_col, linewidth = 0.2),
      #                  panel.grid.minor = element_line(color = grid_col, linewidth = 0.15))

    if (is.character(grid)) {
      if (!grepl("X", grid)) ret <- ret + theme(panel.grid.major.x = element_blank())
      if (!grepl("Y", grid)) ret <- ret + theme(panel.grid.major.y = element_blank())
      if (!grepl("x", grid)) ret <- ret + theme(panel.grid.minor.x = element_blank())
      if (!grepl("y", grid)) ret <- ret + theme(panel.grid.minor.y = element_blank())
    }
  } else {
    ret <- ret + theme(panel.grid = element_blank())
  }

  # Axis lines settings
  if (isTRUE(axis) || is.character(axis)) {
    ret <- ret + theme(axis.line = element_line(color = "#2b2b2b", linewidth = 0.5))

    if (is.character(axis)) {
      if (!grepl("x", axis, ignore.case = TRUE)) ret <- ret + theme(axis.line.x = element_blank())
      if (!grepl("y", axis, ignore.case = TRUE)) ret <- ret + theme(axis.line.y = element_blank())
    } else {
      ret <- ret + theme(axis.line.x = element_line(color = axis_col, linewidth = 0.5),
                         axis.line.y = element_line(color = axis_col, linewidth = 0.5))
    }
  } else {
    ret <- ret + theme(axis.line = element_blank())
  }

  # Ticks settings
  if (!ticks) {
    ret <- ret + theme(axis.ticks = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
  } else {
    ret <- ret + theme(axis.ticks = element_line(linewidth = 0.3),
                       axis.ticks.x = element_line(linewidth = 0.3),
                       axis.ticks.y = element_line(linewidth = 0.3),
                       axis.ticks.length = grid::unit(4, "pt"))
  }

  # Axis title justification
  xj <- switch(tolower(substr(axis_title_just, 1, 1)), b = 0, l = 0, m = 0.5, c = 0.5, r = 1, t = 1)
  yj <- switch(tolower(substr(axis_title_just, 2, 2)), b = 0, l = 0, m = 0.5, c = 0.5, r = 1, t = 1)

  # Adding all text elements
  ret <- ret + theme(axis.text.x = element_text(size = axis_text_size, margin = margin(t = tick_label_gap)),
                     axis.text.y = element_text(size = axis_text_size, margin = margin(r = tick_label_gap)),
                     axis.title = element_text(size = axis_title_size, family = axis_title_family),
                     axis.title.x = element_text(hjust = xj, size = axis_title_size, family = axis_title_family, face = axis_title_face),
                     axis.title.y = element_text(hjust = yj, size = axis_title_size, family = axis_title_family, face = axis_title_face),
                     axis.title.y.right = element_text(hjust = yj, size = axis_title_size, angle = 90, family = axis_title_family, face = axis_title_face),
                     strip.text = element_text(hjust = 0, size = strip_text_size, face = strip_text_face, family = strip_text_family),
                     panel.spacing = grid::unit(2, "lines"),
                     plot.title = element_text(hjust = 0.5, size = plot_title_size, margin = margin(b = plot_title_margin), family = plot_title_family, face = plot_title_face),
                     plot.subtitle = element_text(hjust = 0, size = subtitle_size, margin = margin(b = subtitle_margin), family = subtitle_family, face = subtitle_face),
                     plot.caption = element_text(hjust = 1, size = caption_size, margin = margin(t = caption_margin), family = caption_family, face = caption_face),
                     plot.margin = plot_margin)

  ret
}


save_png <- function(plot, prefix, width, height, res=300, units="in") {
  print(paste0(prefix,".png"))
  png(filename = paste0(prefix,".png"), width = width, height = height, units = units, res = res)
  print(plot)
  dev.off()
}


save_pdf <- function(plot, prefix, width, height) {
  print(paste0(prefix,".pdf"))
  pdf(file = paste0(prefix,".pdf"), width = width, height = height)
  print(plot)
  dev.off()
}


func="
Loop over a list of experiments, read relevant file, combine to one df. 
makes a line plot showing enrich ratio over time. 
Requires output `/timestamp/allinfo_timebin/summary_allReadTypes_ontarget_readsPerTimebin`.
"
format_enrich_ratio_timebin <- function(exp_list, names) {
    df <- data.frame() # Init 
    for (i in 1:length(exp_list)) {
        exp = exp_list[[i]]
        name = names[[i]]
        print(exp)
        print(name)
        file_path <- paste0(base_path, exp, "/timestamp/allinfo_timebin/summary_allReadTypes_ontarget_readsPerTimebin")
        temp_df <- read.table(file_path, sep = "\t", header = TRUE)

        temp_df$name <- name
        temp_df$target <- strsplit(name," ")[[1]][1]
        temp_df$sample <- strsplit(name," ")[[1]][2]

        temp_df$exp <- exp
        # parts <- strsplit(exp, "_")[[1]]
        # sample <- parts[1]
        # target <- paste(parts[3], parts[4], parts[5], sep = "_")
        # temp_df$sample <- sample
        # temp_df$target <- target
        # temp_df$name <- paste(sample, target, sep = "_")
        df <- rbind(df, temp_df)
    }

    df <- df %>%
      separate(hour_range, into = c("hour_start", "hour_end"), sep = "-") %>%
      mutate(hour_start = as.integer(hour_start)) %>%
      arrange(hour_start) %>%
      mutate(hour_start = factor(hour_start, levels = unique(hour_start), ordered = TRUE)) %>%
      mutate(bcFilter_ratio = accept_bcFilter_splice / control_bcFilter_splice,
        noBCFilter_ratio = accept_noBCFilter_splice / control_noBCFilter_splice,
        allinfo_ratio = accept_bcFilter_splice_allinfo / control_bcFilter_splice_allinfo
        )
}



func="
Loop over a list of experiments, read relevant file, combine to one df. 
Requires output `/timestamp/allinfo_timebin/accept_bcFilter_AllInfo_Ontarget_readsPerGene_Timebin.txt`.
exp_list: list of strings of experiment names. 
returns: a formatted dataframe for input to `cumulative_count_timebin_ridge`.
"
format_readsPerGene_timebin <- function(exp_list, condition="accept", read_type="Target gene") {
    df <- data.frame() # Init 
    for (exp in exp_list) {
      print(exp)
      if (condition=="accept" && read_type=="Target gene") {
        file_path <- paste0(base_path, exp, "/timestamp/allinfo_timebin/accept_bcFilter_AllInfo_Ontarget_readsPerGene_Timebin.txt")
      } else if (condition=="accept" && read_type=="Non-target gene") {
        file_path <- paste0(base_path, exp, "/timestamp/allinfo_timebin/accept_bcFilter_AllInfo_Offtarget_readsPerGene_Timebin.txt")
      } else if (condition=="control" && read_type=="Target gene"){
        file_path <- paste0(base_path, exp, "/timestamp/allinfo_timebin/control_bcFilter_AllInfo_Ontarget_readsPerGene_Timebin.txt")
      } else if (condition=="control" && read_type=="Non-target gene"){
        file_path <- paste0(base_path, exp, "/timestamp/allinfo_timebin/control_bcFilter_AllInfo_Offtarget_readsPerGene_Timebin.txt")
      }

      temp_df <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
      colnames(temp_df) <- c("time_bin","hour_range","gene","count","ident")

      temp_df$exp <- exp
      parts <- strsplit(exp, "_")[[1]]
      sample <- parts[1]
      device <- parts[2]
      target <- paste(parts[3], parts[4], parts[5], sep = "_")
      temp_df$sample <- sample
      temp_df$target <- target
      temp_df$name <- paste(sample, target, sep = "_")
      temp_df$device <- device
      temp_df$condition <- condition
      temp_df$read_type <- read_type
      df <- rbind(df, temp_df)
    }

    df <- df %>%
        separate(hour_range, into = c("hour_start", "hour_end"), sep = "-") %>%
        mutate(hour_start = as.integer(hour_start)) %>%
        arrange(hour_start) %>%
        mutate(hour_start = factor(hour_start, levels = unique(hour_start), ordered = TRUE))
}


func="
make a ridge plot for one `experiment` in the df and return the plot.
  df_all: output of `format_readsPerGene`.
  experiment: str of the name of an experiment.
  color: the fill color for ridge plot.
"
cumulative_count_timebin_ridge <- function(df_all, experiment, color, name, xmax=FALSE) {

    df <- df_all %>% filter(exp == experiment)
    
		type <- unique(df$read_type) # is on-target or non-target 
    device <- unique(df$device)
    
    if (device == "minion") {
      duration = 15
      prob_cutoff = 0.98
    } else if (device == "prom") {
      duration = 45 
      prob_cutoff = 0.95
    }	else {
      print("wrong device")
    }
    
    # Transform to cumulative counts
    df <- df %>%
      arrange(gene, time_bin) %>%
      group_by(gene) %>%
      mutate(cumulative_count = cumsum(count)) %>%
      ungroup()
    
    # Filter the dataframe to include only the first 50 hours
    # Drop unused columns 
    df <- df %>%
      filter(as.integer(as.character(hour_start)) <= duration) %>%
      select(gene, hour_start, cumulative_count, name)

    # Calculate cumulative count at 90th percentile of genes
    if (!(xmax)) {
      cutoff_cumulative_count <- quantile(df$cumulative_count, probs = prob_cutoff)
      print(cutoff_cumulative_count)
    } else {
      cutoff_cumulative_count <- xmax
    }
    p <- ggplot(df, aes(x = cumulative_count, y = hour_start) ) +
        geom_density_ridges(aes(height = after_stat(density)), rel_min_height = 0.005,
                            fill = color, scale = 10, alpha=0.6, linewidth=1.1) +
        labs(title = glue("{name} ({device})"),
              x = glue("Cumulative allinfo reads per {type} gene"), y = "Sequencing time (h)") +
        scale_y_discrete(expand = c(0.01, 0)) + scale_x_continuous(expand = c(0.01, 0)) +  
        xlim( 0, cutoff_cumulative_count+10 ) + 
        scale_fill_brewer(palette = 4) + 
        theme_ridges(font_size=16, grid=FALSE) + 
        theme(
          plot.background = element_rect(fill = "white"), 
          panel.background = element_rect(fill = "white"),
          panel.grid.major.x = element_blank(),  # Hide x-axis grid lines
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(color = "gray", size = 0.2),  # Customize y-axis grid lines
          panel.grid.minor.y = element_blank()  # Hide minor y-axis grid lines (if any)
        )
      
    title = "cumulative_read_over_time"
    prefix = glue("{plotDir}/{title}_{name}")
    save_png(p, prefix, 10, 10)
    
    return (p)
}


cumulative_count_timebin_stacked <- function(df_on, df_off, experiment, name, ymax=FALSE) {
    print(name)
    df_combined <- rbind(df_on, df_off)

    df <- df_combined %>% filter(exp == experiment) # select experiment 
    
    df <- df %>% filter(as.integer(as.character(time_bin)) <= 20)

    if (!ymax) { ymax <- sum(df$count) }
    print(ymax)

    # Calculate total counts for each read_type and hour_start
    df <- df %>%
      select(hour_start, count, read_type) %>%
      group_by(hour_start, read_type) %>%
      summarise(total_count = sum(count), .groups = 'drop')

    # print(df %>% tail(2))
    df <- df %>% 
      arrange(hour_start) %>%
      group_by(read_type) %>%
      mutate(cumulative_count = cumsum(total_count))
    # print(df %>% head(30))

    # Create the stacked bar plot
    p <- ggplot(df, aes(x = hour_start, y = cumulative_count, fill = read_type)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(title = name,
          x = "Time (hours)",
          y = "Cumulative read (K)",
          fill = "Read type") +
      theme_jennie(plot_title_size=16, plot_title_face = "plain", base_size=16, axis_text_size=15) +
      scale_fill_manual(values = c("Target gene" = "#f2a93b", "Non-target gene" = "#1f77b4")) + #b25d91
      scale_y_continuous(expand = c(0, 0), 
                        labels = scales::number_format(big.mark="",scale = 0.001), # y-labels to K unit
                        breaks = seq(0, ymax, length.out = 5)) +  # Create 5 breaks for 4 intervals
      coord_cartesian(ylim = c(0, ymax + 100)) + 
      theme(legend.title = element_blank(),
            axis.title.y = element_text(margin = margin(r = 10)))  # Increase margin

    # # Define the breaks to show on the x-axis
    if (grepl("M", name)) {
      x_breaks <- seq(0, 48, by = 2)  # Shows labels every 4 hours
      p <- p + scale_x_discrete(breaks = x_breaks, labels = x_breaks) 
    }

    return (p)
}

# show only targeted gene count 
log_cumulative_count_timebin_stacked <- function(df_on, df_off, experiment, name, ymax=FALSE) {
    print(name)
    df_combined <- rbind(df_on) #, df_off)

    df <- df_combined %>% filter(exp == experiment) # select experiment 
    # print(tail(df))
    df <- df %>% filter(as.integer(as.character(time_bin)) <= 20)

    # Calculate total counts for each read_type and hour_start
    df <- df %>%
      select(hour_start, count, read_type) %>%
      group_by(hour_start, read_type) %>%
      summarise(total_count = sum(count), .groups = 'drop')

    df <- df %>% 
      arrange(hour_start) %>%
      group_by(read_type) %>%
      mutate(cumulative_count = cumsum(total_count))
    # print(tail(df))

    # # Add a log-transformed cumulative count column
    # df <- df %>%
    #   mutate(cumulative_count_log = log2(cumulative_count + 1))  # Added +1 to avoid log(0) issues
    # # print(tail(df))

    # Determine ymax as the sum of the largest cumulative_count_log for Non-target gene and Target gene
    # if (!ymax) {
    #   ymax_ntg <- max(df %>% filter(read_type == "Non-target gene") %>% pull(cumulative_count_log), na.rm = TRUE)
    #   ymax_tg <- max(df %>% filter(read_type == "Target gene") %>% pull(cumulative_count_log), na.rm = TRUE)
    #   ymax <- ymax_ntg + ymax_tg
    # }
    # print(ymax)

    p <- ggplot(df, aes(x = hour_start, y = cumulative_count, fill = read_type)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(title = name,
          x = "Time (hours)",
          y = "Cumulative read (K)",
          fill = "Read type") +
      theme_jennie(plot_title_size=16, base_size=16, axis_text_size=15)  +
      scale_fill_manual(values = c("Target gene" = "#f2a93b", "Non-target gene" = "#1f77b4")) +
      scale_y_continuous(expand = c(0, 0), 
                        labels = scales::number_format(big.mark="",scale = 0.001), # y-labels to K unit
                        breaks = seq(0, ymax, length.out = 5)) +  # Create 5 breaks for 4 intervals
      coord_cartesian(ylim = c(0, ymax)) + 
      theme(legend.title = element_blank(),
            axis.title.y = element_text(margin = margin(r = 10)))  # Increase margin
    
    # # Define the breaks to show on the x-axis
    if (grepl("M", name)) {
      x_breaks <- seq(0, 48, by = 2)  # Shows labels every 4 hours
      p <- p + scale_x_discrete(breaks = x_breaks, labels = x_breaks) 
    }

    return (p)
}

# Main figure 2: ontarget reads over time + filled bar plot  
ontarget_cumulative_count_accept_control <- function(exp, name, ymax=FALSE) {

    print(exp)
    file_path <- paste0(base_path, exp, "/timestamp/allinfo_timebin/accept_bcFilter_AllInfo_Ontarget_readsPerGene_Timebin.txt")
    accept_target <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    accept_target$condition <- "Accept"
    accept_target$read_type <- "Target"

    file_path <- paste0(base_path, exp, "/timestamp/allinfo_timebin/accept_bcFilter_AllInfo_Offtarget_readsPerGene_Timebin.txt")
    accept_nontarget <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    accept_nontarget$condition <- "Accept"
    accept_nontarget$read_type <- "Non-target"
    
    file_path <- paste0(base_path, exp, "/timestamp/allinfo_timebin/control_bcFilter_AllInfo_Ontarget_readsPerGene_Timebin.txt")
    control_target <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    control_target$condition <- "Control"
    control_target$read_type <- "Target"

    file_path <- paste0(base_path, exp, "/timestamp/allinfo_timebin/control_bcFilter_AllInfo_Offtarget_readsPerGene_Timebin.txt")
    control_nontarget <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    control_nontarget$condition <- "Control"
    control_nontarget$read_type <- "Non-target"

    df <- rbind(accept_target, accept_nontarget, control_target, control_nontarget)
    colnames(df) <- c("time_bin","hour_range","gene","count","ident", "condition", "read_type")

    df <- df %>%
        separate(hour_range, into = c("hour_start", "hour_end"), sep = "-") %>%
        mutate(hour_start = as.integer(hour_start)) %>%
        arrange(hour_start) %>%
        mutate(hour_start = factor(hour_start, levels = unique(hour_start), ordered = TRUE))
    
    # Define the time range to display 
    if (grepl("\\(P\\)", name)) {
      max_timbin = 13 
    } else {
      max_timbin = 20 
    }
    df <- df %>% filter(as.integer(as.character(time_bin)) <= max_timbin)

    df_summary <- df %>%
      dplyr::select(hour_start, count, condition, read_type) %>%
      group_by(hour_start, condition, read_type) %>%
      summarise(total_count = sum(count), .groups = 'drop')

    df_summary <- df_summary %>% 
      arrange(hour_start) %>%
      group_by(condition, read_type) %>%
      mutate(cumulative_count = cumsum(total_count))

    
    df_target <- df_summary[df_summary$read_type=="Target",]


    if (!ymax) { ymax <- max(df_target$cumulative_count) }
    print(ymax)

    df_accept <- df_target %>% filter(condition == "Accept") 
    df_control <- df_target %>% filter(condition == "Control") 

    # accept & control ontarget read count over time, stacked bar plot 
    # p <- ggplot(df_target, aes(x = hour_start, y = cumulative_count, fill = condition)) +
    #   geom_bar(stat = "identity", position = "identity", alpha = 1) +  # Overlay bars with transparency
    #   # geom_bar(stat = "identity", position = "stack") +
    #   labs(title = name,
    #       x = "Time (hours)",
    #       y = "Cumulative read (K)",
    #       fill = "Condition") +
    #   theme_jennie(plot_title_size=16, plot_title_face = "plain", base_size=15, axis_text_size=15) +
    #   scale_fill_manual(values = c("Accept" = "#f2a93b", "Control" = "skyblue")) + 
    #   scale_y_continuous(expand = c(0, 0), 
    #                     labels = scales::number_format(big.mark="",scale = 0.001), # y-labels to K unit
    #                     breaks = seq(0, ymax, length.out = 5)) +  # Create 5 breaks for 4 intervals
    #   coord_cartesian(ylim = c(0, ymax + 100)) + 
    #   theme(legend.title = element_blank(),
    #         legend.position = c(0.01, 0.99),  # Upper left corner (x, y)
    #         legend.justification = c("left", "top") )   # Align the legend to the top left

    # # overlay bars
    # p <- ggplot() +
    #   geom_bar(data = df_accept, aes(x = hour_start, y = cumulative_count, fill = "Accept"), 
    #           width = 0.9, stat = "identity", alpha = 0.6, show.legend = TRUE) +  # Overlay bars with transparency
    #   geom_bar(data = df_control, aes(x = hour_start, y = cumulative_count, fill = "Control"), 
    #           width = 0.5, stat = "identity", alpha = 1, show.legend = TRUE) +  # Overlay bars with transparency
    #   labs(title = name,
    #       x = "Time (hours)",
    #       y = "Cumulative read (K)") +
    #   theme_jennie(plot_title_size = 16, plot_title_face = "plain", base_size = 15, axis_text_size = 15) +
    #   scale_y_continuous(expand = c(0, 0),
    #                     labels = scales::number_format(big.mark = "", scale = 0.001),  # y-labels to K unit
    #                     breaks = seq(0, ymax, length.out = 5)) +  # Create 5 breaks for 4 intervals
    #   scale_fill_manual(values = c("Accept" = "#f2a93b", "Control" = "skyblue")) +  # Define manual fill colors
    #   coord_cartesian(ylim = c(0, ymax + 100)) + 
    #   theme(
    #     legend.title = element_blank(),
    #     legend.position = c(0.01, 0.99),  # Upper left corner (x, y)
    #     legend.justification = c("left", "top"))  # Align the legend to the top left

    # # doged bars
    # p <- ggplot(df_target, aes(x = hour_start, y = cumulative_count, fill = condition)) +
    #   geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 1) +  # Side-by-side bars
    #   labs(title = name,
    #       x = "Time (hours)",
    #       y = "Cumulative read (K)",
    #       fill = "Condition") +
    #   theme_jennie(plot_title_size = 16, plot_title_face = "plain", base_size = 15, axis_text_size = 15) +
    #   scale_fill_manual(values = c("Accept" = "#f2a93b", "Control" = "skyblue")) +  # Define manual fill colors
    #   scale_y_continuous(expand = c(0, 0),
    #                     labels = scales::number_format(big.mark = "", scale = 0.001),  # y-labels to K unit
    #                     breaks = seq(0, ymax, length.out = 5)) +  # Create 5 breaks for 4 intervals
    #   coord_cartesian(ylim = c(0, ymax + 100)) + 
    #   theme(
    #     legend.title = element_blank(),
    #     legend.position = c(0.01, 0.99),  # Upper left corner (x, y)
    #     legend.justification = c("left", "top")   # Align the legend to the top left
    #   )

    # dot and line 
    p <- ggplot(df_target, aes(x = hour_start, y = cumulative_count, color = condition, group = condition)) +
      geom_line(size = 0.5) +  # Line plot
      geom_point(size = 1.5) +    # Add points on the lines for emphasis
      labs(title = name,
          x = "Time (hours)",
          y = "Cumulative read (K)",
          color = "Condition") +
      theme_jennie(plot_title_size = 12, plot_title_face = "plain", base_size = 12, axis_text_size = 12) +
      scale_color_manual(
          values = c("Accept" = "#f2a93b", "Control" = "skyblue"), # Define manual line colors
          labels = c("Accept" = "Accept Target", "Control" = "Control Target")  # Custom labels
          ) +  
      scale_y_continuous(expand = c(0, 0),
                        labels = scales::number_format(big.mark = "", scale = 0.001),  # y-labels to K unit
                        breaks = seq(0, ymax, length.out = 5)) +  # Create 5 breaks for 4 intervals
      coord_cartesian(ylim = c(0, ymax *1.1)) + 
      theme(
        legend.title = element_blank(),
        legend.position = c(0.25, 0.9), 
        legend.text = element_text(size = 11),  # Adjust legend text size
        legend.spacing.y = unit(0.1, "cm"),  # Reduce vertical spacing between legend items
        legend.key.width = unit(0.4, "cm"),  # Reduce the width of the legend keys to minimize the gap
        legend.key.height = unit(0.5, "cm")  # Control height of the legend keys (optional)
  ) 

    # Define the breaks to show on the x-axis
    if (grepl("M|P", name)) {
      x_breaks <- seq(0, 60, by = 2)  # Shows labels every 4 hours
      p <- p + scale_x_discrete(breaks = x_breaks, labels = x_breaks) 
    }

    # pie charts 
    df_piechart <- df_summary %>%
      group_by(condition, read_type) %>%
      summarise(total_count = sum(total_count)) %>% 
      mutate(pct = total_count / sum(total_count)*100)
    print( df_piechart$total_count ) 


    colors <- c("Target" = "#f2a93b", "Non-target" = "#D0D3D6") # 橙 灰
    df_accept <- df_piechart %>% filter(condition == "Accept") 

    pie_accept <- ggplot(df_accept, aes(x = "", y = total_count, fill = read_type)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y") + # Convert to pie chart
      geom_text(aes(label = glue("{round(pct, 1)}%")),
                position = position_stack(vjust = 0.5), size = 3, color="black") +
      labs(x = NULL, y = NULL, fill = "Accept", title = NULL) +
      scale_fill_manual(values = colors) + 
      theme_void(base_family = "ArialMT") +  # Clean theme for pie chart
      theme(
          # plot.title = element_text(hjust = 0.5, size = 16),  # Center the title
          legend.position = "right"                        # Position legend at the bottom
        ) 

    colors <- c("Target" = "skyblue", "Non-target" = "#D0D3D6") # 蓝 灰

    df_control <- df_piechart %>% filter(condition == "Control") 

    pie_control <- ggplot(df_control, aes(x = "", y = total_count, fill = read_type)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y") + # Convert to pie chart
      geom_text(aes(label = glue("{round(pct, 1)}%")),
                position = position_stack(vjust = 0.5), size = 3) +
      labs(x = NULL, y = NULL, fill = "Control", title = NULL) +
      scale_fill_manual(values = colors) + 
      theme_void(base_family = "ArialMT") +  # Clean theme for pie chart
      theme(
          legend.position = "right" 
        ) 
        
    combined_plot = p + (pie_accept / pie_control) + 
                      plot_layout(widths = c(3, 1))

    return (combined_plot)
}



    # # stacked bar chart for ratio 
    # df_stackedbar <- df_summary %>%
    #   group_by(condition, read_type) %>%
    #   summarise(total_count = sum(total_count)) %>% 
    #   mutate(frac = total_count / sum(total_count)) %>%
    #   # mutate(frac=ifelse(frac < 0.05, NA, frac)) %>% # avoid squished labels 
    #   mutate(condition = factor(condition, levels = c("Accept", "Control")))  # Ensure correct order
    # print(df_stackedbar)

    # # Create a stacked bar chart
    # stacked_bar <- ggplot(df_stackedbar, aes(x = condition, y = frac, fill = read_type)) +
    #   geom_bar(stat = "identity") + 
    #   # geom_text(aes(label = total_count), 
    #     geom_text(aes(label = ifelse(frac > 0.05, glue("{total_count} \n{round(frac*100,1)}%"), "")),  # Only show labels for > 5% areas
    #             position = position_stack(vjust = 0.5), 
    #             size = 4, color = "black", family = "ArialMT") +
    #   labs(x = NULL, y = NULL, fill = NULL, title = NULL) +
    #   scale_fill_manual(values = c("Non-target" = "#ff6268", "Target" = "#a9d697")) + 
    #   scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +  
    #   theme_jennie(plot_title_size=16, plot_title_face = "plain", base_size=15, axis_text_size=15) +
    #   theme(aspect.ratio = 1.1/1)  # Control aspect ratio for a more compact plot

    # combined_plot = p + stacked_bar



# Function to calculate the max limit after excluding outliers
calculate_max_excluding_outliers <- function(values) {
    q1 <- quantile(values, 0.25, na.rm = TRUE)
    q3 <- quantile(values, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lower_bound <- q1 - 1.5 * iqr
    upper_bound <- q3 + 1.5 * iqr
    filtered_values <- values[values >= lower_bound & values <= upper_bound]
    return(max(filtered_values, na.rm = TRUE))
}

# Function to exclude outliers and get the max value
calculate_max_exclude_99outliers <- function(values) {
  quantile_99 <- quantile(values, 0.99)
  max(values[values <= quantile_99])
}


func="
Requires output `accept_bcFilter_readPerOntargetExon.tsv`, `control_bcFilter_readPerOntargetExon.tsv`
read the files in an experiment dir, make and save scatter plots for exon coverage.  
"
exon_coverage_scatter_plot <- function(exp, name) {
	  parts <- strsplit(exp, "_")[[1]]
	  device<- parts[2]
	  sample <- parts[1]
	  target <- paste(parts[3], parts[4], parts[5], sep = "_")
    title <- glue("{name}")

########### readPerOntargetExon ###########
		control_df <- read.table(glue('{base_path}/{exp}/exon_coverage/control_bcFilter_readPerOntargetExon.tsv'), sep = '\t', header = TRUE)
		accept_df <- read.table(glue('{base_path}/{exp}/exon_coverage/accept_bcFilter_readPerOntargetExon.tsv'), sep = '\t', header = TRUE)

		# Merge based on gene_id # Replace NA with zero
		merged_df <- merge(control_df, accept_df, by = 'anno_exonblock', all.x = TRUE, all.y = TRUE) %>%
		  replace_na(list(control_bcFilter_unique_readID_count = 0,
		                  accept_bcFilter_unique_readID_count = 0))
    
    # print(merged_df %>%
    #       filter(control_bcFilter_unique_readID_count == 50 & accept_bcFilter_unique_readID_count == 0)
    #     ) # when 100 in control and 0 in accept. 

    # count dots on each side of the diagonal line 
    larger_in_accept = sum(merged_df$accept_bcFilter_unique_readID_count > merged_df$control_bcFilter_unique_readID_count)
    larger_in_control = sum(merged_df$accept_bcFilter_unique_readID_count < merged_df$control_bcFilter_unique_readID_count)
    equal = sum(merged_df$accept_bcFilter_unique_readID_count == merged_df$control_bcFilter_unique_readID_count)
    total <- nrow(merged_df)
    # Calculate correlation and p-value
    cor_test <- cor.test(merged_df$control_bcFilter_unique_readID_count, merged_df$accept_bcFilter_unique_readID_count)
    cor_coeff <- cor_test$estimate
    p_value <- cor_test$p.value
    
    # get dot frequency for coloring  
		frequency_df <- merged_df %>%
		  group_by(control_bcFilter_unique_readID_count, accept_bcFilter_unique_readID_count) %>% 
		  summarise(freq = n(), .groups = 'drop')
	  # Merge the frequency back to the merged DataFrame
		merged_df <- merged_df %>%
		  left_join(frequency_df, by = c("control_bcFilter_unique_readID_count", "accept_bcFilter_unique_readID_count"))
		
    # Define bin labels colors
    if (device == "minion") {
      bin_breaks <- c(-Inf, 10, 20, 50, Inf)
      bin_labels <- c("<=10", "10-20", "20-50", ">50")
    } else {
      bin_breaks <- c(-Inf, 5, 10, 50, Inf)
      bin_labels <- c("<=5", "5-10", "10-50", ">50")
    }
    bin_colors <- c("#66FFFF","#FFCCFF","#FF66FF","#993FFF")
    names(bin_colors) <- names(bin_labels)
    # bin_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")

    # Categorize 'freq' into bins
    merged_df <- merged_df %>%
      mutate(freq_bin = cut(freq, breaks = bin_breaks, labels = bin_labels))
      

    # Use the larger of the two max limits for both x and y limits
    max_limit <- max(merged_df$control_bcFilter_unique_readID_count, merged_df$accept_bcFilter_unique_readID_count)

		# Scatter plot
		p1 <- ggplot(merged_df, aes(x = control_bcFilter_unique_readID_count, 
                                y = accept_bcFilter_unique_readID_count, # color = freq,
                                color=freq_bin) ) +
		  geom_point(size=1) +
		  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') +  # diagonal line
		  labs(
		    title = glue("{title}"),
		    x = 'Read number per exon in control',
		    y = 'Read number per exon in enriched',
        color = 'Freq',
		  ) +
		  xlim(0,max_limit) + ylim(0,max_limit) + 
      scale_color_manual(values = bin_colors) +
      theme_jennie(plot_title_size=12, base_size=12) +
      guides(colour = guide_legend(override.aes = list(size=5)))+ 
      annotate("text", x = 0.99 * max_limit, y = 0.15 * max_limit, 
                label = paste("r =", round(cor_coeff, 2), "\np =", sprintf("%.2e", p_value)), 
                hjust = 1, vjust = 0, size = 4.5) + # format(p_value, scientific = TRUE)
      coord_fixed(ratio = 1)  # ensures the plot area is square 

########### zoomin readPerOntargetExon ###########
    if (device == "minion") {
        max_limit_control <- calculate_max_exclude_99outliers(merged_df$control_bcFilter_unique_readID_count)
        max_limit_accept <- calculate_max_exclude_99outliers(merged_df$accept_bcFilter_unique_readID_count)
		    max_limit <- max(max_limit_control, max_limit_accept)
    }  else {
        max_limit_control <- calculate_max_excluding_outliers(merged_df$control_bcFilter_unique_readID_count)
        max_limit_accept <- calculate_max_excluding_outliers(merged_df$accept_bcFilter_unique_readID_count)
        max_limit <- max(max_limit_control, max_limit_accept)
    }

    p2 <- ggplot(merged_df, aes(x = control_bcFilter_unique_readID_count, 
                                y = accept_bcFilter_unique_readID_count, 
                                color=freq_bin)) +   # color=(freq^coeff) color=freq
      geom_point(size=1) +
      geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') +  # diagonal line
      labs(
        title = glue( "enriched > control: {larger_in_accept} ({sprintf('%.1f', larger_in_accept/total * 100)}%)\nenriched = control: {equal} ({sprintf('%.1f', equal/total * 100)}%)\nenriched < control: {larger_in_control} ({sprintf('%.1f', larger_in_control/total * 100)}%)"),
        x = 'Read number per exon in control',
        y = 'Read number per exon in enriched',
        color = 'Freq',
      ) +
      xlim(0,max_limit) + ylim(0,max_limit) + 
      scale_color_manual(values = bin_colors) +
      theme_jennie(plot_title_size=12, base_size=12) +
      guides(colour = guide_legend(override.aes = list(size=5))) +
      coord_fixed(ratio = 1)  # This ensures the plot area is square


##### A single multi-panel figure
		combined_plot <- (p1 + p2)

		name <- paste(sample, target, device, sep = "_")
		title = glue("combine_exon_coverage_{name}")
		prefix = glue("{plotDir}/{title}")
		# save_png(combined_plot, prefix, 11, 5, res=1200)
		# save_pdf(combined_plot, prefix, 11, 5)
		cat("\n")
    
    return (combined_plot)
}


exon_per_ontarget_gene_bars <- function(exp_list, names, dataset_order, 
                                        metric_order = c("enrich > control", "enrich = control", "enrich < control")
                                        ) {
    combined_data <- data.frame()

    for (i in seq_along(exp_list)) {
      exp <- exp_list[i]
      name <- names[i]
      
      parts <- strsplit(exp, "_")[[1]]
      device <- parts[2]
      sample <- parts[1]
      target <- paste(parts[3], parts[4], parts[5], sep = "_")
      title <- glue("{name} ({device})")
      
      base_path <- "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets"
      
      control_df <- read.table(glue('{base_path}/{exp}/exon_coverage/control_bcFilter_exonPerOntargetGene.tsv'), sep = '\t', header = TRUE)
      accept_df <- read.table(glue('{base_path}/{exp}/exon_coverage/accept_bcFilter_exonPerOntargetGene.tsv'), sep = '\t', header = TRUE)

      merged_df <- merge(control_df, accept_df, by = 'anno_geneID', all.x = TRUE, all.y = TRUE) %>%
        replace_na(list(control_bcFilter_unique_anno_exonblock_count = 0,
                        accept_bcFilter_unique_anno_exonblock_count = 0)) 
      
      larger_in_accept <- sum(merged_df$accept_bcFilter_unique_anno_exonblock_count > merged_df$control_bcFilter_unique_anno_exonblock_count)
      larger_in_control <- sum(merged_df$accept_bcFilter_unique_anno_exonblock_count < merged_df$control_bcFilter_unique_anno_exonblock_count)
      equal <- sum(merged_df$accept_bcFilter_unique_anno_exonblock_count == merged_df$control_bcFilter_unique_anno_exonblock_count)
      total <- nrow(merged_df)
      
      counts_df <- data.frame(
        metric = metric_order,
        count = c(larger_in_accept, equal, larger_in_control),
        ratio = c(larger_in_accept, equal, larger_in_control) / total,
        dataset = name
      )      
      combined_data <- rbind(combined_data, counts_df)
    }

    # Reshape the combined data frame for ggplot
    df_melted <- reshape2::melt(combined_data, id.vars = c("metric", "dataset", "count"))

    # Set the metric order
    df_melted$metric <- factor(df_melted$metric, levels = metric_order)
    # Specify the order of the 'dataset' variable
    df_melted$dataset <- factor(df_melted$dataset, levels = names)

    # Define colors and bar width
    # colors <- c("#fcc5eb", "#FCCCAC", "#98f0fa", "#A3A3A3") 
    # colors <- c("#f9c7c4", "#9acea8", "#c5d3e9") # 粉 绿 蓝
    # colors <- c("#ffe5d4","#efc7c2","#bfd3c1") # 黄 粉 绿
    colors <- c("#FCE1E8","#C3E5C6","#C4EBF7") # 粉 绿 蓝

    bar_width <- 0.99

    # stacked bars 
    # p <- ggplot(df_melted, aes(x = dataset, y = value, fill = metric)) +
    #   geom_bar(stat = "identity", color="white", width = bar_width) +
    #   geom_text(aes(label = count), 
    #             position = position_stack(vjust = 0.5), 
    #             size = 4, color = "black") +
    #   labs(title = "Exon per target gene in accept and control",
    #       x = NULL, y = NULL, fill = "") + 
    #   scale_fill_manual(values = colors) +
    #   coord_cartesian(ylim = c(0, 1)) +  
    #   scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +  
    #   scale_x_discrete(expand = c(0.15, 0)) +
    #   theme_jennie(plot_title_size=14, axis_title_size=14, base_size=14) +
    #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

    ## dodged bars 
    # p <- ggplot(df_melted, aes(x = dataset, y = value, fill = metric)) +
    #   geom_bar(stat = "identity", position = position_dodge(width = bar_width), width = bar_width) +
    #   labs(title = "exon number per target gene between accept and control",
    #       x = NULL, y = "Exon number", fill = "Group") + 
    #   geom_text(aes(label = value), position = position_dodge(width = bar_width), vjust = -0.5, size = 4) +  
    #   scale_fill_manual(values = colors) +
    #   coord_cartesian(ylim = c(0, max(df_melted$value) * 1.1)) +  
    #   scale_y_continuous(expand = c(0, 0)) +  
    #   theme_jennie(plot_title_size=14, axis_title_size=12, base_size=12) +
    #   theme(axis.text.x = element_text(angle = 45, hjust = 1))

    p <- ggplot(df_melted, aes(x = dataset, y = value, fill = metric))+
      geom_bar(stat = "identity") + 
      geom_text(aes(label = count), 
                position = position_stack(vjust = 0.5), 
                size = 3.5, color = "black", family = "ArialMT") +
      labs(title = NULL,
          x = NULL, y = NULL, fill = "Exon per target gene") +
      scale_fill_manual(values = colors) +  
      scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +  
      theme_jennie(plot_title_size=12, axis_title_size=11, base_size=11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

    return (p)
}


ratio_per_ontarget_gene_bars <- function(exp_list, names, dataset_order, 
                                        metric_order = c("enrich > control", "enrich = control", "enrich < control")
                                        ) {
    combined_data <- data.frame()

    for (i in seq_along(exp_list)) {
      exp <- exp_list[i]
      name <- names[i]
      
      parts <- strsplit(exp, "_")[[1]]
      device <- parts[2]
      sample <- parts[1]
      target <- paste(parts[3], parts[4], parts[5], sep = "_")
      title <- glue("{name} ({device})")
      
      # base_path <- "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/all_plots/supplementary/"
      base_path <- "/home/xil4009/store_tilgnerlab/Exome_Enrich_Major_Datasets/genewise_enrich_ratio_tables"
      
      merged_df <- read.table(glue('{base_path}/geneID_acceptReads_cltReads_enrichRatio_{exp}'), sep = '\t', header = FALSE)
      larger_in_accept <- sum(merged_df$V6 > 1)
      larger_in_control <- sum(merged_df$V6 < 1)
      equal <- sum(merged_df$V6 == 1)
      total <- nrow(merged_df)
      
      counts_df <- data.frame(
        metric = metric_order,
        count = c(larger_in_accept, equal, larger_in_control),
        ratio = c(larger_in_accept, equal, larger_in_control) / total,
        dataset = name
      )      
      combined_data <- rbind(combined_data, counts_df)
    }

    # Reshape the combined data frame for ggplot
    df_melted <- reshape2::melt(combined_data, id.vars = c("metric", "dataset", "count"))

    # Set the metric order
    df_melted$metric <- factor(df_melted$metric, levels = metric_order)
    # Specify the order of the 'dataset' variable
    df_melted$dataset <- factor(df_melted$dataset, levels = names)

    # Define colors and bar width
    # colors <- c("#fcc5eb", "#FCCCAC", "#98f0fa", "#A3A3A3") 
    colors <- c("#f9c7c4", "#9acea8", "#c5d3e9")
    bar_width <- 0.99

    p <- ggplot(df_melted, aes(x = dataset, y = value, fill = metric))+
      geom_bar(stat = "identity") + 
      geom_text(aes(label = count), 
                position = position_stack(vjust = 0.5), 
                size = 3.5, color = "black") +
      labs(title = NULL,
          x = NULL, y = NULL, fill = "Enrich ratio per gene") +
      scale_fill_manual(values = colors) +  
      scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +  
      theme_jennie(plot_title_size=12, axis_title_size=11, base_size=11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

    return (p)
}



func="
Requires output `accept_bcFilter_exonPerOntargetGene.tsv`, `control_bcFilter_exonPerOntargetGene.tsv`.
read the files in an experiment dir, make 3D scatter plot.  
"
exon_per_ontarget_gene_3D_scatter_plot <- function(exp, name) {
	  parts <- strsplit(exp, "_")[[1]]
	  device<- parts[2]
	  sample <- parts[1]
	  target <- paste(parts[3], parts[4], parts[5], sep = "_")
    title <- glue("{name} ({device})")

########### exonPerOntargetGene ###########
    control_df <- read.table(glue('{base_path}/{exp}/exon_coverage/control_bcFilter_exonPerOntargetGene.tsv'), sep = '\t', header = TRUE)
    accept_df <- read.table(glue('{base_path}/{exp}/exon_coverage/accept_bcFilter_exonPerOntargetGene.tsv'), sep = '\t', header = TRUE)

    # Merge and process data
    merged_df <- merge(control_df, accept_df, by = 'anno_geneID', all.x = TRUE, all.y = TRUE) %>%
      replace_na(list(control_bcFilter_unique_anno_exonblock_count = 0,
                      accept_bcFilter_unique_anno_exonblock_count = 0)) 

    # Get counts of comparisons
    larger_in_accept <- sum(merged_df$accept_bcFilter_unique_anno_exonblock_count > merged_df$control_bcFilter_unique_anno_exonblock_count)
    larger_in_control <- sum(merged_df$accept_bcFilter_unique_anno_exonblock_count < merged_df$control_bcFilter_unique_anno_exonblock_count)
    equal <- sum(merged_df$accept_bcFilter_unique_anno_exonblock_count == merged_df$control_bcFilter_unique_anno_exonblock_count)
    total <- nrow(merged_df)

    # Get frequency for coloring
    frequency_df <- merged_df %>%
      group_by(control_bcFilter_unique_anno_exonblock_count, accept_bcFilter_unique_anno_exonblock_count) %>% 
      summarise(freq = n(), .groups = 'drop')
    # Merge frequency back
    merged_df <- merged_df %>%
      left_join(frequency_df, by = c("control_bcFilter_unique_anno_exonblock_count", "accept_bcFilter_unique_anno_exonblock_count"))

    count_df <- merged_df %>%
          mutate(comparison = case_when(
            accept_bcFilter_unique_anno_exonblock_count > control_bcFilter_unique_anno_exonblock_count ~ glue("enrich > control: {larger_in_accept}"),
            accept_bcFilter_unique_anno_exonblock_count < control_bcFilter_unique_anno_exonblock_count ~ glue("enrich < control: {larger_in_control}"),
            TRUE ~ glue("enrich = control: {equal}")
          ))

          
    # colors 紫 黄 蓝: Equal, control, accept 
    fig <- plot_ly(count_df, 
                  y = ~accept_bcFilter_unique_anno_exonblock_count,  # the order of x and y matters! 
                  x = ~control_bcFilter_unique_anno_exonblock_count, 
                  z = ~freq, 
                  color = ~comparison, colors = c('#8338ec', '#ffbe0b', '#3a86ff'), 
                  marker = list(size = 5)) %>%
      add_markers() %>%
      layout(title = glue("{title} \nenrich > control: {larger_in_accept} \nenrich < control: {larger_in_control} \nenrich = control: {equal}"),
            scene = list(xaxis = list(title = 'Unique exon number per target gene in control'),
                          yaxis = list(title = 'Unique exon number per target gene in enrich'),
                          zaxis = list(title = 'Frequency')),
      legend = list(
        font = list(size = 16),  # Adjust the size 
        x = 1,
        xanchor = 'left',
        y = 0.5,
        yanchor = 'middle',
        itemsizing = 'constant',  # Scale the legend keys proportionally with the font size
        tracegroupgap = 5
      )
    )
    
    name <- paste(sample, target, device, sep = "_")
    htmlwidgets::saveWidget(fig, glue("{plotDir}/3Dplot_{name}.html"))

}







func="
size-selection -> more exons per gene 
"
size_selection_box_plot <- function(accept_df, control_df, title) {

########## Box plot 
  df <- rbind(accept_df, control_df)
  df$sample <- factor(df$sample)
  ### violin overlay box plot 
  p <- ggplot(df, aes(x = sample, y = exon_per_read, fill=sample)) +
    # geom_violin(trim=TRUE) +
    geom_boxplot(width=0.2) + 
    #scale_fill_manual(values=colors) +
    labs(x = "dataset", y = "Exon per read") +
    ggtitle(glue("Exon per on-target read {title} (biccn vs junction_pc)")) + 
    #ylim(0,6500) +
    theme_jennie(plot_title_size=10) + theme(legend.position="none")
  
  save_png(p, glue("{plotDir}/size_selection_exon_per_read_{title}"), 6, 6)
  save_pdf(p, glue("{plotDir}/size_selection_exon_per_read_{title}"), 6, 6)
}




size_selection_scatter_plot <- function(accept_df, control_df, title) {
########## Scatter plot 
		# Scatter plot: Merge by geneID (Union) # Replace NA with zero
		merged_df <- merge(accept_df, control_df, by = 'anno_geneID', all.x = TRUE, all.y = TRUE) %>%
		  replace_na(list(accept_bcFilter_unique_anno_exonblock_count = 0,
		                  control_bcFilter_unique_anno_exonblock_count = 0)) 
    print(head(merged_df))
    
    # count dots on each side of the diagonal line 
    larger_in_accept = sum(merged_df$accept_bcFilter_unique_anno_exonblock_count > merged_df$control_bcFilter_unique_anno_exonblock_count)
    larger_in_control = sum(merged_df$accept_bcFilter_unique_anno_exonblock_count < merged_df$control_bcFilter_unique_anno_exonblock_count)
    equal = sum(merged_df$accept_bcFilter_unique_anno_exonblock_count == merged_df$control_bcFilter_unique_anno_exonblock_count)
    total <- nrow(merged_df)
    
    # Calculate correlation and p-value
    cor_test <- cor.test(merged_df$control_bcFilter_unique_anno_exonblock_count, merged_df$accept_bcFilter_unique_anno_exonblock_count)
    cor_coeff <- cor_test$estimate
    p_value <- cor_test$p.value

    # get dot frequency for coloring 
		frequency_df <- merged_df %>%
		  group_by(control_bcFilter_unique_anno_exonblock_count, accept_bcFilter_unique_anno_exonblock_count) %>% 
		  summarise(freq = n(), .groups = 'drop')
	  # Merge the frequency back to the merged_df
		merged_df <- merged_df %>%
		  left_join(frequency_df, by = c("control_bcFilter_unique_anno_exonblock_count", "accept_bcFilter_unique_anno_exonblock_count"))
    # print(head(merged_df))

    # # Set x and y limits
		max_limit <- max(merged_df$control_bcFilter_unique_anno_exonblock_count, merged_df$accept_bcFilter_unique_anno_exonblock_count)

    # Rescale the proportions based on the max_freq
    colors = c("#66FFFF","#FFCCFF","#FF66FF","#993FFF")  # 浅蓝"#66FFFF", c("#FFA100", "#0004F1") 
    max_freq = max(merged_df$freq)
    proportions <- c(0, 2.9, 3, 5, 8, 23) # c(0, 0.05, 1, 4, 8, 25)
    scaled_values <- scales::rescale(proportions, to = c(0, max_freq))

		title = glue("exonPerOntargetGene_{title}")
		# Scatter plot
		p1 <- ggplot(merged_df, aes(x = control_bcFilter_unique_anno_exonblock_count, 
                                y = accept_bcFilter_unique_anno_exonblock_count, 
                                color = freq)) +
		  geom_point(size=1) +
		  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') +  # diagonal line
		  labs(
		    title = title,
        subtitle = glue("above = {larger_in_accept}, equal = {equal}, total = {total}"),
		    x = 'Unique exons per target gene in control_bcFilter',
		    y = 'Unique exons per target gene in accept_bcFilter',
		    color = 'Freq'
		  ) +
		  xlim(0, max_limit) + ylim(0, max_limit) + 
		  scale_color_gradientn(colors = colors, 
                        values = scales::rescale(scaled_values) ) + 
		  theme_jennie(plot_title_size=10) +
      annotate("text", x = 0.01 * max_limit, y = 0.99 * max_limit, label = paste("Above diagonal:", larger_in_accept), hjust = 0, vjust = 1, color = "black") +
      annotate("text", x = 0.99 * max_limit, y = 0.01 * max_limit, label = paste("Below diagonal:", larger_in_control), hjust = 1, vjust = 0, color = "black") + 
      annotate("text", x = 0.99 * max_limit, y = 0.15 * max_limit, label = paste("r =", round(cor_coeff, 2), "\np =", sprintf("%.2e", p_value)), hjust = 1, vjust = 0)

    save_png(p1, glue("{plotDir}/size_selection_{title}"), 6, 6)
    # save_pdf(p1, glue("{plotDir}/size_selection_{title}"), 6, 6)

}


# Function to create ridge plot for read lengths over time
plot_ridge_read_length <- function(exp, name, x_range=c(0, 4000)) {
  
  samples = c("accept", "control")
  plot_list <- list() 
  
  for (sample in samples) {
	  timestamp = glue("{base_path}/{exp}/timestamp/{sample}_readID_timestamp_readLength_5hTimeBin")
	  df <- read.csv(timestamp, sep="\t")
	  
	  # Set the order of 'hour_range' factor based on its unique values
	  hour_range_order <- unique(df$hour_range)
	  df$hour_range <- factor(df$hour_range, levels = hour_range_order)
    df$start <- as.numeric(str_extract(df$hour_range, "^[0-9]+"))  # Extract the start hour
    df$start <- factor(df$start, levels = sort(unique(df$start)))
    df$end <- as.numeric(str_extract(df$hour_range, "(?<=-)[0-9]+"))  # Extract the end hour
    df$end <- factor(df$end, levels = sort(unique(df$end)))
	  
	  # Define the color palette
	  #colors <- c('#F0E6E4', '#E8DBDA', '#E1D1D1', '#DAC6C8', '#D2BCBF', '#CBB5BA', '#C4B7BF', '#BEBAC3', '#B7BCC8', '#B0BECC', '#AABAC9', '#A6B4C3', '#A1ADBD', '#9CA7B7', '#98A1B1') # serene
	  colors <- c('#B5CAA0', '#ADC395', '#A5BD8B', '#9DB781', '#96B176', '#8CA86C', '#7E9B62', '#718D58', '#63804E', '#557244', '#5F7A55', '#758D74', '#8AA093', '#A0B3B2', '#B6C7D1') # pond
	  
    # Ridge plot
	  # p <- ggplot(df, aes(x = seq_len, y = hour_range, fill = hour_range)) +
		# 	    geom_density_ridges(scale = 1, alpha = 0.8, color = "black", size = 0.8) +
		# 	    scale_fill_manual(values = colors) +
		# 	    scale_x_continuous(limits = x_range) +
		# 	    theme_ridges(font_size = 12) +
		# 	    theme(axis.title.y = element_text(angle = 90, vjust = 0.5)) +
		# 	    labs(
		# 	      title = glue("Read length over time ({sample})"),
		# 	      subtitle = exp,
		# 	      x = "Read length (bp)",
		# 	      y = "Time range (hours)"
		# 	    ) +
		# 	    theme_jennie(plot_title_size=12, base_size=12) 

    p <- ggplot(df, aes(x = seq_len, y = start, fill = start) ) +
        geom_density_ridges(aes(height = after_stat(density)), rel_min_height = 0.005,
                            scale = 5, alpha=0.6, linewidth=1.1) +
        scale_fill_manual(values = colors) +
        scale_x_continuous(limits = x_range) +
        labs(title = glue("{name} ({sample})"), fill = "Hours", 
              x = "Read length (bp)", y = "Sequencing time (hours)") +
        scale_y_discrete(expand = c(0.01, 0)) + 
        theme_ridges(font_size=16, font_family = "ArialMT", grid=FALSE, center_axis_labels=TRUE) + 
        theme(
          plot.title = element_text(hjust = 0.5, face = "plain"),
          plot.background = element_rect(fill = "white"), 
          panel.background = element_rect(fill = "white"),
          panel.grid.major.x = element_blank(),  # Hide x-axis grid lines
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_line(color = "gray", size = 0.2),  # Customize y-axis grid lines
          panel.grid.minor.y = element_blank()  # Hide minor y-axis grid lines (if any)
        )

	  plot_list[[sample]] <- p
	} 
	combined_plot <- wrap_plots(plot_list, ncol = 2)
  return(combined_plot) 
}






############ NOT USED 

draft_exon_coverage_scatter_plot <- function(exp, name) {
	  parts <- strsplit(exp, "_")[[1]]
	  device<- parts[2]
	  sample <- parts[1]
	  target <- paste(parts[3], parts[4], parts[5], sep = "_")
    title <- glue("{name}")

########### exonPerOntargetGene ###########
		control_df <- read.table(glue('{base_path}/{exp}/exon_coverage/control_bcFilter_exonPerOntargetGene.tsv'), sep = '\t', header = TRUE)
		accept_df <- read.table(glue('{base_path}/{exp}/exon_coverage/accept_bcFilter_exonPerOntargetGene.tsv'), sep = '\t', header = TRUE)

		# Merge by geneID (Union) # Replace NA with zero
		merged_df <- merge(control_df, accept_df, by = 'anno_geneID', all.x = TRUE, all.y = TRUE) %>%
		  replace_na(list(control_bcFilter_unique_anno_exonblock_count = 0,
		                  accept_bcFilter_unique_anno_exonblock_count = 0)) 

    # count dots on each side of the diagonal line 
    larger_in_accept = sum(merged_df$accept_bcFilter_unique_anno_exonblock_count > merged_df$control_bcFilter_unique_anno_exonblock_count)
    larger_in_control = sum(merged_df$accept_bcFilter_unique_anno_exonblock_count < merged_df$control_bcFilter_unique_anno_exonblock_count)
    equal = sum(merged_df$accept_bcFilter_unique_anno_exonblock_count == merged_df$control_bcFilter_unique_anno_exonblock_count)
    total <- nrow(merged_df)
    
    # Calculate correlation and p-value
    cor_test <- cor.test(merged_df$control_bcFilter_unique_anno_exonblock_count, merged_df$accept_bcFilter_unique_anno_exonblock_count)
    cor_coeff <- cor_test$estimate
    p_value <- cor_test$p.value

    # get dot frequency for coloring 
		frequency_df <- merged_df %>%
		  group_by(control_bcFilter_unique_anno_exonblock_count, accept_bcFilter_unique_anno_exonblock_count) %>% 
		  summarise(freq = n(), .groups = 'drop')
	  # Merge the frequency back to the merged_df
		merged_df <- merged_df %>%
		  left_join(frequency_df, by = c("control_bcFilter_unique_anno_exonblock_count", "accept_bcFilter_unique_anno_exonblock_count"))

    # Set x and y limits
		max_limit <- max(merged_df$control_bcFilter_unique_anno_exonblock_count, merged_df$accept_bcFilter_unique_anno_exonblock_count)

    # Rescale the proportions based on the max_freq
    colors = c("#66FFFF","#FFCCFF","#FF66FF","#993FFF")  # 浅蓝"#66FFFF", c("#FFA100", "#0004F1") 
    max_freq = max(merged_df$freq)
    proportions <- c(0, 2.9, 3, 5, 8, 23) # c(0, 0.05, 1, 4, 8, 25)
    scaled_values <- scales::rescale(proportions, to = c(0, max_freq))

		# Scatter plot
		p1 <- ggplot(merged_df, aes(x = control_bcFilter_unique_anno_exonblock_count, 
                                y = accept_bcFilter_unique_anno_exonblock_count, 
                                color = freq)) +
		  geom_point(size=1) +
		  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') +  # diagonal line
		  labs(
		    title = glue("{title} \nExon per target gene"),
        subtitle = glue( "accept > control: {larger_in_accept} \n
                          accept < control: {larger_in_control} \n
                          accept = control: {equal}\n
                          total: {total}"),
		    x = 'Unique exon number per target gene in control',
		    y = 'Unique exon number per target gene in accept',
		    color = 'Freq'
		  ) +
		  xlim(0, max_limit) + ylim(0, max_limit) + 
		  scale_color_gradientn(colors = colors, 
                        values = scales::rescale(scaled_values) ) + 
		  theme_jennie(plot_title_size=14) +
      annotate("text", x = 0.99 * max_limit, y = 0.15 * max_limit, label = paste("r =", round(cor_coeff, 2), "\np =", sprintf("%.2e", p_value)), hjust = 1, vjust = 0) + 
      coord_fixed(ratio = 1) 


########### readPerOntargetExon ###########
		control_df <- read.table(glue('{base_path}/{exp}/exon_coverage/control_bcFilter_readPerOntargetExon.tsv'), sep = '\t', header = TRUE)
		accept_df <- read.table(glue('{base_path}/{exp}/exon_coverage/accept_bcFilter_readPerOntargetExon.tsv'), sep = '\t', header = TRUE)

		# Merge based on gene_id # Replace NA with zero
		merged_df <- merge(control_df, accept_df, by = 'anno_exonblock', all.x = TRUE, all.y = TRUE) %>%
		  replace_na(list(control_bcFilter_unique_readID_count = 0,
		                  accept_bcFilter_unique_readID_count = 0))
    
    print(merged_df %>%
          filter(control_bcFilter_unique_readID_count == 50 & accept_bcFilter_unique_readID_count == 0)
        ) # when 100 in control and 0 in accept. 
    # count dots on each side of the diagonal line 
    larger_in_accept = sum(merged_df$accept_bcFilter_unique_readID_count > merged_df$control_bcFilter_unique_readID_count)
    larger_in_control = sum(merged_df$accept_bcFilter_unique_readID_count < merged_df$control_bcFilter_unique_readID_count)
    equal = sum(merged_df$accept_bcFilter_unique_readID_count == merged_df$control_bcFilter_unique_readID_count)
    total <- nrow(merged_df)
    # Calculate correlation and p-value
    cor_test <- cor.test(merged_df$control_bcFilter_unique_readID_count, merged_df$accept_bcFilter_unique_readID_count)
    cor_coeff <- cor_test$estimate
    p_value <- cor_test$p.value
    
    # get dot frequency for coloring  
		frequency_df <- merged_df %>%
		  group_by(control_bcFilter_unique_readID_count, accept_bcFilter_unique_readID_count) %>% 
		  summarise(freq = n(), .groups = 'drop')
	  # Merge the frequency back to the merged DataFrame
		merged_df <- merged_df %>%
		  left_join(frequency_df, by = c("control_bcFilter_unique_readID_count", "accept_bcFilter_unique_readID_count"))
		
    # Define bin labels colors
    if (device == "minion") {
      bin_breaks <- c(-Inf, 10, 20, 50, Inf)
      bin_labels <- c("<=10", "10-20", "20-50", ">50")
    } else {
      bin_breaks <- c(-Inf, 5, 10, 50, Inf)
      bin_labels <- c("<=5", "5-10", "10-50", ">50")
    }
    bin_colors <- c("#66FFFF","#FFCCFF","#FF66FF","#993FFF")
    names(bin_colors) <- names(bin_labels)
    # bin_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")

    # Categorize 'freq' into bins
    merged_df <- merged_df %>%
      mutate(freq_bin = cut(freq, breaks = bin_breaks, labels = bin_labels))
      

    # Use the larger of the two max limits for both x and y limits
    max_limit <- max(merged_df$control_bcFilter_unique_readID_count, merged_df$accept_bcFilter_unique_readID_count)
   
    # # Rescale the proportions based on the max_freq
    # colors = c("#66FFFF","#FFCCFF","#FF66FF","#993FFF")  # 浅蓝"#66FFFF", c("#FFA100", "#0004F1") 
    # max_freq = max(merged_df$freq)
    # proportions <- c(0, 0.05, 1, 4, 8, 25) # c(0, 2.9, 3, 5, 8, 23)
    # scaled_values <- scales::rescale(proportions, to = c(0, max_freq))
    # # coeff = 0.2313782 
    # # break_label <- c(0,2,10,50,500,2000)

    # 1. Calculate the linear regression model
    fit <- lm(accept_bcFilter_unique_readID_count ~ control_bcFilter_unique_readID_count, data = merged_df)
    # 2. Extract the slope and intercept
    slope <- coef(fit)[2]
    intercept <- coef(fit)[1]
    # 3. Create the equation text
    equation_text <- glue::glue("y = {round(slope, 2)}x + {round(intercept, 2)}\nR² = {round(summary(fit)$r.squared, 2)}")


		# Scatter plot
		p2 <- ggplot(merged_df, aes(x = control_bcFilter_unique_readID_count, 
                                y = accept_bcFilter_unique_readID_count, # color = freq,
                                color=freq_bin) ) +
		  geom_point(size=1) +
		  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') +  # diagonal line
      geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "solid", size =0.5) +  # Fitted line
		  labs(
		    title = glue("{title}"),
		    x = 'Read number per exon in control',
		    y = 'Read number per exon in accepted',
        color = 'Freq',
		  ) +
		  xlim(0,max_limit) + ylim(0,max_limit) + 
      scale_color_manual(values = bin_colors) +
 		  # scale_color_gradientn(colors = colors,
                  # values = scales::rescale(scaled_values), 
                  # breaks = scales::rescale(break_label^coeff, to = c(0, max_freq)),
                  # labels = round(break_label)
                  # ) + 
      theme_jennie(plot_title_size=12, base_size=12) +
      guides(colour = guide_legend(override.aes = list(size=5)))+ 
      # annotate("text", x = 0.99 * max_limit, y = 0.85 * max_limit, 
      #      label = equation_text, 
      #      hjust = 1, vjust = 0, size = 5, color = "blue") +  # Equation text annotation
      annotate("text", x = 0.99 * max_limit, y = 0.15 * max_limit, 
                label = paste("r =", round(cor_coeff, 2), "\np =", sprintf("%.2e", p_value)), 
                hjust = 1, vjust = 0, size = 5) + # format(p_value, scientific = TRUE)
      coord_fixed(ratio = 1)  # ensures the plot area is square 

########### zoomin readPerOntargetExon ###########
    if (device == "minion") {
        max_limit_control <- calculate_max_exclude_99outliers(merged_df$control_bcFilter_unique_readID_count)
        max_limit_accept <- calculate_max_exclude_99outliers(merged_df$accept_bcFilter_unique_readID_count)
		    max_limit <- max(max_limit_control, max_limit_accept)
    }  else {
        max_limit_control <- calculate_max_excluding_outliers(merged_df$control_bcFilter_unique_readID_count)
        max_limit_accept <- calculate_max_excluding_outliers(merged_df$accept_bcFilter_unique_readID_count)
        max_limit <- max(max_limit_control, max_limit_accept)
    }

    p3 <- ggplot(merged_df, aes(x = control_bcFilter_unique_readID_count, 
                                y = accept_bcFilter_unique_readID_count, 
                                color=freq_bin)) +   # color=(freq^coeff) color=freq
      geom_point(size=1) +
      geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') +  # diagonal line
      # geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "solid", size =0.5) +  # Fitted line
      labs(
        title = glue( "accept > control: {larger_in_accept} ({sprintf('%.1f', larger_in_accept/total * 100)}%)\naccept < control: {larger_in_control} ({sprintf('%.1f', larger_in_control/total * 100)}%)\naccept = control: {equal} ({sprintf('%.1f', equal/total * 100)}%)"),
        x = 'Read number per exon in control',
        y = 'Read number per exon in accepted',
        color = 'Freq',
      ) +
      xlim(0,max_limit) + ylim(0,max_limit) + 
      scale_color_manual(values = bin_colors) +
      # scale_color_gradient(low = "lightblue", high = "darkblue") +  # Gradient from light blue to dark blue
      theme_jennie(plot_title_size=12, base_size=12) +
      guides(colour = guide_legend(override.aes = list(size=5))) +
      # annotate("text", x = 0.99 * max_limit, y = 0.1 * max_limit, 
      #      label = equation_text, 
      #      hjust = 1, vjust = 0, size = 5, color = "blue") +  # Equation text annotation
      coord_fixed(ratio = 1)  # This ensures the plot area is square
      # scale_color_gradientn(colors = colors,
      #             values = scales::rescale(scaled_values),
                  # # values = scales::rescale(break_label)^coeff,
                  # breaks = scales::rescale(break_label^coeff, to = c(0, max_freq)),
                  # labels = round(break_label)
                  # ) + 


##### A single multi-panel figure
		combined_plot <- (p2 + p3)

		name <- paste(sample, target, device, sep = "_")
		title = glue("combine_exon_coverage_{name}")
		prefix = glue("{plotDir}/{title}")
		save_png(combined_plot, prefix, 12, 6)
		save_pdf(combined_plot, prefix, 12, 6)
		cat("\n")
}
