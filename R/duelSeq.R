require(tidyverse)
require(plotly)
require(htmlwidgets)

#' Analyze and visualize dual-genome sequencing alignment data
#' 
#' This function analyzes dual-genome sequencing alignment data from MultiQC samtools flagstat output.
#' It generates visualizations showing the distribution of reads between two genomes (e.g., host and pathogen)
#' and calculates enrichment relative to positive controls. The function produces two plots:
#' 1. Alignment distribution showing fraction of reads mapping to each genome
#' 2. Enrichment analysis relative to positive control
#' 
#' @param flagstat_file Path to the MultiQC samtools flagstat file
#' @param primary_genome Name of the primary genome (e.g., "human")
#' @param secondary_genome Name of the secondary genome (e.g., "mtb")
#' @param control_sample Name or pattern of the control sample (e.g., "posControl")
#' @param exclude_samples Vector of sample names or patterns to exclude from analysis (e.g., c("negControl", "blank"))
#' @param output_dir Directory to save the plots (optional)
#' @param save_plots Whether to save the plots to files (default: FALSE)
#' @return A list containing:
#'   \item{alignment_summary}{Plotly object showing alignment distribution}
#'   \item{enrichment}{Plotly object showing enrichment analysis}
#'   \item{combined_plot}{Plotly object with both plots combined}
#'   \item{processed_data}{List containing processed alignment and enrichment data}
#' @export
#' @importFrom tidyverse %>%
#' @importFrom plotly plot_ly layout subplot
#' @importFrom htmlwidgets saveWidget
analyze_dual_seq <- function(flagstat_file, 
                           primary_genome = "human",
                           secondary_genome = "mtb",
                           control_sample = "posControl",
                           exclude_samples = NULL,
                           output_dir = NULL,
                           save_plots = FALSE) {
  
  # Input validation
  if (!file.exists(flagstat_file)) {
    stop("Flagstat file does not exist: ", flagstat_file)
  }
  
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Read and process data
  dat <- read_tsv(flagstat_file)
  
  dat <- transform(dat, 
                  Sample_name = str_extract(Sample, ".*(?=_[^_]*$)"),
                  Genome = str_extract(Sample, "(?<=_)[^_]*$")
  )
  
  # Exclude samples if specified
  if (!is.null(exclude_samples)) {
    exclude_pattern <- paste(exclude_samples, collapse = "|")
    dat <- dat %>%
      filter(!grepl(exclude_pattern, Sample_name, ignore.case = TRUE))
  }
  
  dat <- dat %>% 
    filter(Genome %in% c(primary_genome, secondary_genome)) %>%
    select(Sample_name, mapped_passed, Genome) %>% 
    group_by(Sample_name) %>% 
    mutate(total_aligned = sum(mapped_passed))
  
  dat <- dat %>%
    mutate(Grp = str_remove(Sample_name, "_rep.*"))
  
  dat <- dat %>% 
    group_by(Grp, Genome) %>% 
    summarize(
      mean_passed = mean(mapped_passed, na.rm = T),
      sd_passed = sd(mapped_passed, na.rm = T),
      .groups = "drop"
    ) %>% 
    group_by(Grp) %>% 
    mutate(total_passed = sum(mean_passed),
           fraction = mean_passed/total_passed) %>%
    ungroup()
  
  dat <- dat %>%
    mutate(Genome = factor(Genome, levels = c(secondary_genome, primary_genome)))
  
  # Create alignment summary plot
  p1 <- dat %>% 
    plot_ly(x = ~Grp, y = ~fraction, type = 'bar', color = ~Genome) %>% 
    layout(
      barmode = 'stack',
      title = list(
        text = paste(secondary_genome, "vs", primary_genome, "Alignment Distribution"),
        font = list(size = 16)
      ),
      xaxis = list(
        title = "Sample Group",
        tickangle = 45
      ),
      yaxis = list(
        title = list(
          text = "Fraction of Aligned Reads",
          standoff = 5
        ),
        range = c(0, 1)
      ),
      showlegend = TRUE,
      margin = list(l = 100)
    )
  
  # Process secondary genome data
  dat2 <- dat %>% 
    filter(Genome == secondary_genome) %>% 
    filter(!grepl(paste0(secondary_genome, "Contro"), Grp, ignore.case = TRUE)) %>% 
    select(Grp, fraction)
  
  # Check if control sample exists in the data
  if (!any(grepl(control_sample, dat2$Grp, ignore.case = TRUE))) {
    stop(sprintf("Control sample '%s' not found in the data", control_sample))
  }
  
  dat2 <- dat2 %>% 
    mutate(is_control = grepl(control_sample, Grp, ignore.case = TRUE))
  
  dat2 <- dat2 %>% 
    mutate(Inc = fraction / (first(fraction[is_control == TRUE])))
  
  # Create enrichment plot
  p2 <- dat2 %>% 
    arrange(desc(Inc)) %>%
    plot_ly(x = ~reorder(Grp, -Inc), y = ~Inc, type = 'bar') %>% 
    layout(
      title = list(
        text = paste(secondary_genome, "Enrichment Relative to", control_sample),
        font = list(size = 16)
      ),
      xaxis = list(
        title = "Sample Group",
        tickangle = 45
      ),
      yaxis = list(
        title = list(
          text = "Enrichment (Fold Change)",
          standoff = 5
        ),
        rangemode = "tozero"
      ),
      showlegend = FALSE,
      margin = list(l = 100)
    )
  
  # Create combined plot
  p_combined <- subplot(
    p1 %>% layout(
      xaxis = list(title = ""),
      yaxis = list(
        title = list(
          text = "Fraction of Aligned Reads",
          standoff = 5
        ),
        range = c(0, 1)
      ),
      margin = list(l = 100)
    ),
    p2 %>% layout(
      xaxis = list(title = "Sample Group"),
      yaxis = list(
        title = list(
          text = "Enrichment (Fold Change)",
          standoff = 5
        ),
        rangemode = "tozero"
      ),
      margin = list(l = 100)
    ),
    nrows = 2,
    heights = c(0.5, 0.5),
    shareX = TRUE,
    titleY = TRUE
  ) %>%
    layout(
      title = list(
        text = paste("Dual-Genome Sequencing Analysis:", secondary_genome, "vs", primary_genome),
        font = list(size = 18)
      ),
      showlegend = TRUE,
      margin = list(t = 100, l = 100),
      grid = list(
        rows = 2,
        columns = 1,
        pattern = "independent"
      )
    )
  
  # Save plots if requested
  if (save_plots && !is.null(output_dir)) {
    htmlwidgets::saveWidget(p1, file.path(output_dir, "alignment_summary.html"))
    htmlwidgets::saveWidget(p2, file.path(output_dir, "enrichment.html"))
    htmlwidgets::saveWidget(p_combined, file.path(output_dir, "combined_plots.html"))
  }
  
  # Return results
  results <- list(
    alignment_summary = p1,
    enrichment = p2,
    combined_plot = p_combined,
    processed_data = list(
      alignment_data = dat,
      enrichment_data = dat2
    )
  )
  
  return(results)
}

