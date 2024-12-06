#' Create an Interactive Data Table with Expandable Columns
#'
#' This function creates an interactive datatable with expandable text columns,
#' filtering capabilities, and Excel export functionality.
#'
#' @param df A data frame containing the data to display
#' @param caption_text Character string for the table caption
#' @param expandable_cols Numeric vector specifying which columns should be expandable
#'        (default: c(8, 9, 10))
#' @param round_cols Character vector of column names to round
#'        (default: c("log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
#' @param digits Integer specifying the number of decimal places for rounding (default: 3)
#'
#' @return A DT::datatable object with interactive features
#'
#' @details
#' The function implements several features:
#' \itemize{
#'   \item Expandable text columns with "More/Less" buttons
#'   \item Excel export functionality
#'   \item Column filtering
#'   \item Virtual scrolling for performance
#'   \item Automatic number formatting
#' }
#'
#' @import DT
#' @import htmlwidgets
#'
#' @examples
#' df <- data.frame(
#'   gene = c("BRCA1", "TP53"),
#'   log2FoldChange = c(1.23, -0.85),
#'   pvalue = c(0.001, 0.003)
#' )
#' create_interactive_datatable(df, "Gene Expression Results")
#'
#' @export

create_interactive_datatable <- function(df, caption_text, expandable_cols = c(8, 9, 10),
                                         only_sig = FALSE,
                                         round_cols = c("log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
                                         digits = 3) {
  # Limit initial rows shown
  initial_rows <- 5

  if (only_sig) {
    df <- df %>% filter(!is.na(padj) & !is.na(log2FoldChange))
  }
  
  # Create column definitions for expandable text
  column_defs <- lapply(expandable_cols - 1, function(col) {
    list(
      targets = col,
      render = DT::JS(
        "function(data, type, row, meta) {
          if (type === 'display' && data != null && data.length > 20) {
            return '<div class=\"cell-content\" style=\"max-height:20px; overflow:hidden;\">' +
                   data.substring(0, 20) + '...</div><button class=\"btn btn-default btn-xs show-more\">(More)</button>';
          }
          return data;
        }"
      )
    )
  })

  # Create the base datatable with optimizations
  dt <- DT::datatable(
    df,
    caption = caption_text,
    filter = "top",
    extensions = c('Buttons', 'Scroller'),  # Add Scroller extension
    options = list(
      pageLength = initial_rows,
      scrollX = TRUE,
      scrollY = 300,          # Fixed height
      scroller = TRUE,        # Enable virtual scrolling
      deferRender = TRUE,     # Defer rendering until needed
      columnDefs = column_defs,
      dom = 'Bfrtip',
      buttons = list(
        list(
          extend = 'excel',
          text = 'Download Excel',
          filename = paste0(caption_text, '_data')
        )
      ),
      searchDelay = 1000     # Delay search to reduce server load
    )
  )

  # Add number formatting if round_cols are present
  present_cols <- round_cols[round_cols %in% names(df)]
  if (length(present_cols) > 0) {
    dt <- dt %>% formatRound(columns = present_cols, digits = digits)
  }

  # Add the show more/less functionality
  dt <- dt %>%
    htmlwidgets::onRender("
      function(el) {
        el.querySelectorAll('.show-more').forEach(btn => {
          btn.addEventListener('click', function(e) {
            const content = this.previousElementSibling;
            if (content.style.maxHeight) {
              content.style.maxHeight = '';
              this.textContent = '(Less)';
            } else {
              content.style.maxHeight = '20px';
              this.textContent = '(More)';
            }
            e.stopPropagation();
          });
        });
      }")

  return(dt)
}