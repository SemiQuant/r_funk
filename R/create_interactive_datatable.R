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
#' @param show_table Logical indicating whether to display the table (default: TRUE).
#'        If FALSE, only shows the Excel download button.
#'
#' @return A DT::datatable object with interactive features
#'
#' @details
#' The function implements several features:
#' \itemize{
#'   \item Expandable text columns with "More/Less" buttons
#'   \item Excel export functionality
#'   \item Column filtering (when show_table = TRUE)
#'   \item Virtual scrolling for performance (when show_table = TRUE)
#'   \item Automatic number formatting
#' }
#'
#' @import DT
#' @import htmlwidgets
#' @import writexl
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
                                         digits = 3,
                                         show_table = TRUE) {
  if (only_sig) {
    df <- df %>% filter(!is.na(padj) & !is.na(log2FoldChange))
  }

  if (show_table) {
    # Limit initial rows shown
    initial_rows <- 5
    
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
      extensions = c('Buttons', 'Scroller'),
      options = list(
        pageLength = initial_rows,
        scrollX = TRUE,
        scrollY = 300,
        scroller = TRUE,
        deferRender = TRUE,
        columnDefs = column_defs,
        dom = 'Bfrtip',
        buttons = list(
          list(
            extend = 'excel',
            text = 'Download Excel',
            filename = paste0(caption_text, '_data')
          )
        ),
        searchDelay = 1000
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
    
  } else {
    # First create and save the Excel file with the actual data
    filename <- paste0(caption_text, '_data.xlsx')
    writexl::write_xlsx(df, filename)
    
    # Then create minimal table with only download button
    dt <- DT::datatable(
      data.frame("Click 'Download Excel' to download the data" = character(0)),  # empty dataframe with better message
      extensions = 'Buttons',
      options = list(
        dom = 'B',  # Show only buttons
        buttons = list(
          list(
            extend = 'excel',
            text = 'Download Excel',
            filename = paste0(caption_text, '_data'),
            action = htmlwidgets::JS(sprintf(
              "function(e, dt, button, config) {
                window.location.href = '%s';
              }", filename))
          )
        ),
        paging = FALSE,     # Disable paging
        info = FALSE,       # Disable info
        searching = FALSE,  # Disable search
        language = list(
          emptyTable = "Table not shown - use download button above to download the complete dataset"
        )
      ),
      selection = 'none',
      filter = 'none',
      rownames = FALSE,
      colnames = ''
    )
    
    return(dt)
  }
}
