#' Internal utilities and namespace imports
#'
#' These imports ensure that functions used non-explicitly (e.g., in ggplot2) are available during package checks.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_step geom_hline ggtitle xlab ylab xlim ylim theme theme_minimal element_blank element_rect element_text scale_fill_discrete ggplotGrob geom_abline scale_color_manual
#' @importFrom grid grid.newpage grid.draw gpar
#' @importFrom gridExtra arrangeGrob
#' @importFrom DSI datashield.connections_find
#' @importFrom stats approx as.formula var
#' @importFrom rlang sym
#' @importFrom methods is
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  "Std.Mean.Diff", "type", "var", "distance", "domain", "ecdf",
  "qq_control", "qq_treated"
))

