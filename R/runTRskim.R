#' Launch Shiny App for TRskim
#'
#' This functions launches the Shiny app for TRskim. The app is intended for
#' those with less coding experience and performs the entire TRskim workflow.
#'
#' @return launches the shiny app
#'
#' @examples
#' \dontrun{
#' TRskim::runTRskim()
#' }
#'
#' @export runTRskim
#' @importFrom shiny runApp

runTRskim <- function(){
  app_loc <- system.file("shiny-scripts", package = "TRskim")
  shiny_launch <- shiny::runApp(app_loc, display.mode = "normal")

  return(shiny_launch)
}

#[END]
