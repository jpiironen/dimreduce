.onAttach <- function(...) {
  ver <- utils::packageVersion("dimreduce")
  packageStartupMessage("This is dimreduce version ", ver)
}