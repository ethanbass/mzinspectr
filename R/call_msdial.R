#' Call MS-DIAL console app
#' For help configuring the MSDIAL console app on mac OSX or linux, please see
#' the \href{https://github.com/Jiung-Wen/msdial}{instructions} helpfully
#' compiled by \href{https://github.com/Jiung-Wen/}{Jiung-Wen Chen}.
#' @param system Either \code{gcms}, \code{lcmsdda}, \code{lcsdia},
#' \code{lcmsdia}, \code{lcimmsdda}, or \code{lcimmsdia}.
#' @param path_in Path to files.
#' @param path_out Path to output directory
#' @param method A method file.
#' @param settings Settings in lieu of a method file.
#' @param p Logical.
#' @param mce Logical
#' @return Returns MSDIAL alignment.
#' @export

ms_call_msdial <- function(system, path_in, path_out, method, settings, p=FALSE, mce=FALSE){
  system <- match.arg(system, c("gcms", "lcmsdda", "lcsdia", "lcmsdia", "lcimmsdda", "lcimmsdia"))
  exists_msdial <- configure_msdial(check=TRUE)
  ms_check_files(path_in)
  if (!exists_msdial){
    stop()
  }
  msdial_path <- readLines(system.file("shell/msdial_path.txt", package = 'msdialreadr'))
  if (missing(path_in)){
    stop("Path to directory of path_in files must be provided.")
  }
  if (missing(path_out)){
    stop("Path to path_out directory must be provided.")
  }
  if (missing(method) & missing(settings)){
    stop("Please provide an argument to method or settings")
  }
  if(!dir.exists(path_in)){
    stop("Input directory not found.")
  }
  if(!fs::dir_exists(path_out)){
    warning("Output directory not found.", immediate.=TRUE)
    ans <- readline("Create directory? (y/n)")
    if (ans == "y"){
      fs::dir_create(path_out)
      if (!fs::dir_exists(path_out)){
        stop("Directory could not be created.")
      }
    }
  }
  p <- ifelse(p, "-p", "")
  mce <- ifelse(mce, "-mce", "")
  command <- paste(msdial_path, system, "-i", path_in, "-o", path_out, "-m", method, p, mce)
  system(command)
  fout <- list.files(path_out,full.names = TRUE)
  alignments <- grep("AlignResult", fout)
  if (length(alignments) > 1){
    idx<-order(as.POSIXct(file.info(fout[alignments])$ctime), decreasing = TRUE)[1]
    alignments <- alignments[idx]
  }
  ms_read_alignment(fout[alignments])
}


#' Configure shell script to call MS-DIAL console app
#'
#' @name configure_msdial
#' @param check Logical.
#' @return No return value.
#' @author Ethan Bass
#' @noRd
configure_msdial <- function(check = FALSE){
  path_msdial <- readLines(system.file("shell/msdial_path.txt", package = 'msdialreadr'))
  exists <- file.exists(path_msdial)
  if (!exists){
    path_msdial <- readline(prompt="Please provide path to `MsdialConsoleApp`):")
    if (.Platform$OS.type == "windows"){
      path_msdial <- gsub("/", "\\\\", path_msdial)
    }
    writeLines(path_msdial, con=system.file('shell/msdial_path.txt', package='msdialreadr'))
  }
  if (check){
    exists
  }
}

#' @noRd
ms_check_files <- function(path_in){
  pattern="abf|cdf"
  files <- list.files(path_in, pattern = pattern, ignore.case = TRUE, full.names = TRUE)
  info <- file.info(files)
  idx <- which(info$size < 1000)
  if (length(idx > 1)){
    warning("Some of the files are less than 1 kb suggesting that they do not contain any valid data.",immediate.=TRUE)
  }
  ans <- readline(paste("Delete", length(idx), "empty files? (y/n)? (This will permanently remove the files from the path_in directory)."))
  if (ans == "y"){
    fs::file_delete(files[idx])
  }
}
