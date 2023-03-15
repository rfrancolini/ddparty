#' retreive example type quantalife edited file name
#'
#' @export
#' @return filename
ManuallyEditedData_example <- function(){
  system.file("exampledata/SampleDataManuallyEdited.csv",
              package="ddparty")
}

#' retreive example type quantalife raw file name
#'
#' @export
#' @return filename
RawData_example <- function(){
  system.file("exampledata/SampleDataRaw.csv",
              package="ddparty")
}



#' retreive example metadata
#'
#' @export
#' @return filename
metadata_example <- function(){
  system.file("exampledata/SampleMetadata.csv",
              package="ddparty")
}


#' read manually edited ddpcr data file, process it, write as a new file
#'
#' @export
#' @param filename character, the name of the file
#' @param metadata character, the name of the sample metadata file
#' @param target character, the name of the target of interest if only want to work on one
#' @param calltype character, if "auto", indicates the data was not manipulated in QuantaLife program. This is the default setting. If "manual", indicates some samples were manually called and csv exported
#' @param output character, the name for the outputted QAQC file
#' @return tibble
read_ddpcr <- function(filename = ManuallyEditedData_example(),
                       metadata = metadata_example(),
                       target = NA,
                       calltype = c("auto", "manual")[2],
                       output = NA){

  stopifnot(inherits(filename, "character"))
  stopifnot(file.exists(filename[1]))

  #read in file
  x <- suppressMessages(readr::read_csv(filename[1]))

  #filter to only work with target of interest
  if (!is.na(target)) {
    x <- x %>% dplyr::filter(.data$Target == target)
  }

  ####Split into case situation - raw vs. man data

  #only select columns of interest
  #mu needs to be used for two column names, "u" is incorrect
  if(tolower(calltype [1]) == "manual") {
  x <- x %>% dplyr::select("Well",
                            "Sample",
                            "Status",
                            "Conc(copies/\u00b5L)",
                            "Copies/20\u00b5LWell",
                            "PoissonConfMax",
                            "PoissonConfMin",
                            "Positives",
                            "Negatives",
                            "Accepted Droplets",
                            "Threshold1",
                            "MeanAmplitudeOfPositives",
                            "MeanAmplitudeOfNegatives",
                            "MeanAmplitudeTotal",
                            "PoissonConfidenceMax68",
                            "PoissonConfidenceMin68")

  #fix names of columns so they can match to rawdata output
  x <- x %>% dplyr::rename(Concentration = "Conc(copies/\u00b5L)") %>%
             dplyr::rename(CopiesPer20uLWell = "Copies/20\u00b5LWell") %>%
             dplyr::rename(AcceptedDroplets = "Accepted Droplets") %>%
             dplyr::rename(Threshold = "Threshold1") %>%
             dplyr::rename(PoissonConfMax68 = "PoissonConfidenceMax68") %>%
             dplyr::rename(PoissonConfMin68 = "PoissonConfidenceMin68") %>%
             dplyr::rename(MeanAmplitudeofPositives = "MeanAmplitudeOfPositives") %>%
             dplyr::rename(MeanAmplitudeofNegatives = "MeanAmplitudeOfNegatives")

  } #end of "manual" column processing
  else if (tolower(calltype [1]) == "auto") {
    x <- x %>% dplyr::select("Well",
                             "Sample",
                             "Status",
                             "Concentration",
                             "CopiesPer20uLWell",
                             "PoissonConfMax",
                             "PoissonConfMin",
                             "Positives",
                             "Negatives",
                             "AcceptedDroplets",
                             "Threshold",
                             "MeanAmplitudeofPositives",
                             "MeanAmplitudeofNegatives",
                             "MeanAmplitudeTotal",
                             "PoissonConfMax68",
                             "PoissonConfMin68")
  } #end of auto column processing
  else {
    stop("options for calltype are auto or manual. what is ", calltype, "?")
  }

  ##Make Concentration column numeric - turn "No Call" into NA then to numeric
  x <- x %>% dplyr::mutate(ConcentrationNum = .data$Concentration)
  x$ConcentrationNum <- dplyr::na_if(x$ConcentrationNum,"No Call")
  x$ConcentrationNum <- as.numeric(x$ConcentrationNum)


  ##Calculate copies/uL based on our input of 3.3ul into 22ul rxn(also equal to conc*3??)
  x <- x %>% dplyr::mutate(CalcCopies = (.data$ConcentrationNum*22)/3.3)


  #Add metadata columns to Data for Viewing
  metadata <- suppressMessages(readr::read_csv(metadata[1]))

  x_meta <- dplyr::left_join(x,metadata, by = "Sample")

  #fulldataman.noNA <- na.omit(fulldataman)

  if (!is.na(output)) {
    readr::write_csv(x_meta, file = output) }

  return(x_meta)

}





#' Take in processed csv file and calculate mean gene copy number per site per visit
#'
#' @export
#' @param inputdata tibble, ddpcr data processed and formatted with this package
#' @param n numeric, number of samples per site, default 6
#' @param negID character, partial string that identifies samples as negative control, these will be removed in this case, default "G"
#' @param output character, the name for the outputted file with mean values calculated
#' @return tibble
mean_ddpcr <- function(inputdata,
                       n = 6,
                       negID = "G",
                       output = NA){


  #only look at calculated gene copy number
  x <- inputdata %>%
        dplyr::select("Sample", "CalcCopies", "season", "site") %>%
        tidyr::drop_na("site") %>%
        dplyr::mutate_all(~replace_na(.,0)) #replace with 0 so std/serr can be calculated


  ##Remove negative samples
  ##currently using partial string as given by user
  ##create smaller dataframe with all samples that are negatives
  negs <- x[stringr::str_detect(x$Sample, negID), ]
  x <- dplyr::anti_join(x, negs, by = "Sample")

  #mean of the samples per site
  x_mean <- x %>%
    dplyr::group_by(.data$season, .data$site) %>%
    dplyr::mutate(Sum = sum(.data$CalcCopies),
                  Mean = .data$Sum/n,
                  Sdev = sd(.data$CalcCopies),
                  Serr = sd(.data$CalcCopies)/sqrt(n))%>%
    dplyr::select(-"CalcCopies", -"Sample") %>%
    dplyr::distinct() %>%
    dplyr::ungroup()

  if (!is.na(output)) {
          readr::write_csv(x_mean, file = output) }

   return(x_mean)
}
