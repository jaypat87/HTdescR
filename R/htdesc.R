#' HTdesc
#'
#' Get Hamett-Taft descriptor value for a single chemical fragment
#'
#' Get values of Hammett-Taft descriptors for a given chemical fragment in SMILES string format by iterating through a lookup table.
#' Returns values from an exact matching fragment from the library or uses a substituted fragment.
#' In case an exact match isnt found, this function uses a mismatch-tolerant maximum common substructure (fMCS) based fragment substitution library to get the HT descriptors with highest tanimoto coefficient.
#'
#'
#' @param smile SMILES string for a chemical fragment in character, factor, or SMIset datatype.
#' @param HT.type The type of Hammett-Taft (HT) descriptor; valid inputs include "taft", "meta", "para", "ortho", "inductive", "es" and "user" for a user defined sigma descriptor.
#' @param sigma.selection \itemize{The type of sigma to be returned; valid inputs include "A", "B", "C", "D", "E", "F", "G", and "H":
#' \item A = Avgerage value
#' \item B = Hansch preferred value if available then EPIsuite value if available then Average value
#' \item C = EPIsuite value if available then Hansch preferred value if available then Average of distinct values
#' \item D = EPIsuite value if available then Average of distinct value
#' \item E = Hansch preferred value if available then Average of distinct value
#' \item F = Hansch preferred value if available then mode value if available then median value
#' \item G = EPIsuite value if available then Hansch preferred value if available then median value
#' \item H = Mode if available then Average of distinct value
#' \item Avg.dist = Average of Distinct value
#' \item Med = Median Value
#' }
#' @param ... arguments from the fmcsR::fmcsBatch such as al, au, bl, and bu
#'
#' @return A list containing the tanimoto coefficient for the closest matching MCS, SMILES string of the maximum common substructure (MCS), and index of the matched fragment from the library.
#'
#' @import ChemmineR
#' @import methods
#' @import fmcsR
#'
#' @export
#'
#' @examples
#' htdesc (smile = "CCC", HT.type = "taft", sigma.selection = "A")
#'
htdesc <- function(smile, HT.type = "taft", sigma.selection = "A", ...) {
  if (class(smile) == "SMIset") {

    #Convert SMILES object to SDF
    sampleSDF <- ChemmineR::smiles2sdf(smile)

  } else if (class(smile) == "character") {

      #smile <- methods::as (smile, "SMIset")

      #Convert SMILES character string to SDF
      sampleSDF <- ChemmineR::smiles2sdf(smile)
      # add escape character fix
      # possibly can use stringr and use '\\\\' for regex match and change, but note that this is not a real problem, since the inputcsvfiles will have smiles which are already forced with double \\ characters by R itself when it read it

  } else if (class(smile) == "factor") {
      #Convert SMILES as factor to SDF
      smile <- as.character (smile)
      smile <- methods::as (smile, "SMIset")
      sampleSDF <- ChemmineR::smiles2sdf(smile)

  } else {
    stop ("Please enter a valid SMIset object, character string, or factor as SMILES.")
  }

  if (HT.type == "meta") {

    #Initialize SDF libraries and lookup tables for searching
    sigmalibrarySDF <- metaSDF
    siglookuptable <- metasigma
  } else if (HT.type == "para")  {
      sigmalibrarySDF <- paraSDF
      siglookuptable <- parasigma
  } else if (HT.type == "user")  {
      sigmalibrarySDF <- userSDF
      siglookuptable <- usersigma
  } else if (HT.type == "ortho") {
      sigmalibrarySDF <- orthoSDF
      siglookuptable <- orthosigma
  } else if (HT.type == "taft")  {
      sigmalibrarySDF <- taftSDF
      siglookuptable <- taftsigma
  } else if (HT.type == "es")    {
      sigmalibrarySDF <- esSDF
      siglookuptable <- es
  } else if (HT.type == "inductive")    {
      sigmalibrarySDF <- indSDF
      siglookuptable <- indsigma
  } else {
      stop("Specify valid HT.type")
  }

  #Conduct search
  fmcsoutput <- fmcsR::fmcsBatch(sampleSDF[1], sigmalibrarySDF)
  fmcsoutputframe <- data.frame(fmcsoutput, siglookuptable)

  #Order output to find highest Tanimoto coefficient
  fmcsoutputframe <- fmcsoutputframe[order(-fmcsoutputframe$Tanimoto_Coefficient),]

  #Return sigma value based on user's selection of preferred sigma
  if (sigma.selection =="A") {
    # A: reg.avg
    returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$reg.avg[1])
    return(returnlist)

  } else if (sigma.selection == "B") {
      #B: priority order: hansch preferred > epi.value > reg.avg
      if (is.na(fmcsoutputframe$hansch.pref[1]) == FALSE) {
        returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$hansch.pref[1])
        return(returnlist)
      } else if (is.na(fmcsoutputframe$epi.value[1]) == FALSE) {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$epi.value[1])
          return(returnlist)
      } else {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$reg.avg[1])
          return(returnlist)
      }
  } else if (sigma.selection == "C") {
      #C: priority order: epi.value > hansch preferred > avg.dist > reg.avg
      if (is.na(fmcsoutputframe$epi.value[1]) == FALSE) {
        returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$epi.value[1])
        return(returnlist)
      } else if (is.na(fmcsoutputframe$hansch.pref[1]) == FALSE) {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$hansch.pref[1])
          return(returnlist)
      } else if (is.na(fmcsoutputframe$avg.dist[1]) == FALSE) {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$avg.dist[1])
          return(returnlist)
      } else {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$reg.avg[1])
          return(returnlist)
      }
  } else if (sigma.selection == "D") {
      #D: priority order: epi.value > avg.dist > reg.avg
      if (is.na(fmcsoutputframe$epi.value[1]) == FALSE) {
        returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$epi.value[1])
        return(returnlist)
      } else {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$avg.dist[1])
          return(returnlist)
      }
  } else if (sigma.selection == "E") {
    #E: priority order: hansch preffered > avg.dist
      if (is.na(fmcsoutputframe$hansch.pref[1]) == FALSE) {
        returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$hansch.pref[1])
        return(returnlist)
      } else {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$avg.dist[1])
          return(returnlist)
      }
  } else if (sigma.selection == "F") {
    #F: priority order: Hansch preferred > regular mode > regular median
      if (is.na(fmcsoutputframe$hansch.pref[1]) == FALSE) {
        returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$hansch.pref[1])
        return(returnlist)
      } else if (is.na(fmcsoutputframe$reg.mode[1]) == FALSE) {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$reg.mode[1])
          return(returnlist)
      } else {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$reg.median[1])
          return(returnlist)
      }
  } else if (sigma.selection =="G")  {
    #G: priority order: epi.value > hansch preffered > median
      if (is.na(fmcsoutputframe$epi.value[1]) == FALSE) {
        returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$epi.value[1])
        return(returnlist)
      } else if (is.na(fmcsoutputframe$hansch.pref[1]) == FALSE) {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$hansch.pref[1])
          return(returnlist)
      } else {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$reg.median[1])
          return(returnlist)
      }
  } else if (sigma.selection == "H") {
    #H:  mode > avg.dist
      if (is.na(fmcsoutputframe$reg.mode[1]) == FALSE) {
        returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$reg.mode[1])
        return(returnlist)
      } else {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$avg.dist[1])
          return(returnlist)
      }
  } else if (sigma.selection == "Avg.dist") {
    #Avg.dist:  average of distinct
    if (is.na(fmcsoutputframe$reg.mode[1]) == FALSE) {
      returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$avg.dist[1])
      return(returnlist)
    }
  } else if (sigma.selection == "Med") {
    #Med:  Median
    if (is.na(fmcsoutputframe$reg.mode[1]) == FALSE) {
      returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = as.character(fmcsoutputframe$fragments[1]), value = fmcsoutputframe$avg.dist[1])
      return(returnlist)
    }
  }
    else {
      stop ("Specify valid sigma.selection")
    }
}
