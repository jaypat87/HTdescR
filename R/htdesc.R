# A R function to return sigma value from an exact matching fragment from the library or using substituted fragment.
# devtools::use_package("fmcsR")

#' Get Hamett-Taft descriptor for a single chemical fragment
#'
#' Get values of Hammett-Taft descriptors for a given chemical fragment in SMILES string format by iterating through a lookup table. In case an exact match isnt found, this function uses a mismatch tolerant
#' maximum common substructure (fMCS) based fragment substitution library to get the HT dexcriptors with highest tanimoto coefficient.
#'
#'
#' @param smile SMILES string for a chemical fragment in character, factor, or SMIset datatype.
#' @param HT.type The type of Hammett-Taft (HT) descriptor; valid inputs include "taft", "meta", "para", "ortho", "induction", "es" and "user" for a user defined sigma descriptor.
#' @param ... arguments from the fmcsR::fmcbatch such as a and b
#' @return A list containing tanimoto coefficient for the closest matching MCS, SMILES string of the MCS, and index nuber of the matched fragment from the library.
#'
#' @examples
#' htdesc (smile = "CCC", HT.type = "taft", sigma.selection = "A")
#'
htdesc <- function(smile, HT.type = "taft", sigma.selection = "A", ...) {
  if (class(smile) == "SMIset") {

    sampleSDF <- ChemmineR::smiles2sdf(smile)

  } else if (class(smile) == "character") {

    #smile <- methods::as (smile, "SMIset")
      sampleSDF <- ChemmineR::smiles2sdf(smile)
# add escape character fix
# possibly can use stringr and use '\\\\' for regex match and change, but note that this is not a real problem, since the inputcsvfiles will have smiles which are already forced with double \\ characters by R itself when it read it
  } else if (class(smile) == "factor") {
      smile <- as.character (smile)
      smile <- methods::as (smile, "SMIset")
      sampleSDF <- ChemmineR::smiles2sdf(smile)
  } else {
    stop ("please enter a valid SMIset object, characters or factors as SMILES")
  }

  if (HT.type == "meta") {

    sigmalibrarySDF <- HTdescR::metaSDF
    siglookuptable <- HTdescR::metasigma
  } else if (HT.type == "para")  {
      sigmalibrarySDF <- HTdescR::paraSDF
      siglookuptable <- HTdescR::parasigma
  } else if (HT.type == "user")  {
      sigmalibrarySDF <- HTdescR::userSDF
      siglookuptable <- HTdescR::usersigma
  } else if (HT.type == "ortho") {
      sigmalibrarySDF <- HTdescR::orthoSDF
      siglookuptable <- HTdescR::orthosigma
  } else if (HT.type == "taft")  {
      sigmalibrarySDF <- HTdescR::taftSDF
      siglookuptable <- HTdescR::taftsigma
  } else if (HT.type == "es")    {
      sigmalibrarySDF <- HTdescR::esSDF
      siglookuptable <- HTdescR::es
  } else if (HT.type == "induction")    {
      sigmalibrarySDF <- HTdescR::indSDF
      siglookuptable <- HTdescR::indsigma
  } else {
      stop("Specify valid HT.type")
  }


  fmcsoutput <- fmcsR::fmcsBatch(sampleSDF[1], sigmalibrarySDF)
  fmcsoutputframe <- data.frame(fmcsoutput, siglookuptable)
  fmcsoutputframe <- fmcsoutputframe[order(-fmcsoutputframe$Tanimoto_Coefficient),]
  if (sigma.selection =="A") {
    # A: reg.avg
    returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = fmcsoutputframe$fragments[1], value = fmcsoutputframe$reg.avg[1])
    return(returnlist)

  } else if (sigma.selection == "B") {
      #B: priority order: hansch preferred > epi.value > reg.avg
      if (is.na(fmcsoutputframe$hansch.pref[1]) == FALSE) {
        returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = fmcsoutputframe$fragments[1], value = fmcsoutputframe$hansch.pref[1])
        return(returnlist)
      } else if (is.na(fmcsoutputframe$epi.value[1]) == FALSE) {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = fmcsoutputframe$fragments[1], value = fmcsoutputframe$epi.value[1])
          return(returnlist)
      } else {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = fmcsoutputframe$fragments[1], value = fmcsoutputframe$reg.avg[1])
          return(returnlist)
      }
  } else if (sigma.selection == "C") {
      #C: priority order: epi.value > hansch preferred > avg.dist > reg.avg
      if (is.na(fmcsoutputframe$epi.value[1]) == FALSE) {
        returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = fmcsoutputframe$fragments[1], value = fmcsoutputframe$epi.value[1])
        return(returnlist)
      } else if (is.na(fmcsoutputframe$hansch.pref[1]) == FALSE) {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = fmcsoutputframe$fragments[1], value = fmcsoutputframe$hansch.pref[1])
          return(returnlist)
      } else if (is.na(fmcsoutputframe$avg.dist[1]) == FALSE) {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = fmcsoutputframe$fragments[1], value = fmcsoutputframe$avg.dist[1])
          return(returnlist)
      } else {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = fmcsoutputframe$fragments[1], value = fmcsoutputframe$reg.avg[1])
          return(returnlist)
      }
  } else if (sigma.selection == "D") {
      #D: priority order: epi.value > avg.dist
      if (is.na(fmcsoutputframe$epi.value[1]) == FALSE) {
        returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = fmcsoutputframe$fragments[1], value = fmcsoutputframe$epi.value[1])
        return(returnlist)
      } else {
          returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = fmcsoutputframe$fragments[1], value = fmcsoutputframe$avg.dist[1])
          return(returnlist)
      }
  } else if (sigma.selection == "E") {
    #E: Hansch preffered first plus distinct avg for when hydrowin not available
    returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = fmcsoutputframe$fragments[1], value = fmcsoutputframe$E[1])
    return(returnlist)
  } else if (sigma.selection == "F") {
    #F: Hansch preffered first plus duplicate value of highest occurance aka mode plus single values
    returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = fmcsoutputframe$fragments[1], value = fmcsoutputframe$F[1])
    return(returnlist)
  } else if (sigma.selection =="G")  {
    #G: hydrowin plus hansch preffered
    returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = fmcsoutputframe$fragments[1], value = fmcsoutputframe$G[1])
    return(returnlist)
  } else if (sigma.selection == "H") {
    #H:  median of distinct values ??
    returnlist <- list (tanimoto = fmcsoutputframe$Tanimoto_Coefficient[1], index = fmcsoutputframe$index[1], sub = fmcsoutputframe$fragments[1], value = fmcsoutputframe$H[1])
    return(returnlist)
  } else if (sigma.selection == "U") {
    # user created sigma selection
    stop ("user has not created a custom sigma values option")
    } else {
      stop ("Specify valid sigma.selection")
    }
}
