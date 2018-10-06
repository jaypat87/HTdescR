#' Batch search and match of fragment sigma values
#'
#'
#' Get values of Hammett-Taft descriptors for a given chemical fragment in SMILES string format by iterating through a lookup table. In case an exact match isnt found, this function uses a mismatch tolerant
#' maximum common substructure (fMCS) based fragment substitution library to get the HT dexcriptors with highest tanimoto coefficient. This function iterates through a loop to complete a batch file of sigma
#' values.
#'
#' @usage htbatch_hash(file, sigma.selection = "A", ...)
#'
#'
#' @param file path to csv file
#' @param sigma.selection The type of sigma to be returned; valid inputs include "A", "B", "C", "D", "E", "F", "G", "H", and "U"
#' @param ... inherit arguments
#'
#' @return Filled dataframe columns resulting from similarity search and value extraction from
#' esSDF, indSDF, metaSDF, paraSDF, orthoSDF, taftSDF, userSDF
#'
#' @export
#'
#' @examples ## NOT RUN: htbatch("./data/dataacidester.csv", sigma.selection ="A")
htbatch_hash <- function (file, sigma.selection = "A", ...) {

  # create an empty hash table

  hash_taft <- hash::hash()
  hash_ind <- hash::hash()
  hash_es <- hash::hash()
  hash_meta <- hash::hash()
  hash_para <- hash::hash()
  hash_ortho <- hash::hash()

  #reading the csv file as a dataframe

  qsardataframe <- utils::read.csv(file, stringsAsFactors = FALSE,na.strings = "", encoding = "UTF-8")


  # colnames(qsardataframe)[colnames(qsardataframe)=="ï..no"] <- "no"

  # initializing the iterator

  i = 1
  n <- nrow (qsardataframe)

  for (i in 1:n) {
    # create value for hash table to save instead of "r1.taft.smiles[i]" character string
      taftsmiles1 <- qsardataframe$r1.taft.smiles[[i]]
      indsmiles1 <- qsardataframe$r1.ind.smiles[[i]]
      essmiles1 <- qsardataframe$r1.es.smiles[[i]]
      meta1smiles1 <- qsardataframe$r1.meta1.smiles[[i]]
      meta1smiles2 <- qsardataframe$r1.meta2.smiles[[i]]
      parasmiles1 <- qsardataframe$r1.para.smiles[[i]]
      ortho1smiles1 <- qsardataframe$r1.ortho1.smiles[[i]]
      ortho1smiles2 <- qsardataframe$r1.ortho2.smiles[[i]]
      taftsmiles2 <- qsardataframe$r2.taft.smiles[[i]]
      indsmiles2 <- qsardataframe$r2.ind.smiles[[i]]
      essmiles2 <- qsardataframe$r2.es.smiles[[i]]
      meta2smiles1 <- qsardataframe$r2.meta1.smiles[[i]]
      meta2smiles2 <- qsardataframe$r2.meta2.smiles[[i]]
      parasmiles2 <- qsardataframe$r2.para1.smiles[[i]]
      ortho2smiles1 <- qsardataframe$r2.ortho1.smiles[[i]]
      ortho2smiles2 <- qsardataframe$r2.ortho2.smiles[[i]]

    if (is.na(qsardataframe$r1.meta1.smiles[i]) & is.na(qsardataframe$r1.ortho1.smiles[i]) & is.na(qsardataframe$r1.para1.smiles[i]) == TRUE) {

      if (hash::has.key(taftsmiles1, hash_taft) == TRUE) {

        t <- hash_taft[[taftsmiles1]]
        qsardataframe$r1.taft.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r1.taft.mcs.index[i] <- t$index
        qsardataframe$r1.taft.value[i] <- t$value
        rm (t)

      } else {

          t <- htdesc (smile = qsardataframe$r1.taft.smiles[i], HT.type = "taft", sigma.selection)

          hash_taft[[taftsmiles1]] = t

          qsardataframe$r1.taft.sub.smiles[i] <- as.character (t$sub)
          qsardataframe$r1.taft.mcs.index[i] <- t$tanimoto
          qsardataframe$r1.taft.value[i] <- t$value
          rm (t)
      }

      if (hash::has.key(indsmiles1, hash_ind) == TRUE) {

        t <- hash_ind[[indsmiles1]]
        qsardataframe$r1.ind.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r1.ind.mcs.index[i] <- t$index
        qsardataframe$r1.ind.value[i] <- t$value
        rm (t)


      } else {

        t <- htdesc (smile = qsardataframe$r1.ind.smiles[i], HT.type = "inductive", sigma.selection)

        hash_ind[[indsmiles1]] = t

        t <- htdesc (smile = qsardataframe$r1.ind.smiles[i], HT.type = "inductive", sigma.selection)
        qsardataframe$r1.ind.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r1.ind.mcs.index[i] <- t$tanimoto
        qsardataframe$r1.ind.value[i] <- t$value
        rm (t)

      }

      if (hash::has.key(essmiles1, hash_es) == TRUE) {

        t <- hash_es[[essmiles1]]
        qsardataframe$r1.es.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r1.es.mcs.index[i] <- t$index
        qsardataframe$r1.es.value[i] <- t$value
        rm (t)

      } else {

        t <- htdesc (smile = qsardataframe$r1.es.smiles[i], HT.type = "es", sigma.selection)

        hash_es[[essmiles1]] = t

        qsardataframe$r1.es.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r1.es.mcs.index[i] <- t$tanimoto
        qsardataframe$r1.es.value[i] <- t$value
        rm (t)
      }


    } else {

      # if the structure is aromatic, than we will set these to default values

      qsardataframe$r1.taft.smiles[i] <- "*C1=CC=CC=C1"
      qsardataframe$r1.es.smiles[i] <- "*C1=CC=CC=C1"
      qsardataframe$r1.ind.smiles[i] <- "*C1=CC=CC=C1"

      # calling htdesc helper function to fill substitute mcs values
      t <- helper (type = "taft", sigma.select = sigma.selection)
      qsardataframe$r1.taft.sub.smiles[i] <- as.character (t$sub.smiles)
      qsardataframe$r1.taft.mcs.index[i] <- t$tanimoto
      qsardataframe$r1.taft.value[i] <- t$value
      rm (t)

      t <- helper (type = "inductive", sigma.select = sigma.selection)
      qsardataframe$r1.ind.sub.smiles[i] <- as.character (t$sub.smiles)
      qsardataframe$r1.ind.mcs.index[i] <- t$tanimoto
      qsardataframe$r1.ind.value[i] <- t$value
      rm (t)

      t <- helper (type = "es", sigma.select = sigma.selection)
      qsardataframe$r1.es.sub.smiles[i] <- as.character (t$sub.smiles)
      qsardataframe$r1.es.mcs.index[i] <- t$tanimoto
      qsardataframe$r1.es.value[i] <- t$value
      rm (t)

    }

    if (is.na(qsardataframe$r1.meta1.smiles[i]) == FALSE) {

      if (hash::has.key(meta1smiles1, hash_meta) == TRUE) {

        t <- hash_meta[[meta1smiles1]]
        qsardataframe$r1.meta1.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r1.meta1.mcs.index[i] <- t$tanimoto
        qsardataframe$r1.meta1.value[i] <- t$value
        rm (t)

      } else {

          t <- htdesc (smile = qsardataframe$r1.meta1.smiles[i], HT.type = "meta", sigma.selection)

          hash_meta[[meta1smiles1]] = t

          qsardataframe$r1.meta1.sub.smiles[i] <- as.character (t$sub)
          qsardataframe$r1.meta1.mcs.index[i] <- t$tanimoto
          qsardataframe$r1.meta1.value[i] <- t$value
          rm (t)
      }
    }

    if (is.na(qsardataframe$r1.meta2.smiles[i]) == FALSE) {

      if (hash::has.key(meta1smiles2, hash_meta) == TRUE) {

        t <- hash_meta[[meta1smiles2]]
        qsardataframe$r1.meta2.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r1.meta2.mcs.index[i] <- t$tanimoto
        qsardataframe$r1.meta2.value[i] <- t$value
        rm (t)

      } else {

          t <- htdesc (smile = qsardataframe$r1.meta2.smiles[i], HT.type = "meta", sigma.selection)

          hash_meta[[meta1smiles2]] = t

          qsardataframe$r1.meta2.sub.smiles[i] <- as.character (t$sub)
          qsardataframe$r1.meta2.mcs.index[i] <- t$tanimoto
          qsardataframe$r1.meta2.value[i] <- t$value
          rm (t)
      }
    }

    if (is.na(qsardataframe$r1.ortho1.smiles[i]) == FALSE) {

      if (hash::has.key(ortho1smiles1, hash_ortho) == FALSE) {

        t <- hash_ortho[[ortho1smiles1]]
        qsardataframe$r1.ortho1.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r1.ortho1.mcs.index[i] <- t$tanimoto
        qsardataframe$r1.ortho1.value[i] <- t$value
        rm(t)

      } else {

          t <- htdesc (smile = qsardataframe$r1.ortho1.smiles[i], HT.type = "ortho", sigma.selection)

          hash_ortho[[ortho1smiles1]] = t

          qsardataframe$r1.ortho1.sub.smiles[i] <- as.character (t$sub)
          qsardataframe$r1.ortho1.mcs.index[i] <- t$tanimoto
          qsardataframe$r1.ortho1.value[i] <- t$value
          rm (t)
      }
    }

    if (is.na(qsardataframe$r1.ortho2.smiles[i]) == FALSE) {

      if (hash::has.key(ortho1smiles2, hash_ortho) == FALSE) {

        t <- hash_ortho[[ortho1smiles2]]
        qsardataframe$r1.ortho2.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r1.ortho2.mcs.index[i] <- t$tanimoto
        qsardataframe$r1.ortho2.value[i] <- t$value
        rm(t)

      } else {

          t <- htdesc (smile = qsardataframe$r1.ortho2.smiles[i], HT.type = "ortho", sigma.selection)

          hash_ortho[[ortho1smiles2]] = t

          qsardataframe$r1.ortho2.sub.smiles[i] <- as.character (t$sub)
          qsardataframe$r1.ortho2.mcs.index[i] <- t$tanimoto
          qsardataframe$r1.ortho2.value[i] <- t$value
          rm (t)
      }
    }

    if (is.na(qsardataframe$r1.ortho2.smiles[i]) == FALSE) {
      if (hash::has.key(parasmiles1, hash_para) == FALSE) {

        t <- hash_para[[parasmiles1]] = t
        qsardataframe$r1.para1.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r1.para1.mcs.index[i] <- t$tanimoto
        qsardataframe$r1.para1.value[i] <- t$value
        rm (t)

      } else {

          t <- htdesc (smile = qsardataframe$r1.para1.smiles[i], HT.type = "para", sigma.selection)

          hash_para[[parasmiles1]] = t

          qsardataframe$r1.para1.sub.smiles[i] <- as.character (t$sub)
          qsardataframe$r1.para1.mcs.index[i] <- t$tanimoto
          qsardataframe$r1.para1.value[i] <- t$value
          rm (t)
      }
    }

    #For R2

    if (is.na(qsardataframe$r2.meta1.smiles[i]) & is.na(qsardataframe$r2.ortho1.smiles[i]) & is.na(qsardataframe$r2.para1.smiles[i]) == TRUE) {

      if (hash::has.key(taftsmiles2, hash_taft) == TRUE) {

        t <- hash_taft[[taftsmiles2]]
        qsardataframe$r2.taft.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r2.taft.mcs.index[i] <- t$tanimoto
        qsardataframe$r2.taft.value[i] <- t$value
        rm (t)

      } else {

        # calling htdesc to fill substitute mcs values

        t <- htdesc (smile = qsardataframe$r2.taft.smiles[i], HT.type = "taft", sigma.selection)

        hash_taft[[taftsmiles2]] = t

        qsardataframe$r2.taft.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r2.taft.mcs.index[i] <- t$tanimoto
        qsardataframe$r2.taft.value[i] <- t$value
        rm(t)
      }

      if (hash::has.key(indsmiles2, hash_ind) == TRUE) {

        t <- hash_ind[[indsmiles2]]
        qsardataframe$r2.ind.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r2.ind.mcs.index[i] <- t$tanimoto
        qsardataframe$r2.ind.value[i] <- t$value
        rm (t)

      } else {

          t <- htdesc (smile = qsardataframe$r2.ind.smiles[i], HT.type = "inductive", sigma.selection)

          hash_ind[[indsmiles2]] = t

          qsardataframe$r2.ind.sub.smiles[i] <- as.character (t$sub)
          qsardataframe$r2.ind.mcs.index[i] <- t$tanimoto
          qsardataframe$r2.ind.value[i] <- t$value
          rm (t)

      }

      if (hash::has.key(essmiles2, hash_es) == TRUE) {

        t <- hash_es[[essmiles2]]
        qsardataframe$r2.es.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r2.es.mcs.index[i] <- t$tanimoto
        qsardataframe$r2.es.value[i] <- t$value
        rm (t)

      } else {

          t <- htdesc (smile = qsardataframe$r2.es.smiles[i], HT.type = "es", sigma.selection)

          hash_es[[essmiles2]] = t

          qsardataframe$r2.es.sub.smiles[i] <- as.character (t$sub)
          qsardataframe$r2.es.mcs.index[i] <- t$tanimoto
          qsardataframe$r2.es.value[i] <- t$value
          rm (t)

      }


    } else {

      # if the structure is aromatic, than we will set these to default values

      qsardataframe$r2.taft.smiles[i] <- "*C1=CC=CC=C1"
      qsardataframe$r2.es.smiles[i] <- "*C1=CC=CC=C1"
      qsardataframe$r2.ind.smiles[i] <- "*C1=CC=CC=C1"

      t <- helper (type = "taft", sigma.select = sigma.selection)
      qsardataframe$r2.taft.sub.smiles[i] <- as.character (t$sub.smiles)
      qsardataframe$r2.taft.mcs.index[i] <- t$tanimoto
      qsardataframe$r2.taft.value[i] <- t$value
      rm (t)

      t <- helper (type = "inductive", sigma.select = sigma.selection)
      qsardataframe$r2.ind.sub.smiles[i] <- as.character (t$sub.smiles)
      qsardataframe$r2.ind.mcs.index[i] <- t$tanimoto
      qsardataframe$r2.ind.value[i] <- t$value
      rm (t)

      t <- helper (type = "es", sigma.select = sigma.selection)
      qsardataframe$r2.es.sub.smiles[i] <- as.character (t$sub.smiles)
      qsardataframe$r2.es.mcs.index[i] <- t$tanimoto
      qsardataframe$r2.es.value[i] <- t$value
      rm (t)

    }

    if (is.na(qsardataframe$r2.meta1.smiles[i]) == FALSE) {

      if (hash::has.key(meta2smiles1, hash_meta) == TRUE) {

        t <- hash_meta[[meta2smiles1]]
        qsardataframe$r2.meta1.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r2.meta1.mcs.index[i] <- t$tanimoto
        qsardataframe$r2.meta1.value[i] <- t$value
        rm (t)

      } else {

        t <- htdesc (smile = qsardataframe$r2.meta1.smiles[i], HT.type = "meta", sigma.selection)

        hash_meta[[meta2smiles1]] = t

        qsardataframe$r2.meta1.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r2.meta1.mcs.index[i] <- t$tanimoto
        qsardataframe$r2.meta1.value[i] <- t$value
        rm (t)

      }
    }

    if (is.na(qsardataframe$r2.meta2.smiles[i]) == FALSE) {

      if (hash::has.key(meta2smiles2, hash_meta) == TRUE) {

        t <- hash_meta[[meta2smiles2]]
        qsardataframe$r2.meta2.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r2.meta2.mcs.index[i] <- t$tanimoto
        qsardataframe$r2.meta2.value[i] <- t$value
        rm (t)

      } else {

        t <- htdesc (smile = qsardataframe$r2.meta2.smiles[i], HT.type = "meta", sigma.selection)

        hash_meta[[meta2smiles2]] = t

        qsardataframe$r2.meta2.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r2.meta2.mcs.index[i] <- t$tanimoto
        qsardataframe$r2.meta2.value[i] <- t$value
        rm (t)

      }

    }

    if (is.na(qsardataframe$r2.ortho1.smiles[i]) == FALSE) {

      if (hash::has.key(ortho2smiles1, hash_ortho) == TRUE) {

        t <- hash_ortho[[ortho2smiles1]]
        qsardataframe$r2.ortho1.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r2.ortho1.mcs.index[i] <- t$tanimoto
        qsardataframe$r2.ortho1.value[i] <- t$value
        rm (t)

      } else {

        t <- htdesc (smile = qsardataframe$r2.ortho1.smiles[i], HT.type = "ortho", sigma.selection)

        hash_ortho[[ortho2smiles1]] = t

        qsardataframe$r2.ortho1.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r2.ortho1.mcs.index[i] <- t$tanimoto
        qsardataframe$r2.ortho1.value[i] <- t$value
        rm (t)

      }

    }

    if (is.na(qsardataframe$r2.ortho2.smiles[i]) == FALSE) {

      if (hash::has.key(ortho2smiles2, hash_ortho) == TRUE) {

        t <- hash_ortho[[ortho2smiles2]]
        qsardataframe$r2.ortho2.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r2.ortho2.mcs.index[i] <- t$tanimoto
        qsardataframe$r2.ortho2.value[i] <- t$value
        rm (t)

      } else {

        t <- htdesc (smile = qsardataframe$r2.ortho2.smiles[i], HT.type = "ortho", sigma.selection)

        hash_ortho[[ortho2smiles2]] = t

        qsardataframe$r2.ortho2.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r2.ortho2.mcs.index[i] <- t$tanimoto
        qsardataframe$r2.ortho2.value[i] <- t$value
        rm (t)

      }

    }

    if (is.na(qsardataframe$r2.para1.smiles[i]) == FALSE) {

      if (hash::has.key(parasmiles2, hash_para) == TRUE) {

        t <- hash_para[[parasmiles2]]
        qsardataframe$r2.para1.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r2.para1.mcs.index[i] <- t$tanimoto
        qsardataframe$r2.para1.value[i] <- t$value
        rm (t)

      } else {

        t <- htdesc (smile = qsardataframe$r2.para1.smiles[i], HT.type = "para", sigma.selection)

        hash_para[[parasmiles2]] = t

        qsardataframe$r2.para1.sub.smiles[i] <- as.character (t$sub)
        qsardataframe$r2.para1.mcs.index[i] <- t$tanimoto
        qsardataframe$r2.para1.value[i] <- t$value
        rm (t)

      }

    }
  }

  closeAllConnections()
  # I am not sure why this is here, but dont remove it!!
  return (qsardataframe)

  #work still left
  # Low priority
  # insert output file format as a function attribute

  # medium priority

  # use appropriate HT.type after creating the data files
  # use appropriate sigma.selection method

  # high priority

  # insert if statement for es and ind
  # replace NA with 0 for .values cells
  # you can leave NA for .mcs.index and .sub.smiles for ones where we didnt trigger getsmiles function
  # insert a summation method which adds up values for hammetts

}
