#initialization scripts

# inputting training and test set for carboxylic acid esters

acidester <- read.csv(file = "./inst/extdata/carboxylicacidesterinputfile-fixed.csv", stringsAsFactors = FALSE,na.strings = "", encoding = "UTF-8")
colnames(acidester)[colnames(acidester)=="X.U.FEFF.no"] <- "index"
save (acidester, file = "./data/acidester.rda")



#####

# script to create metasigma

metasigma <- read.csv(file = "./inst/extdata/sigmameta.csv", encoding = "UTF-8")

colnames(metasigma)[colnames(metasigma)=="X.U.FEFF.no"] <- "index"

colnames(metasigma)[colnames(metasigma)=="fragment.daylight"] <- "fragments"

metasigma$std.dev <- as.numeric(format(round(as.numeric(metasigma$std.dev), 2), nsmall = 2))
metasigma$reg.avg <- as.numeric(format(round(metasigma$reg.avg, 2), nsmall = 2))
metasigma$reg.median <- as.numeric(format(round(metasigma$reg.median, 2), nsmall = 2))
metasigma$reg.mode <- as.numeric(format(round(metasigma$reg.mode, 2), nsmall = 2))
metasigma$avg.dist <- as.numeric(format(round(metasigma$avg.dist, 2), nsmall = 2))
metasigma$hansch.pref <- as.numeric(format(round(metasigma$hansch.pref, 2), nsmall = 2))

save (metasigma, file = "./data/metasigma.rda")


# create metaSDF from SDF file

metaSDF <- ChemmineR::read.SDFset("./inst/extdata/sigmameta.SDF")
save (metaSDF, file = "./data/metaSDF.rda")



#####

# type out similar scripts to create taftsigma, parasigma, orthosigma, essigma, inductionsigma

# type out similar scripts to create orthoSDF.....



# script to create indsigma

indsigma <- read.csv(file = "./inst/extdata/sigmaind.csv", encoding = "UTF-8")

colnames(indsigma)[colnames(indsigma)=="X.U.FEFF.no"] <- "index"

colnames(indsigma)[colnames(indsigma)=="fragment.daylight"] <- "fragments"

indsigma$std.dev <- as.numeric(format(round(indsigma$std.dev, 2), nsmall = 2))
indsigma$reg.avg <- as.numeric(format(round(indsigma$reg.avg, 2), nsmall = 2))
indsigma$reg.median <- as.numeric(format(round(indsigma$reg.median, 2), nsmall = 2))
indsigma$reg.mode <- as.numeric(format(round(indsigma$reg.mode, 2), nsmall = 2))
indsigma$avg.dist <- as.numeric(format(round(indsigma$avg.dist, 2), nsmall = 2))
indsigma$hansch.pref <- as.numeric(format(round(indsigma$hansch.pref, 2), nsmall = 2))

save (indsigma, file = "./data/indsigma.rda")


# create indSDF from SDF file

indSDF <- ChemmineR::read.SDFset("./inst/extdata/sigmaind.SDF")
save (indSDF, file = "./data/indSDF.rda")



#####

# script to create parasigma

parasigma <- read.csv(file = "./inst/extdata/sigmapara.csv", encoding = "UTF-8")

colnames(parasigma)[colnames(parasigma)=="X.U.FEFF.no"] <- "index"

colnames(parasigma)[colnames(parasigma)=="fragment.daylight"] <- "fragments"

parasigma$std.dev <- as.numeric(format(round(parasigma$std.dev, 2), nsmall = 2))
parasigma$reg.avg <- as.numeric(format(round(parasigma$reg.avg, 2), nsmall = 2))
parasigma$reg.median <- as.numeric(format(round(parasigma$reg.median, 2), nsmall = 2))
parasigma$reg.mode <- as.numeric(format(round(parasigma$reg.mode, 2), nsmall = 2))
parasigma$avg.dist <- as.numeric(format(round(parasigma$avg.dist, 2), nsmall = 2))
parasigma$hansch.pref <- as.numeric(format(round(parasigma$hansch.pref, 2), nsmall = 2))

save (parasigma, file = "./data/parasigma.rda")


# create paraSDF from SDF file

paraSDF <- ChemmineR::read.SDFset("./inst/extdata/sigmapara.SDF")
save (paraSDF, file = "./data/paraSDF.rda")



#####

# script to create orthosigma

orthosigma <- read.csv(file = "./inst/extdata/sigmaortho.csv", encoding = "UTF-8")

colnames(orthosigma)[colnames(orthosigma)=="X.U.FEFF.no"] <- "index"

colnames(orthosigma)[colnames(orthosigma)=="fragment.daylight"] <- "fragments"

orthosigma$std.dev <- as.numeric(format(round(orthosigma$std.dev, 2), nsmall = 2))
orthosigma$reg.avg <- as.numeric(format(round(orthosigma$reg.avg, 2), nsmall = 2))
orthosigma$reg.median <- as.numeric(format(round(orthosigma$reg.median, 2), nsmall = 2))
orthosigma$reg.mode <- as.numeric(format(round(orthosigma$reg.mode, 2), nsmall = 2))
orthosigma$avg.dist <- as.numeric(format(round(orthosigma$avg.dist, 2), nsmall = 2))
orthosigma$hansch.pref <- as.numeric(format(round(orthosigma$hansch.pref, 2), nsmall = 2))

save (orthosigma, file = "./data/orthosigma.rda")


# create orthoSDF from SDF file

orthoSDF <- ChemmineR::read.SDFset("./inst/extdata/sigmaortho.SDF")
save (orthoSDF, file = "./data/orthoSDF.rda")



#####

# script to create taftsigma

taftsigma <- read.csv(file = "./inst/extdata/sigmataft.csv", encoding = "UTF-8")

colnames(taftsigma)[colnames(taftsigma)=="X.U.FEFF.no"] <- "index"

colnames(taftsigma)[colnames(taftsigma)=="fragment.daylight"] <- "fragments"

taftsigma$std.dev <- as.numeric(format(round(taftsigma$std.dev, 2), nsmall = 2))
taftsigma$reg.avg <- as.numeric(format(round(taftsigma$reg.avg, 2), nsmall = 2))
taftsigma$reg.median <- as.numeric(format(round(taftsigma$reg.median, 2), nsmall = 2))
taftsigma$reg.mode <- as.numeric(format(round(taftsigma$reg.mode, 2), nsmall = 2))
taftsigma$avg.dist <- as.numeric(format(round(taftsigma$avg.dist, 2), nsmall = 2))
taftsigma$hansch.pref <- as.numeric(format(round(taftsigma$hansch.pref, 2), nsmall = 2))

save (taftsigma, file = "./data/taftsigma.rda")


# create taftSDF from SDF file

taftSDF <- ChemmineR::read.SDFset("./inst/extdata/sigmataft.sdf")
save (taftSDF, file = "./data/taftSDF.rda")



#####

# script to create es

es <- read.csv(file = "./inst/extdata/es.csv", encoding = "UTF-8")

colnames(es)[colnames(es)=="X.U.FEFF.no"] <- "index"

colnames(es)[colnames(es)=="fragment.daylight"] <- "fragments"

es$std.dev <- as.numeric(format(round(es$std.dev, 2), nsmall = 2))
es$reg.avg <- as.numeric(format(round(es$reg.avg, 2), nsmall = 2))
es$reg.median <- as.numeric(format(round(es$reg.median, 2), nsmall = 2))
es$reg.mode <- as.numeric(format(round(es$reg.mode, 2), nsmall = 2))
es$avg.dist <- as.numeric(format(round(es$avg.dist, 2), nsmall = 2))
es$hansch.pref <- as.numeric(format(round(es$hansch.pref, 2), nsmall = 2))

save (es, file = "./data/es.rda")


# create esSDF from SDF file

esSDF <- ChemmineR::read.SDFset("./inst/extdata/es.SDF")
save (esSDF, file = "./data/esSDF.rda")
