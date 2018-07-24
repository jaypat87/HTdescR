#Helper function for htdesc
#Grabs values for phenyl fragments without running htdesc

helper <- function(type, sigma.select){
  hfilter <- filter(htdescHelper, HT.type == type, sigma.selection == sigma.select)
  return(hfilter)
}

