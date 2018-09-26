# simple illustration of using hash function from package hash to use create key-value pairs with user imputted SMILES as key and a list as value


   initiallist <- list (tanimoto = 1, index = 10, sub = "C", value = 100)
   #x <- hash("*C", initiallist)
   #creating an empty hash

   x <- hash::hash()

# setting a new key

  x[["*CC"]] <- initiallist

  #looking up its value from the hash

  y <- x [["smile"]]
   y$tanimoto

   #checking if a key is present in the hashtable
   hash::has.key("*CC", x )

   # getting all the keys

   hash::keys(x)


