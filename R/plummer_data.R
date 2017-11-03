#'@rdname plummer_data
#'@title Plummer Model - Data
#'@description dataset from Plummer's article
#'@export

plummer_data <- function(){
  nhpv <- c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4)
  Npart <- c(111, 71, 162, 188, 145, 215, 166, 37, 173,
             143, 229, 696, 93)
  J <- length(nhpv)
  # Data for module 2
  ncases <- c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194)
  Npop <- c(26983, 250930, 829348, 157775, 150467, 352445, 553066,
            26751, 75815, 150302, 354993, 3683043, 507218)
  D <- log(10**(-3) * Npop)
  return(list(nhpv = nhpv, Npart = Npart, J = J, ncases = ncases, Npop = Npop, D = D))
}
