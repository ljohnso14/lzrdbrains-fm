# 3 functions to import and process the brain, trait, and phylogeny data

#' Make tree data
#'
#' This function imports a phylogeny (.nex), converts species ID's to all lowercase, and saves the file as .Rda
#' @param data Raw data set of .nex filetype. Defaults to 'lizard_tree.nex' from the raw_data folder
#' @return returns .Rda file named lzrd_tree in the data folder
#' @export
#'

make_tree <- function(data = 'lizard_tree.nex'){

  # Read data file
  lzrd_tree <- ape::read.nexus(here::here(paste('data_raw', data, sep = '/')))

  # Convert spieces ID to all lowercase since impute function requires it
  lzrd_tree$tip.label <- tolower(lzrd_tree$tip.label)

  # Save data
  save(lzrd_tree, file = here::here('data/lizard_tree.Rda')) #hopefuly changing to Rda is okay!
}



#' Make trait data
#'
#' This function imports lizard trait data, converts species ID's to all lowercase, imputes missing continuous data for species, and saves the file as .Rda
#' @param data Raw data set of .csv filetype. Defaults to 'lizard_trait_data.csv' from the raw_data folder
#' @return returns .Rda file named lizard_trait_data in the data folder
#' @export

make_trait <- function(data = 'lizard_trait_data.csv'){

  # Read data file
  trait <- read.csv(here::here(paste('data_raw', data, sep = '/')))

  # Load lizard_tree for imputation
  load('data/lizard_tree.Rda')

  # Convert species ID to all lowercase since impute function requires it
  trait$Id <- tolower(trait$Id)
  trait$foraging_mode <- tolower(trait$foraging_mode)


  # Convert rownames from numbers to the names of the species
  rownames(trait) <- trait$Id

  ### MULTIPLE STEPS TO IMPUTE MISSING CONTINUOUS TRAIT DATA

  # Create new df with just the continous traits
  cont.trait <- trait[,c("Id", "maxSVL", "fSVL", "hatchlingSVL", "clutch_size",
                         "precip_seasonality", "temp_seasonality")]

  # The impute function (phylopars) requires the column with teh species names to be called 'species'
  names(cont.trait)[1] <- "species"

  # Imputed values for missing continuous trait data by using the relationships in tree
  impute.cont.trait <- Rphylopars::phylopars(trait_data = cont.trait, tree = lzrd_tree)

  # Create new df with the imputed results
  values.impute.cont.trait <- impute.cont.trait$anc_recon

  # The impute function creates empty rows at the end, so removing those
  values.impute.cont.trait <- values.impute.cont.trait[c(1:29),]

  # Manually add imputed values into the NA position in trait
  trait["amphisbaena_scutigerum","fSVL"] = 325.891
  trait["amphisbaena_scutigerum", "hatchlingSVL"] <- 98.704
  trait["amphisbaena_scutigerum","clutch_size"] <- 5.507
  trait["amphisbaena_scutigerum","precip_seasonality"] = 54.69
  trait["amphisbaena_scutigerum","temp_seasonality"] = 6485.59



  trait["plestiodon_marginatus","hatchlingSVL"] <- 33.62
  trait["plestiodon_marginatus","clutch_size"] <- 4.78
  trait["plestiodon_marginatus","precip_seasonality"] <- 70.22
  trait["plestiodon_marginatus","temp_seasonality"] <- 3369.46

  trait["rieppeleon_brevicaudatus","hatchlingSVL"] <- 25.52

  # For analyses we want to use absolute latitude
  # not latitude since we care about difference between species
  # near the equator and species near the poles
  # whether they be south or north of the equator
  trait$abslatitude <- abs(trait$latitude)

  # Save data
  save(trait, file = here::here('data/lizard_trait_data.Rda'))

}



#' Make brain data
#'
#' This function imports lizard and snake brain coordinate data (x,y,z), removes snakes, converts species ID's to all lowercase, converts 2D array into 3D array, creates separate data arrays for each brain region of interest, and saves the file as .Rda
#' @param data Raw data set of .csv filetype. Defaults to 'whole_brain_coords_S1.csv' from the raw_data folder
#' @return returns six .Rda files named whole_brain_coords_S1, tel_coords_S1, dien_coords_S1, mes_coords_S1, cere_coords_S1, medob_coords_S1 in the data folder
#' @export
#' @details the 3D array is [a, b, c]
#' @details where
#' @details a = landmark id (e.g., landmark 1, landmark 2, landmark 3, etc.) with a max of 61 values
#' @details b = the 3D coordinates (e.g., x coord, y coord, z coord) with a max of 3 values
#' @details c = species id coded by position in dataframe (e.g., Actonias_meleagris is 1)
#' @details species id       #    , , 1
#' @details
#' @details 3D coordinates   #    [,1]     [,2]     [,3]
#' @details landmark 1       #    [1,] 9104.361 14396.24 7085.975
#' @details landmark 2       #    [2,] 9740.752 13144.82 6818.917
#' @details
#' @details
#' @details Which landmarks ids (1 through 61) correspond to what major brain region
#' @details *3 and 45 are shared between Diencephalon, Mesencephalon, and Medulla Oblongata
#' @details (even though 3 isn't duplicated between them in their data - weird)
#' @details
#' @details Telencephalon
#' @details 1:5, 10:14, 19:26
#' @details Anterior-most extent of the olfactory bulb 1-2
#' @details Lateral-most extent of the olfactory bulb 3-4
#' @details
#' @details Diencephalon
#' @details 6, 15, 43:45, 54:55, 61
#' @details Optic chiasm - mid-sagittal plane 20
#' @details
#' @details Mesencephalon
#' @details 7:8, 16:17, 27:28, 45, 48:49, 52:53, 58:59
#' @details
#' @details Cerebellum
#' @details 30:42, 50:51
#' @details
#' @details Medulla oblongota
#' @details 9, 18, 29, 45:47, 56:57, 60


make_brain <- function(data = 'whole_brain_coords_S1.csv'){

  # Read data file
  brain <- read.csv(here::here(paste('data_raw', data, sep = '/')))


  # Remove snakes from the data since we are only interested in lizards
  snakes <- c("Xerotyphlops_vermicularis",
              "Python_regius",
              "Epicrates_cenchria",
              "Eryx_jaculus",
              "Cerastes_cerastes",
              "Hydrophis_platurus",
              "Boaedon_fuliginosus",
              "Chrysopelea_ornata",
              "Dendrelaphis_pictus",
              "Dasypeltis_gansi",
              "Pantherophis_guttatus")

  lzrd_brain <- brain[!((brain$Id) %in% snakes),]

  # Convert species ID to all lowercase
  lzrd_brain$Id <- tolower(lzrd_brain$Id)

  # Arrange species in alphabetical order, so they line up with tree and trait row
  lzrd_brain <- lzrd_brain[order(lzrd_brain$Id),]

  # Convert rownames from numbers to the names of the species
  rownames(lzrd_brain) <- lzrd_brain$Id

  # Remove column with species names, so that I can convert coordinates from 2D to 3D array
  lzrd_brain_3D <- as.matrix(lzrd_brain[,-(1)])

  # Converting coordinates from 2D to 3D array
  lzrd_brain_3D <- geomorph::arrayspecs(lzrd_brain_3D, 61, 3)

  # Creating individual 3D arrays for each brain region of interest

  wholebrain <- lzrd_brain_3D
  tel <- lzrd_brain_3D[c(1:5, 10:14, 19:26), 1:3, 1:29]
  dien <- lzrd_brain_3D[c(6, 15, 43:45, 54:55, 61), 1:3, 1:29]
  mes <- lzrd_brain_3D[c(7:8, 16:17, 27:28, 45, 48:49, 52:53, 58:59), 1:3, 1:29]
  cere <- lzrd_brain_3D[c(30:42, 50:51), 1:3, 1:29]
  medob <- lzrd_brain_3D[c(9, 18, 29, 45:47, 56:57, 60), 1:3, 1:29]

  # Save data
  save(wholebrain, file = here::here('data/whole_brain_coords_S1.Rda'))
  save(tel, file = here::here('data/tel_coords_S1.Rda'))
  save(dien, file = here::here('data/dien_coords_S1.Rda'))
  save(mes, file = here::here('data/mes_coords_S1.Rda'))
  save(cere, file = here::here('data/cere_coords_S1.Rda'))
  save(medob, file = here::here('data/medob_coords_S1.Rda'))
}



### the 3D array is [a, b, c]
### where
### a = landmark id (e.g., landmark 1, landmark 2, landmark 3, etc.) with a max of 61 values
### b = the 3D coordinates (e.g., x coord, y coord, z coord) with a max of 3 values
### c = species id coded by position in dataframe (e.g., Actonias_meleagris is 1)
###  species id       #    , , 1

###  3D coordinates   #    [,1]     [,2]     [,3]
###  landmark 1       #    [1,] 9104.361 14396.24 7085.975
###  landmark 2       #    [2,] 9740.752 13144.82 6818.917


### Create individual dataframes for each major brain region
### source for subsetting 3D arrays: http://adv-r.had.co.nz/Subsetting.html

### Which landmarks ids (1 through 61) correspond to what major brain region
### *3 and 45 are shared between Diencephalon, Mesencephalon, and Medulla Oblongata
### (even though 3 isn't duplicated between them in their data - weird)

# Telencephalon
# 1:5, 10:14, 19:26
# Anterior-most extent of the olfactory bulb 1-2
# Lateral-most extent of the olfactory bulb 3-4

# Diencephalon
# 6, 15, 43:45, 54:55, 61
# Optic chiasm - mid-sagittal plane 20

# Mesencephalon
# 7:8, 16:17, 27:28, 45, 48:49, 52:53, 58:59

# Cerebellum
# 30:42, 50:51

# Medulla oblongota
# 9, 18, 29, 45:47, 56:57, 60

