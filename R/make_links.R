#' Make links
#'
#' This isn't a function to run. But I'm creating the links between the different 3D brain coordinates so that we can use them to visualize the brain with wire frames. I originally had this in the R markdown, but that doesn't work because the links function is a active function where I'd have to recreate the links each time.
#' @param data whole_brain_coords_S1.Rda from the data folder
#' @return returns .Rda file named whole_brain_links in the data folder
#' @export
#'
#'
make_links <- function(){
load(here::here('data/whole_brain_coords_S1.Rda'))
load(here::here('data/cere_coords_S1.Rda'))
load(here::here('data/dien_coords_S1.Rda'))
load(here::here('data/medob_coords_S1.Rda'))
load(here::here('data/mes_coords_S1.Rda'))
load(here::here('data/tel_coords_S1.Rda'))

library(geomorph)

# GPA to get procrustes coordinates
wholebrain_gpa <- gpagen(wholebrain)
cere_gpa <- gpagen(cere)
dien_gpa <- gpagen(dien)
medob_gpa <- gpagen(medob)
mes_gpa <- gpagen(mes)
tel_gpa <- gpagen(tel)

# PCA on procrustes coordinates
wholebrain_pca <- gm.prcomp(wholebrain_gpa$coords)

# creating reference brain shape based on procrustes coordinates
wholebrain_ref <- mshape(wholebrain_gpa$coords)
cere_ref <- mshape(cere_gpa$coords)
dien_ref <- mshape(dien_gpa$coords)
medob_ref <- mshape(medob_gpa$coords)
mes_ref <- mshape(mes_gpa$coords)
tel_ref <- mshape(tel_gpa$coords)

# Note: I actually think it may be easier to do this by doing hte indivual subregions first...
# Then maybe I could either bind the different files or use the images to help when drawing whole brain links

cere_links <- define.links(spec = cere_ref, ptsize = 5)
dien_links <- define.links(spec = dien_ref, ptsize = 5)
medob_links <- define.links(spec = medob_ref, ptsize = 5)
mes_links <- define.links(spec = mes_ref, ptsize = 5)
tel_links <- define.links(spec = tel_ref, ptsize = 5)
wholebrain_links <- define.links(spec = wholebrain_ref, ptsize = 5)


# Save data
save(cere_links, file = here::here('data/cere_links.Rda'))
save(dien_links, file = here::here('data/dien_links.Rda'))
save(medob_links, file = here::here('data/medob_links.Rda'))
save(mes_links, file = here::here('data/mes_links.Rda'))
save(tel_links, file = here::here('data/tel_links.Rda'))
save(wholebrain_links, file = here::here('data/wholebrain_links.Rda'))
}





