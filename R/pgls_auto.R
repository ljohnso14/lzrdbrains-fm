# functions to automatically run gpa and pgls and produce a pretty table output

#' Phylogenetic ANOVA/Regression for Procrustes Shape Variables - SHAPE
#'
#'
#' @param brain 3D array of brain coordinates
#' @param tree tree in nexus format
#' @param traitdata dataframe with traits of interest
#' @param region string identifying the region of the brain coordinates (e.g., 'Whole Brain')
#' @param outputname string indicating the name of the file in a .rtf format (e.g., 'wholebrainshape.rft')
#' @return returns .rtf file containing a gt table into the vignettes folder
#' @export
#'
#' @details brain, tree, and traitdata need to be organized with the species in the same order and have the same exact name format (for example the tree can't have upper case while the trait has lower case)


pgls_shape <- function(brain, tree, traitdata, region, outputname){
  library(gt)
  library(geomorph)
  library(knitr)
  gpa <- gpagen(brain, ProcD = T, verbose = T)
  gdf <- geomorph.data.frame(gpa, traitdata)

  write(region, file = outputname, append = T)

  pgls1 <- procD.pgls(coords ~ Csize + abslatitude, phy = tree, SS.type = "II", data = gdf)

  pgls2 <- procD.pgls(coords ~ Csize + temp_seasonality, phy = tree, SS.type = "II", data = gdf)

  pgls3 <- procD.pgls(coords ~ Csize + precip_seasonality, phy = tree, SS.type = "II", data = gdf)

  pgls4 <- procD.pgls(coords ~ Csize + substrate, phy = tree, SS.type = "II", data = gdf)

  pgls5 <- procD.pgls(coords ~ Csize + microhabitat, phy = tree, SS.type = "II", data = gdf)

  pgls6 <- procD.pgls(coords ~ Csize + activity_time, phy = tree, SS.type = "II", data = gdf)

  pgls7 <- procD.pgls(coords ~ Csize + foraging_mode, phy = tree, SS.type = "II", data = gdf)


  pgls8 <- procD.pgls(coords ~ Csize + fSVL, phy = tree, SS.type = "II", data = gdf)

  pgls9 <- procD.pgls(coords ~ Csize + hatchlingSVL, phy = tree, SS.type = "II", data = gdf)

  pgls10 <- procD.pgls(coords ~ Csize + reproductive_mode, phy = tree, SS.type = "II", data = gdf)

  pgls11 <- procD.pgls(coords ~ Csize + clutch_size, phy = tree, SS.type = "II", data = gdf)

  pgls12 <- procD.pgls(coords ~ Csize + clutch_size + reproductive_mode, phy = tree, SS.type = "II", data = gdf)



  t1 <- broom::tidy(pgls1$aov.table)
  t2 <- broom::tidy(pgls2$aov.table)
  t3 <- broom::tidy(pgls3$aov.table)
  t4 <- broom::tidy(pgls4$aov.table)
  t5 <- broom::tidy(pgls5$aov.table)
  t6 <- broom::tidy(pgls6$aov.table)
  t7 <- broom::tidy(pgls7$aov.table)
  t8 <- broom::tidy(pgls8$aov.table)
  t9 <- broom::tidy(pgls9$aov.table)
  t10 <- broom::tidy(pgls10$aov.table)
  t11 <- broom::tidy(pgls11$aov.table)
  t12 <- broom::tidy(pgls12$aov.table)



  mergedtable <- rbind(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12)


  mergedtablefancy <- gt(mergedtable) %>%
    cols_label(term = 'Source', statistic = 'F', p.value = 'p-value') %>%
    fmt_number(columns = c('SS', 'MS', 'Rsq', 'statistic', 'Z', 'p.value'), decimals = 4) %>%
    fmt_missing(columns = 2:8, missing_text = " ") %>% tab_header(title = region) %>%
    tab_footnote(footnote = "df, degree of freedom.", locations = cells_column_labels(columns = df)) %>%
    tab_footnote(footnote = "SS, sum of squares.", locations = cells_column_labels(columns = SS)) %>%
    tab_footnote(footnote = "MS, mean squares.", locations = cells_column_labels(columns = MS)) %>%
    tab_footnote(footnote = "Significant values (p-value < 0.05) from permutation tests (1,000 permutation rounds) are bolded.",
                 locations = cells_column_labels(columns = p.value)) %>%
    tab_row_group(label = 'Latitude', rows = 1:4) %>%
    tab_row_group(label = 'Temperature Seasonality', rows = 5:8) %>%
    tab_row_group(label = 'Precipitation Seasonality', rows = 9:12) %>%
    tab_row_group(label = 'Substrate', rows = 13:16) %>%
    tab_row_group(label = 'Microhabitat', rows = 17:20) %>%
    tab_row_group(label = 'Activity Time', rows = 21:24) %>%
    tab_row_group(label = 'Foraging Mode', rows = 25:28) %>%
    tab_row_group(label = 'Female SVL', rows = 29:32) %>%
    tab_row_group(label = 'Hatchling SVL', rows = 33:36) %>%
    tab_row_group(label = 'Reproductive Mode', rows = 37:40) %>%
    tab_row_group(label = 'Clutch Size', rows = 41:44) %>%
    tab_row_group(label = 'Clutch Size and Reproductive Mode', rows = 45:49)


  gtsave(mergedtablefancy, filename = outputname, exapnd = 10)

}



#' Phylogenetic ANOVA/Regression for Procrustes Shape Variables - SIZE
#'
#'
#' @param brain 3D array of brain coordinates
#' @param tree tree in nexus format
#' @param traitdata dataframe with traits of interest
#' @param region string identifying the region of the brain coordinates (e.g., 'Whole Brain')
#' @param outputname string indicating the name of the file in a .rtf format (e.g., 'wholebrainshape.rft')
#' @return returns .rtf file containing a gt table into the vignettes folder
#' @export
#'
#' @details brain, tree, and traitdata need to be organized with the species in the same order and have the same exact name format (for example the tree can't have upper case while the trait has lower case)


pgls_size <- function(brain, tree, traitdata, region, outputname){
  library(gt)
  library(geomorph)
  library(knitr)
  gpa <- gpagen(brain, ProcD = T, verbose = T)
  gdf <- geomorph.data.frame(gpa, traitdata)

  write(region, file = outputname, append = T)

  pgls1 <- procD.pgls(Csize ~ abslatitude, phy = tree, SS.type = "II", data = gdf)

  pgls2 <- procD.pgls(Csize ~ temp_seasonality, phy = tree, SS.type = "II", data = gdf)

  pgls3 <- procD.pgls(Csize ~ precip_seasonality, phy = tree, SS.type = "II", data = gdf)

  pgls4 <- procD.pgls(Csize ~ substrate, phy = tree, SS.type = "II", data = gdf)

  pgls5 <- procD.pgls(Csize ~ microhabitat, phy = tree, SS.type = "II", data = gdf)

  pgls6 <- procD.pgls(Csize ~ activity_time, phy = tree, SS.type = "II", data = gdf)

  pgls7 <- procD.pgls(Csize ~ foraging_mode, phy = tree, SS.type = "II", data = gdf)


  pgls8 <- procD.pgls(Csize ~ fSVL, phy = tree, SS.type = "II", data = gdf)

  pgls9 <- procD.pgls(Csize ~ hatchlingSVL, phy = tree, SS.type = "II", data = gdf)

  pgls10 <- procD.pgls(Csize ~ reproductive_mode, phy = tree, SS.type = "II", data = gdf)

  pgls11 <- procD.pgls(Csize ~ clutch_size, phy = tree, SS.type = "II", data = gdf)

  pgls12 <- procD.pgls(Csize ~ clutch_size + reproductive_mode, phy = tree, SS.type = "II", data = gdf)




  t1 <- broom::tidy(pgls1$aov.table)
  t2 <- broom::tidy(pgls2$aov.table)
  t3 <- broom::tidy(pgls3$aov.table)
  t4 <- broom::tidy(pgls4$aov.table)
  t5 <- broom::tidy(pgls5$aov.table)
  t6 <- broom::tidy(pgls6$aov.table)
  t7 <- broom::tidy(pgls7$aov.table)
  t8 <- broom::tidy(pgls8$aov.table)
  t9 <- broom::tidy(pgls9$aov.table)
  t10 <- broom::tidy(pgls10$aov.table)
  t11 <- broom::tidy(pgls11$aov.table)
  t12 <- broom::tidy(pgls12$aov.table)



  mergedtable <- rbind(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12)

  library(gt)
  mergedtablefancy <- gt(mergedtable) %>%
    cols_label(term = 'Source', statistic = 'F', p.value = 'p-value') %>%
    fmt_number(columns = c('SS', 'MS'), decimals = 0) %>%
    fmt_number(columns = c('Rsq', 'statistic', 'Z', 'p.value'), decimals = 4) %>%
    fmt_missing(columns = 2:8, missing_text = " ") %>% tab_header(title = region) %>%
    tab_footnote(footnote = "df, degree of freedom.", locations = cells_column_labels(columns = df)) %>%
    tab_footnote(footnote = "SS, sum of squares.", locations = cells_column_labels(columns = SS)) %>%
    tab_footnote(footnote = "MS, mean squares.", locations = cells_column_labels(columns = MS)) %>%
    tab_footnote(footnote = "Significant values (p-value < 0.05) from permutation tests (1,000 permutation rounds) are bolded.",
                 locations = cells_column_labels(columns = p.value)) %>%
    tab_row_group(label = 'Latitude', rows = 1:3) %>%
    tab_row_group(label = 'Temperature Seasonality', rows = 4:6) %>%
    tab_row_group(label = 'Precipitation Seasonality', rows = 7:9) %>%
    tab_row_group(label = 'Substrate', rows = 10:12) %>%
    tab_row_group(label = 'Microhabitat', rows = 13:15) %>%
    tab_row_group(label = 'Activity Time', rows = 16:18) %>%
    tab_row_group(label = 'Foraging Mode', rows = 19:21) %>%
    tab_row_group(label = 'Female SVL', rows = 22:24) %>%
    tab_row_group(label = 'Hatchling SVL', rows = 25:27) %>%
    tab_row_group(label = 'Reproductive Mode', rows = 28:30) %>%
    tab_row_group(label = 'Clutch Size', rows = 31:33) %>%
    tab_row_group(label = 'Clutch Size and Reproductive Mode', rows = 34:37)




  gtsave(mergedtablefancy, filename = outputname, exapnd = 10)

}















pgls_shape_logsize <- function(brain, tree, traitdata, region, outputname){
  library(gt)
  library(geomorph)
  library(knitr)
  gpa <- gpagen(brain, ProcD = T, verbose = T)
  gdf <- geomorph.data.frame(gpa, traitdata)

  write(region, file = outputname, append = T)

  pgls1 <- procD.pgls(coords ~ log10(Csize) + abslatitude, phy = tree, SS.type = "II", data = gdf)

  pgls2 <- procD.pgls(coords ~ log10(Csize) + temp_seasonality, phy = tree, SS.type = "II", data = gdf)

  pgls3 <- procD.pgls(coords ~ log10(Csize) + precip_seasonality, phy = tree, SS.type = "II", data = gdf)

  pgls4 <- procD.pgls(coords ~ log10(Csize) + substrate, phy = tree, SS.type = "II", data = gdf)

  pgls5 <- procD.pgls(coords ~ log10(Csize) + microhabitat, phy = tree, SS.type = "II", data = gdf)

  pgls6 <- procD.pgls(coords ~ log10(Csize) + activity_time, phy = tree, SS.type = "II", data = gdf)

  pgls7 <- procD.pgls(coords ~ log10(Csize) + foraging_mode, phy = tree, SS.type = "II", data = gdf)


  pgls8 <- procD.pgls(coords ~ log10(Csize) + fSVL, phy = tree, SS.type = "II", data = gdf)

  pgls9 <- procD.pgls(coords ~ log10(Csize) + hatchlingSVL, phy = tree, SS.type = "II", data = gdf)

  pgls10 <- procD.pgls(coords ~ log10(Csize) + reproductive_mode, phy = tree, SS.type = "II", data = gdf)

  pgls11 <- procD.pgls(coords ~ log10(Csize) + clutch_size, phy = tree, SS.type = "II", data = gdf)

  pgls12 <- procD.pgls(coords ~ log10(Csize) + clutch_size + reproductive_mode, phy = tree, SS.type = "II", data = gdf)



  t1 <- broom::tidy(pgls1$aov.table)
  t2 <- broom::tidy(pgls2$aov.table)
  t3 <- broom::tidy(pgls3$aov.table)
  t4 <- broom::tidy(pgls4$aov.table)
  t5 <- broom::tidy(pgls5$aov.table)
  t6 <- broom::tidy(pgls6$aov.table)
  t7 <- broom::tidy(pgls7$aov.table)
  t8 <- broom::tidy(pgls8$aov.table)
  t9 <- broom::tidy(pgls9$aov.table)
  t10 <- broom::tidy(pgls10$aov.table)
  t11 <- broom::tidy(pgls11$aov.table)
  t12 <- broom::tidy(pgls12$aov.table)



  mergedtable <- rbind(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12)


  mergedtablefancy <- gt(mergedtable) %>%
    cols_label(term = 'Source', statistic = 'F', p.value = 'p-value') %>%
    fmt_number(columns = c('SS', 'MS', 'Rsq', 'statistic', 'Z', 'p.value'), decimals = 4) %>%
    fmt_missing(columns = 2:8, missing_text = " ") %>% tab_header(title = region) %>%
    tab_footnote(footnote = "df, degree of freedom.", locations = cells_column_labels(columns = df)) %>%
    tab_footnote(footnote = "SS, sum of squares.", locations = cells_column_labels(columns = SS)) %>%
    tab_footnote(footnote = "MS, mean squares.", locations = cells_column_labels(columns = MS)) %>%
    tab_footnote(footnote = "Significant values (p-value < 0.05) from permutation tests (1,000 permutation rounds) are bolded.",
                 locations = cells_column_labels(columns = p.value)) %>%
    tab_row_group(label = 'Latitude', rows = 1:4) %>%
    tab_row_group(label = 'Temperature Seasonality', rows = 5:8) %>%
    tab_row_group(label = 'Precipitation Seasonality', rows = 9:12) %>%
    tab_row_group(label = 'Substrate', rows = 13:16) %>%
    tab_row_group(label = 'Microhabitat', rows = 17:20) %>%
    tab_row_group(label = 'Activity Time', rows = 21:24) %>%
    tab_row_group(label = 'Foraging Mode', rows = 25:28) %>%
    tab_row_group(label = 'Female SVL', rows = 29:32) %>%
    tab_row_group(label = 'Hatchling SVL', rows = 33:36) %>%
    tab_row_group(label = 'Reproductive Mode', rows = 37:40) %>%
    tab_row_group(label = 'Clutch Size', rows = 41:44) %>%
    tab_row_group(label = 'Clutch Size and Reproductive Mode', rows = 45:49)


  gtsave(mergedtablefancy, filename = outputname, exapnd = 10)

}

pgls_size_logsize <- function(brain, tree, traitdata, region, outputname){
  library(gt)
  library(geomorph)
  library(knitr)
  gpa <- gpagen(brain, ProcD = T, verbose = T)
  gdf <- geomorph.data.frame(gpa, traitdata)

  write(region, file = outputname, append = T)

  pgls1 <- procD.pgls(log10(Csize) ~ abslatitude, phy = tree, SS.type = "II", data = gdf)

  pgls2 <- procD.pgls(log10(Csize) ~ temp_seasonality, phy = tree, SS.type = "II", data = gdf)

  pgls3 <- procD.pgls(log10(Csize) ~ precip_seasonality, phy = tree, SS.type = "II", data = gdf)

  pgls4 <- procD.pgls(log10(Csize) ~ substrate, phy = tree, SS.type = "II", data = gdf)

  pgls5 <- procD.pgls(log10(Csize) ~ microhabitat, phy = tree, SS.type = "II", data = gdf)

  pgls6 <- procD.pgls(log10(Csize) ~ activity_time, phy = tree, SS.type = "II", data = gdf)

  pgls7 <- procD.pgls(log10(Csize) ~ foraging_mode, phy = tree, SS.type = "II", data = gdf)


  pgls8 <- procD.pgls(log10(Csize) ~ fSVL, phy = tree, SS.type = "II", data = gdf)

  pgls9 <- procD.pgls(log10(Csize) ~ hatchlingSVL, phy = tree, SS.type = "II", data = gdf)

  pgls10 <- procD.pgls(log10(Csize) ~ reproductive_mode, phy = tree, SS.type = "II", data = gdf)

  pgls11 <- procD.pgls(log10(Csize) ~ clutch_size, phy = tree, SS.type = "II", data = gdf)

  pgls12 <- procD.pgls(log10(Csize) ~ clutch_size + reproductive_mode, phy = tree, SS.type = "II", data = gdf)




  t1 <- broom::tidy(pgls1$aov.table)
  t2 <- broom::tidy(pgls2$aov.table)
  t3 <- broom::tidy(pgls3$aov.table)
  t4 <- broom::tidy(pgls4$aov.table)
  t5 <- broom::tidy(pgls5$aov.table)
  t6 <- broom::tidy(pgls6$aov.table)
  t7 <- broom::tidy(pgls7$aov.table)
  t8 <- broom::tidy(pgls8$aov.table)
  t9 <- broom::tidy(pgls9$aov.table)
  t10 <- broom::tidy(pgls10$aov.table)
  t11 <- broom::tidy(pgls11$aov.table)
  t12 <- broom::tidy(pgls12$aov.table)



  mergedtable <- rbind(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12)

  library(gt)
  mergedtablefancy <- gt(mergedtable) %>%
    cols_label(term = 'Source', statistic = 'F', p.value = 'p-value') %>%
    fmt_number(columns = c('SS', 'MS'), decimals = 0) %>%
    fmt_number(columns = c('Rsq', 'statistic', 'Z', 'p.value'), decimals = 4) %>%
    fmt_missing(columns = 2:8, missing_text = " ") %>% tab_header(title = region) %>%
    tab_footnote(footnote = "df, degree of freedom.", locations = cells_column_labels(columns = df)) %>%
    tab_footnote(footnote = "SS, sum of squares.", locations = cells_column_labels(columns = SS)) %>%
    tab_footnote(footnote = "MS, mean squares.", locations = cells_column_labels(columns = MS)) %>%
    tab_footnote(footnote = "Significant values (p-value < 0.05) from permutation tests (1,000 permutation rounds) are bolded.",
                 locations = cells_column_labels(columns = p.value)) %>%
    tab_row_group(label = 'Latitude', rows = 1:3) %>%
    tab_row_group(label = 'Temperature Seasonality', rows = 4:6) %>%
    tab_row_group(label = 'Precipitation Seasonality', rows = 7:9) %>%
    tab_row_group(label = 'Substrate', rows = 10:12) %>%
    tab_row_group(label = 'Microhabitat', rows = 13:15) %>%
    tab_row_group(label = 'Activity Time', rows = 16:18) %>%
    tab_row_group(label = 'Foraging Mode', rows = 19:21) %>%
    tab_row_group(label = 'Female SVL', rows = 22:24) %>%
    tab_row_group(label = 'Hatchling SVL', rows = 25:27) %>%
    tab_row_group(label = 'Reproductive Mode', rows = 28:30) %>%
    tab_row_group(label = 'Clutch Size', rows = 31:33) %>%
    tab_row_group(label = 'Clutch Size and Reproductive Mode', rows = 34:37)




  gtsave(mergedtablefancy, filename = outputname, exapnd = 10)

}

