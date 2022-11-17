

allometry_plot <- function(brain, tree, traitdata, region, outputname){
library(gt)
library(geomorph)
library(knitr)
gpa <- gpagen(brain, ProcD = T, verbose = T)
gdf <- geomorph.data.frame(gpa, traitdata)


fit <- procD.pgls(coords ~ log10(Csize), phy = lzrd_tree, data=gdf, iter=999,
                  print.progress = FALSE)

pdf(file = outputname, width = 4, height = 4)
plot(fit, type="regression", reg.type="RegScore",
     predictor = log10(gdf$Csize))
dev.off()

# pdf(file = outputname, width = 4, height = 4)
# plotAllometry(fit, size = gdf$Csize, logsz = T, method = 'RegScore', pch = 19)
# dev.off()

}




# https://cran.r-project.org/web/packages/geomorph/vignettes/geomorph.assistance.html

allometry_test <- function(brain, tree, traitdata, region, outputname) {

  gpa <- gpagen(brain, ProcD = T, verbose = T)
  gdf <- geomorph.data.frame(gpa, traitdata)

  fit <- procD.pgls(coords ~ log10(Csize), phy = tree, data=gdf, iter=999,
                    print.progress = FALSE)

  t <- broom::tidy(fit$aov.table)

  mergedtablefancy <- gt(t) %>%
    cols_label(term = 'Source', statistic = 'F', p.value = 'p-value') %>%
    fmt_number(columns = c('SS', 'MS', 'Rsq', 'statistic', 'Z', 'p.value'), decimals = 4) %>%
    fmt_missing(columns = 2:8, missing_text = " ") %>% tab_header(title = region) %>%
    tab_footnote(footnote = "df, degree of freedom.", locations = cells_column_labels(columns = df)) %>%
    tab_footnote(footnote = "SS, sum of squares.", locations = cells_column_labels(columns = SS)) %>%
    tab_footnote(footnote = "MS, mean squares.", locations = cells_column_labels(columns = MS)) %>%
    tab_footnote(footnote = "Significant values (p-value < 0.05) from permutation tests (1,000 permutation rounds) are bolded.",
                 locations = cells_column_labels(columns = p.value))

  gtsave(mergedtablefancy, filename = outputname, exapnd = 10)

}


# not liking this because we can't take into account phylogentic signal and remove it
two.block.plot <- function(brain, tree, traitdata, region, outputname) {

  gpa <- gpagen(brain, ProcD = T, verbose = T)
  gdf <- geomorph.data.frame(gpa, traitdata)

  fit <- procD.pgls(coords ~ log10(Csize), phy = tree, data=gdf, iter=999,
                    print.progress = FALSE)


  PLS <- two.b.pls(log(gdf$Csize), gdf$coords, print.progress = FALSE)
  summary(PLS)
  pdf(file = outputname, width = 4, height = 4)
  plot(PLS)
  dev.off()

}
