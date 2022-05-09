## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, results = "hide", message = FALSE)
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))

## -----------------------------------------------------------------------------
# load the package
library(fossilbrush)
# load the example datasets
data(brachios)
data(sepkoski)
# trim the Sepkoski Compendium to the relevant entries
sepkoski <- sepkoski[which(sepkoski$PHYLUM == "Brachiopoda"),]

## -----------------------------------------------------------------------------
# update chronostratigraphy
brachios <- GTS2020_scale(brachios, srt = "early_interval", end = "late_interval",
                          max_ma = "max_ma", min_ma = "min_ma", verbose = FALSE)
# brachios <- get_pbdb("Brachiopoda", interval = "Palaeozoic, tscale = "GTS2020")

## ---- results='markup', message=TRUE------------------------------------------
# combine the datasets
occs <- cbind.data.frame(phylum = c(brachios$phylum, sepkoski$PHYLUM),
                         class = c(brachios$class, sepkoski$CLASS),
                         order = c(brachios$order, sepkoski$ORDER), 
                         family = c(brachios$family, rep(NA, nrow(sepkoski))),
                         genus = c(brachios$genus, sepkoski$GENUS),
                         max_ma = c(brachios$GTS_FAD, sepkoski$RANGE_BASE),
                         min_ma = c(brachios$GTS_LAD, sepkoski$RANGE_TOP),
                         coll_no = c(brachios$collection_no, rep(NA, nrow(sepkoski))))
# define the taxonomic ranks used in the dataset (re-used elsewhere)
b_ranks <- c("phylum", "class", "order", "family", "genus")
# define a list of suffixes to be used at each taxonomic level when scanning for synonyms
b_suff = list(NULL, NULL, NULL, NULL, c("ina", "ella", "etta"))
# scan for errors
occs_c <- check_taxonomy(occs, suff_set = b_suff, ranks = b_ranks)

## -----------------------------------------------------------------------------
# manually resolve the clear synonymous spellings
occs$family[which(occs$family == "Disciniidae")] <- "Discinidae"
occs$genus[which(occs$genus == "Ptychomalotoechia")] <- "Ptychomaletoechia"
occs$genus[which(occs$genus == "Hipparionix")] <- "Hipparionyx"
occs$genus[which(occs$genus == "Sphaenospira")] <- "Sphenospira"
# strip out the PBDB 'missing taxon' format
for(i in c("phylum", "class", "order", "family", "genus")) {
  occs[grep("^NO_", occs[,i]),i] <- NA
}

## -----------------------------------------------------------------------------
# clean the data, this time resolving classifications
occs_c <- check_taxonomy(occs, suff_set = b_suff, ranks = b_ranks, verbose = FALSE,
                         clean_name = TRUE, resolve_duplicates = TRUE, jump = 5)
# plot a taxon before and after cleaning to confirm that it is correct
par(mfrow = c(1, 2))
plot_taxa(occs, "Atrypa", trank = "genus", ranks = b_ranks, mode = "parent")
plot_taxa(occs_c$data, "Atrypa", trank = "genus", ranks = b_ranks, mode = "parent")

## -----------------------------------------------------------------------------
# extract PBDB
sepkoski_c <- occs_c$data[(nrow(brachios) + 1):nrow(occs_c$data),]
# extract Sepkoski
brachios_c <- occs_c$data[1:nrow(brachios),]

## -----------------------------------------------------------------------------
# drop occurrences with older LADs than FADs
brachios_c <- brachios_c[brachios_c$max_ma > brachios_c$min_ma,]
# chunk to a small size for better run time
set.seed(1)
samp <- sample(1:nrow(brachios_c), 1000)
# flag and resolve against the Sepkoski Compendium, collection-wise
revrng <- revise_ranges(x = brachios_c[samp,], y = sepkoski_c, do.flag = TRUE, verbose = F,
                        taxon = "genus", assemblage = "coll_no",
                        srt = "max_ma", end = "min_ma")
# append the revised occurrence ages and error codes to the dataset
brachios_c$newfad <- brachios_c$newlad <- brachios_c$errcode <- NA
brachios_c$newfad[samp] <- revrng$occurrence$FAD
brachios_c$newlad[samp] <- revrng$occurrence$LAD
brachios_c$errcode[samp] <- revrng$occurrence$tax_flag

## -----------------------------------------------------------------------------
# densify ranges
dens <- densify(brachios_c)
# plot an example taxon
plot_dprofile(dens, "Atrypa")

## -----------------------------------------------------------------------------
# pacmacro trimming
pacm <- pacmacro_ranges(brachios_c, tail.flag = c(0.3, 0.35, 0.4),
                        rank = "genus", srt = "max_ma", end = "min_ma")
# replot the taxon and mark its truncated ranges
plot_dprofile(dens, "Atrypa")
abline(v = pacm$kdensity["Atrypa","FAD95"], col = "blue")
abline(v = pacm$kdensity["Atrypa","LAD95"], col = "blue")

## -----------------------------------------------------------------------------
# extract the truncated ranges
pranges <- pacm$kdensity
# for those identified as anomalous (tflag = 1), update the range values to the 95% trim
pranges$FAD[which(pranges$tflag0.35 == 1)] <- pranges$FAD95[which(pranges$tflag0.35 == 1)]
pranges$LAD[which(pranges$tflag0.35 == 1)] <- pranges$LAD95[which(pranges$tflag0.35 == 1)]
# format the range table
pranges <- cbind.data.frame(genus = rownames(pranges),
                            max_ma = pranges$FAD,
                            min_ma = pranges$LAD)
# perform the flagging and append to the dataset
pflags <- flag_ranges(brachios_c, pranges, verbose = FALSE)
brachios_c$pflag <- pflags$occurrence$status

## -----------------------------------------------------------------------------
# interpeak thresholding
itp <- threshold_ranges(brachios_c, win = 8, thresh = 10,
                        rank = "genus", srt = "max_ma", end = "min_ma")
# append the stratigraphically thresholded taxon names to the dataset
brachios_c$newgen <- itp$data
# plot the taxon, now identifying the peaks
plot_dprofile(dens, "Atrypa")
add_itp(itp, "Atrypa")

