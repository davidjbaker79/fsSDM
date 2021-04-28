#' Create multi-species occupancy distribution data
#'
#' The function create a data frame of species occupancy for a defined number of
#' species by defining the probability of occurrence as a function of Gaussian
#' responses to the axes of the principle components of the environmental data.
#'
#' @param enviro_dat Environmental data (requires \code{id}, \code{X}, \code{Y},
#' plus some environmental variables).
#' @param enviro_vars A named vector of environmental variables to be included.
#' @param n_spp The number of species to simulate.
#' @param pca_axes The number of pca axes to use in creating a species.
#' @param region_sf An sf object defining the spatial region.
#' @param plotRichness Create richness plot.
#' @param plotSpecies Create example species plot (plots species number 1).
#'
#' @import sf
#' @import ggplot2
#' @import stats
#'
#' @return A data frame
#'
#' @export
occupancy.multi <- function(enviro_dat, enviro_vars, n_spp = 50,
                            pca_axes = 2, region_sf,
                            plotRichness = FALSE, plotSpecies = FALSE) {

  # Coordinates and id
  idXY <- enviro_dat[,c("id", "X", "Y")]

  # For each species independently create a occupancy distribution based on
  # unique responses to climate and habitat gradients.
  # These gradients will be defined based on pca axes defined using a different
  # random subset of variables. n_spp = 50; sp = 1
  spp_sim <- lapply(1:n_spp, function(sp) {

    # Randomly choose environmental variable and run PCA
    pca_dat <- enviro_dat[ ,names(enviro_dat) %in% sample(enviro_vars, 12)]
    pca_dat[is.na(pca_dat)] <- 0
    pca_dat_res <- stats::prcomp(pca_dat, center = TRUE, scale = TRUE)
    pca_dat <- pca_dat_res$x[, 1:pca_axes]

    # Gaussian responses to PCA
    gaus_m <- rnorm(pca_axes, 0, 0.5)
    gaus_sd <- abs(rnorm(pca_axes, 0, 2))
    for (i in 1:pca_axes) {
      pca_dat[, i] <- dnorm(pca_dat[, i], mean = gaus_m[i], sd = gaus_sd[i])
    }

    # Calculate suitability
    p_occ <- apply(pca_dat, 1, function(x) prod(x) )
    p_occ <- round((p_occ - min(p_occ)) / (max(p_occ) - min(p_occ)), 4)

    # p1 <- ggplot() +
    #   geom_point(aes(x = pca_dat_res$x[,1], y = pca_dat_res$x[,2],colour = p_occ, alpha = p_occ), size = 1) +
    #   scale_colour_scico(palette = "lajolla", direction = -1, end = 0.8) +
    #   xlab("PC1") + ylab("PC2") +
    #   theme_minimal(base_size = 14) %+replace% theme(legend.position = 'right')
    # PC1 <- ggplot() +
    #   geom_line(aes(x = seq(-5,5,0.1), y = dnorm(seq(-5,5,0.1), gaus_m[1],  gaus_sd[1]))) +
    #   theme_minimal() %+replace% theme(axis.text = element_blank(),
    #                                    axis.title = element_blank(),
    #                                    panel.grid = element_blank())
    # PC2 <- ggplot() +
    #   geom_line(aes(x = seq(-5,5,0.1), y = dnorm(seq(-5,5,0.1), gaus_m[2],  gaus_sd[2]))) +
    #   theme_minimal() %+replace% theme(axis.text = element_blank(),
    #                                    axis.title = element_blank(),
    #                                    panel.grid = element_blank()) + coord_flip()
    # blank_p <- ggplot() + geom_blank(fill = "white")
    # p2 <- plot_grid(PC2, blank_p, ncol = 1, rel_heights = c(4.5,1.5))
    # p3 <- plot_grid(p1, PC1, ncol = 1, rel_heights = c(5,1))
    # p_all <- plot_grid(p2, p3, ncol = 2, rel_widths = c(1,5))
    # p_all
    # # ggsave(filename = "PCA_occ_gaus_plot.png", plot = p_all, width = 7, height = 7)

    # Convert to PA
    pa <- rbinom(n = length(p_occ), size = 1, prob = p_occ)

    # Join with XY coordinates
    pca_dat <- data.frame(id = enviro_dat$id, p_occ, pa)
    colnames(pca_dat)[2:3] <-  c(paste0("po_", sp), paste0("pa_", sp))

    #hex1 <- fsSDM::makeHexGrid(region_sf, 1)
    # spp_ex <- left_join(hex1, pca_dat)
    # sp_ex_p <- ggplot() +
    #   geom_sf(aes(fill = pa_1), alpha = 0.75, colour = "transparent", data = spp_ex) +
    #   scale_fill_scico(name = "POC", begin = 0.2, end = 1) +
    #   ggsn::north(region_sf, "topleft", 0.3, 16) +
    #   theme_minimal() %+replace% theme(axis.title = element_blank(),legend.position = 'none')
    # sp_ex_p
    # ggsave(filename = "PA_map_plot.png", plot = sp_ex_p, width = 7, height = 7)

    pca_dat

  })
  spp_sim <- Reduce(function(...) merge(..., by = "id", all = TRUE), spp_sim)
  spp_sim <- merge(idXY, spp_sim, by = "id", all = TRUE)

  # Create plots if required
  if (plotRichness) {

    # hex1 <- fsSDM::makeHexGrid(region_sf, 1)
    # spp_ex <- left_join(hex1, pca_dat)
    sim_sr <- spp_sim[,c(1:3, grep(names(spp_sim), pattern = "pa"))]
    sim_sr$sr <- rowSums(sim_sr[,-c(1:3)])
    sim_sr <- merge(hex1, sim_sr)
    sr_plot <-
      ggplot() +
        geom_sf(aes(fill = sr), alpha = 0.75, colour = "transparent", data = sim_sr) +
        #scale_fill_scico(name = "POC", begin = 0.2, end = 1) +
        #geom_sf(data = region_sf, fill = "grey95", colour = "grey95") +
        #geom_point(aes_string(x = 'X', y = 'Y', colour = 'sr'), data = sim_sr) +
        scico::scale_fill_scico(name = "#Spp", begin = 0.2, end = 1) +
        ggsn::north(region_sf, "topleft", 0.3, 16) +
        theme_minimal() %+replace% theme(axis.title = element_blank())
    ggsave(filename = "SR_map_plot.png", plot = sr_plot, width = 7, height = 7)
    print(sr_plot)

  }

  if (plotSpecies) {

    sim_sp1 <- spp_sim[spp_sim$pa_1 == 1, c("X", "Y", "pa_1")]
    sp_plot <-
      ggplot() +
        geom_sf(data = region_sf, fill = "grey95", colour = "grey95") +
        geom_point(aes_string(x = 'X', y = 'Y'),
                   colour = "Black", data = sim_sp1) +
        theme_minimal() %+replace% theme(axis.title = element_blank())

    print(sp_plot)

  }

  # Check calibration
  # calib <- spp_sim[spp_sim$pa_1 == 1,]
  # calib$po_1 <- round(calib$po_1, 1)
  # round(table(round( calib$po_1 , 1)) / table(round(spp_sim$po_1 , 1)),1)

  spp_sim

}

#' Virtual sampling function
#'
#' Virtual citizen scientists visit sites:
#' * Ideally, visits are spatially and temporally random, they record everything
#' they detect, record the area surveyed.
#' * and time spend in that area, and some index of their skill-level.
#'
#' @param multi_occ The output from \code{occupancy.multi}.
#' @param p_sites The proportion of sites available for sampling.
#' @param spat_bias Whether or not sampling is spatially biased.
#' @param n_yrs The number of years to simulate.
#' @param spat_bias The sampling spatially biased on random.
#' @param samp_weights A pre-generated data.frame of sample weights.
#' @param bias_level A number specifying the degree of spatial biased. The
#' higher the number the more spatially biased the sampling scheme.
#' @param rep_visit A number (exponent) controlling the distribution of repeat
#' visits to sites.
#' @param maxVisit The maximum number of visits to a site in a year.
#' @param reli_score How likely is an observer to record all species detected on
#'  a visit.
#'
#' @import gstat
#' @import stats
#'
#' @return A data.frame of visit x species detections.
#'
#' @export
virtual.sampling <- function(multi_occ, p_sites = 0.10, n_yrs = 10,
                            spat_bias = TRUE, samp_weights = NULL,
                            bias_level = 3, rep_visit = 2.5, maxVisit = 10,
                            reli_score = c("low", "medium", "high")) {

  # Remove unneeded columns from multi_occ
  multi_occ <- multi_occ[,c(1:3, grep(names(multi_occ), pattern = "pa_"))]

  # Number of species
  n_spp <- length(grep(names(multi_occ), pattern = "pa_"))

  # All available survey sites
  sites <- gtools::mixedsort(as.numeric(unique(multi_occ$id)))

  # How many potential sites are available
  n_sites <- length(sites)

  # Add a spatial bias to the spatial allocation of survey effort by creating a spatially correlated
  # random field (following Begueria) and then converting to probability weighting for sampling survey sites
  if (spat_bias) {

    if (is.null(samp_weights)) {

      # Generate a sampling probability surface with bias a function of x and y coordinates
      g.dummy <- gstat(formula = z ~ 1 + X + Y, locations = ~ X + Y, dummy = TRUE, beta = 1, model = vgm(psill = 0.05, range = 2000, model = "Exp"), nmax = 25)
      samp_weights <- scales::rescale(predict(g.dummy, newdata = multi_occ[, c("X", "Y")], nsim = 1)[, 3], to = c(0, 1))
      samp_weights <- samp_weights^bias_level #
      # bias <- ggplot() +
      #   geom_point(aes(x = multi_occ$X, y = multi_occ$Y, colour = samp_weights)) +
      #   scico::scale_colour_scico(name = "#Spp", begin = 0, end = 1) +
      #   ggsn::north(region_sf, "topleft", 0.3, 16) +
      #   theme_minimal() %+replace% theme(axis.title = element_blank(), legend.position = 'none')
      # ggsave(filename = "Bias_plot.png", plot = bias, width = 7, height = 7)


    } else {

      samp_weights <- samp_weights^bias_level

    }

    multi_occ <- cbind(samp_weights, multi_occ)
    #ggplot(multi_occ) + geom_point(aes(x = X, y = Y, colour = samp_weights)) + scale_colour_scico()

  } else {

    samp_weights <- NULL

  }

  # For each year sample the survey sites yr = 1
  sample_yrs <- lapply(1:n_yrs, function(yr) {

    # Sample sites either randomly or with spatial bias
    surveySites <- sample(sites, (n_sites * p_sites), replace = FALSE, prob = samp_weights)

    # Take random proportion of sites, defined by power function scaled by rep_visit, and given them a repeat visit each year
    # Which sites visited more than once -
    probVisits <- seq(1, maxVisit)^-rep_visit # This bit gets the probability of sites getting 1 - maxVisit
    surveySitesVisits <- vector("list", maxVisit)

    if (maxVisit > 1) {

      for (i in 2:maxVisit) {

        n <- floor(n_sites * p_sites * probVisits[i])
        surveySitesAddVis <- sample(surveySites, n, replace = FALSE)
        surveySitesVisits[[i]] <- data.frame(surveySites = gtools::mixedsort(rep(surveySitesAddVis, i)), visitNo = rep(1:i, n))
        surveySites <- surveySites[!(surveySites %in% surveySitesAddVis)]

      }

    }
    surveySitesVisits[[1]] <- data.frame(surveySites = gtools::mixedsort(surveySites), visitNo = rep(1, length(surveySites)))
    surveySitesVisits <- do.call(rbind, surveySitesVisits)

    # Extract survey data for each visit
    multi_occ_s <- multi_occ[surveySitesVisits$surveySites, ] #
    # ggplot() + geom_point(aes(x = X, y = Y), colour = 'grey95', data = multi_occ) + geom_point(aes(x = X, y = Y), size=3 , data = multi_occ_s) +
    #  scale_colour_scico() + theme_minimal()
    multi_occ_s <- cbind(surveySitesVisits, multi_occ_s)
    multi_occ_s <- multi_occ_s[gtools::mixedorder(multi_occ_s$surveySites), ]

    # Assign detection failures based on detection probability
    # First, draw values of detectability from a beta distribution
    det_prob_spp <- round(rbeta(n_spp, 1, 0.4), 1)

    for (sp in 1:n_spp) {

      # Species occupancy
      sp_occ <- multi_occ_s[, paste0("pa_", sp)]
      # Create detection column for species
      multi_occ_s[paste0("det_", sp)] <- sp_occ
      # Identify occupied survey sites
      occ <- which(sp_occ == 1)
      # Assign detection failure based on detection probability
      detectFail <- occ[rbinom(length(occ), 1, det_prob_spp[sp]) == 0]
      # Assign detection failure
      multi_occ_s[[paste0("det_", sp)]][detectFail] <- 0

    }

    # Add year and order, remove pa columns as these exist in the input data
    multi_occ_s$yr <- yr
    multi_occ_s <-
      multi_occ_s[,c("surveySites", "visitNo", "samp_weights",
                     "id", "X", "Y", "yr",
                     names(multi_occ_s)[grep(names(multi_occ_s),
                                             pattern = "det_")])]

    # Surveyor have differing probabilities of reporting a complete species list
    if (reli_score == "low") {
      a_r <- 0.1
      b_r <- 1
    } else if (reli_score == "medium") {
      a_r <- 0.5
      b_r <- 0.4
    } else if (reli_score == "high") {
      a_r <- 1
      b_r <- 0.1
    }
    reli_score <- round(rbeta(nrow(multi_occ_s), a_r, b_r), 1) # probability of reporting complete list #

    # Assign charisma scores to species, which affects the likelihood of them being recorded
    charima_scores_spp <- round(rbeta(n_spp, 1.25, 3), 1) # hist(round(rbeta(nrow(multi_occ_s), 1.25, 3), 1))

    det_v <- grep(names(multi_occ_s), pattern = "det_") # cols holding det data
    for (v in 1:nrow(multi_occ_s)) {

      if (reli_score[v] < runif(1, 0, 1)) {

        # Which of the detected species is recorded
        rep_p <- rbinom(length(charima_scores_spp), 1, charima_scores_spp)

        # Modify detected species to give reported species
        multi_occ_s[v, det_v] <- multi_occ_s[v, det_v] * rep_p

      }

    }

    return(multi_occ_s)

  })

  multi_occ_yrs <- do.call(rbind, sample_yrs)

}

#' Plot species and sampling
#'
#' @param multi_occ The output from \code{occupancy.multi}.
#' @param multi_det The output from \code{virual.sampling}.
#' @param sp A numeric value indicating the species to plot.
#' @param year A numeric value indicating the year to plot.
#'
#' @import ggplot2
#'
#' @return A map of the 'real' distribution and the sample, with sampling
#' weighting as shown as the background.
#'
#' @export
# multi_occ <- spp50
# multi_det <- surv50
plot.virtual.species <- function(multi_occ, multi_det, sp = 1, year = 1:10) {

  sp_n <- paste0("det_", sp)
  det <- multi_det[multi_det[[sp_n]] == 1 & multi_det$yr %in% year, ]
  det <- det[, c("id", "X", "Y", sp_n)]
  names(det)[4] <- "det"
  det <- ddply(det, .(id, X, Y), summarise, det = max(det))

  p1 <-
    ggplot() +
      geom_point(aes_string(x = "X", y = "Y", colour = "samp_weights"),
                 data = multi_det) +
      scico::scale_colour_scico(name = "Sample\nweights") +
      geom_point(aes_string(x = "X", y = "Y"), colour = "red", size = 2,
                 data = det) +
      geom_point(aes_string(x = "X", y = "Y"), colour = "blue", size = 0.5,
                 data = multi_occ[multi_occ[paste0("pa_", sp)] == 1, ]) +
      theme_minimal() %+replace% theme(axis.title = element_blank())

  p1

}
#'
#'
#' #~#~# #~#~# TEST APPLICATION #~#~# #~#~#
#'
#' #~# Generate virtual species
#' library(tidyverse)
#' library(sf)
#' library(rgdal)
#' library(gstat)
#' library(ggsn)
#' library(scico)
#' library(cowplot)
#'
#' #' Shared data path
#' shareDir <- "C:/Users/DB625/OneDrive - University of Exeter/Fellowship_data/"
#'
#' #' UK shapefile
#' set_transform_wkt_comment(TRUE)
#' cornwall <- readOGR(paste0(shareDir,
#'                            "EnglishCountiesCornwall_2011"),"england_ct_2011")
#' cornwall_sf <- st_as_sf(cornwall)
#'
#' #' Load CEH land use data at 1 km2 with buffer 0 and 1
#' land_cover_c19 <- get(load(paste0(shareDir,"/cornwall_ceh2019_landcover_1km.Rdata")))
#'
#' #' Load bioclimate variables
#' bioClim <- get(load(paste0(shareDir,"/cornwall_climvars_1km_means.Rdata")))
#'
#' #' Merge into single enviro data frame
#' enviro_dat <- left_join(land_cover_c19, bioClim, by = c('id', 'X', 'Y', 'grain'))
#'
#' #' Vector of environmental variables
#' enviro_vars <- c("LC_1_0", "LC_3_0", "LC_4_0", "LC_7_0", "LC_9_0", "LC_10_0",
#'                  "LC_11_0", "LC_14_0", "LC_15_0","LC_16_0", "LC_17_0", "LC_18_0",
#'                  "LC_19_0", "LC_20_0", "LC_21_0", "temp_ann_m", "temp_gs_m",
#'                  "mtcm_m", "mtcq_m", "frost_m", "precip_gs_m", "soilm_gs_m")
#' multi_occ = enviro_dat
#' n_spp = 50
#' pca_axes = 3
#' region_sf = cornwall_sf
#'
#' #' Generate multi-species occupancy
#' spp_occ <- occupancy.multi(enviro_dat, enviro_vars, n_spp = 45, pca_axes = 3, cornwall_sf, plotRichness = T, plotSpecies = F)
#'
#' #~# Generate virtual sampling
#' s1 <- virtual.sampling(spp_occ, p_sites = 0.10,  n_yrs = 30, spat_bias = TRUE, samp_weights = NULL, bias_level = 3,
#'                       rep_visit = 2.5, maxVisit = 10, reli_score = 'medium')
#'
#' #'
#' sample_sr <- s1 %>%
#'   dplyr::select(c('id','X','Y',grep(names(.), pattern = 'det_'))) %>%
#'   dplyr::rowwise() %>%
#'   dplyr::mutate(sr = sum(across(starts_with("det_")))) %>%
#'   dplyr::select(c('id','X','Y','sr')) %>%
#'   ggplot() +
#'   geom_sf(data=cornwall_sf) +
#'   geom_point(aes(x = X, y =Y, colour = sr)) +
#'   scale_colour_scico(name = 'SR', begin = 0.2) +
#'   theme_minimal(base_size = 14) %+replace% theme(axis.title = element_blank()) +
#'   north(cornwall_sf, "topleft", 0.3, 16)
#' print(sample_sr)
#'
#'
#' real_sr %>% filter(sr == 16)
#'
#' #' Test rarefaction methods
#' library(rareNMtests)
#' #'
#' tmp <- s1 %>% filter(id == 1626) %>% dplyr::select(grep(names(.),pattern = 'det_'))
#' rarefaction.sample(tmp, q = 0, method = "coverage")
#'
#'
#'
#'
#' #' Example species plots
#' plot.virtual.species(spp_occ, s1, 1, 20)
#' plot.virtual.species(spp_occ, s1, 40, 4)
#'
#' #' Summarise distribution of effort over cells
#' s1 %>% dplyr::group_by(id) %>% dplyr::summarise(N = n()) %>%
#'   ggplot() + geom_histogram(aes(x = N), colour = 'White', binwidth = 1) +
#'   theme_bw(base_size = 14)
#'
#' #' Summarise distribution of visits
#' s1 %>% as_tibble %>% dplyr::group_by(id, yr) %>% dplyr::summarise(n_visits = max(visitNo))%>%
#'   ggplot() + geom_histogram(aes(x = n_visits)) +
#'   theme_bw() + facet_wrap(~yr)
#'
#' #' Target group bias
#' sr_sim <- s1 %>%
#'   dplyr::select(c('id', 'samp_weights','X','Y',grep(names(.), pattern = 'det_'))) %>%
#'   dplyr::rowwise() %>%
#'   dplyr::mutate(sr = sum(across(starts_with("det_")))) %>%
#'   dplyr::group_by(id, X, Y, samp_weights) %>%
#'   dplyr::summarise(sr = sum(sr))
#' sr_sim$sr <- sr_sim$sr / max(sr_sim$sr, na.rm = TRUE)
#'   wt_plot <- ggplot(sr_sim) +
#'     geom_point(aes(x = X, y =Y, colour = samp_weights)) +
#'     scale_colour_scico(name = 'SR') +
#'     theme_minimal(base_size = 14) %+replace% theme(axis.title = element_blank()) +
#'     north(region_sf, "topleft", 0.3, 16)
#'   print(wt_plot)
#'
#' ggplot() + geom_point(aes(x = sr_sim$sr, y = sr_sim$samp_weights))
#'
#'
#'
#'
