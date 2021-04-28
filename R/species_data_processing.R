#' Load ERCCIS data
#'
#' Create a connection with the ERCCIS SQLite db and
#' filter missing coordinates and create required date formats
#'
#' @param sqlite_file A path name to .sqlite file.
#' @param table_name The name of .sqlite table to load (if multiple)
#' @param sppGrps Species group required.
#' @param sppNameLookUp A table with species names can be provided for filtering. A column can also
#' be added to add better names (e.g. for plotting).
#'
#' @import sf
#'
#' @return An tibble data.frame
#'
#' @export
erccis.db.tbl <- function(sqlite_file, table_name, sppGrps = NULL, sppNameLookUp = NULL) {

  # Create db connection
  sq_con <- RSQLite::dbConnect(drv = RSQLite::SQLite(), dbname = sqlite_file)

  # Select species group
  sq_dat <- RSQLite::dbGetQuery(sq_con, 'SELECT * FROM erccis_oct2020 WHERE "SpeciesGroup" == :x', params = list(x = sppGrps))
  if(is.null(sppGrps) | nrow(sq_dat) == 0) {
    sppGrps <- RSQLite::dbGetQuery(sq_con, 'SELECT DISTINCT SpeciesGroup FROM erccis_oct2020')
    cat('Please enter one or more of the following species group names')
    print(sppGrps)
  }

  # Remove records with no spatial coordinates (typically a small percentage)
  names(sq_dat)[which(names(sq_dat) == "SpatialLong")] <- 'X'
  names(sq_dat)[which(names(sq_dat) == "SpatialLat")] <- 'Y'
  sq_dat <- sq_dat[!is.na("X") & !is.na("Y"), ]
  sq_dat <- sq_dat[sq_dat$SpeciesTaxonLevel == "Species",]
  sq_dat <- sq_dat[!is.na(sq_dat$SpeciesGroup),]
  sq_dat$StartDate <- as.Date(sq_dat$StartDate, format = "%Y-%m-%d")
  sq_dat$EndDate <- as.Date(sq_dat$EndDate, format = "%Y-%m-%d")

  # Add in better names if required
  if (!is.null(sppNameLookUp)) sq_dat <- merge(sq_dat, sppNameLookUp, by = "SpeciesGroup", all.x = TRUE)
  sq_dat

}

#' Re-project point data
#'
#' This function uses a reference spatial polygon to transform the point data to a new projection
#' and returns an sf spatial object.
#'
#' @param bioDat A data frame with SpatialLong given as \code{X} and SpatialLat given as \code{Y}.
#' @param crs_sf_ref An sf object with the desired crs.
#' @param region A region defined within  \code{crs_sf_ref}.
#'
#' @import sf
#'
#' @return An sf object.
#'
#' @export
#' bioDat <- grp; crs_sf_ref <- cornwall_sf; region = 'Cornwall'; spatial = FALSE
db.transform.proj <- function(bioDat, crs_sf_ref, region = NULL, df = FALSE) {

  # Convert to sf
  dat_sf <- st_as_sf(
    x = bioDat,
    coords = c("X", "Y"),
    crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  )

  # Spatial transform and mask out records off the coast and Isle of Scilly (also few in Devon)
  dat_sf <- st_transform(dat_sf, crs = st_crs(crs_sf_ref))
  dat_sf <- st_join(dat_sf, crs_sf_ref, join = st_intersects)

  if (!is.null(region)) dat_sf <- dat_sf[dat_sf$name %in% region, ]

  if(df) {
    dat_df <- st_drop_geometry(dat_sf)
    dat_df$lon <- st_coordinates(dat_sf)[,1]
    dat_df$lat <- st_coordinates(dat_sf)[,2]
    dat_sf <- dat_df
  }
  dat_sf

}

#' Obtain species group data
#'
#' Wrapper function for filtering by species groups and projecting point data
#'
#' @param bioDat A data frame with SpatialLong given as \code{X} and SpatialLat given as \code{Y}.
#' @param sppGrp Name of species group (unwise to run on all species!!!).
#' @param region_sf An sf object defining the spatial region.
#' @param regionSelect A region defined within  \code{crs_sf_ref}.
#' @param sppSelect A SpeciesVenacular name or vector of names if just a subset of species is wanted.
#'
#' @return A filtered and transformed data frame.
#'
#' @export
species.group.select <- function(bioDat, sppGrp, region_sf, regionSelect = NULL, sppSelect = NULL) {

  # Subset to just species group
  bioDat <- bioDat[bioDat$SpeciesGroup %in% sppGrp, ]

  # If we only want a subset of species from a group
  if (!is.null(sppSelect)) bioDat <- bioDat[bioDat$SpeciesVenacular %in% sppSelect, ]

  # Create spatial object from db subset and transform to BNG
  bioDat <- db.transform.proj(bioDat, region_sf, regionSelect)

}

#' Plots for evaluating species distributions
#'
#' @param spp_dat An export from \code{db.transform.proj}
#' @param region_sf An sf object defining the spatial region.
#' @param focal_sp The SpeciesVenacular name.
#' @param start_yr A numeric year defining the start of the focal period.
#' @param end_yr A numeric year defining the end of the focal period.
#' @param plotType One of Map, Hist, or listLength.
#'
#' @import data.table
#' @import sf
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @return A ggplot figure.
#'
#' @export
#spp_dat = grpPrj; focal_sp =  "Hazel Dormouse"; start_yr = 2000; end_yr = 2019; region_sf = cornwall_sf; plotType = 'ListLength'
species.evaluation <- function(spp_dat, focal_sp, start_yr = 2001, end_yr = NULL,
                               region_sf, plotType = c("Map", "Hist", "ListLength")) {

  # Subset time period
  if (!is.null(start_yr) & !is.null(end_yr)) spp_dat <- spp_dat[spp_dat$Yr >= start_yr & spp_dat$Yr <= end_yr, ]
  if ( plotType !=  "Map" ) spp_dat <-  data.table(spp_dat)
  fspp <- spp_dat[spp_dat$SpeciesVenacular == focal_sp, ]

  # Plot target species occurrence data
  if (plotType == "Map") {

    plotOut <-
      ggplot() +
        geom_sf(fill = "grey75", data = region_sf) +
        geom_sf(aes_string(colour = "Yr"), alpha = 0.5, size = 1, data = fspp) +
        scico::scale_colour_scico(palette = "bamako", direction = -1, name = "Year") +
        theme_minimal(base_size = 14) %+replace% theme(axis.title = element_blank()) +
        ggsn::north(region_sf, "topleft", 0.3, 16) +
        ggtitle(focal_sp)

  }

  # Temporal trend in regional record
  if (plotType == "Hist") {

    # Convert to data.table
    fspp_n <- fspp[, .N, by = list(Yr)]

    plotOut <-
      ggplot(aes_string(x = "Yr", y = "N"), data = fspp_n) +
        geom_col(fill = "grey75") +
        theme_bw(base_size = 14) +
        xlab("Year") +
        ylab("Number of Records") +
        ggtitle(focal_sp)

  }

  # List length per year
  if (plotType == "ListLength") {

    sp_grp <- spp_dat[, .N, by = list(SpatialReference, StartDate, Yr)]
    ll_n <- spGrp[, .N, by = list(Yr, N)]
    names(ll_n) <- c('Yr', 'll', 'N')


    sp_foc <- fspp[, .N, by = list(SpatialReference, StartDate, Yr)]
    foc_n <- sp_foc[, .N, by = list(Yr, N)]
    names(foc_n) <- c('Yr', 'll', 'N')

    plotOut <-
      ggplot() +
        geom_point(aes_string(x = "ll", y = "log10(N)", group = "Yr"), colour = "black", data = ll_n) +
        geom_line(aes_string(x = "ll", y = "log10(N)", group = "Yr"), colour = "grey90", data = ll_n) +
        geom_point(aes_string(x = "ll", y = "log10(N)", group = "Yr"), colour = "red", data = foc_n) +
        theme_bw(base_size = 14) +
        xlab("List length") + ylab("Number of lists (log10)") +
        ggtitle(focal_sp)

  }

  plotOut

}


#' Make hexagonal grid
#'
#' Analysis is conducted on hexagonal grids and this function creates the grid with a specified area.
#'
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea A numeric value specifying the area of the hexagonal cell (in km2).
#'
#' @import sf
#'
#' @return An sf object with 'id' used as a unique cell id.
#'
#' @export
makeHexGrid <- function(region_sf, hexCellArea) {

  hexGrid <- st_make_grid(region_sf, cellsize = (1074.57 * sqrt(hexCellArea)), square = FALSE)
  hexGrid <- st_as_sf(hexGrid)
  hexGrid$id <- rownames(hexGrid)
  hexGrid$grain <- hexCellArea
  hexGrid <- st_intersection(hexGrid, region_sf)
  hexGrid$area <- round(st_area(hexGrid) / 10^6, 2)

  return(hexGrid)

}

#' Assign to hexagonal grid by any time slice
#'
#' Function assigns records to grid cells based on spatial overlap.
#'
#' @param bioDat A data frame with SpatialLong given as \code{X} and SpatialLat given as \code{Y}.
#' @param timeSlice A numeric value specifying the length on each timeslice (in years).
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea A numeric value specifying the area of the hexagonal cell (in km2).
#' @param hexGrid A pre-made hexagonal grid, as producted by \code{makeHexGrid}.
#'
#' @import sf
#'
#' @return A tibble data frame with hexagonal cell id and X, Y coordinates.
#'
#' @export
recToHexSetUp <- function(bioDat, timeSlice, region_sf, hexCellArea, hexGrid = NULL) {

  # Create time slices based on intervals
  start_date <- seq.Date(as.Date("1960/01/01"), as.Date("2020/01/01"), "year")[seq(1, 61, timeSlice)]
  start_date <- start_date[!is.na(start_date)]

  # Create TS labels
  TSLab <- sub("-01-01", "", start_date)
  TSLab <- paste(TSLab[-length(start_date)], TSLab[-1], sep = "_")

  # For each time slice filter bioDat into time slices using start and end dates
  bioDat[["TS"]] <-  cut(bioDat$StartDate, breaks = start_date, right = TRUE, labels = TSLab)
  bioDat <- bioDat[!is.na(bioDat$TS),]
  bioDat <- bioDat[,!(names(bioDat) %in% c("label", "name", "code"))]

  # Create hexagonal grid with id attribute
  if (is.null(hexGrid)) hexGrid <- makeHexGrid(region_sf, hexCellArea)

  # Join biological records to hexGrid
  bioHex <- st_join(hexGrid, bioDat, join = st_intersects)

  # Add in XY
  idXY <- data.frame(st_coordinates(st_centroid(hexGrid)), id = hexGrid$id)
  bioHex <- merge(bioHex, idXY, by = 'id')

}

#' Number of records \code{x} hexagonal grid \code{x} time slice
#'
#' This function summarises number of records per cell and time slice.
#' This is useful for exploring trends in record collection for identifying periods when
#' records may be more or less representative.
#'
#' @param bioDat A data frame with SpatialLong given as X and SpatialLat given as Y.
#' @param sppGrp Name of species group (unwise to run on all species!!!).
#' @param timeSlice A numeric value specifying the length on each timeslice (in years).
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea A numeric value specifying the area of the hexagonal cell (in km2).
#' @param bioHex An export from \code{recToHexSetUp}.
#' @param hexGrid A pre-made hexagonal grid, as producted by \code{makeHexGrid}.
#'
#' @import sf
#' @import data.table
#' @importFrom rlang .data
#'
#' @return A data frame with hexagonal cell id and X, Y coordinates.
#'
#' @export
recXgridXtime <- function(bioDat, sppGrp, timeSlice, region_sf, hexCellArea, bioHex = NULL, hexGrid = NULL) {

  if (is.null(bioHex)) bioHex <- recToHexSetUp(bioDat, timeSlice, region_sf, hexCellArea, hexGrid)

  # Summarise by id and TS
  nRecHex <- st_drop_geometry(bioHex)
  nRecHex <- data.table(nRecHex)
  nRecHex <- nRecHex[, .N, by = list(id, X, Y, TS)]
  nRecHex$sppGrp <- sppGrp
  nRecHex <- as.data.frame(nRecHex)

  return(nRecHex)

}

#' Site \code{x} Species \code{x} Time data frame
#'
#' This function creates a site \code{x} species \code{x} time matrix
#' (presence/absence), which can then be subset into specific time slices.
#'
#' @param bioDat A data frame with SpatialLong given as X and SpatialLat given as Y.
#' @param timeSlice A numeric value specifying the length on each timeslice (in years).
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea Area of each hexagonal grid cell.
#' @param bioHex An export from \code{recToHexSetUp}.
#' @param hexGrid A pre-made hexagonal grid, as produced by \code{makeHexGrid}.
#' @param randTS If TRUE the date of surveys are randomised within each cell.
#'
#' @import data.table
#' @import sf
#' @importFrom rlang .data
#'
#' @return A binary (0,1) presence absence data frame is returned
#'
#' @export
siteXSpeciesXTime <- function(bioDat, timeSlice, region_sf, hexCellArea,
                              bioHex = NULL, hexGrid = NULL) {

  # Call setup function to create TS and hexGrid association
  if (is.null(bioHex)) bioHex <- recToHexSetUp(bioDat, timeSlice, region_sf, hexCellArea, hexGrid)

  # Summarise
  sumHex <- st_drop_geometry(bioHex)
  sumHex <- data.table(sumHex)
  sumHex <- sumHex[, .N, by = list(id, TS, SpeciesScientific, SpatialReference)]
  sumHex$siteTime <- paste0(sumHex$id, "_", sumHex$TS, "_", sumHex$SpatialReference)
  sumHex <- sumHex[,c("siteTime", "SpeciesScientific", "N")]
  sumHex <- as.data.frame(sumHex)

  # Convert to site x species df
  siteXSpp <- tidyr::pivot_wider(sumHex, names_from = SpeciesScientific, values_from = N)
  siteXSpp[is.na(siteXSpp)] <- 0
  siteXSpp <- as.data.frame(siteXSpp)

}

#' Focal species select
#'
#' @param bioDat A data frame with SpatialLong given as X and SpatialLat given as Y.
#' @param focal_sp A focal species to select.
#' @param region_sf An sf object defining the spatial region.
#' @param grain Area of each hexagonal grid cell.
#' @param hexGrid A pre-made hexagonal grid, as produced by \code{makeHexGrid}.
#' @param collapse_time Should all year's data be collapsed into a single presences/absences data frame
#'
#' @import sf
#' @import plyr
#'
#' @return A tibble data frame holding data for \code{focal_sp}.
#'
#' @export
focal.species.select <- function(bioDat, focal_sp, region_sf, grain, hexGrid = NULL, collapse_time = TRUE) {

  # Create hexagonal grid with id attribute
  if (is.null(hexGrid)) hexGrid <- makeHexGrid(region_sf, hexCellArea = grain)

  # Subset to just species group
  bioDat <- bioDat[bioDat$SpeciesVenacular == focal_sp, ]

  # Grid species data
  sppHexGrid <- recToHexSetUp(bioDat, 1, region_sf, grain, hexGrid)

  # Create smaller output dataframe
  bioDat <- st_drop_geometry(sppHexGrid)
  bioDat <- bioDat[,c("id", "X", "Y", "grain", "SpeciesVenacular", "Mth", "Yr", "Precision_")]
  names(bioDat)[which(names(bioDat) == "SpeciesVenacular")] <- "species"
  names(bioDat)[which(names(bioDat) == "Mth")] <- "mth"
  names(bioDat)[which(names(bioDat) == "Yr")] <- "yr"
  names(bioDat)[which(names(bioDat) == "Precision_")] <- "precision"
  bioDat <- bioDat[!is.na(bioDat$species),]

  if (collapse_time) {

    bioDat <- ddply(bioDat, c("id", "X", "Y", "grain", "species"), summarise,
                    mth = stats::median(mth),
                    yr = stats::median(yr),
                    precision = min(precision))

  }

}

#' Target group weightings for modelling
#'
#' Create target group for evaluating bias in data collection
#' and for sampling background data in presence-only modelling
#'
#' @param bioDat A data frame with SpatialLong given as X and SpatialLat given as Y, filtered for target group species.
#' @param sppGrp Name of species group (Just for plot labels.
#' @param region_sf An sf object defining the spatial region.
#' @param grain Area of each hexagonal grid cell.
#' @param hexGrid A pre-made hexagonal grid, as produced by \code{makeHexGrid}.
#' @param showPlot Produce a plot.
#'
#' @import plyr
#' @import sf
#' @import ggplot2
#'
#' @return A tibble data frame holding data with target group weights.
#'
#' @export
target.group.weights <- function(bioDat, sppGrp, region_sf, grain, hexGrid = NULL, showPlot = TRUE) {

  # Create hexagonal grid with id attribute
  if (is.null(hexGrid)) hexGrid <- makeHexGrid(region_sf, hexCellArea = grain)

  # Grid species data
  sppHexGrid <- recToHexSetUp(bioDat, 1, region_sf, grain, hexGrid)

  # Target group layer
  targetGrp <- ddply(sppHexGrid, .(id), summarise, N = length(id))
  targetGrp$wt <-  round(targetGrp$N / max(targetGrp$N), 4)
  targetGrp <- merge(sppHexGrid, targetGrp, by = 'id')

  if (showPlot) {

    wt_plot <- ggplot() +
        geom_sf(aes_string(fill = "wt"), data = targetGrp) +
        scico::scale_fill_scico(name = "Weight") +
        theme_minimal(base_size = 14) %+replace% theme(axis.title = element_blank()) +
        ggsn::north(region_sf, "topleft", 0.3, 16) +
        ggtitle(sppGrp)
    print(wt_plot)

  }

  # Drop geometry from output
  targetGrp <- st_drop_geometry(targetGrp)

}

#' Mapping inventory completeness
#'
#' This function calculates inventory completeness using species accumulation curves
#'
#' @param bioDat A data frame with SpatialLong given as X and SpatialLat given as Y, filtered for target group species.
#' @param sppGrp Name of species group (Just for plot labels.
#' @param start_yr A numeric year defining the start of the focal period.
#' @param end_yr A numeric year defining the end of the focal period.
#' @param timeSlice A numeric value specifying the length on each timeslice (in years).
#' @param region_sf An sf object defining the spatial region.
#' @param grain Area of each hexagonal grid cell (in km2).
#' @param hexGrid A pre-made hexagonal grid, as produced by \code{makeHexGrid}.
#'
#' @import vegan
#' @import gtools
#'
#' @export
inventory.completeness <- function(bioDat, sppGrp, start_yr = 2000, end_yr = 2019, timeSlice, region_sf, grain, hexGrid = NULL) {

  #~# Select all years
  siteXSppYr <- siteXSpeciesXTime(bioDat, timeSlice, region_sf, grain, NULL, hexGrid = NULL)
  siteXSppYr$site <- sapply(siteXSppYr$siteTime, function(x) strsplit(x, split = "_")[[1]][1])
  siteXSppYr$year <- sapply(siteXSppYr$siteTime, function(x) strsplit(x, split = "_")[[1]][2])

  #~# Results list
  nyrs <- end_yr - start_yr - 1
  resOut <- vector(mode = "list", length = nyrs)

  #~# Species coverage estimation yr <- 2005
  for(yr in seq((start_yr + timeSlice), end_yr, timeSlice)) {

    #~# Sample increasing amounts of time
    ts_dat <- siteXSppYr[siteXSppYr$year <= yr, ]

    #~# Calculate site based sampling coverage i = 13
    covRes <- data.frame(site = mixedsort(unique(siteXSppYr$site)),
                         year = yr, n_samp = NA, obs_spp = NA,
                         chao = NA, chao.se = NA, jack1 = NA, jack1.se = NA,
                         jack2 = NA, boot = NA, boot.se = NA)
    covRes$site <- as.numeric(as.character(covRes$site))

    #~# Run through all sites
    for(i in 1:nrow(covRes)) {

      # Subset to time period
      ts_dat_tmp <- ts_dat[ts_dat$site == covRes[i,1], ]
      ts_dat_tmp <- ts_dat_tmp[, !(names(ts_dat_tmp) %in% c('siteTime','site','year'))]

      # N species
      covRes[i,4] <- length(which(colSums(ts_dat_tmp) > 0))
      covRes[i,3] <- nrow(ts_dat_tmp)

      if( covRes[i,3] >= 2 & covRes[i,4] > 1) {

        specRes <- specpool(ts_dat_tmp, smallsample = T)

        # Output
        covRes[i,3] <- specRes$n
        covRes[i,4] <- specRes$Species
        covRes[i,5] <- specRes$chao
        covRes[i,6] <- specRes$chao.se
        covRes[i,7] <- specRes$jack1
        covRes[i,8] <- specRes$jack1.se
        covRes[i,9] <- specRes$jack2
        covRes[i,10] <- specRes$boot
        covRes[i,11] <- specRes$boot.se

      }

    }

    resOut[[ (yr - start_yr) ]] <- covRes

  }

  resOut <- do.call(rbind, resOut)

}


