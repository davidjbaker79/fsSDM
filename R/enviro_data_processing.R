# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' Process spatial data from OS and other sources
#'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Distance decay to environmental features
#'
#' @param enviro_sf An environmental sf object
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea A numeric value specifying the area of the hexagonal cell (in km2).
#' @param hexGrid A pre-made hexagonal grid, as producted by \code{makeHexGrid}.
#' @param var_name Name of environmental variable from \code{enviro_sf}.
#'
#' @import sf
#'
#' @return An sf object.
#' @export
distance.decay <- function(enviro_sf, region_sf, hexCellArea, hexGrid = NULL, var_name) {

  # ~# Create hexagonal grid with id attribute
  if (is.null(hexGrid)) hexGrid <- makeHexGrid(region_sf, hexCellArea)

  # ~# Clip to region
  enviro_clip <- st_crop(enviro_sf, region_sf)

  # ~# Calculate distance from each grid cell to poly / point
  enviro_dist <- st_distance(hexGrid, enviro_clip)
  enviro_hex_dist <- cbind(hexGrid, value = apply(enviro_dist, 1, min) / 1000)
  enviro_hex_dist$grain <- hexCellArea

  # ~# Drop geometry and return
  enviro_hex_dist <- st_drop_geometry(enviro_hex_dist)
  enviro_hex_dist <- enviro_hex_dist[,c("id", "grain", "value")]
  enviro_hex_dist$var_name <- var_name

  return(enviro_hex_dist)
}

#' Count feature occurrence within cells
#'
#' @param enviro_sf An environmental sf object
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea A numeric value specifying the area of the hexagonal cell (in km2).
#' @param hexGrid A pre-made hexagonal grid, as producted by \code{makeHexGrid}.
#' @param var_name Name of environmental variable from \code{enviro_sf}.
#'
#' @import sf
#'
#' @return An sf object.
#'
#' @export
feature.occurence <- function(enviro_sf, region_sf, hexCellArea, hexGrid = NULL, var_name) {

  # ~# Create hexagonal grid with id attribute
  if (is.null(hexGrid)) hexGrid <- makeHexGrid(region_sf, hexCellArea)

  # ~# Count features within cells
  featureCount <- st_join(hexGrid, enviro_sf, join = st_intersects)
  featureCount <- ddply(featureCount, c("id", "grain"), summarise, value = as.numeric(length()))
  featureCount$grain <- hexCellArea

  # ~# Prep and return
  featureCount <- st_drop_geometry(featureCount)
  featureCount <- featureCount[, c("id", "grain", "value")]
  featureCount$var_name <- var_name

  return(featureCount)
}

#' Count feature occurrence within cells
#'
#' @param enviro_sf An environmental sf object
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea A numeric value specifying the area of the hexagonal cell (in km2).
#' @param hexGrid A pre-made hexagonal grid, as producted by \code{makeHexGrid}.
#' @param var_name Name of environmental variable from \code{enviro_sf}.
#'
#' @import  sf
#'
#' @return An sf object.
#'
#' @export
feature.length <- function(enviro_sf, region_sf, hexCellArea, hexGrid = NULL, var_name) {

  # ~# Create hexagonal grid with id attribute
  if (is.null(hexGrid)) hexGrid <- makeHexGrid(region_sf, hexCellArea)

  # ~# Count features within cells
  fl <- st_join(enviro_sf, hexGrid, join = st_intersects)

  # ~# Calculate road/path lengths
  fl.res <- data.frame(id = na.omit(unique(hexGrid$id)), value = NA)
  for (i in fl.res$id) {
    flg <- fl["id" == i,]
    lKm <- as.numeric(sum(st_length(st_cast(flg, "LINESTRING"))) / 1000)
    fl.res[which(fl.res[, 1] == i), 2] <- lKm
  }
  fl.res$grain <- hexCellArea
  fl.res <- fl.res[,c("id", "grain", "value")]
  fl.res$var_name <- var_name

  # # Test plot
  # demo_plot <- left_join(hexGrid, fl.res, by='id')
  # print(ggplot() + geom_sf(aes(fill=value),data=demo_plot) + scale_fill_distiller(palette = 'Greys', direction = 1))

  return(fl.res)
}

#' Calculate feature cover within cells
#'
#' @param enviro_sf An environmental sf object
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea A numeric value specifying the area of the hexagonal cell (in km2).
#' @param hexGrid A pre-made hexagonal grid, as producted by \code{makeHexGrid}.
#' @param var_name Name of environmental variable from \code{enviro_sf}.
#'
#' @import sf
#'
#' @return An sf object.
#'
#' @export
feature.cover <- function(enviro_sf, region_sf, hexCellArea, hexGrid = NULL, var_name) {

  # ~# Create hexagonal grid with id attribute
  if (is.null(hexGrid)) hexGrid <- makeHexGrid(region_sf, hexCellArea)

  # ~# Count features within cells
  fa <- st_intersection(hexGrid, enviro_sf)
  fa.out <- data.frame(id = unique(hexGrid$id), grain = hexCellArea, value = NA)

  # ~# run through cells
  for (i in unique(hexGrid$id)) {
    print(i)
    fai <- fa["id" == i, ]
    if (nrow(fai) != 0) {
      fai <- st_cast(fai)
      fai <- st_area(fai)
      fa.out[which(fa.out[, 1] == i), 3] <- round(as.numeric(fai / 10^6), 3)
    }
  }
  fa.out$var_name <- var_name

  # ~# Test plot
  # demo_plot <- left_join(hexGrid, fa.out, by='id')
  # print(ggplot() + geom_sf(aes(fill=value),data=demo_plot) +
  #         geom_sf(aes(),data=enviro_sf) +
  #         scale_fill_distiller(palette = 'Greys', direction = 1,na.value='Yellow'))

  return(fa.out)
}

#' Function for extracting CEH landcover classess to grid cells
#'
#' LC_1 = Deciduous woodland (wd_dec)
#' LC_2 = Coniferous woodland (wd_con)
#' LC_3 = Arable (arab)
#' LC_4 = Improved grassland (grImp)
#' LC_5 = Neutral grassland (grNeu)
#' LC_6 = Calcareous grassland (grCal)
#' LC_7 = Acid grassland (grAci)
#' LC_8 = Fen (fen)
#' LC_9 = Heather (heath)
#' LC_10 = Heather grassland (grHeath)
#' LC_11 = Bog (bog)
#' LC_12 = Inland rock (rkInld)
#' LC_13 = Saltwater (saltwat)
#' LC_14 = Freshwater (freshwat)
#' LC_15 = Supralittoral rock (rkSup)
#' LC_16 = Supralittoral sediment (sediSup)
#' LC_17 = Littoral rock (rkLit)
#' LC_18 = Littoral sediment (sediLit)
#' LC_19 = Saltmarsh (saltmar)
#' LC_20 = Urban (urban)
#' LC_21 = Suburban (surban)
#'
#' @param ceh_dat A CEH Land cover dataset (e.g. 2019) as an sf object.
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea A numeric value specifying the area of the hexagonal cell (in km2).
#' @param hexGrid A pre-made hexagonal grid, as produced by \code{makeHexGrid}.
#' @param buffer The buffer distance in kilometers.
#'
#' @import sf
#'
#' @return A data frame with the area of cover from each land class category (LC_1 to LC_21).
#'
#' @export
ceh.landcov.extract <- function(ceh_dat, region_sf, hexCellArea, hexGrid = NULL, buffer = 1) {

  # Create hexagonal grid with id attribute
  if (is.null(hexGrid)) hexGrid <- makeHexGrid(region_sf, hexCellArea)

  # Transform to sort coord issues with CEH
  hexGridProj <- st_transform(hexGrid, crs = raster::crs(ceh_dat))

  # Extract cover values for each land cover class
  res_out <- data.frame(matrix(NA, length(unique(hexGrid$id)), 22))
  names(res_out) <- c("id", paste0("LC_", 1:21, "_", buffer))
  res_out$id <- unique(hexGrid$id)
  res_out$area <- hexGrid$area
  res_out$grain <- hexCellArea

  # Loop through cell
  for (i in 1:nrow(res_out)) {

    hexGrainCell <- hexGridProj[id == res_out$id[i],]

    if (buffer > 0) {

      hexGrainCell <- st_buffer(hexGrainCell, dist = (buffer * 1000))

    }

    lc_cell <- raster::mask(raster::crop(ceh_dat, hexGrainCell), hexGrainCell)

    # plot(lc_cell)

    # Loop through habitat and calculate area - j <- 15; plot(lc_cell); plot(hexGrainCell, add = TRUE)
    for (j in unique(na.omit(getValues(lc_cell)))) {

      lc_tmp <- lc_cell
      lc_tmp[lc_tmp != j] <- NA
      lc_tmp <- st_as_sf(rasterToPolygons(lc_tmp, dissolve = TRUE))
      lc_tmp <-  st_intersection(hexGridProj, lc_tmp)

      # Extract and return
      try(res_out[which(res_out$id == i), (j + 1)] <- as.numeric(round(rgeos::gArea(sf::as(lc_tmp, "Spatial")) / 10^6, 3)))

    }
  }

  # Return
  return(res_out)

}

#' Bioclimate data is from 1983 to 2017 - need to aggregate
#'
#' add a period argument to select years
#'
#' @param dir_name The path to the climate files
#'
#' @import raster
#'
#' @return A list of data frames
#'
#' @export
extract.bioclim <- function(dir_name) {

  fyrs <- list.files(dir_name, full.names = TRUE)
  r_tmp <- lapply(fyrs, function(x) tmp <- raster(x))
  r_tmp <- do.call(stack, r_tmp)
  r_m <- stackApply(r_tmp, 1, mean)
  r_sd <- stackApply(r_tmp, 1, sd)

  return(list(r_m, r_sd))

}

#' Calculate the mean and standard deviation of bioclimate variables across each cell
#'
#' @param clim_grd A grid of 100 m x 100 m micro-bioclimate variables
#' @param region_sf An sf object defining the spatial region.
#' @param hexCellArea A numeric value specifying the area of the hexagonal cell (in km2).
#' @param hexGrid A pre-made hexagonal grid, as produced by \code{makeHexGrid}.
#' @param var_name Name of climate variable.
#' @param buffer_dist Extract values from buffer if presence in focal cell is likely to be
#' affected by availability of surrounding climate
#'
#' @import sf
#'
#' @return A tibble data frame.
#'
#' @export
climate.hex.variation <- function(clim_grd, region_sf, hexCellArea, hexGrid = NULL, var_name, buffer_dist) {

  # Create hexagonal grid with id attribute
  if (is.null(hexGrid)) hexGrid <- makeHexGrid(region_sf, hexCellArea)

  # Transform to sort coordinates issues with CEH
  hexGridProj <- st_transform(hexGrid, crs = crs(clim_grd))

  # Extract values and calculate mean se within cells
  tmp_m <- st_as_sf(raster::extract(clim_grd, hexGridProj, fun = mean, buffer = buffer_dist, na.rm = TRUE, sp = TRUE))
  tmp_sd <- st_as_sf(raster::extract(clim_grd, hexGridProj, fun = sd, na.rm = TRUE, sp = TRUE))

  # # Test plot
  # print(ggplot() +
  #   geom_sf(aes(fill = index_1), data = tmp_m) +
  #   scale_fill_scico())
  # print(ggplot() +
  #   geom_sf(aes(fill = index_1), data = tmp_sd) +
  #   scale_fill_scico())

  # Name variables
  tmp_m <- tmp_m["index_1" == paste0(var_name, "_m"), ]
  tmp_sd <- tmp_m["index_1" == paste0(var_name, "_sd"), ]

  # Drop geometry and join mean and sd
  tmp_m <- st_drop_geometry(tmp_m)
  tmp_sd <- st_drop_geometry(tmp_sd)
  tmp <- merge(tmp_m, tmp_sd, by = c("id", "grain", "label", "name", "code", "area"), all.x = TRUE)
  tmp <- tmp[, !(colnames(tmp) %in% c("label", "name", "code", "area"))]

  return(tmp)

}
