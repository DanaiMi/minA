#' area_calc – R function for estimating the minA for each starting point
#'
#' Description
#' An R function that estimates the minA for each specified starting point using the output of dist_calc function.
#' `geosphere´, `sf´, `sp´ packages should have been installed before using the functions.
#'
#' Usage
#' area_calc (cb)
#'
#' @param cb A dataframe as produced by the dist_calc function.
#'
#' @return A list with two elements. The first is a dataframe called data with 5 columns: pntid is the unique id of the starting point, spnum is the species richness found in the corresponding area, area is the estimated minA in km2, plon and plat are the coordinates of the starting point. The second element is a list containing another list for each starting point with the coordinates of the minimum convex polygon’s vertices.

area_calc <- function(cb) {
  # cb is the output of dist_calc
  requireNamespace(sf)
  requireNamespace(geosphere)
  requireNamespace(sp)
  requireNamespace(plyr)

  # set the coordinates for the species present at the starting point equal to the starting point's coordinates
  cb[cb$inp == 1, 2] <- cb[cb$inp == 1, 7]
  cb[cb$inp == 1, 3] <- cb[cb$inp == 1, 8]
  cb <- cb[complete.cases(cb), ]
  # convert starting points to sf points
  cbp <- st_as_sf(cb, coords = c("lon", "lat"))
  st_crs(cbp) <- 4326
  cb$Row.names <- NULL

  pt <- unique(cb$pntid) # the starting points' ids
  spp <- 1:max(cb$spnum) # the number of species
  cb <- cb[order(cb$pntid, cb$spnum), ]


  res <- data.frame(pntid = numeric(), spnum = numeric(), area = double(), plon = double(), plat = double())

  cvlist <- vector("list", length(pt))
  cvlist <- setNames(cvlist, paste("pt", pt,sep = "_"))

  for (p in seq_along(pt)) { # repeat for each starting point
    xy <- cb[cb$pntid == pt[p], ] # the dataframe from dist_calc filtered for the starting point
    if (sum(xy$inp) == 0) {
      xy <- rbind(
        data.frame(distance = 0, lon = xy[1, "plon"], lat = xy[1, "plat"], inp = 1, polyid = NA, id = NA, plon = xy[1, 'plon'],
                   plat = xy[1, "plat"], pntid = pt[p], spnum = 0), xy)
    }
    n1 <- paste("pt", pt[p], sep = "_")
    tmp <- data.frame(pntid = numeric(), spnum = numeric(), area = double(),  plon = double(), plat = double())
    splist <- vector("list", max(spp))
    splist <- setNames(splist, paste("sp", spp, sep = "_"))


    for (sp in seq_along(spp)) { # repeat for each species
      x <- xy[xy$spnum <= sp, ] # filter the dataframe according to the number of species
      tm <- x[ ,c(2,3)]
      tm <- round(tm, 5)
      ind <- duplicated(tm[, 1:2], fromLast = TRUE)
      tm <- tm[!ind, ] # get available unique point coordinates

      if (nrow(tm) <= 2) { # if less than 3 points available, the polygon cannot be formed and area is set to zero
        tm <- data.frame(pntid = pt[p], spnum = sp, area = 0, plon = xy[1, 7], plat = xy[1, 8])
        tmp <- rbind(tmp, tm)
      }
      else {
        hull <- chull(tm) # create minimum convex polygon from the points set
        chulls <- tm[hull, ]
        n2 <- paste("sp", sp, sep = "_")
        splist[[n2]] <- chulls # store the vertices of the minimum convex polygon
        area <- areaPolygon(chulls)/10^6 # calculate the area of the minimum convex polygon
        mkpl <- makePoly(chulls, 1000000, sp = TRUE)
        hull <- SpatialPoints(tm, crs(mkpl))
        tab <- over(hull, mkpl)

        tm <- data.frame(pntid = pt[p], spnum = sp, area = area, plon = xy[1, 7], plat = xy[1, 8])
        tmp <- rbind(tmp, tm)
      }
    }
    cvlist[[n1]] <- splist
    res <- plyr::rbind.fill(res, tmp)
 }
  res <- list(data = res, chulls = cvlist) #' export result
  return(res)
}
