#' dist_calc – R function for estimating the minimum distance between starting point and species

#' Description
#' An R function that estimates the minimum distance between a specified starting point and the species under study based on the distVincentyEllipsoid function of `geosphere´ package (Hijmans, 2017).
#' There is an option to run the function in parallel. 
#' `geosphere´, `sp´, `foreach´, `doParallel´ packages should have been installed before using the functions (Bivand, Pebesma, & Gómez-Rubio, 2013; Hijmans, 2017; Microsoft corporation & Weston, 2014; Microsoft corporation & Weston, 2020; Pebesma & Bivand, 2005).
#' Save the file minA_functions.R in your working directory. 
#' From R menu select file - > Source R code… and select the file minA_functions.R
#' Usage
#' dist_calc (poly, pnt, parallel = FALSE, cores)

#' @param poly	A SpatialPolygonsDataFrame object of the species distribution following EPSG:4326 geodetic coordinate system with two columns in data slot. The first column is named id and refers to the taxa unique id (integer or character). The second column is named polyid, a unique id that refers to the polygon feature (integer). 
#' @param pnt	A SpatialPointsDataFrame object of the starting points following EPSG:4326 geodetic coordinate system with three columns. The first is named pntid and refers to the starting point’s unique id (integer), the other two are called plon, plat and are the starting point’s longitude and latitude.
#' @param parallel	Logical. Should the function run in parallel?
#' @param cores	Integer. The number of cores to use when parallel = TRUE. If missing, the cores are automatically set equal to the available cores.     


#' @return A dataframe with ten columns: distance is the minimum distance between the starting point and the corresponding species closest occurrence measured in km, lon and lat are the coordinates of the closest occurrence of that species, inp is binary for species already present at the starting point, polyid is the unique id of the species polygon feature, id is the unique identifier for the species, plon and plat are the coordinates of the starting point, pntid is the unique id of the starting point, spnum is the species richness found in the corresponding distance from the starting point sorted according to distance and per starting point.
#' @export

dist_calc <- function(poly, pnt, parallel = FALSE, cores) { 
	requireNamespace(sf)
	requireNamespace(geosphere)
	requireNamespace(sp)
	requireNamespace(plyr)
	requireNamespace(foreach)
	
  
	`%fun%` <- `%do%`
	s <- unique(poly$id) # get species' unique ids
	
	# check if the function should run in parallel
	if (parallel == TRUE) { 
		requireNamespace(doParallel)
		if (missing(cores)) { cores <- detectCores() }
		cl <- makePSOCKcluster(cores)
		registerDoParallel(cl)
		`%fun%` <- `%dopar%`
	}
	 
res <- foreach (i = 1:nrow(pnt), .packages = c("geosphere", "sp", "sf"), .combine = rbind, .errorhandling = 'remove', .inorder = FALSE) %fun% { # repeat this process for each starting point
		
		length(poly$polyid)
		ce <- pnt[i, ]
		inptab <- over(ce,poly) # identify which species are present at the starting point
		l <- unique(inptab[!is.na(inptab$id), 1]) # count the species present at the starting point
		if (length(l) == length(s)) { stop("All species present at starting point") }
		
		if (length(l) == 0) {
			nnot <- 0
			indx <- s
			tmp <- data.frame(matrix(ncol = 10, nrow = nnot))
		}
		else {
			not <- l # species present at centroid
			nnot <- length(l)
			indx <- setdiff(s, not) # get the rest of the species
			tmp <- data.frame(matrix(ncol = 10, nrow = nnot)) #' create dataframe to store temporary values
			tmp[ ,1:3] <- 0 # set coordinates of closest occurrence and distance to zero
			tmp[ ,4] <- 1 # set inp to 1
			tmp[ ,5] <- unique(inptab$polyid) # fill in with polygons ids
			tmp[ ,6] <- not # fill in with species ids
			}
		
	inp <- data.frame(matrix(ncol = 10, nrow = length(s) - nnot)) # create dataframe to store temporary values

# calculate the shortest distance between each species distribution polygon and the starting point and finds the closest polygon
		dists <- sapply(1:length(indx), function(b) tryCatch(dist2Line(ce, poly[poly@data$id == indx[b], ], distVincentyEllipsoid), error = function(e) c(NA, NA, NA, NA)), simplify = FALSE) 

		dists <- as.data.frame(matrix(ncol = 4, byrow = TRUE, unlist(dists))) 
		polyid <- sapply(1:length(indx), function(b) unique(poly@data[poly@data$id == indx[b], 2]))
		
		# fill the temporary dataframe
		inp[ ,1] <- dists[ ,1]/1000
		inp[ ,2] <- dists[ ,2]
		inp[ ,3] <- dists[ ,3]
		inp[ ,4] <- 0
		pntid <- dists[ ,4]
		inp[ ,5] <- sapply(1:length(pntid), function(i) polyid[[i]][pntid[i]], simplify = TRUE)
		inp[ ,6] <- indx
		inp <- rbind(inp, tmp)
		inp[ ,7] <- ce@data$plon
		inp[ ,8] <- ce@data$plat
		inp[ ,9] <- ce@data$pntid
		inp[ ,10] <- NA
		
		res <- inp # attach to the final dataframe
		
	}
	
		colnames(res) <- c("distance", "lon", "lat", "inp", "polyid", "id", "plon", "plat", "pntid", "spnum")
		res <- with(res, res[order(pntid,distance), ]) #' sort species according to distance from starting point in ascending order
		res$spnum <- ave(res$distance, res$pntid, FUN = seq_along) #' count species richness with every point added
	if (parallel == TRUE) { 
stopCluster(cl) 
}
	return(res)
}
