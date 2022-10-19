help_count <- function(img, watershed = TRUE, filter = FALSE, invert = FALSE, fill_hull = FALSE,
                       object_size = "medium", index = "NB", object_index = NULL, 
                       threshold = "Otsu", tolerance = NULL, extension = NULL, lower_size = NULL, 
                       upper_size = NULL, topn_lower = NULL, topn_upper = NULL, 
                       lower_eccent = NULL, upper_eccent = NULL, lower_circ = NULL, 
                       upper_circ = NULL) {
  if (filter != FALSE) {
    if (!is.numeric(filter)) {
      stop("Argument `filter` must be numeric.", call. = FALSE)
    }
    img <- image_filter(img, size = filter)
  }
  
  img2 <- image_binary(img, index = index, 
                       invert = invert, fill_hull = fill_hull, threshold = threshold, 
                       resize = FALSE, show_image = FALSE)[[1]]
  if (isTRUE(watershed)) {
    parms <- read.csv(file = system.file("parameters.csv", 
                                         package = "pliman", mustWork = TRUE), header = T, 
                      sep = ";")
    res <- length(img2)
    parms2 <- parms[parms$object_size == object_size, 
    ]
    rowid <- which(sapply(as.character(parms2$resolution), 
                          function(x) {
                            eval(parse(text = x))
                          }))
    ext <- ifelse(is.null(extension), parms2[rowid, 
                                             3], extension)
    tol <- ifelse(is.null(tolerance), parms2[rowid, 
                                             4], tolerance)
    nmask <- EBImage::watershed(EBImage::distmap(img2), 
                                tolerance = tol, ext = ext)
  }
  else {
    nmask <- EBImage::bwlabel(img2)
  }
  ID <- which(img2 == 1)
  ID2 <- which(img2 == 0)
  
  shape <- cbind(data.frame(EBImage::computeFeatures.shape(nmask)), 
                 data.frame(EBImage::computeFeatures.moment(nmask)))
  object_contour <- EBImage::ocontour(nmask)
  ch <- conv_hull(object_contour)
  area_ch <- trunc(as.numeric(unlist(poly_area(ch))))
  shape <- transform(shape, id = 1:nrow(shape), radius_ratio = s.radius.max/s.radius.min, 
                     diam_mean = s.radius.mean * 2, diam_min = s.radius.min * 
                       2, diam_max = s.radius.max * 2, area_ch = area_ch, 
                     solidity = s.area/area_ch, circularity = 4 * pi * 
                       (s.area/s.perimeter^2), minor_axis = m.majoraxis * 
                       sqrt(1 - m.eccentricity^2))
  shape <- shape[, c("id", "m.cx", "m.cy", "s.area", "area_ch", 
                     "s.perimeter", "s.radius.mean", "s.radius.min", "s.radius.max", 
                     "s.radius.sd", "radius_ratio", "diam_mean", "diam_min", 
                     "diam_max", "m.majoraxis", "minor_axis", "m.eccentricity", 
                     "m.theta", "solidity", "circularity")]
  colnames(shape) <- c("id", "x", "y", "area", "area_ch", 
                       "perimeter", "radius_mean", "radius_min", "radius_max", 
                       "radius_sd", "radius_ratio", "diam_mean", "diam_min", 
                       "diam_max", "major_axis", "minor_axis", "eccentricity", 
                       "theta", "solidity", "circularity")
  if (!is.null(lower_size) & !is.null(topn_lower) | !is.null(upper_size) & 
      !is.null(topn_upper)) {
    stop("Only one of 'lower_*' or 'topn_*' can be used.")
  }
  ifelse(!is.null(lower_size), shape <- shape[shape$area > 
                                                lower_size, ], shape <- shape[shape$area > mean(shape$area) * 
                                                                                0.1, ])
  if (!is.null(upper_size)) {
    shape <- shape[shape$area < upper_size, ]
  }
  if (!is.null(topn_lower)) {
    shape <- shape[order(shape$area), ][1:topn_lower, 
    ]
  }
  if (!is.null(topn_upper)) {
    shape <- shape[order(shape$area, decreasing = TRUE), 
    ][1:topn_upper, ]
  }
  if (!is.null(lower_eccent)) {
    shape <- shape[shape$eccentricity > lower_eccent, 
    ]
  }
  if (!is.null(upper_eccent)) {
    shape <- shape[shape$eccentricity < upper_eccent, 
    ]
  }
  if (!is.null(lower_circ)) {
    shape <- shape[shape$circularity > lower_circ, ]
  }
  if (!is.null(upper_circ)) {
    shape <- shape[shape$circularity < upper_circ, ]
  }
  object_contour <- object_contour[shape$id]
  
  invisible(
    list("obj_ctr" = object_contour,
         "results" = shape)
  )
}












measure_mite_gall <- function(img.path, qc.check.path){
  img <- image_import(img.path)
  img.seg <- image_segment_iter(img, 
                                c("HUE","GR"),
                                nseg = 2,
                                invert = c(FALSE,FALSE),
                                threshold = c("Otsu","Otsu"), 
                                show_image = FALSE,
                                verbose = FALSE)
  img.an <- img.seg$images$seg2 %>% 
    help_count(watershed = TRUE,
               filter = 4,
               lower_size = 300,
               upper_eccent = 0.9,
               object_size = "small")
  
  if(!is.null(qc.check.path)){
    png(qc.check.path)
    plot(img)
    plot_contour(img.an$obj_ctr, col = "red", 
                 lwd = 1)
    dev.off()
  }
  return(sum(img.an$results$area))
}
