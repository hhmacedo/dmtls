#' Graduate grouped data
#'
#' @description
#' A wrapper function for \code{"sprague"} and \code{"uniform"} graduation
#' methods
#'
#' @param Value numeric vector, presumably counts in grouped ages
#' @param Age integer vector, lower bounds of age groups
#' @param AgeInt integer vector, age interval widths
#' @param OAG logical, default = \code{TRUE} is the final age group open?
#' @param OAnew integer, optional new open age, higher than \code{max(Age)}. See
#'   details.
#' @param method character, either \code{"sprague"} or \code{"uniform")}
#' @param keep0 logical. Default \code{FALSE}. If available, should the value in
#'   the infant age group be maintained, and ages 1-4 constrained?
#' @param constrain logical. Default \code{FALSE}. Should output be constrained
#'   to sum within the input age groups?
#' @param ... extra arguments
#'
#' @seealso \code{\link{graduate_sprague}}, \code{\link{graduate_uniform}}
#'
#' @export

graduate <-
  function(Value,
           Age,
           AgeInt = age2int(Age),
           OAG = TRUE,
           OAnew = max(Age),
           method = c("sprague", "uniform"),
           keep0 = FALSE,
           constrain = FALSE,
           ...) {

    method <- tolower(method)

    # validate method choice
    method <- match.arg(method)

    # handle infant age for sprague and beers, if required.
    if (keep0) {
      if (Age[1] == 0 & AgeInt[1] == 1) {
        V0 <- Value[1]
        V5 <- sum(Value[Age < 5])
        V4 <- V5 - V0
        a1 <- names2age(out)
        ind <- a1 < 5 & a1 > 0
        out[ind] <- rescale_vector(out[ind], scale = V4)
        out[1] <- V0
      }
    }

    n  <- length(out)
    a1 <- min(Age):(min(Age) + n - 1)

    # TR: those methods that require AgeInt may depend on final value not being
    # NA even if it's an open age, in essence, keep all value inside this same
    # "single" age. NA coding isn't the best choice here, but we anticipate this
    # and assign 1 to AgeInt[length(AgeInt)] IFF OAG & max(Age) == OAnew
    # This is inconsequential for those methods that don't use AgeInt
    N <- length(AgeInt)
    if (OAG & is.na(AgeInt[N])) {
      nlast <- OAnew - max(Age) + 1
      AgeInt[N] <- nlast
    }

    # Sprague in strict 5-year age groups
    if (method == "sprague") {
      out <- graduate_sprague(Value, Age = Age, OAG = OAG)
    }

    # Uniform respects irregular intervals
    if (method == "uniform") {
      OAvalue <- OAnew - max(Age) + 1
      out <- graduate_uniform(
        Value = Value,
        Age = Age,
        AgeInt = AgeInt,
        OAG = OAG,
        OAvalue = OAvalue
      ) # OAvalue only if OAG and
      # extrapolation desired?
    }

    # option to contrain to sum to original age groups
    if (constrain) {
      out <- rescaleAgeGroups(Value1 = out,
                              AgeInt1 = rep(1, n),
                              Value2 = Value,
                              AgeInt2 = AgeInt,
                              splitfun = graduate_uniform,
                              recursive = FALSE)
      out[is.nan(out)] <- 0
    }

    # last min names assure
    names(out) <- a1

    out

  }

#' Convert arbitrary age groupings into single years of age.
#'
#' @description
#' Uniformly splits aggregate counts in age groups into single year age groups.
#'
#' @details
#' Assumes that the population is uniformly distributed across each age
#' interval, and that initial age intervals are integers greater than or equal
#' to 1. If \code{AgeInt} is given, its final value is used as the interval for
#' the final age group. If \code{AgeInt} is missing, then \code{Age} must be
#' given, and the open age group is by default preserved \code{OAvalue} rather
#' than split. To instead split the final age group into, e.g., a 5-year age
#' class, either give \code{AgeInt}, *or* give \code{Age}, \code{OAG = TRUE},
#' and \code{OAvalue = 5}.  `Age` be any age range, it does not need to start at
#' 0.
#'
#' @inheritParams graduate
#' @param OAvalue Desired width of open age group. See details.
#'
#' @return Numeric vector of counts for single year age groups.
#'
#' @export
#'
#' @examples
#' MalePop <- c(9544406,7471790,11590109,11881844,11872503,12968350,
#' 		11993151,10033918,14312222,8111523,15311047,6861510,13305117,7454575,
#' 		9015381,10325432,9055588,5519173)
#' Ages <- seq(0, 85, by = 5)
#' graduate_uniform(MalePop, Age = Ages)

graduate_uniform <-
  function(Value,
           Age,
           AgeInt,
           OAG = TRUE,
           OAvalue = 1) {

    if (missing(Age) & missing(AgeInt)) {
      Age <- names2age(Value)
    }
    if (missing(AgeInt)) {
      # give 1 to final interval to preserve
      AgeInt <- age2int(Age, OAG = OAG, OAvalue = OAvalue)
    }
    if (missing(Age)) {
      Age <- cumsum(AgeInt) - AgeInt
    }
    # discount for single
    out <- rep(Value / AgeInt, times = AgeInt)
    names(out) <- min(Age):(min(Age) + length(out) - 1)
    out
  }

#' The basic Sprague age-splitting method.
#'
#' @description
#' This method is used to interpolate counts based on the Sprague formula. It is
#' based on the first stage of the Sprague R script prepared by Thomas Buettner
#' and Patrick Gerland, itself based on the description in Siegel and Swanson,
#' 2004, p. 727.
#'
#' @details Ages should refer to lower age bounds, ending in the open age group
#' in the last row (not a closed terminal age). Dimension labeling is necessary.
#' There must be at least six age groups (including the open group). One year of
#' data will work as well, as long as it's given as or coercible to a
#' single-column matrix. This method may produce negative values, most likely in
#' the youngest or oldest ages. This case is dealt with in the \code{graduate()}
#' wrapper function but not in this function.
#'
#' If the highest age does not end in a 0 or 5, and \code{OAG == TRUE}, then the
#' open age will be grouped down to the next highest age ending in 0 or 5. If
#' the highest age does not end in a 0 or 5, and \code{OAG == FALSE}, then
#' results extend to single ages covering the entire 5-year age group.
#'
#' @inheritParams graduate
#'
#' @return Numeric vector of counts split into single ages.
#'
#' @export

graduate_sprague <- function(Value,
                             Age,
                             AgeInt,
                             OAG = TRUE) {

  if (missing(Age) & missing(AgeInt)) {
    Age <- names2age(Value)
  }
  if (missing(AgeInt)) {
    # give 1 to final interval to preserve
    AgeInt <- age2int(Age, OAG = OAG, OAvalue = 1)
  }
  if (missing(Age)) {
    Age <- int2age(AgeInt)
  }

  punif1 <- graduate_uniform(
    Value = Value,
    AgeInt = AgeInt,
    Age = Age,
    OAG = OAG)
  # this is innocuous if ages are already grouped
  a1 <- as.integer(names(punif1))
  pop5 <- groupAges(
    punif1,
    Age = a1,
    N = 5,
    shiftdown = 0
  )
  # depending on OAG, highest age may shift down.
  a5 <- as.integer(names(pop5))

  # generate coefficient matrix
  scm <- graduate_sprague_expand(
    Value = pop5,
    Age = a5,
    OAG = OAG)

  # redistribute
  pop1 <- scm %*% pop5

  dim(pop1) <- NULL
  # label and return
  names(pop1) <- min(Age):(length(pop1) - 1)

  pop1
}

#' Create the Sprague coefficient matrix.
#'
#' @description
#' The resulting coefficient matrix is based on the number of rows in
#' \code{popmat} where is assumed that each row of data is a 5-year age group.
#' The final row may be an open or closed age group, as indicated by the
#' \code{OAG} argument.
#'
#' @details
#' The \code{popmat} matrix is really just a placeholder in this case. This
#' function is a utility called by the Sprague family of functions, where it is
#' most convenient to just pass in the same matrix being used in those
#' calculations to determine the layout of the coefficient matrix.
#'
#' @inheritParams graduate

graduate_sprague_expand <-
  function(Value, Age, OAG = TRUE) {
  popmat <- as.matrix(Value)

  # figure out ages and years
  Age5 <- Age
  Age1 <- min(Age5):max(Age5)

  # nr 5-year age groups
  m <- length(Value)
  # nr rows in coef mat.
  n <- m * 5 - ifelse(OAG, 4, 0)
  # number of middle blocks
  MP <- m - ifelse(OAG, 5, 4)

  # get the split coefficients

  # block for ages 0-9
  # TR: 5-5-2021, this assumes ages start at 0...
  g1g2 <-
    matrix(
      c(
         0.3616, -0.2768,  0.1488, -0.0336,  0.0000,
         0.2640, -0.0960,  0.0400, -0.0080,  0.0000,
         0.1840,  0.0400, -0.0320,  0.0080,  0.0000,
         0.1200,  0.1360, -0.0720,  0.0160,  0.0000,
         0.0704,  0.1968, -0.0848,  0.0176,  0.0000,
         0.0336,  0.2272, -0.0752,  0.0144,  0.0000,
         0.0080,  0.2320, -0.0480,  0.0080,  0.0000,
        -0.0080,  0.2160, -0.0080,  0.0000,  0.0000,
        -0.0160,  0.1840,  0.0400, -0.0080,  0.0000,
        -0.0176,  0.1408,  0.0912, -0.0144,  0.0000
      ),
      nrow = 10,
      ncol = 5,
      byrow = TRUE
    )

  # block for middle ages
  g3 <-
    matrix(
      c(
        -0.0128,  0.0848,  0.1504, -0.0240,  0.0016,
        -0.0016,  0.0144,  0.2224, -0.0416,  0.0064,
         0.0064, -0.0336,  0.2544, -0.0336,  0.0064,
         0.0064, -0.0416,  0.2224,  0.0144, -0.0016,
         0.0016, -0.0240,  0.1504,  0.0848, -0.0128
      ),
      nrow = 5,
      ncol = 5,
      byrow = TRUE
    )

  # block prior to closeout
  g4g5 <-
    matrix(
      c(
        0.0000, -0.0144,  0.0912,  0.1408, -0.0176,
        0.0000, -0.0080,  0.0400,  0.1840, -0.0160,
        0.0000,  0.0000, -0.0080,  0.2160, -0.0080,
        0.0000,  0.0080, -0.0480,  0.2320,  0.0080,
        0.0000,  0.0144, -0.0752,  0.2272,  0.0336,
        0.0000,  0.0176, -0.0848,  0.1968,  0.0704,
        0.0000,  0.0160, -0.0720,  0.1360,  0.1200,
        0.0000,  0.0080, -0.0320,  0.0400,  0.1840,
        0.0000, -0.0080,  0.0400, -0.0960,  0.2640,
        0.0000, -0.0336,  0.1488, -0.2768,  0.3616
      ),
      nrow = 10,
      ncol = 5,
      byrow = TRUE
    )

  ## create a Sprague coefficient matrix for 5-year age groups
  bm <- matrix(0, nrow = n, ncol =  m)
  ## insert upper left block
  bm[1:10, 1:5] <- g1g2

  # determine positions of middle blocks
  rowpos <- matrix(11:((MP * 5) + 10), ncol = 5, byrow = TRUE)
  colpos <- row(rowpos) + col(rowpos) - 1
  for (i in (1:MP)) {
    # calculate the slices and add middle panels accordingly
    bm[rowpos[i,], colpos[i,]] <- g3
  }

  ## insert last two panels

  fr <- nrow(bm) - ifelse(OAG, 10, 9)
  lr <- fr + 9
  fc <- ncol(bm) - ifelse(OAG, 5, 4)
  lc <- fc + 4
  bm[fr:lr, fc:lc] <- g4g5

  if (OAG) {
    # preserve open ended age group
    bm[nrow(bm), ncol(bm)] <- 1
  }

  bm
}
