# AGEINT.R ----------------------------------------------------------------

#' Interpolate between two population age distributions.
#'
#' @description
#' The interpolation is done by age using a linear, exponential, or power
#' function. This comes from the PAS spreadsheet called \code{AGEINT}. Be aware
#' that this is not cohort-component interpolation.
#'
#' @param popmat numeric. An age-period matrix of data to interpolate over. Age
#'   in rows, time in columns.
#' @param datesIn vector of dates. The exact dates in each column. See details
#'   for ways to express it.
#' @param datesOut vector of dates. The desired dates to interpolate to. See
#'   details for ways to express it.
#' @param method string. The method to use for the interpolation, either
#'   \code{"linear"}, \code{"exponential"}, or \code{"power"}. Default
#'   \code{"linear"}.
#' @param power numeric power to interpolate by, if \code{method = "power"}.
#'   Default 2.
#' @param extrap logical. In case \code{datesOut} is out of range of `datesIn`,
#'   do extrapolation using slope in extreme pairwise. Default \code{FALSE}.
#' @param negatives logical. In case negative output are accepted, set to
#'   \code{TRUE}. Default \code{FALSE}.
#' @param ... arguments passed to \code{stats::approx}. For example,
#'   \code{rule}, which controls extrapolation behavior.
#'
#' @details
#' The age group structure of the output is the same as that of the input.
#' Ideally, \code{datesOut} should be within the range of \code{datesIn}. If
#' not, the left-side and right-side output are held constant outside the range
#' if \code{rule = 2} is passed in, otherwise \code{NA} is returned (see
#' examples). Dates can be given in three ways 1) a \code{Date} class object, 2)
#' an unambiguous character string in the format \code{"YYYY-MM-DD"}, or 3) as a
#' decimal date consisting in the year plus the fraction of the year passed as
#' of the given date.
#'
#' For methods \code{"exponential"} and \code{"power"}, calculations are carried
#' out using linear interpolation through \code{log(popmat)}, or
#' \code{popmat^(1/power)} respectively, then back-transformed. If the data
#' contain 0s, \code{method = "exponential"} will fail, but \code{"power"} would
#' still generally work.
#'
#' @return numeric matrix (age-period) (or vector if \code{length(datesOut) ==
#'   1} of the interpolated data for the requested dates.
#'
#' @importFrom stats approx

interp <- function(popmat,
                   datesIn,
                   datesOut,
                   method = c("linear", "exponential", "power"),
                   power = 2,
                   extrap = FALSE,
                   negatives = FALSE,
                   ...) {
  # ... args passed to stats::approx . Can give control over extrap assumptions
  # IW: extrap=T for extrapolate following each slope in extreme pairwise.
  # If not is explicit extrap=T, returns NA at those points

  # a basic check
  stopifnot(ncol(popmat) == length(datesIn))

  # no sense documenting this wrapper ...
  .approxwrap <- function(x, y, xout, extrap, ...) {

    # interp
    yout = stats::approx(x = x,
                         y = y,
                         xout = xout,
                         ...)$y

    if (extrap){
      # extrap (each side)
      rg <- range(x)
      xext <- xout < rg[1]
      if(any(xext))
        yout[xext] <- (y[2]-y[1])/(x[2]-x[1])*(xout[xext]-x[1])+y[1]
      xext <- xout > rg[2]
      n <- length(y)
      if(any(xext))
        yout[xext] <- (y[n]-y[n-1])/(x[n]-x[n-1])*(xout[xext]-x[n-1])+y[n-1]
    }

    return(yout)
  }


  # clean method declaration
  # match.arg does partial matching and it's safer:
  # match.arg("lin", c("linear", "exponential", "power"))
  method <- tolower(match.arg(method,
                              choices = c("linear", "exponential", "power")))

  # coerce dates to decimal if necessary
  datesIn <- dec.date(datesIn)
  datesOut <- dec.date(datesOut)

  # carry out transform 1
  if (method == "exponential") {
    if (any(popmat == 0)) {
      stop(
        "popmat contains 0s, so exponential interpolation won't work.\n
							Either handle 0s then do this, or\n
							maybe try method = 'power'."
      )
    }
    popmat <- log(popmat)
  }
  if (method == "power") {
    popmat <- popmat ^ (1 / power)
  }

  int <- apply(popmat,
               1,
               .approxwrap,
               x = datesIn,
               xout = datesOut,
               extrap = extrap,
               ...)
  dims            <- dim(int)
  if (!is.null(dims)) {
    int           <- t(int)
    rownames(int) <- rownames(popmat)
    colnames(int) <- datesOut
  } else {
    names(int)    <- rownames(popmat)
  }

  # transform back
  if (method == "exponential") {
    int           <- exp(int)
  }
  if (method == "power") {
    int           <- int ^ power
  }

  # IW: no negatives when extrapolate. Thinking in pop and lt expressions.
  # Inactive when explicitly are accepted negatives
  if(!negatives & all(!is.na(int)) & any(int < 0)){
    cat(paste0("Negative values have been replaced with 0s.\n",
               "Negatives not accepted in population counts,\n",
               "fertility rates or life table functions.\n","
               You can allow negatives (e.g. interpolating net migration)\n",
               "by specifying negatives = TRUE"))
    int[int < 0] <- 0
  }

  int
}

# utils.R -----------------------------------------------------------------

#' Rescale a vector proportionally to a new sum.
#'
#' @description
#' Function used to rescale a vector to a given value. This is a frequently
#' needed operation.
#'
#' @details
#' For a distribution, use \code{scale = 1}. For percentages, use
#' \code{scale = 100}, etc.
#'
#' @param x numeric vector.
#' @param scale numeric. Value the vector should sum to. Default 1.
#'
#' @return The vector rescaled.

rescale_vector <- function(x, scale = 1) {
  scale * x / sum(x, na.rm = TRUE)
}

#' Convert date to decimal year fraction.
#'
#' @description
#' Convert a character or date class to decimal, taking into account leap years.
#'
#' @details
#' This makes use of the \code{lubridate::decimal_date} to compute the
#' proportion of the year that has passed. If the date is numeric, it is
#' returned as such. If it is \code{"character"}, we try to coerce to date
#' through \code{lubridate::ymd}, ergo, it is best to specify a character string
#' in an unambiguous \code{"YYYY-MM-DD"} format.  If \code{date} is given in a
#' \code{"Date"} class it is dealt with accordingly.
#'
#' @param date Either a \code{Date} class object or an unambiguous character
#'   string in the format \code{"YYYY-MM-DD"}.
#'
#' @return Numeric expression of the date, year plus the fraction of the year
#'   passed as of the date.

dec.date  <- function(date) {

  if (inherits(date, "numeric")) {
    return(date)
  }

  res <- lubridate::decimal_date(lubridate::ymd(date))
  res
}

#' Rescale counts in age groups to match counts in different age groups
#'
#' @description
#' This method rescales a vector of counts in arbitrary (integer) age groups to
#' approximate a vector of counts in a potentially different age grouping.
#' Common use cases will be to scale single ages (whose age pattern we wish to
#' roughly maintain) to sum to abridged or 5-year age groups from another
#' source. The counts to be rescaled could potentially be in any grouping (see
#' example).
#'
#' @details If the final age group is open, define its age interval as 1.
#'
#' Presently the intermediate splitting function assumes that counts inside the
#' age groups of population 1 are uniformly distributed, although this may be
#' relaxed if other methods become available whose behavior matches that of
#' \code{splitUniform()}. \code{splitMono()} will be modified soon to be
#' applicable here.
#'
#' The method is an original contribution. It works by first splitting the
#' counts of \code{Value1} to single ages using the assumptions of
#' \code{splitfun()}, which presently only works for \code{splitUniform()}.
#' \code{Value1} is then rescaled such that were it re-grouped to match the age
#' classes of \code{Value2} they would be identical. If
#' \code{recursive = FALSE}, the single-age rescaled \code{Value1} data are
#' returned regrouped to their original ages. If \code{recursive = TRUE}, the
#' process is repeated until \code{Value1} is rescaled such that it could be
#' split and regrouped to \code{Value2} using the same process a single time
#' with no need for further rescaling. If age groups in \code{Value1} are very
#' irregular, \code{recursive = TRUE} can induce noise (see example). If the age
#' groups of \code{Value1} nest cleanly within the age groups of \code{Value2}
#' then recursion is unnecessary. This is the case, for example, whenever
#' \code{Value1} is in single ages and \code{Value2} is in grouped ages, which
#' is likely the most common usage scenario.
#'
#' @param Value1 numeric vector. A vector of demographic counts for population
#'   1.
#' @param AgeInt1 integer vector. Age interval widths for population 1.
#' @param Value2 numeric vector. A vector of demographic counts for population
#'   2.
#' @param AgeInt2 integer vector. Age interval widths for population 2.
#' @param splitfun function to use for splitting \code{pop1}. Presently on
#'   \code{splitUniform()} works.
#' @param recursive logical. Shall we repeat the split/regroup/rescale process
#'   until stable? See details. Default \code{FALSE}.
#' @param tol numeric. Default 1e-3. The numerical tolerance for the residual.
#'   Used to detect stability if \code{recursive = TRUE}.
#'

rescaleAgeGroups <- function(
  Value1,
  AgeInt1,
  Value2,
  AgeInt2,
  splitfun = splitUniform,
  recursive = FALSE,
  tol = 1e-3) {
  N1 <- length(Value1)
  # ages must cover same span
  stopifnot(sum(AgeInt1) == sum(AgeInt2))

  Age1 <- int2age(AgeInt1)
  Age2 <- int2age(AgeInt2)

  stopifnot(N1 == length(Age1))

  AgeN <- rep(Age2, times = AgeInt2)

  # scale from single.
  ValueS <- splitfun(Value1, AgeInt = AgeInt1)
  # right now splitMono() doesn't have AgeInt, so does not create the right spread.
  # comparison forthcoming.
  AgeS <- names2age(ValueS)

  #splitMono(Value = Value1, Age1, F)
  #plot(splitMono(Value = Value1, Age1, T))

  #
  AgeN2 <- rep(Age2, times = AgeInt2)
  beforeN <- groupAges(ValueS, AgeS, AgeN = AgeN2)

  beforeNint <- rep(beforeN, times = AgeInt2)
  afterNint <- rep(Value2, times = AgeInt2)
  ratio <- afterNint / beforeNint
  SRescale <- ValueS * ratio



  # group back to original, in case these weren't single
  AgeN1 <- rep(Age1, times = AgeInt1)
  out <- groupAges(SRescale, AgeS, AgeN = AgeN1)

  # check for recursion
  newN <- splitfun(out, AgeInt = AgeInt1)
  check <- groupAges(newN, AgeS, AgeN = AgeN2)
  if (max(abs(check-pop2)) < tol | !recursive){
    return(out)
  } else {
    rescaleAgeGroups(
      Value1 = out,
      AgeInt1 = AgeInt1,
      Value2 = Value2,
      AgeInt2 = AgeInt2,
      splitfun = splitfun,
      # res = res,
      recursive = recursive)
  }
}

# utilsAge.R --------------------------------------------------------------

#' Detect ages from names of vector(s)
#'
#' @description
#' Often as a shorthand we pull lower age bounds from the names of a vector.
#' This modulates that action, and allows for giving several vectors to check
#' for names.
#'
#' @details
#' If more than one vector is given, names are pulled from the first available
#' named vector. All vectors must be of the same length. If no names are
#' available on any vector, then `NA`s are returned. This clearly won't work if
#' the names on the vector are of something else.
#'
#' @param ... one or more vectors of any class, which has/have a names
#'   attribute.
#'
#' @return integer vector of ages (presumably).

names2age <- function(...) {
  XL <- list(...)
  # if multiple vectors given must be same lengths
  if (length(XL) > 1) {
    LL <-
      unlist(
        lapply(XL,
               function(x) {
                 length(x)
               }
        )
      )
    stopifnot(all(diff(LL) == 0))
  }
  # which have names
  TF <-
    unlist(
      lapply(
        XL,
        function(x) {
          name_fun <-
            ifelse(is.null(dim(x)),
                   names,
                   function(y) {
                     dimnames(y)[[1]]
                   }
            )
          !is.null(name_fun(x))
        }
      )
    )

  if (any(TF)) {
    # takes names from first available
    x <- XL[[which(TF)[1]]]
    name_fun <-
      ifelse(
        is.null(dim(x)),
        names,
        function(y) {
          dimnames(y)[[1]]
        }
      )
    Age <- as.integer(name_fun(x))
  } else {
    x <- XL[[1]]
    Age <- rep(NA, length(x))
  }
  Age
}

#' Infer age class intervals from lower age bounds.
#'
#' @description
#' Determine age class intervals based on a vector of age class lower bounds.
#'
#' @details
#' If the final age group is open, it is given a value of \code{NA} by default,
#' or else a user-determined value. If the final age group is closed, it is
#' assumed to be equal to the next-lower interval. If the final age interval is
#' known and not equal to the next lowest interval, specify \code{OAG = TRUE}
#' and assign its value to \code{OAvalue}.
#'
#' @param Age integer or numeric. Vector of lower age group bounds.
#' @param OAG logical. Whether or not the final age group is open.
#'   Default \code{TRUE}.
#' @param OAvalue numeric or integer. The value to use for the final age
#'   interval if \code{OAG = TRUE}. Default \code{NA}.
#'
#' @return Age interval vector, of same length as \code{Age}.

age2int <- function(Age, OAG = TRUE, OAvalue = NA) {
  fd <- diff(Age)
  c(fd, ifelse(OAG, OAvalue, fd[length(fd)]))
}

#' Infer lower age bounds from age class intervals.
#'
#' @description
#' Determine lower bounds of age classes based on a vector of age intervals and
#' a starting age.
#'
#' @param AgeInt integer or numeric. Vector of age intervals.
#' @param ageMin integer. The lowest age, default 0.
#'
#' @return Age vector of same length as \code{AgeInt}.

int2age <- function(AgeInt, ageMin = 0) {
  n <- length(AgeInt)
  # if final AgeInt is NA, then assume it's OAG,
  # count as zero for this calc
  if (is.na(AgeInt[n])) {
    AgeInt[n] <- 0
  }
  cumsum(AgeInt) - AgeInt + ageMin
}

#' Calculate which large age group single ages belong to.
#'
#' @description
#' Assign single ages to age groups of equal and arbitrary width, and also
#' optionally shifted.
#'
#' @param Age integer. Vector of single ages (lower bound).
#' @param N integer. Desired width of resulting age groups.
#' @param shiftdown integer. Move the grouping down by one or more single ages. Optional argument.
#'
#' @details
#' If you shift the groupings, then the first age groups may have a negative
#' lower bound (for example of -5). These counts would be discarded for the
#' oscillatory version of Sprague smoothing, for example, but they are preserved
#' in this function. The important thing to know is that if you shift the
#' groups, the first and last groups won't be N years wide. For example if
#' \code{shiftdown} is 1, the first age group is 4-ages wide.
#'
#' @return An integer vector of \code{length(Age)} indicating the age group that
#'   each single age belongs to.

calcAgeN <- function(Age, N = 5, shiftdown = 0) {
  shift <- abs(shiftdown)
  stopifnot(shift < N)
  Ngroups <- (Age + shift) - (Age + shift) %% N
  l <- rle(Ngroups)$lengths
  inds <- cumsum(l) - l + 1
  rep(Age[inds], times = l)
}

#' Group down to a new open age class.
#'
#' @description
#' This simple utility lowers the open age group. It only returns the input
#' value vector, not the age vector.
#'
#' @param Value numeric. Vector of counts.
#' @param Age integer. Vector of age classes.
#' @param OAnew The desired open age group.
#'
#' @return Value vector potentially of reduced length up to OAG.

groupOAG <- function(Value, Age, OAnew) {
  stopifnot(OAnew %in% Age)
  N        <- length(Value[Age <= OAnew])
  Value[N] <- sum(Value[Age >= OAnew])
  Value    <- Value[1:N]
  names(Value) <- Age[1:N]
  Value
}

#' Group single ages into equal age groups of arbitrary width
#'
#' @description
#' This can be useful to check constrained sums, or as an intermediate step for
#' smoothing.
#'
#' @details
#' If you shift the groupings, then the first age groups may have a negative
#' lower bound (for example of -5). These counts would be discarded for the
#' oscillatory version of Sprague smoothing, for example, but they are preserved
#' in this function. The important thing to know is that if you shift the
#' groups, the first and last groups will not be N years wide. For example if
#' \code{shiftdown} is 1, the first age group is 4-ages wide. The ages
#' themselves are not returned but they are the name attribute of the output
#' count vector. Note this will also correctly group abridged ages into equal
#' 5-year age groups if the \code{Age} argument is explicitly given.
#'
#' \code{OAnew} (optional) must be less than or equal to \code{max(Age)} to have
#' any effect.
#'
#' @param Value numeric. Vector of single age counts.
#' @param Age integer. Vector of lower bounds of single age groups.
#' @param N integer. The desired width of resulting age groups. Default 5.
#' @param shiftdown integer. Optionally shift age groupings down by single ages.
#'   Default 0.
#' @param AgeN integer vector, otherwise calculated using \code{calcAgeN()}.
#'   Optional argument.
#' @param OAnew integer. Value of lower bound of new open age group.
#'
#' @return Vector of counts in N-year age groups.

groupAges <- function(Value,
                      Age = 1:length(Value) - 1,
                      N = 5,
                      shiftdown = 0,
                      AgeN,
                      OAnew = max(Age)) {
  if (missing(AgeN)) {
    AgeN <- calcAgeN(Age = Age,
                     N = N,
                     shiftdown = shiftdown)
  }
  out <- tapply(Value, AgeN, sum)

  # group down open age
  if (OAnew < max(AgeN)) {
    AgeA <- sort(unique(AgeN))
    out  <- groupOAG(Value = out,
                     Age = AgeA,
                     OAnew = OAnew)
  }
  out
}
