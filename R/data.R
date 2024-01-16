#' Glacier weather data – Sonnblick observatory
#'
#' Weather data from Austria's highest weather station, situated in the Austrian Central Alps on the glaciated mountain
#' "Hoher Sonnblick", standing 3106 meters above sea level.
#'
#' @usage data(weather)
#'
#' @format
#' An array of dimension \eqn{(p,q,n)}, comprising \eqn{n = 136} observations,
#' each represented by a \eqn{p = 5} times \eqn{q = 12} dimensional matrix.
#' Observed parameters are monthly averages of
#' \itemize{
#'   \item{air pressure (AP)}
#'   \item{precipitation (P)}
#'   \item{sunshine hours (SH)}
#'   \item{temperature (T)}
#'   \item{proportion of solid precipitation (SP)}
#' }
#' from 1891 to 2022.
#'
#' @source Datasource: GeoSphere Austria - \url{https://data.hub.geosphere.at}
"weather"


#' DARWIN (Diagnosis AlzheimeR WIth haNdwriting)
#'
#' The DARWIN (Diagnosis AlzheimeR WIth haNdwriting) dataset
#' comprises handwriting samples from 174 individuals.
#' Among them, 89 have been diagnosed with Alzheimer's disease (AD), while the remaining 85 are considered healthy subjects (H).
#' Each participant completed 25 handwriting tasks on paper, and their pen movements were recorded using a graphic tablet.
#' From the raw handwriting data, a set of 18 features was extracted.
#'
#' @usage data(darwin)
#'
#' @format
#' An array of dimension \eqn{(p,q,n)}, comprising \eqn{n = 174} observations,
#' each represented by a \eqn{p = 18} times \eqn{q = 25} dimensional matrix.
#' The observed parameters are:
#' \itemize{
#'   \item{Total Time}
#'   \item{Air Time}
#'   \item{Paper Time}
#'   \item{Mean Speed on paper}
#'   \item{Mean Acceleration on paper}
#'   \item{Mean Acceleration in air}
#'   \item{Mean Jerk on paper}
#'   \item{Pressure Mean}
#'   \item{Pressure Variance}
#'   \item{Generalization of the Mean Relative Tremor (GMRT) on paper}
#'   \item{GMTR in air}
#'   \item{Mean GMRT}
#'   \item{Pendowns Number}
#'   \item{Max X Extension}
#'   \item{Max Y Extension}
#'   \item{Dispersion Index}
#' }
#'
#' @references
#' \insertRef{cilia2018experimental}{robustmatrix} \cr
#' \insertRef{cilia2022diagnosing}{robustmatrix}
#'
#' @source UC Irvine Machine Learning Repository - DARWIN - \doi{https://doi.org/10.24432/C55D0K}
"darwin"
