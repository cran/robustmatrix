% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{weather}
\alias{weather}
\title{Glacier weather data – Sonnblick observatory}
\format{
An array of dimension \eqn{(p,q,n)}, comprising \eqn{n = 136} observations,
each represented by a \eqn{p = 5} times \eqn{q = 12} dimensional matrix.
Observed parameters are monthly averages of
\itemize{
  \item{air pressure (AP)}
  \item{precipitation (P)}
  \item{sunshine hours (SH)}
  \item{temperature (T)}
  \item{proportion of solid precipitation (SP)}
}
from 1891 to 2022.
}
\source{
Datasource: GeoSphere Austria - \url{https://data.hub.geosphere.at}
}
\usage{
data(weather)
}
\description{
Weather data from Austria's highest weather station, situated in the Austrian Central Alps on the glaciated mountain
"Hoher Sonnblick", standing 3106 meters above sea level.
}
\keyword{datasets}
