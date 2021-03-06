% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matlab.R
\name{readMATmod}
\alias{readMATmod}
\title{Read matlab model}
\usage{
readMATmod(file)
}
\arguments{
\item{file}{Full path to matlab model file}
}
\description{
The generic function \code{readMATmod} imports matlab cobra models into sybil model files
}
\details{
Returns sybil model object (time needed: bacterial model ~ 10s, recon2 ~ 60s)
}
