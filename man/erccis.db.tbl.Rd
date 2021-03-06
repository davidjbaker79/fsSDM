% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/species_data_processing.R
\name{erccis.db.tbl}
\alias{erccis.db.tbl}
\title{Load ERCCIS data}
\usage{
erccis.db.tbl(sqlite_file, table_name, sppGrps = NULL, sppNameLookUp = NULL)
}
\arguments{
\item{sqlite_file}{A path name to .sqlite file.}

\item{table_name}{The name of .sqlite table to load (if multiple)}

\item{sppGrps}{Species group required.}

\item{sppNameLookUp}{A table with species names can be provided for filtering. A column can also
be added to add better names (e.g. for plotting).}
}
\value{
An tibble data.frame
}
\description{
Create a connection with the ERCCIS SQLite db and
filter missing coordinates and create required date formats
}
