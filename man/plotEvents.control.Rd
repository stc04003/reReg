% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plotEvents.control}
\alias{plotEvents.control}
\title{Plot options for plotEvents}
\usage{
plotEvents.control(
  xlab = NULL,
  ylab = NULL,
  main = NULL,
  terminal.name = NULL,
  recurrent.name = NULL,
  recurrent.type = NULL,
  legend.position = "top",
  base_size = 12,
  cex = NULL,
  alpha = 0.7
)
}
\arguments{
\item{xlab}{a character string indicating the label for the x axis.
The default value is "Time".}

\item{ylab}{a character string indicating the label for the y axis.
The default value is "Subject".}

\item{terminal.name}{a character string indicating the label for the terminal event
displayed in the legend. The default value is "Terminal event".}

\item{recurrent.name}{a character string indicating the label for the recurrent event
displayed in the legend. The default value is "Recurrent events".}

\item{recurrent.type}{a factor indicating the labels for the different recurrent event types.
This option is only available when there are more than one types of recurrent events.
The default value is "Recurrent events 1", "Recurrent events 2", ....}

\item{legend.position}{a character string specifies the position of the legend.
The available options are "none", "left", "right", "bottom", "top",
or a two-element numeric vector specifies the coordinate of the legend.
This argument is passed to the \code{ggplot} theme environment.
The default value is "top".}

\item{base_size}{a numerical value to specify the base font size, given in pts.
This argument is passed to the \code{ggplot} theme environment.
The default value is 12.}

\item{cex}{a numerical value specifies the size of the points.}

\item{alpha}{a numerical value specifies the transparency of the points.}
}
\description{
This function provides the plotting options for the \code{plotEvents()} function.
}
\seealso{
\code{\link{plotEvents}}
}
