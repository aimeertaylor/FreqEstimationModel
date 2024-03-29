% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/draw_moi_proposed.R
\name{draw_moi_proposed}
\alias{draw_moi_proposed}
\title{Propose MOI}
\usage{
draw_moi_proposed(moi_current, moi_current_equal_moi_min, moi_max)
}
\description{
Function for drawing new MOI vector. Note that there is no option to stay still, and that Proposal density does not have equal probability over state space [m_min,..,m_max] and is not symmetric: Pr(m*|m)=Pr(m|m*)
}
