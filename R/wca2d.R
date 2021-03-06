#' Runs wca2d for a single time step.
#'
#' Function to run the weighted cellular automata hydrological model for a single time step (dt).
#'
#' The parameter \code{vcheck} is used to update the time step.
#' The new value for dt is returned as \code{newdt}
#'
#' Updated 2018-04-16 to output in and out flows
#'
#' @param dem A matrix with elevation values in each cell (m)
#' @param ppt A matrix with precipitation values in each cell (mm)
#' @param wse A matrix with the water surface elevation values from the previous iteration(m)
#' @param deltat The iteration time step (s)
#' @param deltatu The update time step (s)
#' @param currt The current simulation time (from 0, s)
#' @param nullval Value for any grid cells to be ignored
#' @param mannN Manning's coefficent
#' @param cellem The cell dimensions (m). Note that this currently assumes a constant and regular grid
#' @param vcheck Run velocity check?
#' @param tolwd Water depth tolerance (m)
#' @param tolslope Slope tolerance (%)
#' @return \code{simcf}: a list of output from the Fortran routine
#' @examples
#' wca2d_1t(dem, ppt, wse, 1, 60, 0)
wca2d_1t <- function(dem, ppt, wse,
                     deltat, deltatu, currt,
                     nullvall=1e6, mannN=0.05,
                     cellem=50, vcheck=0,
                     tolwd=0.0001, tolslope=0.1) {
  simcf = .Fortran("wca2d_1t",
                   m = as.integer(gridx), n = as.integer(gridy),
                   dem = as.double(dem), ppt = as.double(ppt),
                   wse = as.double(wse),
                   otot = as.double(dem), itot = as.double(dem),
                   dt = as.double(deltat), dtu = as.double(deltatu),
                   currt = as.double(currt), nullcell = as.double(nullval),
                   mannn = as.double(mannN), cellx = as.double(cellem),
                   cellem = as.double(cellem), cella = as.double(cellem*cellem),
                   vcheck = as.integer(vcheck),
                   newdt = as.double(0), vamax = as.double(0),
                   tolwd = as.double(tolwd), tolslope = as.double(tolslope))
  return(simcf)

}

#' Runs wca2d for a set number of time steps
#'
#' Function to run the weighted cellular automata hydrological model for multiple time steps.
#' TBD.
#'
#' @param dem A matrix with elevation values in each cell (m)
#' @param ppt A matrix with precipitation values in each cell (mm)
#' @param wse A matrix with the water surface elevation values from the previous iteration(m)
#' @param deltat The iteration time step (s)
#' @param deltatu The update time step (s)
#' @param currt The current simulation time (from 0, s)
#' @param nullval Value for any grid cells to be ignored
#' @param mannN Manning's coefficent
#' @param cellem The cell dimensions (m). Note that this currently assumes a constant and regular grid
#' @return \code{simcf}: a list of output from the Fortran routine
#' @examples
#' wca2d_1t(dem, ppt, wse, 1, 60, 0)
wca2d <- function(dem, ppt, wse,
                     deltat, deltatu, currt,
                     nullvall=1e6, mannN=0.05,
                     cellem=50) {
  # simcf = .Fortran("wca2d_1t",
  #                  m = as.integer(gridx), n = as.integer(gridy),
  #                  dem = as.double(dem), ppt = as.double(ppt),
  #                  wse = as.double(wse),
  #                  dt = as.double(deltat), dtu = as.double(deltatu),
  #                  currt = as.double(currt), nullcell = as.double(nullval),
  #                  mannn = as.double(mannN), cellx = as.double(cellem),
  #                  cellem = as.double(cellem), cella = as.double(cellem*cellem),
  #                  vcheck = as.integer(vcheck),
  #                  newdt = as.double(0), vamax = as.double(0))
  # return(simcf)
  return(NULL) ## Return null for now

}
