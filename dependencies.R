install.packages('pacman', repos = 'https://cran.stat.unipd.it/')

pacman::p_load('rgl', 'Matrix', 'plot3D', 'RcppEigen', 'Rcpp', 'MASS', 'testthat',
               'rcmdcheck', 'devtools', 'webshot')
webshot::install_phantomjs(force = TRUE)
Sys.setenv(OPENSSL_CONF="/dev/null")