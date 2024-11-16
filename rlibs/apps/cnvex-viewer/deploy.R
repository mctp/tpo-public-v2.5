if('cnvex' %in% (.packages())){detach("package:cnvex",unload=T)}
# if the local repo is not in options(), add it
if(!("file:///mctp/share/cnvex/cran" %in% options('repos')$repos)) {
    options(repos = c(options('repos')$repos, "cnvex-fake-cran"="file:///mctp/share/cnvex/cran"))
}

devtools::build(pkg='./rlibs/tpolib', path='/mctp/share/cnvex/cran/src/contrib')
devtools::build(pkg='./rlibs/cnvex', path='/mctp/share/cnvex/cran/src/contrib')
tools::write_PACKAGES('/mctp/share/cnvex/cran/src/contrib', type="source")

remotes::install_version('tpolib',version='0.9.9')
remotes::install_version('cnvex',version='1.1.8')
