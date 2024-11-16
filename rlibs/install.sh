#!/usr/bin/env bash
R CMD INSTALL -l /code/rlibs/install /code/rlibs/tpolib
R CMD INSTALL -l /code/rlibs/install /code/rlibs/cnvex &
R CMD INSTALL -l /code/rlibs/install /code/rlibs/carat &
R CMD INSTALL -l /code/rlibs/install /code/rlibs/codac &
wait
R CMD INSTALL -l /code/rlibs/install /code/rlibs/cargo &
R CMD INSTALL -l /code/rlibs/install /code/rlibs/curve &
wait
chmod -R 777 /code/rlibs/install || :
