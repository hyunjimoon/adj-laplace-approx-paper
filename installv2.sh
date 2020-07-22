# install
# v 2.0
# Run this bash file the directory where you want to install cmdStan
# with prototype functions for the embedded Laplace approximation.
#
# Not fully tested, report if you find a bug.

#! /bin/bash
# git clone https://github.com/stan-dev/cmdStan
# cd cmdStan
# git checkout laplace_approximation_stanc
# make stan-update
# cd stan/lib/stan_math
# git checkout try-laplace_approximation


git clone https://github.com/stan-dev/stanc3.git
cd stanc3
git checkout try-laplace_approximation
cd ..
git clone https://github.com/stan-dev/cmdstan.git
cd cmdstan
make stan-update
cd stan
# git checkout feature/issue-2522-laplace
cd lib/stan_math
# git checkout feature/issue-755-laplace
git checkout try-laplace_approximation
cd ..
cd ..
cd ..
# make build
STANC3=../stanc3  
