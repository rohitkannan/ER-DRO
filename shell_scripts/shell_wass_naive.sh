#!/bin/bash

# Script for running Wasserstein ER-DRO instances with Algorithm 1
julia-0.6.4/bin/julia main_wass_naive.jl $1 $2 $3 $4
tar -czvf case1_wass_naive_$1_$2_$3_$4.tar.gz case1_wass_naive
