#!/bin/bash

# Script for running Hellinger ER-DRO instances with Algorithm 2
julia-0.6.4/bin/julia main_hellinger_scrambled.jl $1 $2 $3 $4
tar -czvf case1_hellinger_scrambled_$1_$2_$3_$4.tar.gz case1_hellinger_scrambled
