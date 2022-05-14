#!/bin/bash

# Script for running FI-SAA instances
julia-0.6.4/bin/julia main_saa.jl $1 $2 $3
tar -czvf case1_saa_$1_$2_$3.tar.gz case1_saa
