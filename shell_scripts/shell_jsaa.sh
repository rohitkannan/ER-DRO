#!/bin/bash

# Script for running J-SAA and J+-SAA instances
julia-0.6.4/bin/julia main_jsaa.jl $1 $2 $3 $4
tar -czvf case3_jsaa_$1_$2_$3_$4.tar.gz case3_jsaa
