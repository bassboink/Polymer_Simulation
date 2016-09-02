#!/bin/bash

play vis_state_ionomer.vmd

 if { 1 } {
	display stereo "RowInterleaved"
	axes location "LowerLeft"
        display rendermode GLSL
    }

sleep 8

imd connect localhost 4000

pbc get -now

pbc set {17.5 17.5 17.5} -all

lhall

play falcon.tcl
