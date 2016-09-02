#!/bin/bash

cd ~/Desktop/TEK8/diblock

play vis_state_diblock.vmd

 if { 1 } {
	display stereo "RowInterleaved"
	axes location "LowerLeft"
        display rendermode GLSL
    }

sleep 25

play falcon.tcl

imd connect localhost 5682
