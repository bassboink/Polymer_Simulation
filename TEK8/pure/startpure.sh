#!/bin/bash

cd ~/Desktop/TEK8/pure

play vis_state_pure.vmd

 if { 1 } {
	display stereo "RowInterleaved"
	axes location "LowerLeft"
        display rendermode GLSL
    }


sleep 25

play falcon.tcl

imd connect localhost 3691
