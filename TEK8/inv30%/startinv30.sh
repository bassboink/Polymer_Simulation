#!/bin/bash

cd ~/Desktop/TEK8/inv30%

play vis_state_inv30.vmd

 if { 1 } {
	display stereo "RowInterleaved"
	axes location "LowerLeft"
        display rendermode GLSL
    }

sleep 20

play falcon.tcl

imd connect localhost 5683
