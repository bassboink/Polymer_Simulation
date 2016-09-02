#!/bin/bash

cd ~/Desktop/TEK8/homo

play vis_state_homo.vmd

 if { 1 } {
	display stereo "RowInterleaved"
	axes location "LowerLeft"
        display rendermode GLSL
    }


sleep 8

play falcon.tcl

imd connect localhost 3690
