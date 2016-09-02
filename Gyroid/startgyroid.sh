#!/bin/bash

play vis_state_gyroid.vmd

 if { 1 } {
	display stereo "RowInterleaved"
	axes location "LowerLeft"
        display rendermode GLSL
    }

sleep 40

play falcon.tcl

imd connect localhost 5000
