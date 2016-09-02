#!/bin/bash

play vis_state_30.vmd

 if { 1 } {
	display stereo "RowInterleaved"
	axes location "LowerLeft"
        display rendermode GLSL
    }

sleep 20

play falcon.tcl

imd connect localhost 5678
