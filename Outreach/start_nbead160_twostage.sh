#!/bin/bash

cd ~/Desktop/Outreach/Nbead160_twostage

play outreach160.vmd

 if { 1 } {
	display stereo "RowInterleaved"
	axes location "LowerLeft"
    }

sleep 8

play imd.tcl

play falcon.tcl






