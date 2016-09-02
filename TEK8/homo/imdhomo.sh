#!/bin/bash

cd ~/Desktop/TEK8/homo

gnome-terminal -x ./falcon_calib.sh &

gnome-terminal -x ./homoalone.sh &

gnome-terminal -x  ./homoalonevmd.sh

