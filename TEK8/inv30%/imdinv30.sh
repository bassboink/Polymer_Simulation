#!/bin/bash

cd ~/Desktop/TEK8/inv30%

gnome-terminal -x ./falcon_calib.sh &

gnome-terminal -x ./inv30alone.sh &

gnome-terminal -x  ./inv30alonevmd.sh

