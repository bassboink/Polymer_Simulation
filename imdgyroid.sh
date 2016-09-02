#!/bin/bash

cd ~/Desktop/Gyroid

gnome-terminal -x ./falcon_calib.sh &

gnome-terminal -x ./gyroidalone.sh &

gnome-terminal -x  ./gyroidalonevmd.sh

