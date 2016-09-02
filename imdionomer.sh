#!/bin/bash

cd ~/Desktop/Ionomers

gnome-terminal -x ./falcon_calib.sh &

gnome-terminal -x ./ionomeralone.sh &

gnome-terminal -x  ./ionomeralonevmd.sh

