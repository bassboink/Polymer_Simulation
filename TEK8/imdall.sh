#!/bin/bash

cd ~/Desktop/TEK8/pure

gnome-terminal -x ./falcon_calib.sh &

gnome-terminal -x ./purealone.sh &

gnome-terminal -x  ./purealonevmd.sh

wait

cd ~/Desktop/TEK8/diblock

gnome-terminal -x ./falcon_calib.sh &

gnome-terminal -x ./diblockalone.sh &

gnome-terminal -x  ./diblockalonevmd.sh





