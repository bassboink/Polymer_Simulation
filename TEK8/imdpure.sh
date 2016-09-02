#!/bin/bash

cd pure

gnome-terminal -x ./falcon_calib.sh &

gnome-terminal -x ./purealone.sh &

gnome-terminal -x  ./purealonevmd.sh

wait

