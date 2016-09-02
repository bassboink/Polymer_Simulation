#!/bin/bash

cd homo

gnome-terminal -x ./falcon_calib.sh &

gnome-terminal -x ./homoalone.sh &

gnome-terminal -x  ./homoalonevmd.sh

wait

