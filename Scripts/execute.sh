#!/bin/bash
clear
cwd=$(pwd)

cd "/cygdrive/c/Users/Amirhossein/Google Drive/Research/Exciton/Data/Environmental Effect"

echo "Simulation started in background!"
$cwd/main.exe &

