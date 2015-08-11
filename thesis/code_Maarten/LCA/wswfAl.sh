#!/usr/bin/sh

echo 0 0 1 -1
./wswf . Al 27 13 0 0 1 -1 -50 > cout
echo 0 1 3 -1
./wswf . Al 27 13 0 1 3 -1 -35 > cout
echo 0 1 1 -1
./wswf . Al 27 13 0 1 1 -1 -20 > cout
echo 0 2 5 -1
./wswf . Al 27 13 0 2 5 -1 -10 > cout

echo
echo 0 0 1 1
./wswf . Al 27 13 0 0 1 1 -50 > cout
echo 0 1 3 1
./wswf . Al 27 13 0 1 3 1 -30 > cout
echo 0 1 1 1
./wswf . Al 27 13 0 1 1 1 -20 > cout
echo 0 2 5 1
./wswf . Al 27 13 0 2 5 1 -10 > cout
