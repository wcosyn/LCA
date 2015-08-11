#!/usr/bin/sh

echo 0 0 1 -1
./wswf . Fe 56 26 0 0 1 -1 -60 > cout
echo 0 1 3 -1
./wswf . Fe 56 26 0 1 3 -1 -45 > cout
echo 0 1 1 -1
./wswf . Fe 56 26 0 1 1 -1 -40 > cout
echo 0 2 5 -1
./wswf . Fe 56 26 0 2 5 -1 -30 > cout
echo 1 0 1 -1
./wswf . Fe 56 26 1 0 1 -1 -20 > cout
echo 0 2 3 -1
./wswf . Fe 56 26 0 2 3 -1 -20 > cout
echo 0 3 7 -1
./wswf . Fe 56 26 0 3 7 -1 -15 > cout
echo 1 1 3 -1
./wswf . Fe 56 26 1 1 3 -1 -10 > cout

echo 0 0 1 1
./wswf . Fe 56 26 0 0 1 1 -51 > cout
echo 0 1 3 1
./wswf . Fe 56 26 0 1 3 1 -35 > cout
echo 0 1 1 1
./wswf . Fe 56 26 0 1 1 1 -35 > cout
echo 0 2 5 1
./wswf . Fe 56 26 0 2 5 1 -20 > cout
echo 1 0 1 1
./wswf . Fe 56 26 1 0 1 1 -15 > cout
echo 0 2 3 1
./wswf . Fe 56 26 0 2 3 1 -15 > cout
echo 0 3 7 1
./wswf . Fe 56 26 0 3 7 1 -10 > cout
