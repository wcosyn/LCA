#!/bin/bash
./dens ../input/ ./result_02/overlap/ $1 $2 $3 -1 -1;
./dens ../input/ ./result_02/overlap/ $1 $2 $3 0 0;
./dens ../input/ ./result_02/overlap/ $1 $2 $3 0 1;
./dens ../input/ ./result_02/overlap/ $1 $2 $3 0 2;
./dens ../input/ ./result_02/overlap/ $1 $2 $3 1 0;
./dens ../input/ ./result_02/overlap/ $1 $2 $3 1 1;
./dens ../input/ ./result_02/overlap/ $1 $2 $3 1 2;
