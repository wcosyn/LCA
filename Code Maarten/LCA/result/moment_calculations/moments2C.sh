#!/bin/bash

echo "# C"
echo "# nn"
echo "0" `./moments.py hocommom_nn_0-1.C | grep std_dev | cut -f2 -d" "`
echo "1" `./moments.py hocommom_nn_1-1.C | grep std_dev | cut -f2 -d" "`
echo "2" `./moments.py hocommom_nn_2-1.C | grep std_dev | cut -f2 -d" "`

echo -e "\n\n"
echo "# nn all"
echo "-1" `./moments.py hocommom_nn_-1-1.C | grep std_dev | cut -f2 -d" "`

echo -e "\n\n"
echo "# np S=0"
echo "0" `./moments.py hocommom_np_00.C | grep std_dev | cut -f2 -d" "`
echo "1" `./moments.py hocommom_np_10.C | grep std_dev | cut -f2 -d" "`
echo "2" `./moments.py hocommom_np_20.C | grep std_dev | cut -f2 -d" "`
echo -e "\n\n"
echo "# np S=0 all"
echo "-1" `./moments.py hocommom_np_-10.C | grep std_dev | cut -f2 -d" "`

echo -e "\n\n"
echo "# np S=1"
echo "0" `./moments.py hocommom_np_01.C | grep std_dev | cut -f2 -d" "`
echo "1" `./moments.py hocommom_np_11.C | grep std_dev | cut -f2 -d" "`
echo "2" `./moments.py hocommom_np_21.C | grep std_dev | cut -f2 -d" "`
echo -e "\n\n"
echo "# np S=1 all"
echo "-1" `./moments.py hocommom_np_-11.C | grep std_dev | cut -f2 -d" "`

echo -e "\n\n"
echo "# pp"
echo "0" `./moments.py hocommom_pp_0-1.C | grep std_dev | cut -f2 -d" "`
echo "1" `./moments.py hocommom_pp_1-1.C | grep std_dev | cut -f2 -d" "`
echo "2" `./moments.py hocommom_pp_2-1.C | grep std_dev | cut -f2 -d" "`
echo -e "\n\n"
echo "# pp all"
echo "-1" `./moments.py hocommom_pp_-1-1.C | grep std_dev | cut -f2 -d" "`
