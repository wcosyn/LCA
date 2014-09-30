#!/bin/bash

echo "# Fe"
echo "# nn"
echo "0" `./moments.py hocommom_nn_0-1.Fe | grep std_dev | cut -f2 -d" "`
echo "1" `./moments.py hocommom_nn_1-1.Fe | grep std_dev | cut -f2 -d" "`
echo "2" `./moments.py hocommom_nn_2-1.Fe | grep std_dev | cut -f2 -d" "`
echo "3" `./moments.py hocommom_nn_3-1.Fe | grep std_dev | cut -f2 -d" "`
echo "4" `./moments.py hocommom_nn_4-1.Fe | grep std_dev | cut -f2 -d" "`
echo "5" `./moments.py hocommom_nn_5-1.Fe | grep std_dev | cut -f2 -d" "`
echo "6" `./moments.py hocommom_nn_6-1.Fe | grep std_dev | cut -f2 -d" "`

echo -e "\n\n"
echo "# nn all"
echo "-1" `./moments.py hocommom_nn_-1-1.Fe | grep std_dev | cut -f2 -d" "`

echo -e "\n\n"
echo "# np S=0"
echo "0" `./moments.py hocommom_np_00.Fe | grep std_dev | cut -f2 -d" "`
echo "1" `./moments.py hocommom_np_10.Fe | grep std_dev | cut -f2 -d" "`
echo "2" `./moments.py hocommom_np_20.Fe | grep std_dev | cut -f2 -d" "`
echo "3" `./moments.py hocommom_np_30.Fe | grep std_dev | cut -f2 -d" "`
echo "4" `./moments.py hocommom_np_40.Fe | grep std_dev | cut -f2 -d" "`
echo "5" `./moments.py hocommom_np_50.Fe | grep std_dev | cut -f2 -d" "`
echo "6" `./moments.py hocommom_np_60.Fe | grep std_dev | cut -f2 -d" "`

echo -e "\n\n"
echo "# np S=0 all"
echo "-1" `./moments.py hocommom_np_-10.Fe | grep std_dev | cut -f2 -d" "`

echo -e "\n\n"
echo "# np S=1"
echo "0" `./moments.py hocommom_np_01.Fe | grep std_dev | cut -f2 -d" "`
echo "1" `./moments.py hocommom_np_11.Fe | grep std_dev | cut -f2 -d" "`
echo "2" `./moments.py hocommom_np_21.Fe | grep std_dev | cut -f2 -d" "`
echo "3" `./moments.py hocommom_np_31.Fe | grep std_dev | cut -f2 -d" "`
echo "4" `./moments.py hocommom_np_41.Fe | grep std_dev | cut -f2 -d" "`
echo "5" `./moments.py hocommom_np_51.Fe | grep std_dev | cut -f2 -d" "`
echo "6" `./moments.py hocommom_np_61.Fe | grep std_dev | cut -f2 -d" "`
echo -e "\n\n"
echo "# np S=1 all"
echo "-1" `./moments.py hocommom_np_-11.Fe | grep std_dev | cut -f2 -d" "`

echo -e "\n\n"
echo "# pp"
echo "0" `./moments.py hocommom_pp_0-1.Fe | grep std_dev | cut -f2 -d" "`
echo "1" `./moments.py hocommom_pp_1-1.Fe | grep std_dev | cut -f2 -d" "`
echo "2" `./moments.py hocommom_pp_2-1.Fe | grep std_dev | cut -f2 -d" "`
echo "3" `./moments.py hocommom_pp_3-1.Fe | grep std_dev | cut -f2 -d" "`
echo "4" `./moments.py hocommom_pp_4-1.Fe | grep std_dev | cut -f2 -d" "`
echo "5" `./moments.py hocommom_pp_5-1.Fe | grep std_dev | cut -f2 -d" "`
echo "6" `./moments.py hocommom_pp_6-1.Fe | grep std_dev | cut -f2 -d" "`
echo -e "\n\n"
echo "# pp all"
echo "-1" `./moments.py hocommom_pp_-1-1.Fe | grep std_dev | cut -f2 -d" "`
