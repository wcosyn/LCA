#!/bin/bash

echo -e "# Al\n# nn" > momentsAl.p
echo -e "# Al\n# nn" > momentsAl.x
echo -e "# Al\n# nn" > momentsAl.y
echo -e "# Al\n# nn" > momentsAl.z
echo '0' `./simulate.py hocommom_nn_0-1.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsAl.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsAl.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsAl.z

echo "1" `./simulate.py hocommom_nn_1-1.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsAl.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsAl.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsAl.z

echo "2" `./simulate.py hocommom_nn_2-1.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsAl.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsAl.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsAl.z

echo "3" `./simulate.py hocommom_nn_3-1.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsAl.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsAl.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsAl.z

echo "4" `./simulate.py hocommom_nn_4-1.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsAl.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsAl.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsAl.z

echo -e "\n\n# nn all" >> momentsAl.p
echo -e "\n\n# nn all" >> momentsAl.x
echo -e "\n\n# nn all" >> momentsAl.y
echo -e "\n\n# nn all" >> momentsAl.z
echo "-1" `./simulate.py hocommom_nn_-1-1.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsAl.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsAl.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsAl.z

echo "NN DONE..."

echo -e "\n\n# np S=0" >> momentsAl.p
echo -e "\n\n# np S=0" >> momentsAl.x
echo -e "\n\n# np S=0" >> momentsAl.y
echo -e "\n\n# np S=0" >> momentsAl.z
echo "0" `./simulate.py hocommom_np_00.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsAl.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsAl.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsAl.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsAl.z

echo "1" `./simulate.py hocommom_np_10.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsAl.p
echo '1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsAl.x
echo '1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsAl.y
echo '1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsAl.z

echo "2" `./simulate.py hocommom_np_20.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsAl.p
echo '2' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsAl.x
echo '2' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsAl.y
echo '2' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsAl.z

echo "3" `./simulate.py hocommom_np_30.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsAl.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsAl.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsAl.z

echo "4" `./simulate.py hocommom_np_40.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsAl.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsAl.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsAl.z

echo -e "\n\n" >> momentsAl.p
echo -e "\n\n" >> momentsAl.x
echo -e "\n\n" >> momentsAl.y
echo -e "\n\n" >> momentsAl.z
echo "# np S=0 all" >> momentsAl.p
echo "# np S=0 all" >> momentsAl.x
echo "# np S=0 all" >> momentsAl.y
echo "# np S=0 all" >> momentsAl.z
echo "-1" `./simulate.py hocommom_np_-10.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsAl.p
echo '-1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsAl.x
echo '-1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsAl.y
echo '-1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsAl.z

echo "NP S=0 DONE..."

echo -e "\n\n" >> momentsAl.p
echo -e "\n\n" >> momentsAl.x
echo -e "\n\n" >> momentsAl.y
echo -e "\n\n" >> momentsAl.z
echo "# np S=1" >> momentsAl.p
echo "# np S=1" >> momentsAl.x
echo "# np S=1" >> momentsAl.y
echo "# np S=1" >> momentsAl.z

echo "0" `./simulate.py hocommom_np_01.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsAl.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsAl.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsAl.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsAl.z

echo "1" `./simulate.py hocommom_np_11.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsAl.p
echo '1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsAl.x
echo '1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsAl.y
echo '1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsAl.z

echo "2" `./simulate.py hocommom_np_21.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsAl.p
echo '2' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsAl.x
echo '2' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsAl.y
echo '2' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsAl.z

echo "3" `./simulate.py hocommom_np_31.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsAl.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsAl.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsAl.z

echo "4" `./simulate.py hocommom_np_41.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsAl.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsAl.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsAl.z

echo -e "\n\n" >> momentsAl.p
echo -e "\n\n" >> momentsAl.x
echo -e "\n\n" >> momentsAl.y
echo -e "\n\n" >> momentsAl.z
echo "# np S=1 all" >> momentsAl.p
echo "# np S=1 all" >> momentsAl.x
echo "# np S=1 all" >> momentsAl.y
echo "# np S=1 all" >> momentsAl.z
echo "-1" `./simulate.py hocommom_np_-11.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsAl.p
echo '-1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsAl.x
echo '-1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsAl.y
echo '-1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsAl.z

echo "NP S=1 DONE..."

echo -e "\n\n# Al\n# pp" >> momentsAl.p
echo -e "\n\n# Al\n# pp" >> momentsAl.x
echo -e "\n\n# Al\n# pp" >> momentsAl.y
echo -e "\n\n# Al\n# pp" >> momentsAl.z
echo '0' `./simulate.py hocommom_pp_0-1.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsAl.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsAl.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsAl.z

echo "1" `./simulate.py hocommom_pp_1-1.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsAl.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsAl.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsAl.z


echo "2" `./simulate.py hocommom_pp_2-1.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsAl.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsAl.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsAl.z

echo "3" `./simulate.py hocommom_pp_3-1.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsAl.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsAl.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsAl.z

echo "4" `./simulate.py hocommom_pp_4-1.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsAl.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsAl.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsAl.z

echo -e "\n\n# pp all" >> momentsAl.p
echo -e "\n\n# pp all" >> momentsAl.x
echo -e "\n\n# pp all" >> momentsAl.y
echo -e "\n\n# pp all" >> momentsAl.z
echo "-1" `./simulate.py hocommom_pp_-1-1.Al | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsAl.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsAl.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsAl.z

echo "PP S=1 DONE..."
