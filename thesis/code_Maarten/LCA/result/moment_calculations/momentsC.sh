#!/bin/bash

echo -e "# C\n# nn" > momentsC.p
echo -e "# C\n# nn" > momentsC.x
echo -e "# C\n# nn" > momentsC.y
echo -e "# C\n# nn" > momentsC.z
echo '0' `./simulate.py hocommom_nn_0-1.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsC.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsC.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsC.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsC.z

echo "1" `./simulate.py hocommom_nn_1-1.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsC.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsC.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsC.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsC.z

echo "2" `./simulate.py hocommom_nn_2-1.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsC.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsC.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsC.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsC.z


echo -e "\n\n# nn all" >> momentsC.p
echo -e "\n\n# nn all" >> momentsC.x
echo -e "\n\n# nn all" >> momentsC.y
echo -e "\n\n# nn all" >> momentsC.z
echo "-1" `./simulate.py hocommom_nn_-1-1.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsC.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsC.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsC.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsC.z

echo "NN DONE..."

echo -e "\n\n# np S=0" >> momentsC.p
echo -e "\n\n# np S=0" >> momentsC.x
echo -e "\n\n# np S=0" >> momentsC.y
echo -e "\n\n# np S=0" >> momentsC.z
echo "0" `./simulate.py hocommom_np_00.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsC.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsC.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsC.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsC.z

echo "1" `./simulate.py hocommom_np_10.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsC.p
echo '1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsC.x
echo '1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsC.y
echo '1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsC.z

echo "2" `./simulate.py hocommom_np_20.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsC.p
echo '2' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsC.x
echo '2' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsC.y
echo '2' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsC.z


echo -e "\n\n" >> momentsC.p
echo -e "\n\n" >> momentsC.x
echo -e "\n\n" >> momentsC.y
echo -e "\n\n" >> momentsC.z
echo "# np S=0 all" >> momentsC.p
echo "# np S=0 all" >> momentsC.x
echo "# np S=0 all" >> momentsC.y
echo "# np S=0 all" >> momentsC.z
echo "-1" `./simulate.py hocommom_np_-10.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsC.p
echo '-1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsC.x
echo '-1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsC.y
echo '-1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsC.z

echo "NP S=0 DONE..."

echo -e "\n\n" >> momentsC.p
echo -e "\n\n" >> momentsC.x
echo -e "\n\n" >> momentsC.y
echo -e "\n\n" >> momentsC.z
echo "# np S=1" >> momentsC.p
echo "# np S=1" >> momentsC.x
echo "# np S=1" >> momentsC.y
echo "# np S=1" >> momentsC.z

echo "0" `./simulate.py hocommom_np_01.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsC.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsC.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsC.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsC.z

echo "1" `./simulate.py hocommom_np_11.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsC.p
echo '1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsC.x
echo '1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsC.y
echo '1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsC.z

echo "2" `./simulate.py hocommom_np_21.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsC.p
echo '2' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsC.x
echo '2' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsC.y
echo '2' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsC.z


echo -e "\n\n" >> momentsC.p
echo -e "\n\n" >> momentsC.x
echo -e "\n\n" >> momentsC.y
echo -e "\n\n" >> momentsC.z
echo "# np S=1 all" >> momentsC.p
echo "# np S=1 all" >> momentsC.x
echo "# np S=1 all" >> momentsC.y
echo "# np S=1 all" >> momentsC.z
echo "-1" `./simulate.py hocommom_np_-11.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsC.p
echo '-1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsC.x
echo '-1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsC.y
echo '-1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsC.z

echo "NP S=1 DONE..."

echo -e "\n\n# C\n# pp" >> momentsC.p
echo -e "\n\n# C\n# pp" >> momentsC.x
echo -e "\n\n# C\n# pp" >> momentsC.y
echo -e "\n\n# C\n# pp" >> momentsC.z
echo '0' `./simulate.py hocommom_pp_0-1.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsC.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsC.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsC.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsC.z

echo "1" `./simulate.py hocommom_pp_1-1.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsC.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsC.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsC.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsC.z


echo "2" `./simulate.py hocommom_pp_2-1.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsC.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsC.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsC.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsC.z

echo "3" `./simulate.py hocommom_pp_3-1.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsC.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsC.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsC.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsC.z

echo "4" `./simulate.py hocommom_pp_4-1.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsC.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsC.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsC.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsC.z

echo -e "\n\n# pp all" >> momentsC.p
echo -e "\n\n# pp all" >> momentsC.x
echo -e "\n\n# pp all" >> momentsC.y
echo -e "\n\n# pp all" >> momentsC.z
echo "-1" `./simulate.py hocommom_pp_-1-1.C | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsC.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsC.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsC.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsC.z

echo "PP S=1 DONE..."
