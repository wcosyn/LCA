#!/bin/bash

echo -e "# Fe\n# nn" > momentsFe.p
echo -e "# Fe\n# nn" > momentsFe.x
echo -e "# Fe\n# nn" > momentsFe.y
echo -e "# Fe\n# nn" > momentsFe.z
echo '0' `./simulate.py hocommom_nn_00-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "1" `./simulate.py hocommom_nn_10-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "2" `./simulate.py hocommom_nn_20-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "3" `./simulate.py hocommom_nn_30-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo -e "\n\n# nn all" >> momentsFe.p
echo -e "\n\n# nn all" >> momentsFe.x
echo -e "\n\n# nn all" >> momentsFe.y
echo -e "\n\n# nn all" >> momentsFe.z
echo "-1" `./simulate.py hocommom_nn_0-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "NN DONE..."

echo -e "\n\n# np S=0" >> momentsFe.p
echo -e "\n\n# np S=0" >> momentsFe.x
echo -e "\n\n# np S=0" >> momentsFe.y
echo -e "\n\n# np S=0" >> momentsFe.z
echo "0" `./simulate.py hocommom_np_000.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "1" `./simulate.py hocommom_np_100.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
echo '1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "2" `./simulate.py hocommom_np_200.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
echo '2' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '2' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '2' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo -e "\n\n" >> momentsFe.p
echo -e "\n\n" >> momentsFe.x
echo -e "\n\n" >> momentsFe.y
echo -e "\n\n" >> momentsFe.z
echo "# np S=0 all" >> momentsFe.p
echo "# np S=0 all" >> momentsFe.x
echo "# np S=0 all" >> momentsFe.y
echo "# np S=0 all" >> momentsFe.z
echo "-1" `./simulate.py hocommom_np_00.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
echo '-1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '-1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '-1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "NP S=0 DONE..."

echo -e "\n\n" >> momentsFe.p
echo -e "\n\n" >> momentsFe.x
echo -e "\n\n" >> momentsFe.y
echo -e "\n\n" >> momentsFe.z
echo "# np S=1" >> momentsFe.p
echo "# np S=1" >> momentsFe.x
echo "# np S=1" >> momentsFe.y
echo "# np S=1" >> momentsFe.z

echo "0" `./simulate.py hocommom_np_001.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "1" `./simulate.py hocommom_np_101.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
echo '1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "2" `./simulate.py hocommom_np_201.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
echo '2' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '2' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '2' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo -e "\n\n" >> momentsFe.p
echo -e "\n\n" >> momentsFe.x
echo -e "\n\n" >> momentsFe.y
echo -e "\n\n" >> momentsFe.z
echo "# np S=1 all" >> momentsFe.p
echo "# np S=1 all" >> momentsFe.x
echo "# np S=1 all" >> momentsFe.y
echo "# np S=1 all" >> momentsFe.z
echo "-1" `./simulate.py hocommom_np_01.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
echo '-1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '-1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '-1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "NP S=1 DONE..."

echo -e "\n\n# Fe\n# pp" >> momentsFe.p
echo -e "\n\n# Fe\n# pp" >> momentsFe.x
echo -e "\n\n# Fe\n# pp" >> momentsFe.y
echo -e "\n\n# Fe\n# pp" >> momentsFe.z
echo '0' `./simulate.py hocommom_pp_00-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "1" `./simulate.py hocommom_pp_10-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "2" `./simulate.py hocommom_pp_20-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo -e "\n\n# pp all" >> momentsFe.p
echo -e "\n\n# pp all" >> momentsFe.x
echo -e "\n\n# pp all" >> momentsFe.y
echo -e "\n\n# pp all" >> momentsFe.z
echo "-1" `./simulate.py hocommom_pp_0-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "PP S=1 DONE..."
