#!/bin/bash

echo -e "# Fe\n# nn" > momentsFe.p
echo -e "# Fe\n# nn" > momentsFe.x
echo -e "# Fe\n# nn" > momentsFe.y
echo -e "# Fe\n# nn" > momentsFe.z
echo '0' `./simulate.py hocommom_nn_0-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "1" `./simulate.py hocommom_nn_1-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "2" `./simulate.py hocommom_nn_2-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "3" `./simulate.py hocommom_nn_3-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "4" `./simulate.py hocommom_nn_4-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "5" `./simulate.py hocommom_nn_5-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "5" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "5" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "5" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "6" `./simulate.py hocommom_nn_6-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "6" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "6" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "6" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo -e "\n\n# nn all" >> momentsFe.p
echo -e "\n\n# nn all" >> momentsFe.x
echo -e "\n\n# nn all" >> momentsFe.y
echo -e "\n\n# nn all" >> momentsFe.z
echo "-1" `./simulate.py hocommom_nn_-1-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "NN DONE..."

echo -e "\n\n# np S=0" >> momentsFe.p
echo -e "\n\n# np S=0" >> momentsFe.x
echo -e "\n\n# np S=0" >> momentsFe.y
echo -e "\n\n# np S=0" >> momentsFe.z
echo "0" `./simulate.py hocommom_np_00.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "1" `./simulate.py hocommom_np_10.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
echo '1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "2" `./simulate.py hocommom_np_20.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
echo '2' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '2' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '2' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "3" `./simulate.py hocommom_np_30.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "4" `./simulate.py hocommom_np_40.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "5" `./simulate.py hocommom_np_50.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "5" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "5" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "5" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "6" `./simulate.py hocommom_np_60.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "6" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "6" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "6" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo -e "\n\n" >> momentsFe.p
echo -e "\n\n" >> momentsFe.x
echo -e "\n\n" >> momentsFe.y
echo -e "\n\n" >> momentsFe.z
echo "# np S=0 all" >> momentsFe.p
echo "# np S=0 all" >> momentsFe.x
echo "# np S=0 all" >> momentsFe.y
echo "# np S=0 all" >> momentsFe.z
echo "-1" `./simulate.py hocommom_np_-10.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
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

echo "0" `./simulate.py hocommom_np_01.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "1" `./simulate.py hocommom_np_11.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
echo '1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "2" `./simulate.py hocommom_np_21.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
echo '2' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '2' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '2' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "3" `./simulate.py hocommom_np_31.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "4" `./simulate.py hocommom_np_41.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "5" `./simulate.py hocommom_np_51.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "5" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "5" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "5" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "6" `./simulate.py hocommom_np_61.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "6" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "6" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "6" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo -e "\n\n" >> momentsFe.p
echo -e "\n\n" >> momentsFe.x
echo -e "\n\n" >> momentsFe.y
echo -e "\n\n" >> momentsFe.z
echo "# np S=1 all" >> momentsFe.p
echo "# np S=1 all" >> momentsFe.x
echo "# np S=1 all" >> momentsFe.y
echo "# np S=1 all" >> momentsFe.z
echo "-1" `./simulate.py hocommom_np_-11.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsFe.p
echo '-1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '-1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '-1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "NP S=1 DONE..."

echo -e "\n\n# Fe\n# pp" >> momentsFe.p
echo -e "\n\n# Fe\n# pp" >> momentsFe.x
echo -e "\n\n# Fe\n# pp" >> momentsFe.y
echo -e "\n\n# Fe\n# pp" >> momentsFe.z
echo '0' `./simulate.py hocommom_pp_0-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo "1" `./simulate.py hocommom_pp_1-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z


echo "2" `./simulate.py hocommom_pp_2-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "3" `./simulate.py hocommom_pp_3-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "4" `./simulate.py hocommom_pp_4-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "5" `./simulate.py hocommom_pp_5-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "5" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "5" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "5" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "6" `./simulate.py hocommom_pp_6-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "6" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "6" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "6" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo -e "\n\n# pp all" >> momentsFe.p
echo -e "\n\n# pp all" >> momentsFe.x
echo -e "\n\n# pp all" >> momentsFe.y
echo -e "\n\n# pp all" >> momentsFe.z
echo "-1" `./simulate.py hocommom_pp_-1-1.Fe | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "PP S=1 DONE..."
