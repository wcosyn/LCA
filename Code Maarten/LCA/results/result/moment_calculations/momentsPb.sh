#!/bin/bash

echo -e "# Pb\n# nn" > momentsPb.p
echo -e "# Pb\n# nn" > momentsPb.x
echo -e "# Pb\n# nn" > momentsPb.y
echo -e "# Pb\n# nn" > momentsPb.z
echo '0' `./simulate.py hocommom_nn_0-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "1" `./simulate.py hocommom_nn_1-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "2" `./simulate.py hocommom_nn_2-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "3" `./simulate.py hocommom_nn_3-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "4" `./simulate.py hocommom_nn_4-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "5" `./simulate.py hocommom_nn_5-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "5" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "5" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "5" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "6" `./simulate.py hocommom_nn_6-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "6" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "6" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "6" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo -e "\n\n# nn all" >> momentsPb.p
echo -e "\n\n# nn all" >> momentsPb.x
echo -e "\n\n# nn all" >> momentsPb.y
echo -e "\n\n# nn all" >> momentsPb.z
echo "-1" `./simulate.py hocommom_nn_-1-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "NN DONE..."

echo -e "\n\n# np S=0" >> momentsPb.p
echo -e "\n\n# np S=0" >> momentsPb.x
echo -e "\n\n# np S=0" >> momentsPb.y
echo -e "\n\n# np S=0" >> momentsPb.z
echo "0" `./simulate.py hocommom_np_00.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "1" `./simulate.py hocommom_np_10.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "2" `./simulate.py hocommom_np_20.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '2' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '2' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '2' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "3" `./simulate.py hocommom_np_30.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "4" `./simulate.py hocommom_np_40.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "5" `./simulate.py hocommom_np_50.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "5" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "5" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "5" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "6" `./simulate.py hocommom_np_60.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "6" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "6" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "6" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo -e "\n\n" >> momentsPb.p
echo -e "\n\n" >> momentsPb.x
echo -e "\n\n" >> momentsPb.y
echo -e "\n\n" >> momentsPb.z
echo "# np S=0 all" >> momentsPb.p
echo "# np S=0 all" >> momentsPb.x
echo "# np S=0 all" >> momentsPb.y
echo "# np S=0 all" >> momentsPb.z
echo "-1" `./simulate.py hocommom_np_-10.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '-1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '-1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '-1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "NP S=0 DONE..."

echo -e "\n\n" >> momentsPb.p
echo -e "\n\n" >> momentsPb.x
echo -e "\n\n" >> momentsPb.y
echo -e "\n\n" >> momentsPb.z
echo "# np S=1" >> momentsPb.p
echo "# np S=1" >> momentsPb.x
echo "# np S=1" >> momentsPb.y
echo "# np S=1" >> momentsPb.z

echo "0" `./simulate.py hocommom_np_01.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "1" `./simulate.py hocommom_np_11.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "2" `./simulate.py hocommom_np_21.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '2' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '2' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '2' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "3" `./simulate.py hocommom_np_31.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "4" `./simulate.py hocommom_np_41.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "5" `./simulate.py hocommom_np_51.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "5" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "5" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "5" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "6" `./simulate.py hocommom_np_61.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "6" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "6" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "6" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo -e "\n\n" >> momentsPb.p
echo -e "\n\n" >> momentsPb.x
echo -e "\n\n" >> momentsPb.y
echo -e "\n\n" >> momentsPb.z
echo "# np S=1 all" >> momentsPb.p
echo "# np S=1 all" >> momentsPb.x
echo "# np S=1 all" >> momentsPb.y
echo "# np S=1 all" >> momentsPb.z
echo "-1" `./simulate.py hocommom_np_-11.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '-1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '-1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '-1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "NP S=1 DONE..."

echo -e "\n\n# Pb\n# pp" >> momentsPb.p
echo -e "\n\n# Pb\n# pp" >> momentsPb.x
echo -e "\n\n# Pb\n# pp" >> momentsPb.y
echo -e "\n\n# Pb\n# pp" >> momentsPb.z
echo '0' `./simulate.py hocommom_pp_0-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "1" `./simulate.py hocommom_pp_1-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z


echo "2" `./simulate.py hocommom_pp_2-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "3" `./simulate.py hocommom_pp_3-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "4" `./simulate.py hocommom_pp_4-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "5" `./simulate.py hocommom_pp_5-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "5" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "5" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "5" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "6" `./simulate.py hocommom_pp_6-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "6" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "6" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "6" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo -e "\n\n# pp all" >> momentsPb.p
echo -e "\n\n# pp all" >> momentsPb.x
echo -e "\n\n# pp all" >> momentsPb.y
echo -e "\n\n# pp all" >> momentsPb.z
echo "-1" `./simulate.py hocommom_pp_-1-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "PP S=1 DONE..."
