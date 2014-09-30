#!/bin/bash

echo -e "# Pb\n# nn" > momentsPb.p
echo -e "# Pb\n# nn" > momentsPb.x
echo -e "# Pb\n# nn" > momentsPb.y
echo -e "# Pb\n# nn" > momentsPb.z
echo '0' `./simulate.py hocommom_nn_00-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "1" `./simulate.py hocommom_nn_10-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "2" `./simulate.py hocommom_nn_20-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "3" `./simulate.py hocommom_nn_30-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "4" `./simulate.py hocommom_nn_40-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo -e "\n\n# nn all" >> momentsPb.p
echo -e "\n\n# nn all" >> momentsPb.x
echo -e "\n\n# nn all" >> momentsPb.y
echo -e "\n\n# nn all" >> momentsPb.z
echo "-1" `./simulate.py hocommom_nn_0-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "NN DONE..."

echo -e "\n\n# np S=0" >> momentsPb.p
echo -e "\n\n# np S=0" >> momentsPb.x
echo -e "\n\n# np S=0" >> momentsPb.y
echo -e "\n\n# np S=0" >> momentsPb.z
echo "0" `./simulate.py hocommom_np_000.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "1" `./simulate.py hocommom_np_100.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "2" `./simulate.py hocommom_np_200.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '2' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '2' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '2' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "3" `./simulate.py hocommom_np_300.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '3' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '3' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '3' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "4" `./simulate.py hocommom_np_400.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '4' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '4' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '4' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo -e "\n\n" >> momentsPb.p
echo -e "\n\n" >> momentsPb.x
echo -e "\n\n" >> momentsPb.y
echo -e "\n\n" >> momentsPb.z
echo "# np S=0 all" >> momentsPb.p
echo "# np S=0 all" >> momentsPb.x
echo "# np S=0 all" >> momentsPb.y
echo "# np S=0 all" >> momentsPb.z
echo "-1" `./simulate.py hocommom_np_00.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
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

echo "0" `./simulate.py hocommom_np_001.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "1" `./simulate.py hocommom_np_101.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "2" `./simulate.py hocommom_np_201.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '2' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '2' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '2' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "3" `./simulate.py hocommom_np_301.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '3' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '3' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '3' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "4" `./simulate.py hocommom_np_401.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '4' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '4' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '4' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo -e "\n\n" >> momentsPb.p
echo -e "\n\n" >> momentsPb.x
echo -e "\n\n" >> momentsPb.y
echo -e "\n\n" >> momentsPb.z
echo "# np S=1 all" >> momentsPb.p
echo "# np S=1 all" >> momentsPb.x
echo "# np S=1 all" >> momentsPb.y
echo "# np S=1 all" >> momentsPb.z
echo "-1" `./simulate.py hocommom_np_01.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsPb.p
echo '-1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '-1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '-1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "NP S=1 DONE..."

echo -e "\n\n# Pb\n# pp" >> momentsPb.p
echo -e "\n\n# Pb\n# pp" >> momentsPb.x
echo -e "\n\n# Pb\n# pp" >> momentsPb.y
echo -e "\n\n# Pb\n# pp" >> momentsPb.z
echo '0' `./simulate.py hocommom_pp_00-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo "1" `./simulate.py hocommom_pp_10-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "2" `./simulate.py hocommom_pp_20-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "3" `./simulate.py hocommom_pp_30-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "4" `./simulate.py hocommom_pp_40-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo -e "\n\n# pp all" >> momentsPb.p
echo -e "\n\n# pp all" >> momentsPb.x
echo -e "\n\n# pp all" >> momentsPb.y
echo -e "\n\n# pp all" >> momentsPb.z
echo "-1" `./simulate.py hocommom_pp_0-1.Pb | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "PP S=1 DONE..."
