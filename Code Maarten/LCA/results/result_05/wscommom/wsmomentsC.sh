#!/bin/bash


echo -e "# C\n# pp" > wsmomentsC.p
echo -e "# C\n# pp" > wsmomentsC.x
echo -e "# C\n# pp" > wsmomentsC.y
echo -e "# C\n# pp" > wsmomentsC.z
echo '0' `./wssimulate.py wscommom_0_4.C 2 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsC.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> wsmomentsC.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> wsmomentsC.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> wsmomentsC.z

echo "1" `./wssimulate.py wscommom_0_4.C 3 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsC.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsC.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsC.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsC.z


echo "2" `./wssimulate.py wscommom_0_4.C  4 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsC.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsC.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsC.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsC.z

echo "3" `./wssimulate.py wscommom_0_4.C 5 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsC.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsC.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsC.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsC.z

echo "4" `./wssimulate.py wscommom_0_4.C 6 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsC.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsC.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsC.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsC.z

echo -e "\n\n# pp all" >> wsmomentsC.p
echo -e "\n\n# pp all" >> wsmomentsC.x
echo -e "\n\n# pp all" >> wsmomentsC.y
echo -e "\n\n# pp all" >> wsmomentsC.z
echo "-1" `./wssimulate.py wscommom_0_4.C 1 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsC.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsC.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsC.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsC.z

echo "PP DONE..."
