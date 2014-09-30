#!/bin/bash


echo -e "# Al\n# pp" > wsmomentsAl.p
echo -e "# Al\n# pp" > wsmomentsAl.x
echo -e "# Al\n# pp" > wsmomentsAl.y
echo -e "# Al\n# pp" > wsmomentsAl.z
echo '0' `./wssimulate.py wscommom_0_4.Al 2 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsAl.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> wsmomentsAl.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> wsmomentsAl.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> wsmomentsAl.z

echo "1" `./wssimulate.py wscommom_0_4.Al 3 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsAl.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsAl.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsAl.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsAl.z


echo "2" `./wssimulate.py wscommom_0_4.Al  4 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsAl.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsAl.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsAl.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsAl.z

echo "3" `./wssimulate.py wscommom_0_4.Al 5 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsAl.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsAl.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsAl.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsAl.z

echo "4" `./wssimulate.py wscommom_0_4.Al 6 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsAl.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsAl.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsAl.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsAl.z

echo -e "\n\n# pp all" >> wsmomentsAl.p
echo -e "\n\n# pp all" >> wsmomentsAl.x
echo -e "\n\n# pp all" >> wsmomentsAl.y
echo -e "\n\n# pp all" >> wsmomentsAl.z
echo "-1" `./wssimulate.py wscommom_0_4.Al 1 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsAl.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsAl.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsAl.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsAl.z

echo "PP DONE..."
