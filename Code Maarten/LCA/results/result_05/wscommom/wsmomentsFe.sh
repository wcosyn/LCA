#!/bin/bash


echo -e "# Fe\n# pp" > wsmomentsFe.p
echo -e "# Fe\n# pp" > wsmomentsFe.x
echo -e "# Fe\n# pp" > wsmomentsFe.y
echo -e "# Fe\n# pp" > wsmomentsFe.z
echo '0' `./wssimulate.py wscommom_0_4.Fe 2 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsFe.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> wsmomentsFe.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> wsmomentsFe.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> wsmomentsFe.z

echo "1" `./wssimulate.py wscommom_0_4.Fe 3 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsFe.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsFe.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsFe.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsFe.z


echo "2" `./wssimulate.py wscommom_0_4.Fe  4 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsFe.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsFe.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsFe.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsFe.z

echo "3" `./wssimulate.py wscommom_0_4.Fe 5 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsFe.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsFe.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsFe.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsFe.z

echo "4" `./wssimulate.py wscommom_0_4.Fe 6 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsFe.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsFe.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsFe.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsFe.z

echo "5" `./wssimulate.py wscommom_0_4.Fe 7 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsFe.p
echo "5" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsFe.x
echo "5" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsFe.y
echo "5" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsFe.z

echo "6" `./wssimulate.py wscommom_0_4.Fe 8 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsFe.p
echo "6" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsFe.x
echo "6" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsFe.y
echo "6" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsFe.z

echo -e "\n\n# pp all" >> wsmomentsFe.p
echo -e "\n\n# pp all" >> wsmomentsFe.x
echo -e "\n\n# pp all" >> wsmomentsFe.y
echo -e "\n\n# pp all" >> wsmomentsFe.z
echo "-1" `./wssimulate.py wscommom_0_4.Fe 1 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsFe.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsFe.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsFe.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsFe.z

echo "PP DONE..."
