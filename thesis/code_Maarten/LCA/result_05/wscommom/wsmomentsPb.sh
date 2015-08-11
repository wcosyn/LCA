#!/bin/bash


echo -e "# Pb\n# pp" > wsmomentsPb.p
echo -e "# Pb\n# pp" > wsmomentsPb.x
echo -e "# Pb\n# pp" > wsmomentsPb.y
echo -e "# Pb\n# pp" > wsmomentsPb.z
echo '0' `./wssimulate.py wscommom_0_2.Pb 2 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsPb.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> wsmomentsPb.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> wsmomentsPb.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> wsmomentsPb.z

echo "1" `./wssimulate.py wscommom_0_2.Pb 3 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsPb.p
echo "1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsPb.x
echo "1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsPb.y
echo "1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsPb.z


echo "2" `./wssimulate.py wscommom_0_2.Pb  4 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsPb.p
echo "2" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsPb.x
echo "2" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsPb.y
echo "2" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsPb.z

echo "3" `./wssimulate.py wscommom_0_2.Pb 5 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsPb.p
echo "3" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsPb.x
echo "3" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsPb.y
echo "3" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsPb.z

echo "4" `./wssimulate.py wscommom_0_2.Pb 6 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsPb.p
echo "4" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsPb.x
echo "4" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsPb.y
echo "4" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsPb.z

echo "5" `./wssimulate.py wscommom_0_2.Pb 7 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsPb.p
echo "5" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsPb.x
echo "5" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsPb.y
echo "5" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsPb.z

echo "6" `./wssimulate.py wscommom_0_2.Pb 8 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsPb.p
echo "6" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsPb.x
echo "6" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsPb.y
echo "6" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsPb.z

echo -e "\n\n# pp all" >> wsmomentsPb.p
echo -e "\n\n# pp all" >> wsmomentsPb.x
echo -e "\n\n# pp all" >> wsmomentsPb.y
echo -e "\n\n# pp all" >> wsmomentsPb.z
echo "-1" `./wssimulate.py wscommom_0_2.Pb 1 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> wsmomentsPb.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> wsmomentsPb.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> wsmomentsPb.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> wsmomentsPb.z

echo "PP DONE..."
