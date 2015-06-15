#!/bin/bash

echo -e "# Al\n# pp" > momentsAl.p
echo -e "# Al\n# pp" > momentsAl.x
echo -e "# Al\n# pp" > momentsAl.y
echo -e "# Al\n# pp" > momentsAl.z
echo '00' `./simulate.py dens_com_pp.00-1.Al 0.197327 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo '00' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsAl.x
echo '00' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsAl.y
echo '00' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsAl.z


echo -e "\n\n# pp all" >> momentsAl.p
echo -e "\n\n# pp all" >> momentsAl.x
echo -e "\n\n# pp all" >> momentsAl.y
echo -e "\n\n# pp all" >> momentsAl.z
echo "-1" `./simulate.py hocommom_pp_-1-1.Al 1.41421356237 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsAl.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsAl.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsAl.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsAl.z

echo "PP S=1 DONE..."
