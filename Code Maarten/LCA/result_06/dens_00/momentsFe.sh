#!/bin/bash

echo -e "# Fe\n# pp" > momentsFe.p
echo -e "# Fe\n# pp" > momentsFe.x
echo -e "# Fe\n# pp" > momentsFe.y
echo -e "# Fe\n# pp" > momentsFe.z
echo '00' `./simulate.py dens_com_pp.00-1.Fe 0.197327 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo '00' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsFe.x
echo '00' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsFe.y
echo '00' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsFe.z

echo -e "\n\n# pp all" >> momentsFe.p
echo -e "\n\n# pp all" >> momentsFe.x
echo -e "\n\n# pp all" >> momentsFe.y
echo -e "\n\n# pp all" >> momentsFe.z
echo "-1" `./simulate.py hocommom_pp_-1-1.Fe 1.41421356237 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsFe.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsFe.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsFe.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsFe.z

echo "PP S=1 DONE..."
