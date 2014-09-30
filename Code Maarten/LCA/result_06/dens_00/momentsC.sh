#!/bin/bash

echo -e "# C\n# pp" > momentsC.p
echo -e "# C\n# pp" > momentsC.x
echo -e "# C\n# pp" > momentsC.y
echo -e "# C\n# pp" > momentsC.z
echo '00' `./simulate.py dens_com_pp.00-1.C 0.197327 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsC.p
echo '00' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsC.x
echo '00' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsC.y
echo '00' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsC.z

echo -e "\n\n# pp all" >> momentsC.p
echo -e "\n\n# pp all" >> momentsC.x
echo -e "\n\n# pp all" >> momentsC.y
echo -e "\n\n# pp all" >> momentsC.z
echo "-1" `./simulate.py hocommom_pp_-1-1.C 1.41421356237 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsC.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsC.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsC.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsC.z

echo "PP S=1 DONE..."
