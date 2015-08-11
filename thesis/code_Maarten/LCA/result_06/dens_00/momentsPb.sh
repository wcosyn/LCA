#!/bin/bash

echo -e "# Pb\n# pp" > momentsPb.p
echo -e "# Pb\n# pp" > momentsPb.x
echo -e "# Pb\n# pp" > momentsPb.y
echo -e "# Pb\n# pp" > momentsPb.z
echo '00' `./simulate.py dens_com_pp.00-1.Pb 0.197327 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo '00' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsPb.x
echo '00' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsPb.y
echo '00' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsPb.z

echo -e "\n\n# pp all" >> momentsPb.p
echo -e "\n\n# pp all" >> momentsPb.x
echo -e "\n\n# pp all" >> momentsPb.y
echo -e "\n\n# pp all" >> momentsPb.z
echo "-1" `./simulate.py hocommom_pp_-1-1.Pb 1.41421356237 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsPb.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsPb.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsPb.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsPb.z

echo "PP S=1 DONE..."
