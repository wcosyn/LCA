#!/bin/bash

echo -e "# He\n# np" > momentsHe.p
echo -e "# He\n# np" > momentsHe.x
echo -e "# He\n# np" > momentsHe.y
echo -e "# He\n# np" > momentsHe.z
echo '00' `./simulate.py dens_com_np.00-1.He 0.197327 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsHe.p
echo '00' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsHe.x
echo '00' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsHe.y
echo '00' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsHe.z

echo -e "\n\n# np all" >> momentsHe.p
echo -e "\n\n# np all" >> momentsHe.x
echo -e "\n\n# np all" >> momentsHe.y
echo -e "\n\n# np all" >> momentsHe.z
echo "-1" `./simulate.py hocommom_np_-10.He 1.41421356237 | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsHe.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsHe.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsHe.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsHe.z

echo "np S=1 DONE..."
