#!/bin/bash

echo -e "# He\n# nn" > momentsHe.p
echo -e "# He\n# nn" > momentsHe.x
echo -e "# He\n# nn" > momentsHe.y
echo -e "# He\n# nn" > momentsHe.z
echo '0' `./simulate.py hocommom_nn_0-1.He | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsHe.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsHe.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsHe.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsHe.z


echo -e "\n\n# nn all" >> momentsHe.p
echo -e "\n\n# nn all" >> momentsHe.x
echo -e "\n\n# nn all" >> momentsHe.y
echo -e "\n\n# nn all" >> momentsHe.z
echo "-1" `./simulate.py hocommom_nn_-1-1.He | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsHe.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsHe.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsHe.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsHe.z

echo "NN DONE..."

echo -e "\n\n# np S=0" >> momentsHe.p
echo -e "\n\n# np S=0" >> momentsHe.x
echo -e "\n\n# np S=0" >> momentsHe.y
echo -e "\n\n# np S=0" >> momentsHe.z
echo "0" `./simulate.py hocommom_np_00.He | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsHe.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsHe.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsHe.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsHe.z


echo -e "\n\n" >> momentsHe.p
echo -e "\n\n" >> momentsHe.x
echo -e "\n\n" >> momentsHe.y
echo -e "\n\n" >> momentsHe.z
echo "# np S=0 all" >> momentsHe.p
echo "# np S=0 all" >> momentsHe.x
echo "# np S=0 all" >> momentsHe.y
echo "# np S=0 all" >> momentsHe.z
echo "-1" `./simulate.py hocommom_np_-10.He | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsHe.p
echo '-1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsHe.x
echo '-1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsHe.y
echo '-1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsHe.z

echo "NP S=0 DONE..."

echo -e "\n\n" >> momentsHe.p
echo -e "\n\n" >> momentsHe.x
echo -e "\n\n" >> momentsHe.y
echo -e "\n\n" >> momentsHe.z
echo "# np S=1" >> momentsHe.p
echo "# np S=1" >> momentsHe.x
echo "# np S=1" >> momentsHe.y
echo "# np S=1" >> momentsHe.z

echo "0" `./simulate.py hocommom_np_01.He | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsHe.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsHe.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsHe.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsHe.z


echo -e "\n\n" >> momentsHe.p
echo -e "\n\n" >> momentsHe.x
echo -e "\n\n" >> momentsHe.y
echo -e "\n\n" >> momentsHe.z
echo "# np S=1 all" >> momentsHe.p
echo "# np S=1 all" >> momentsHe.x
echo "# np S=1 all" >> momentsHe.y
echo "# np S=1 all" >> momentsHe.z
echo "-1" `./simulate.py hocommom_np_-11.He | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >>momentsHe.p
echo '-1' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsHe.x
echo '-1' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsHe.y
echo '-1' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsHe.z

echo "NP S=1 DONE..."

echo -e "\n\n# He\n# pp" >> momentsHe.p
echo -e "\n\n# He\n# pp" >> momentsHe.x
echo -e "\n\n# He\n# pp" >> momentsHe.y
echo -e "\n\n# He\n# pp" >> momentsHe.z
echo '0' `./simulate.py hocommom_pp_0-1.He | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsHe.p
echo '0' `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "`  >> momentsHe.x
echo '0' `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "`  >> momentsHe.y
echo '0' `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "`  >> momentsHe.z



echo -e "\n\n# pp all" >> momentsHe.p
echo -e "\n\n# pp all" >> momentsHe.x
echo -e "\n\n# pp all" >> momentsHe.y
echo -e "\n\n# pp all" >> momentsHe.z
echo "-1" `./simulate.py hocommom_pp_-1-1.He | tee tmp.txt | grep "P2 std_dev" | cut -f3 -d" "` >> momentsHe.p
echo "-1" `cat tmp.txt | grep "x std_dev" | cut -f3 -d" "` >> momentsHe.x
echo "-1" `cat tmp.txt | grep "y std_dev" | cut -f3 -d" "` >> momentsHe.y
echo "-1" `cat tmp.txt | grep "z std_dev" | cut -f3 -d" "` >> momentsHe.z

echo "PP S=1 DONE..."
