#!/bin/bash
./plot.py ../figures/He_tb_cm.pdf He_tb_cm.txt He
./plot.py ../figures/Be_tb_cm.pdf Be_tb_cm.txt Be
./plot.py ../figures/C_tb_cm.pdf C_tb_cm.txt C
./plot.py ../figures/O_tb_cm.pdf O_tb_cm.txt O
./plot.py ../figures/Al_tb_cm.pdf Al_tb_cm.txt Al
./plot.py ../figures/Ca40_tb_cm.pdf Ca40_tb_cm.txt Ca40
./plot.py ../figures/Ca48_tb_cm.pdf Ca48_tb_cm.txt Ca48
./plot.py ../figures/Fe_tb_cm.pdf Fe_tb_cm.txt Fe
./plot_multi.py ../figures/multi_tb_cm.pdf Be_tb_cm.txt C_tb_cm.txt O_tb_cm.txt Al_tb_cm.txt Ar_tb_cm.txt Ca40_tb_cm.txt 

