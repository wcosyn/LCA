#!/bin/bash
./plot.py ../figures/He_tb_rel.pdf He_tb_rel.txt He
./plot.py ../figures/Be_tb_rel.pdf Be_tb_rel.txt Be
./plot.py ../figures/C_tb_rel.pdf C_tb_rel.txt C
./plot.py ../figures/O_tb_rel.pdf O_tb_rel.txt O
./plot.py ../figures/Al_tb_rel.pdf Al_tb_rel.txt Al
./plot.py ../figures/Ca40_tb_rel.pdf Ca40_tb_rel.txt Ca40
./plot.py ../figures/Ca48_tb_rel.pdf Ca48_tb_rel.txt Ca48
./plot.py ../figures/Fe_tb_rel.pdf Fe_tb_rel.txt Fe
./plot_multi.py ../figures/multi_tb_rel.pdf Be_tb_rel.txt C_tb_rel.txt O_tb_rel.txt Al_tb_rel.txt Ar_tb_rel.txt Ca40_tb_rel.txt 

