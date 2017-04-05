#include "norm_tb.h"
#include "nucleusall.h"
#include <iostream>
#include "norm_ob.h"

int main()
{
    int A=12;
    int Z=6;
    Nucleusall nuc(".",".",A,Z);
    norm_tb nt(&nuc,true,true,true); // nuc, central,tensor,spin/iso
    norm_tb::norm_tb_params ntp = {-1, -1, -1, -1, -1, -1}; // nA,lA,nB,lB,S,T
    /*
     * norm_tb::sum_me calculates \sum_{a1,a2} \sum_{AB} C_{a1,a2}^{A} C_{a1,a2}^{B} \delta_{AB} * 2/(A*(A-1))
     * */
    std::cout << "[normtb] norm_tb::sum_me               = " << nt.sum_me(&ntp) << std::endl;
    std::cout << "[normtb] norm_tb::sum_me_corr_coefs    = " << nt.sum_me_corr_coefs(&ntp) << std::endl;
    std::cout << "[normtb] norm_tb::sum_me_3b_corr_coefs = " << nt.sum_me_3b_corr_coefs(&ntp) << std::endl;

    norm_ob no(&nuc,true,true,true);
    norm_ob::norm_ob_params nop = {-1,-1,-1,-1,0}; // nA,lA,nB,lB,t
    std::cout << "[normtb] norm_ob::sum_me_corr          = " << no.sum_me_corr(&nop) << std::endl;
    std::cout << "[normtb] norm_ob::sum_me_corr/(A-1)    = " << no.sum_me_corr(&nop)/(A-1.) << std::endl;
    return 0;
}
