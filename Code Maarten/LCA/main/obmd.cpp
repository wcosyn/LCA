#include "density_ob3.h"
#include "nucleusall.h"
#include "nucleusnp.h"
#include "nucleuspp.h"
#include "norm_ob.h"
#include <string>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include "threej.h"


void ob(int A,int Z,std::string name, int isospin){

    // the inputdir and resultdir are only used for storage of the moshinsky brakets, always, everywhere
    //Nucleusall nuc("../data/mosh","../data/mosh",A,Z);   // inputdir,resultdir,A,Z

    NucleusPP  nuc("../data/mosh","../data/mosh",A,Z); // idem
    //NucleusNP  nucnp("../data/mosh","../data/mosh",A,Z); // idem


    norm_ob no(&nuc);
    norm_ob::norm_ob_params nob= {-1, -1, -1, -1, isospin}; // nA,lA,nB,lB,t
    double norm_mf  = no.sum_me_pairs( &nob );
    double norm_corr= no.sum_me_corr( &nob );
    double norm_res = norm_mf+norm_corr;

    density_ob3 dob3(&nuc,true,true,true,norm_res,10); // nuc, central,tensor,isospin,norm,qmax (default=7)
    /* qmax is the maximum value of q in the sum in Eq. D.37 in Maartens thesis
     * note that this equation is incorrect/incomplete, see manual
     */
    double mf,corr;
    int t = isospin; // 1 for proton, -1 for neutron, 0 for both
    dob3.write(".",name.c_str(),-1,-1,-1,-1, t,&mf,&corr); // outputdir, outputname, nA,lA,nB,lB,t,mean field integral,corr integral
}

int main(int argc,char* argv[]){
    if (argc!=5){
        std::cerr << "[Error] expected four arguments:" << std::endl;
        std::cerr << "[executable] [A] [Z] [nucleusname] [proton/neutron/all]"<< std::endl;
        std::cerr << std::endl;
        exit(-1);
    }
    ob( atoi(argv[1]), atoi(argv[2]), argv[3], atoi(argv[4]));
}
