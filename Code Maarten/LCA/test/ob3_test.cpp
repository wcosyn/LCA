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


void ob(int A,int Z,std::string name){

    // the inputdir and resultdir are only used for storage of the moshinsky brakets, always, everywhere
    Nucleusall nuc("../data/mosh","../data/mosh",A,Z);   // inputdir,resultdir,A,Z
    //NucleusPP  nucpp(".",".",40,20); // idem
    //NucleusNP  nucnp(".",".",40,20); // idem


    norm_ob no(&nuc);
    norm_ob::norm_ob_params nob= {-1, -1, -1, -1, 0};
    double norm_mf  = no.sum_me_pairs( &nob );
    double norm_corr= no.sum_me_corr( &nob );
    double norm_res = norm_mf+norm_corr;

    density_ob3 dob3(&nuc,true,true,true,norm_res,10); // nuc, central,tensor,isospin,norm,qmax (default=7)
    /* qmax is the maximum value of q in the sum in Eq. D.37 in Maartens thesis
     * note that this equation is incorrect/incomplete, see manual
     */
    double mf,corr;
    int t = 0; // 1 for proton, -1 for neutron, 0 for both
    dob3.write(".",name.c_str(),-1,-1,-1,-1, t,&mf,&corr); // outputdir, outputname, nA,lA,nB,lB,t,mean field integral,corr integral
}

int main(int argc,char* argv[]){
    /* the following block for He,Be,C,O,Al took:
     * 00.111.-1-1-1-1 : real=7m6s user=26m17s :: original version
     * 00.111.-1-1-1-1 : real=7m8s user=26m17s :: with the new density_ob3::get_me12_factor implemented
     */
    if (argc!=4){
        std::cerr << "[Error] expected three arguments:" << std::endl;
        std::cerr << "[executable] [A] [Z] [nucleusname]"<< std::endl;
        std::cerr << std::endl;
        exit(-1);
    }
    ob( atoi(argv[1]), atoi(argv[2]), argv[3]);

    std::unordered_map< threejobj, double > n = threej::threejsmap;
    unsigned int b = n.bucket_count();
    unsigned int maxb = 0;
    std::cout << "[3jmap] threejsmap has " << b << " buckets" << std::endl;
    for (unsigned int i=0; i<b;i++){
        std::cout << "[3jmap] bucket #" << i << " has " << n.bucket_size(i) << " elements." << std::endl;
        if (n.bucket_size(i) > maxb)
            maxb = n.bucket_size(i);
    }
    std::cout << "[3jmap] map size is       : " << n.size() << std::endl;
    std::cout << "[3jmap] max bucket size is: " << maxb << std::endl;
    std::cout << "[3jmap] map load factor is: " << n.load_factor() << std::endl;

    /*
     * if threejlogging is enabled you can use this cout to inspect the logger
     * */
    /*
    std::unordered_map< threejobj, unsigned int > m = threej::my3jlogger.threejmap;
    unsigned int i=0;
    for( std::unordered_map< threejobj, unsigned int >::iterator it = m.begin(); it != m.end(); ++it){
        std::cout << "[3jlogger]"<< "[" << i << "] ";
        std::cout << " :: " << it->second << " : ";
        std::cout << " (";
        std::cout << it->first.two_j1 << ", ";
        std::cout << it->first.two_j2 << ", ";
        std::cout << it->first.two_j3 << ", ";
        std::cout << it->first.two_m1 << ", ";
        std::cout << it->first.two_m2 << ", ";
        std::cout << it->first.two_m3 << ")" << std::endl;
        i++;
    }*/
}
