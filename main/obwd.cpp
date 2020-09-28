#include "wigner_iso_ob3.h"
#include "nucleus_iso.h"
#include "norm_iso_ob.h"
#include <string>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include "threej.h"


void ob(std::string name, NucleusIso *nuc){

    int N=nuc->getA()-nuc->getZ();
    norm_iso_ob no(nuc,IsoMatrixElement(double(nuc->getZ())*(nuc->getZ()-1)/(nuc->getA()*(nuc->getA()-1)),double(N)*(N-1)/(nuc->getA()*(nuc->getA()-1)),
            double(N)*nuc->getZ()/(nuc->getA()*(nuc->getA()-1)),double(N)*nuc->getZ()/(nuc->getA()*(nuc->getA()-1))));
    norm_iso_ob::norm_ob_params nob= {-1, -1, -1, -1}; // nA,lA,nB,lB
    IsoMatrixElement norm_mf  = no.sum_me_coefs( &nob );
    IsoMatrixElement norm_corr= no.sum_me_corr( &nob );
    IsoMatrixElement norm_all = norm_mf+norm_corr;

    wigner_iso_ob3 dob3(nuc,norm_all,true,true,true,10); // nuc, central,tensor,nucl,norm,qmax (default=7)
    /* qmax is the maximum value of q in the sum in Eq. D.37 in Maartens thesis
     * note that this equation is incorrect/incomplete, see manual
     */
    double mf,corr;
    dob3.write(".",name,mf,corr,-1,-1,-1,-1); // outputdir, outputname, nA,lA,nB,lB,t,mean field integral,corr integral
}

int main(int argc,char* argv[]){
    // if (argc<6){
    //     std::cerr << "[Error] expected five arguments: " << argc << std::endl;
    //     std::cerr << "[executable] [A] [Z] [nucleusname] [p/n/all] [pp/nn/np/pn/all]"<< std::endl;
    //     std::cerr << argv[1] << std::endl;
    //     std::cerr << std::endl;
    //     exit(-1);
    // }

    int A=atoi(argv[1]);
    int Z=atoi(argv[2]);
    std::string nucl_name=argv[3];

    NucleusIso nuc("../data/mosh","../data/mosh",A,Z);   // inputdir,resultdir,A,Z
    ob(nucl_name,&nuc);
}
