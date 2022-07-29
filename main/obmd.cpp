#include "density_ob3.h"
#include "density_iso_ob3.h"
#include "nucleusall.h"
#include "nucleusnp.h"
#include "nucleuspp.h"
#include "nucleusnn.h"
#include "nucleus_iso.h"
#include "norm_ob.h"
#include "norm_iso_ob.h"
#include <string>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include "threej.h"




void ob(int A,int Z,std::string name, int isospin, Nucleus *nuc){


    //NucleusPP  nuc("../data/mosh","../data/mosh",A,Z); // idem
    // NucleusNP  nuc("../data/mosh","../data/mosh",A,Z); // idem
    // NucleusNN  nuc("../data/mosh","../data/mosh",A,Z); // idem

    Nucleusall nucall("../data/mosh","../data/mosh",A,Z);
    norm_ob no(&nucall);
    norm_ob::norm_ob_params nob= {-1, -1, -1, -1, 0}; // nA,lA,nB,lB,t
    double norm_mf  = no.sum_me_pairs( &nob );
    double norm_corr= no.sum_me_corr( &nob );
    double norm_all = norm_mf+norm_corr;

    double me_sum = 0.;
    for (int p=0;p<nucall.get_number_of_pairs();p++){
        Pair* pair = nucall.getPair(p);
        double me = no.get_me(pair, (void*) &nob);
        printf("[Info]: me for pair [% 5d] ",p);
        printf("(n,l,j,mj,t) = (% d,% d,% d,% d,% d)",pair->getn1(),pair->getl1(),pair->gettwo_j1(),pair->gettwo_mj1(),pair->gettwo_t1());
        printf(",(% d,% d,% d,% d,% d)  ",pair->getn2(),pair->getl2(),pair->gettwo_j2(),pair->gettwo_mj2(),pair->gettwo_t2());
        printf(" norm : %.2f  ",pair->getfnorm());
        double normP = pair->getfnorm(); 
        printf("ME is %f \n",me);
        me_sum+= me*normP;
    }
    printf("\n[Info]: ME sum is %f\n",me_sum);
    density_ob3 dob3(nuc,true,true,true,norm_all,10); // nuc, central,tensor,nucl,norm,qmax (default=7)
    /* qmax is the maximum value of q in the sum in Eq. D.37 in Maartens thesis
     * note that this equation is incorrect/incomplete, see manual
     */
    double mf,corr;
    dob3.write(".",name.c_str(),-1,-1,-1,-1, isospin,&mf,&corr); // outputdir, outputname, nA,lA,nB,lB,t,mean field integral,corr integral
}

void ob(std::string name, std::string outputdir ,NucleusIso *nuc, const bool hard, double a, double b,double c){


    //NucleusPP  nuc("../data/mosh","../data/mosh",A,Z); // idem
    // NucleusNP  nuc("../data/mosh","../data/mosh",A,Z); // idem
    // NucleusNN  nuc("../data/mosh","../data/mosh",A,Z); // idem

    // NucleusIso nucall("../data/mosh","../data/mosh",A,Z);
    int N=nuc->getA()-nuc->getZ();
    norm_iso_ob no(nuc,IsoMatrixElement((nuc->getZ()>1? double(nuc->getZ())*(nuc->getZ()-1)/(nuc->getA()*(nuc->getA()-1)):1),(N>1? double(N)*(N-1)/(nuc->getA()*(nuc->getA()-1)):1),
            double(N)*nuc->getZ()/(nuc->getA()*(nuc->getA()-1)),double(N)*nuc->getZ()/(nuc->getA()*(nuc->getA()-1))),a,b,c,hard,true,true,true);
    norm_iso_ob::norm_ob_params nob= {-1, -1, -1, -1}; // nA,lA,nB,lB,t
    IsoMatrixElement norm_mf  = no.sum_me_coefs( &nob );
    IsoMatrixElement norm_corr= no.sum_me_corr( &nob );
    IsoMatrixElement norm_all = norm_mf+norm_corr;

    density_iso_ob3 dob3(nuc,norm_all,a,b,c,hard,true,true,true,10); // nuc, central,tensor,nucl,norm,qmax (default=7)
    /* qmax is the maximum value of q in the sum in Eq. D.37 in Maartens thesis
     * note that this equation is incorrect/incomplete, see manual
     */
    double mf,corr;
    dob3.write(outputdir,name,mf,corr,-1,-1,-1,-1,hard); // outputdir, outputname, nA,lA,nB,lB,t,mean field integral,corr integral
}

int main(int argc,char* argv[]){
    if (argc<8){
        std::cerr << "[Error] expected six arguments: " << argc << std::endl;
        std::cerr << "[executable] [A] [Z] [nucleusname] [p/n/all] [pp/nn/np/pn/all/iso] [hard=1/soft=0]"<< std::endl;
        std::cerr << argv[1] << std::endl;
        std::cerr << std::endl;
        exit(-1);
    }
    std::string isospin = argv[4];
    bool hard = atoi(argv[6]);
    int nucl=0;
    if(!isospin.compare("p")) { nucl=1;}
    else if(!isospin.compare("n")) { nucl=-1;}
    else if(!isospin.compare("all")) {nucl=0;}
    else {std::cerr << "Invalid fourth arguments (isospin): select either p, n or all " << isospin << std::endl; exit(-1); assert(1==0);} 

    std::string pairs = argv[5];
    int A=atoi(argv[1]);
    int Z=atoi(argv[2]);
    std::string nucl_name=argv[3];
    std::string outputdir=argv[7];
    double a=atof(argv[8]);
    double b = atof(argv[9]);
    double c = atof(argv[10]);

    if(!pairs.compare("all")&&!isospin.compare("all")){
        NucleusPP  nucpp("../data/mosh","../data/mosh",A,Z); // idem
        NucleusNP  nucnp("../data/mosh","../data/mosh",A,Z); // idem
        NucleusNN  nucnn("../data/mosh","../data/mosh",A,Z); // idem

        Nucleusall nucall("../data/mosh","../data/mosh",A,Z);

        norm_ob no(&nucall);
        norm_ob::norm_ob_params nob= {-1, -1, -1, -1, 0}; // nA,lA,nB,lB,t
        double norm_mf  = no.sum_me_pairs( &nob );
        double norm_corr= no.sum_me_corr( &nob );
        std::cout << "all all " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
        nob.t=-1;
        norm_mf  = no.sum_me_pairs( &nob );
        norm_corr= no.sum_me_corr( &nob );
        std::cout << "all n " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
        nob.t=1;
        norm_mf  = no.sum_me_pairs( &nob );
        norm_corr= no.sum_me_corr( &nob );
        std::cout << "all p " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
        
        norm_ob nopp(&nucpp);
        nob.t=0;
        norm_mf  = nopp.sum_me_pairs( &nob );
        norm_corr= nopp.sum_me_corr( &nob );
        std::cout << "pp all " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
        nob.t=-1;
        norm_mf  = nopp.sum_me_pairs( &nob );
        norm_corr= nopp.sum_me_corr( &nob );
        std::cout << "pp n " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
        nob.t=1;
        norm_mf  = nopp.sum_me_pairs( &nob );
        norm_corr= nopp.sum_me_corr( &nob );
        std::cout << "pp p " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;

        norm_ob nonn(&nucnn);
        nob.t=0;
        norm_mf  = nonn.sum_me_pairs( &nob );
        norm_corr= nonn.sum_me_corr( &nob );
        std::cout << "nn all " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
        nob.t=-1;
        norm_mf  = nonn.sum_me_pairs( &nob );
        norm_corr= nonn.sum_me_corr( &nob );
        std::cout << "nn n " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
        nob.t=1;
        norm_mf  = nonn.sum_me_pairs( &nob );
        norm_corr= nonn.sum_me_corr( &nob );
        std::cout << "nn p " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;

        norm_ob nonp(&nucnp);
        nob.t=0;
        norm_mf  = nonp.sum_me_pairs( &nob );
        norm_corr= nonp.sum_me_corr( &nob );
        std::cout << "np all " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
        nob.t=-1;
        norm_mf  = nonp.sum_me_pairs( &nob );
        norm_corr= nonp.sum_me_corr( &nob );
        std::cout << "np n " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
        nob.t=1;
        norm_mf  = nonp.sum_me_pairs( &nob );
        norm_corr= nonp.sum_me_corr( &nob );
        std::cout << "np p " << norm_mf << " " << norm_corr << " " << norm_mf+norm_corr <<std::endl;
    }

    // the inputdir and resultdir are only used for storage of the moshinskybrackets, always, everywhere
    if(!pairs.compare("all")){
        Nucleusall nuc("../data/mosh","../data/mosh",A,Z);   // inputdir,resultdir,A,Z
        ob( A, Z, nucl_name,nucl, &nuc);
    }
    else if(!pairs.compare("pp")){
        NucleusPP nuc("../data/mosh","../data/mosh",A,Z);   // inputdir,resultdir,A,Z
        ob( A, Z, nucl_name,nucl, &nuc);
    }
    else if(!pairs.compare("nn")){
        NucleusNN nuc("../data/mosh","../data/mosh",A,Z);   // inputdir,resultdir,A,Z
        ob( A, Z, nucl_name,nucl, &nuc);
    }
    else if(!pairs.compare("np")||!pairs.compare("pn")){
        NucleusNP nuc("../data/mosh","../data/mosh",A,Z);   // inputdir,resultdir,A,Z
        ob( A, Z, nucl_name,nucl, &nuc);
    }
    else if(!pairs.compare("iso")){
        NucleusIso nuc("../data/mosh","../data/mosh",A,Z);   // inputdir,resultdir,A,Z
        ob(nucl_name,outputdir,&nuc,hard,a,b,c);
    }
    else {std::cerr << "Invalid fifth argument (pairs): select either pp, nn, pn(=np) or all or iso " << pairs << std::endl; exit(-1); assert(1==0);} 
}
