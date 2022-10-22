#include "mfShellobmd.h"
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

#include <iostream>
using std::cout;
using std::endl;




void SpecMF(int A,int Z,std::string name, int isospin, Nucleus *nuc,int sh, int ns, int nj){


    //NucleusPP  nuc("../data/mosh","../data/mosh",A,Z); // idem
    // NucleusNP  nuc("../data/mosh","../data/mosh",A,Z); // idem
    // NucleusNN  nuc("../data/mosh","../data/mosh",A,Z); // idem
    WavefunctionP::mapwfp.setA( A, 45, 25);
    WavefunctionP::mapwfp.setpstep( 0.05 );
    Nucleusall nucall("../data/mosh","../data/mosh",A,Z);
    norm_ob no(&nucall);
    norm_ob::norm_ob_params nob= {-1, -1, -1, -1, 0}; // nA,lA,nB,lB,t
    double norm_mf  = no.sum_me_pairs( &nob,sh,ns,nj);
    //no.sum_me_pairs(&nob);
    double me_sum = 0.;
    /*
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
    */
    mfShellobmd mfShellobmd(nuc,true,true,true,norm_mf,10); // nuc, central,tensor,nucl,norm,qmax (default=7)
    /* qmax is the maximum value of q in the sum in Eq. D.37 in Maartens thesis
     * note that this equation is incorrect/incomplete, see manual
     */
    double mf;
    mfShellobmd.write("../../spyder/CosynResults",name.c_str(),-1,-1,-1,-1, isospin,&mf,sh,ns,nj); // outputdir, outputname, nA,lA,nB,lB,t,mean field integral,corr integral
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
    const char* nucl_name=argv[3];
    std::string outputdir=argv[7];
    double a=atof(argv[8]);
    double b = atof(argv[9]);
    double c = atof(argv[10]);
    int sh = atof(argv[11]);
    int ns = atof(argv[12]);
    int nj = atof(argv[13]);

    if(!pairs.compare("all")&&!isospin.compare("all")){
        NucleusPP  nucpp("../data/mosh","../data/mosh",A,Z); // idem
        NucleusNP  nucnp("../data/mosh","../data/mosh",A,Z); // idem
        NucleusNN  nucnn("../data/mosh","../data/mosh",A,Z); // idem

        Nucleusall nucall("../data/mosh","../data/mosh",A,Z);

        norm_ob no(&nucall);
        norm_ob::norm_ob_params nob= {-1, -1, -1, -1, 0}; // nA,lA,nB,lB,t
        double norm_mf  = no.sum_me_pairs( &nob,sh,ns,nj );
        std::cout << "all all " << norm_mf  <<std::endl;
        nob.t=-1;
        norm_mf  = no.sum_me_pairs( &nob,sh,ns,nj );
        std::cout << "all n " << norm_mf  <<std::endl;
        nob.t=1;
        norm_mf  = no.sum_me_pairs( &nob,sh,ns,nj );
        std::cout << "all p " << norm_mf  <<std::endl;
        
        norm_ob nopp(&nucpp);
        nob.t=0;
        norm_mf  = nopp.sum_me_pairs1( &nob );
        std::cout << "pp all " << norm_mf  <<std::endl;
        nob.t=-1;
        norm_mf  = nopp.sum_me_pairs1( &nob );
        std::cout << "pp n " << norm_mf <<std::endl;
        nob.t=1;
        norm_mf  = nopp.sum_me_pairs1( &nob );
        std::cout << "pp p " << norm_mf <<std::endl;

        norm_ob nonn(&nucnn);
        nob.t=0;
        norm_mf  = nonn.sum_me_pairs1( &nob );
        std::cout << "nn all " << norm_mf <<std::endl;
        nob.t=-1;
        norm_mf  = nonn.sum_me_pairs1( &nob );
        std::cout << "nn n " << norm_mf <<std::endl;
        nob.t=1;
        norm_mf  = nonn.sum_me_pairs1( &nob );
        std::cout << "nn p " << norm_mf <<std::endl;

        norm_ob nonp(&nucnp);
        nob.t=0;
        norm_mf  = nonp.sum_me_pairs1( &nob );
        std::cout << "np all " << norm_mf <<std::endl;
        nob.t=-1;
        norm_mf  = nonp.sum_me_pairs1( &nob );
        std::cout << "np n " << norm_mf <<std::endl;
        nob.t=1;
        norm_mf  = nonp.sum_me_pairs1( &nob );
        std::cout << "np p " << norm_mf <<std::endl;
    }

    // the inputdir and resultdir are only used for storage of the moshinskybrackets, always, everywhere
    if(!pairs.compare("all")){
        Nucleusall nuc("../data/mosh","../data/mosh",A,Z);   // inputdir,resultdir,A,Z
        SpecMF( A, Z, nucl_name,nucl, &nuc,sh,ns,nj);
    }
    else if(!pairs.compare("pp")){
        NucleusPP nuc("../data/mosh","../data/mosh",A,Z);   // inputdir,resultdir,A,Z
        SpecMF( A, Z, nucl_name,nucl, &nuc,sh,ns,nj);
    }
    else if(!pairs.compare("nn")){
        NucleusNN nuc("../data/mosh","../data/mosh",A,Z);   // inputdir,resultdir,A,Z
        SpecMF( A, Z, nucl_name,nucl, &nuc,sh,ns,nj);
    }
    else if(!pairs.compare("np")||!pairs.compare("pn")){
        NucleusNP nuc("../data/mosh","../data/mosh",A,Z);   // inputdir,resultdir,A,Z
        SpecMF( A, Z, nucl_name,nucl, &nuc,sh,ns,nj);
    }
    else {std::cerr << "Invalid fifth argument (pairs): select either pp, nn, pn(=np) or all or iso " << pairs << std::endl; exit(-1); assert(1==0);} 
    
}
