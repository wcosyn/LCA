#include "nucleusall.h"
#include "norm_tb.h"
#include "density_rel.h"
#include "wavefunctionp.h"
#include <string>
#include <iostream>

void tbmd_rel(int A,int Z,std::string name )
{
    WavefunctionP::mapwfp.setA( A );
    WavefunctionP::mapwfp.setpstep( 0.05 );
    WavefunctionP::mapwfcentralp_Hard.setA( A );
    WavefunctionP::mapwfcentralp_Hard.setpstep( 0.05 );
    WavefunctionP::mapwfcentralp_VMC.setA( A );
    WavefunctionP::mapwfcentralp_VMC.setpstep( 0.05 );
    WavefunctionP::mapwftensorp.setA( A );
    WavefunctionP::mapwftensorp.setpstep( 0.05 );
    WavefunctionP::mapwfspinisospinp.setA( A );
    WavefunctionP::mapwfspinisospinp.setpstep( 0.05 );


    // the inputdir and resultdir are only used for storage of the moshinskybrackets, always, everywhere
    Nucleusall nucleus("../data/mosh","../data/mosh",A,Z);   // inputdir,resultdir,A,Z

    int nA=-1;
    int lA=-1;
    int nB=-1;
    int lB=-1;
    int S =-1;
    int T =-1;
    bool central = true;
    bool tensor  = true;
    bool isospin = true;
    double norm=1;

    norm_tb::norm_tb_params ntp = {nA, lA, nB, lB, S, T};
    norm_tb* nt = new norm_tb( &nucleus, central, tensor, isospin );
    double norm_mf_sel = nt->sum_me( &ntp )/norm*A*(A-1)*0.5;
    double norm_2b_sel = nt->sum_me_corr_coefs( &ntp )/norm*A*(A-1)*0.5;
    double norm_3b_sel = nt->sum_me_3b_corr_coefs( &ntp )/norm*A*(A-1)*0.5;
//  double norm_3b_sel2 = nt->sum_me_3b_corr( &ntp )/norm*A*(A-1)*0.5;

    /*
    cout << "--------------------------------------------------" << endl;
    cout << norm_3b_sel << "\t" << norm_3b_sel2 << endl;
    cout << "--------------------------------------------------" << endl;
    */

    std::cout << "[tbmd_rel] selected normalisation mf\t" << norm_mf_sel << "\t" << norm_mf_sel*norm << std::endl;
    std::cout << "[tbmd_rel] selected normalisation 2b corr\t" << norm_2b_sel << std::endl;
    std::cout << "[tbmd_rel] selected normalisation 3b corr\t" << norm_3b_sel << std::endl;
    std::cout << "[tbmd_rel] total \t" << norm_mf_sel+ norm_2b_sel+ norm_3b_sel << std::endl;

    delete nt;
//  return;

    double intmf, int2b, int3b;

    density_rel rel( &nucleus, central, tensor, isospin, norm, 10);
    rel.write(".", name.c_str(), nA, lA, nB, lB, S, T, &intmf, &int2b, &int3b );

    std::cout << "[tbmd_rel]   \tnorm\t\tintegral\tintegral/norm" << std::endl;
    std::cout << "[tbmd_rel] mf\t" << norm_mf_sel << "\t\t" << intmf << "\t" << intmf/norm_mf_sel << std::endl;
    std::cout << "[tbmd_rel] 2b\t" << norm_2b_sel << "\t" << int2b << "\t" << int2b/norm_2b_sel << std::endl;
    std::cout << "[tbmd_rel] 3b\t" << norm_3b_sel << "\t" << int3b << "\t" << int3b/norm_3b_sel << std::endl;
}

int main(int argc,char* argv[]){
    if (argc!=4){
        std::cerr << "[Error] expected three arguments:" << std::endl;
        std::cerr << "[executable] [A] [Z] [nucleusname]"<< std::endl;
        std::cerr << std::endl;
        exit(-1);
    }
    tbmd_rel( atoi(argv[1]), atoi(argv[2]), argv[3]);
    return 0;
}
