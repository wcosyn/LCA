#include "density_ob3.h"
#include "nucleusall.h"
#include "nucleusnp.h"
#include "nucleuspp.h"
#include "norm_ob.h"
#include <string>

void ob(int A,int Z,std::string name){
    Nucleusall nuc(".",".",A,Z);   // inputdir,resultdir,A,Z
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

int main(){
    /* the following block for He,Be,C,O,Al took:
     * 00.111.-1-1-1-1 : real=7m6s user=26m17s
     */
    ob(4,2,"He");
    ob(9,4,"Be");
    ob(12,6,"C");
    ob(16,8,"O");
    ob(27,13,"Al");
    return 0;
}
