#include "density_ob3.h"
#include "nucleusall.h"


int main(){
    Nucleusall nuc(".",".",12,6);
    density_ob3 dob3(&nuc,true,true,true,1.,10);
    double mf,corr;
    dob3.write(".","onmd.out",-1,-1,-1,-1,0,&mf,&corr);
    return 0;
}
