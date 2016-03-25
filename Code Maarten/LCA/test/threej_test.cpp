#include "threej.h"
#include <gsl/gsl_sf_coupling.h>
#include <iostream>
#include <cmath>

class Threej_checker {
    /* remember, all j's and m's are actually two j and two m */
public:
    Threej_checker(int j1_max,int j2_max) :
        j1(0),j2(0),j3( abs(j2-j1)),
        m1(0),m2(0),m3( -j3),
        j1_max(j1_max),
        j2_max(j2_max),
        j3_max(j1_max+j2_max) {}
    void runcheck(){
        unsigned pass=0;
        unsigned rej=0;
        // increments of 1 in loop means increment of 1/2 in ang. mom
        for (j1=0;j1<=j1_max;j1++){
            for (j2=0;j2<=j2_max;j2++){
                for (j3=abs(j2-j1);j3<=j1+j2;j3++){
                    for (m1=-j1;m1<=j1;m1++){
                        for (m2=-j2;m2<=j2;m2++){
                            for (m3=-j3;m3<=j3;m3++){
                                check() ? pass++:rej++;
                            }
                        }
                    }
                }
            }
        }
        std::cout << "[CGC] check: " << pass << "/" <<  pass+rej << " passed, " << rej << " rejected cgc's " << std::endl;
    }
    bool check(){
        double cg1 = threej::threejs.get(j1,j2,j3,m1,m2,m3);
        double cg2 = gsl_sf_coupling_3j(j1,j2,j3,m1,m2,m3);
        bool pass = fabs(cg1-cg2) < 1e-10;
        if (!pass){
            std::cerr << "ERROR: MISMATCH ";
            std::cerr << " for quantum numbers (j1,j2,j3,m1,m2,m3): (";
            std::cerr << j1 << ", ";
            std::cerr << j2 << ", ";
            std::cerr << j3 << ", ";
            std::cerr << m1 << ", ";
            std::cerr << m2 << ", ";
            std::cerr << m3 << ") ";
            std::cerr << "\n[3J] threej::" << cg1 << ", GSL::" << cg2 << " diff : " << cg1-cg2 << std::endl;
        }
        return pass;
    }
    int j1,j2,j3,m1,m2,m3;
private:
    int j1_max,j2_max,j3_max;
};

void runcheck(){
    Threej_checker c(24,24);
    c.runcheck();
    // for j1 = 24, j2=24, j3 is in the range [0,48]
    unsigned nonzero=0;
    srand(1234);
    for (int i=0; i<(0x01 << 26);i++){
        int j = rand()%49; // j1 \in [0,48]
        //double wg3j =  threej::threejs.get(24,24,j,j-10,10,-j);
        double wg3j = gsl_sf_coupling_3j(24,24,j,j-10,10,-j);
        fabs(wg3j) > 1e-10? nonzero++ : 0;
    }
    std::cout << " calculated " << nonzero << " nonzero 3j's" << std::endl;
}

void testSignSwitch(){
    double cg1 = threej::threejs.get(2,2,4,0,2,-2);
    double cg2 = threej::threejs.get(2,2,4,0,-2,2);
    std::cout << "c1 : " << cg1 << " , cg2 : " << cg2 << std::endl;
}

int main(){
    //testSignSwitch();
    runcheck();
    return 0;
}
