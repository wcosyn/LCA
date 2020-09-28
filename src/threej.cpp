#include "threej.h"

#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_math.h>
#include <cmath>
#include <utility> // for swap method


threej threej::threejs; // static object threejs
threej::threej_logger threej::my3jlogger;
std::unordered_map< threejobj, double> threej::threejsmap;

double threej::get(int two_j1,int two_j2, int two_j3, int two_m1,int two_m2, int two_m3){
    if ( triangle_selection_fails(two_j1, two_j2, two_j3)
      || m_selection_fails(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3))
        return 0.;
    threejobj c(two_j1,two_j2,two_j3,two_m1,two_m2,two_m3);
    std::unordered_map<threejobj,double>::const_iterator it = threejsmap.find(c);
    if (it != threejsmap.end()){ // key is in the map
/*#pragma omp critical(getthreej) // prevent concurrent writes to the shared threejsmap!
        {
        my3jlogger.add(c);
        }*/
        return it->second;
    } else { // key is not in map
        double t = gsl_sf_coupling_3j( two_j1, two_j2, two_j3, two_m1, two_m2, two_m3);
#pragma omp critical(getthreej) // prevent concurrent writes to the shared threejsmap!
        {
        threejsmap[c] = t;
//        my3jlogger.add(c);
        }
        return t;
    }
}

// this method uses less memory but takes considerably longer
double threej::get_sorted(int two_j1,int two_j2, int two_j3, int two_m1,int two_m2, int two_m3){
    if ( triangle_selection_fails(two_j1, two_j2, two_j3)
      || m_selection_fails(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3))
        return 0.;
    int s = 0b00;

    if (two_j1 > two_j2){ std::swap(two_j1,two_j2); std::swap(two_m1,two_m2); s^= 0b01; }
    if (two_j1 > two_j3){ std::swap(two_j1,two_j3); std::swap(two_m1,two_m3); s^= 0b01; }
    if (two_j2 > two_j3){ std::swap(two_j2,two_j3); std::swap(two_m2,two_m3); s^= 0b01; }

    const int swapphase = s ? (((two_j1 + two_j2 + two_j3) & 0b10) ? -1 : 1) : 1;
    threejobj c(two_j1,two_j2,two_j3,two_m1,two_m2,two_m3);
    std::unordered_map<threejobj,double>::const_iterator it = threejsmap.find(c);
    if (it != threejsmap.end()){ // key is in the map
        return swapphase*(it->second);
    } else { // key is not in map
        double t = gsl_sf_coupling_3j( two_j1, two_j2, two_j3, two_m1, two_m2, two_m3);
#pragma omp critical(getthreej) // prevent concurrent writes to the shared threejsmap!
        {
        threejsmap[c] = t;
        }
        return swapphase*t;
    }
}

int threej::triangle_selection_fails(int two_ja, int two_jb, int two_jc)
{
    /*
     * enough to check the triangle condition for one spin vs. the other two
     */
    return ( (two_jb < std::abs(two_ja - two_jc)) || (two_jb > two_ja + two_jc) ||
             GSL_IS_ODD(two_ja + two_jb + two_jc) );
}

int threej::m_selection_fails(int two_ja, int two_jb, int two_jc,
                              int two_ma, int two_mb, int two_mc)
{
    return (
                  std::abs(two_ma) > two_ja
               || std::abs(two_mb) > two_jb
               || std::abs(two_mc) > two_jc
               || GSL_IS_ODD(two_ja + two_ma)
               || GSL_IS_ODD(two_jb + two_mb)
               || GSL_IS_ODD(two_jc + two_mc)
               || (two_ma + two_mb + two_mc)
           );
}
