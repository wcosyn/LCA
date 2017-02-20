#ifndef THREEJ_H
#define THREEJ_H


#include <unordered_map>

/*
 * Class that calculates and saves 3j- coefficients.
 * This is faster than calculating the same one over and over again
 * Instead of using advanced memory efficient scheme we use a simple
 * unordered (hash) map.
 * This appear to be faster as index calculations taking into account
 * all kinds of symmetry are expensive. The 3j's are called a lot!
 */

struct threejobj{
    threejobj(int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3) :
        two_j1(two_j1), two_j2(two_j2) , two_j3(two_j3), two_m1(two_m1) , two_m2(two_m2) , two_m3(two_m3){}
    bool operator==(const threejobj& rhs) const {
        return (two_j1 == rhs.two_j1 && two_j2 == rhs.two_j2 && two_j3 == rhs.two_j3 && two_m1 == rhs.two_m1 && two_m2 == rhs.two_m2 && two_m3 == rhs.two_m3);
    }
    int two_j1,two_j2,two_j3,two_m1,two_m2,two_m3;
};

namespace std {
    template<>
    struct hash<threejobj>{
    public:
        size_t operator()(const threejobj &t) const{
            //return (t.two_j1 << 24) + (t.two_j2 << 18) + (t.two_j3 << 12) + (t.two_m1 << 6) + t.two_m2;
            return (t.two_j1 << 24) ^ (t.two_j2 << 18) ^ (t.two_j3 << 12) ^ (t.two_m1 << 6) ^ t.two_m2; // this is considerably faster (~15% for Fe)
        }
    };
}
class threej
{
public:
    double get( int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3 );
    double get_sorted( int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3 );
    static threej threejs;
    static std::unordered_map< threejobj, double> threejsmap;

    int triangle_selection_fails(int two_ja, int two_jb, int two_jc);
    int m_selection_fails(int two_ja, int two_jb, int two_jc,
                          int two_ma, int two_mb, int two_mc);
    /** if you want to investigate what 3j's get called and
     *  how many times, you can use something like the logger
     *  below. Remember omp critical block when writing to the map (add)
     */
    /*
    struct threej_logger {
        std::unordered_map< threejobj, unsigned int> threejmap;
        void add(int two_j1,int two_j2,int two_j3,int two_m1, int two_m2,int two_m3){
            threejobj c(two_j1,two_j2,two_j3,two_m1,two_m2,two_m3);
            threejmap[c]++;
        }
    };
    static threej_logger my3jlogger;
    */
};


#endif // THREE_J

