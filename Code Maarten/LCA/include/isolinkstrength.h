#ifndef ISO_LINK_STRENGTH_H
#define ISO_LINK_STRENGTH_H


/**
 * @brief small structure that holds three doubles, holds the linkstrength for the isopaircoefs for pp,nn and np pairs.
 * 
 */
class Isolinkstrength{
    public:
        double pplink, nnlink, nplink;
        Isolinkstrength();
        Isolinkstrength(const double pp, const double nn, const double np);
        ~Isolinkstrength();
        Isolinkstrength operator+(const Isolinkstrength& rhs) const;
        Isolinkstrength operator*(const double sc) const;
        Isolinkstrength& operator+=(const Isolinkstrength& rhs);
        double get_pplink() const{ return pplink;};
        double get_nnlink() const{ return nnlink;};
        double get_nplink() const{ return nplink;};
        
};

#endif //ISO_LINK_STRENGTH_H