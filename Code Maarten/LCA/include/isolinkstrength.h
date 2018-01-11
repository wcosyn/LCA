#ifndef ISO_LINK_STRENGTH_H
#define ISO_LINK_STRENGTH_H


/**
 * @brief small class that holds three doubles, holds the linkstrength for the isopaircoefs for pp,nn and np pairs.
 * 
 */
class Isolinkstrength{
    public:
        double pplink, nnlink, nplink; //link strengths, individual pair contributions
        /**
         * @brief default constructor
         * 
         */
        Isolinkstrength();
        /**
         * @brief constructor
         * 
         * @param pp part of link strength originating from pp pairs
         * @param nn part of link strength originating from nn pairs
         * @param np part of link strength originating from np pairs
         */
        Isolinkstrength(const double pp, const double nn, const double np);
        /**
         * @brief Destructor
         * 
         */
        ~Isolinkstrength();
        /**
         * @brief add two linkstrength objects
         * 
         * @param rhs what we want to add
         * @return Isolinkstrength this + rhs
         */
        Isolinkstrength operator+(const Isolinkstrength& rhs) const;
        /**
         * @brief multiply all linkstrengths with a scalar
         * 
         * @param sc scalar
         * @return Isolinkstrength result of multiplication
         */
        Isolinkstrength operator*(const double sc) const;
        /**
         * @brief add and assign
         * 
         * @param rhs what we want to add
         * @return Isolinkstrength& result of the addition
         */
        Isolinkstrength& operator+=(const Isolinkstrength& rhs);
        double get_pplink() const{ return pplink;};///< returns part of link strength originating from pp pairs
        double get_nnlink() const{ return nnlink;};///< returns part of link strength originating from nn pairs
        double get_nplink() const{ return nplink;};///< returns part of link strength originating from np pairs
        
};

#endif //ISO_LINK_STRENGTH_H