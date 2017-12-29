#ifndef ISO_MATRIX_ELEMENT_H
#define ISO_MATRIX_ELEMENT_H

/**
 * @brief small structure that is used as return value for matrix element calculations.  
 * Holds the results for pp,nn,np(tagged p) and np(tagged n) pairs separately
 * 
 */
class IsoMatrixElement {
    public:
    double pp_res, nn_res, np_p_res, np_n_res;

    IsoMatrixElement();
    IsoMatrixElement(const double pp, const double nn, const double np_p, const double np_n);
    ~IsoMatrixElement();
    IsoMatrixElement operator+(const IsoMatrixElement& rhs) const;
    IsoMatrixElement operator*(const double sc) const;
    IsoMatrixElement& operator+=(const IsoMatrixElement& rhs);

    double norm() const;

    double getValue(const int i) const;
};

#endif