#ifndef ISO_MATRIX_ELEMENT_H
#define ISO_MATRIX_ELEMENT_H

/**
 * @brief small class that is used as return value for matrix element calculations of one-body operators.  
 * Holds the results for pp,nn,np(tagged p) and np(tagged n) pairs separately
 * 
 */
class IsoMatrixElement {
    public:
    double pp_res, nn_res, np_p_res, np_n_res;

    /**
     * @brief Default constructor
     * 
     */
    IsoMatrixElement();
    /**
     * @brief Constructor
     * 
     * @param pp value for pp pairs
     * @param nn value for nn pairs
     * @param np_p value for np pairs tagged on proton
     * @param np_n value for np pairs tagged on neutron
     */
    IsoMatrixElement(const double pp, const double nn, const double np_p, const double np_n);
    /**
     * @brief copy constructor
     * 
     * @param copy 
     */
    IsoMatrixElement(const IsoMatrixElement& copy);
    /**
     * @brief Destructor
     * 
     */
    ~IsoMatrixElement();
    /**
     * @brief adding operator
     * 
     * @param rhs what we want to add
     * @return IsoMatrixElement sum of this + argument
     */
    IsoMatrixElement operator+(const IsoMatrixElement& rhs) const;
    /**
     * @brief right-multiply all elements with a scalar
     * 
     * @param sc what we multiply with
     * @return IsoMatrixElement result of the multiplication
     */
    IsoMatrixElement operator*(const double sc) const;
    /**
     * @brief add another matrixelement to this one
     * 
     * @param rhs what we want to add
     * @return IsoMatrixElement& result of the addition
     */
    IsoMatrixElement& operator+=(const IsoMatrixElement& rhs);

    /**
     * @brief sum of all the individual isospin contributions
     * 
     * @return double norm 
     */
    double norm() const;

    /**
     * @brief returns the different possible contributions
     * 
     * @param i 
     * 0: pp pairs contribution
     * 1: nn pairs contribution
     * 2: np(p) pairs contribution
     * 3: np(n) pairs contribution
     * 4: proton one-body contribution (=0+2)
     * 5: neutron one-body contribution (=1+3)
     * 6: total one-body contribution (=0+1+2+3)
     * @return double Matrix element value
     */
    double getValue(const int i) const;
};

#endif