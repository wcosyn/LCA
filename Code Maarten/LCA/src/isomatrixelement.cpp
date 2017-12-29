#include <iostream>
#include <cmath>

#include "isomatrixelement.h"

IsoMatrixElement::IsoMatrixElement():
pp_res(0.),nn_res(0.),np_p_res(0.),np_n_res(0.){
    
}

IsoMatrixElement::IsoMatrixElement(const double pp, const double nn, const double np_p, const double np_n):
pp_res(pp),nn_res(nn),np_p_res(np_p),np_n_res(np_n)
{
}

IsoMatrixElement::~IsoMatrixElement(){

}

IsoMatrixElement IsoMatrixElement::operator+(const IsoMatrixElement& rhs) const{
    IsoMatrixElement result;
    result.pp_res=this->pp_res+rhs.pp_res;
    result.nn_res=this->nn_res+rhs.nn_res;
    result.np_p_res=this->np_p_res+rhs.np_p_res;
    result.np_n_res=this->np_n_res+rhs.np_n_res;

    return result;
}

IsoMatrixElement IsoMatrixElement::operator*(const double sc) const{
    IsoMatrixElement result;
    result.pp_res=this->pp_res*sc;
    result.nn_res=this->nn_res*sc;
    result.np_p_res=this->np_p_res*sc;
    result.np_n_res=this->np_n_res*sc;

    return result;
}

IsoMatrixElement& IsoMatrixElement::operator+=(const IsoMatrixElement& rhs){
    this->pp_res+=rhs.pp_res;
    this->nn_res+=rhs.nn_res;
    this->np_p_res+=rhs.np_p_res;
    this->np_n_res+=rhs.np_n_res;

    return *this;
}

double IsoMatrixElement::norm() const{
    return fabs(pp_res)+fabs(nn_res)+fabs(np_p_res)+fabs(np_n_res);
}

double IsoMatrixElement::getValue(const int i) const{
    switch(i){
        case 0:
            return pp_res;
            break;
        case 1:
            return nn_res;
            break;
        case 2:
            return np_p_res;
            break;
        case 3:
            return np_n_res;
            break;
        case 4:
            return pp_res+np_p_res;
            break;
        case 5:
            return nn_res+np_n_res;
            break;
        case 6:
            return pp_res+nn_res+np_n_res+np_p_res;
            break;
        default:
            std::cerr << "Invalid choice in IsoMatrixElement::getValue" << std::endl;
            exit(1);
    
    }
}