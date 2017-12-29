#include "isolinkstrength.h"

Isolinkstrength::Isolinkstrength():
pplink(0.),nnlink(0.),nplink(0.){

}

Isolinkstrength::Isolinkstrength(const double pp, const double nn, const double np):
    pplink(pp),nnlink(nn),nplink(np){

}


Isolinkstrength::~Isolinkstrength(){

}

Isolinkstrength Isolinkstrength::operator+(const Isolinkstrength& rhs) const{
    Isolinkstrength result;
    result.pplink=this->pplink+rhs.pplink;
    result.nnlink=this->nnlink+rhs.nnlink;
    result.nplink=this->pplink+rhs.nplink;

    return result;
}

Isolinkstrength Isolinkstrength::operator*(const double sc) const{
    Isolinkstrength result;
    result.pplink=this->pplink*sc;
    result.nnlink=this->nnlink*sc;
    result.nplink=this->nplink*sc;    

    return result;
}

Isolinkstrength& Isolinkstrength::operator+=(const Isolinkstrength &rhs){
    this->pplink+=rhs.pplink;
    this->nnlink+=rhs.nnlink;
    this->nplink+=rhs.nplink;

    return *this;
}