#include "recmosh.h"
#include <iomanip>

using std::map;
using std::fabs;
using std::cout;
using std::endl;
using std::cerr;
using std::ofstream;
using std::ifstream;
using std::stringstream;

MapRecMosh RecMosh::maprecmosh = MapRecMosh(); //initialise static object (empty map)

RecMosh::RecMosh(const int n1, const int l1, const int n2, const int l2, const std::string& path)
    : n1(n1),
      l1(l1),
      n2(n2),
      l2(l2),
      coefficients(),
      path(path)

{
    if( n1*10+l1 <= n2*10+l2 )
        loadFile();

    // if n1*10+l1 > n2*10+l2 wordt in getCoefficient de
    // waarde uit n2l2 n1l2 genomen en de juiste voorfactor berekend

// used is 1, so when it is not used any more, it is still saved, 
//so when a next pair or nucleus needs it, it is still loaded, 
//otherwise it again loads the file from input, and updates written out to output are not used anymore.
// Warning, need to be sure that everything gets deleted and written out at the end !!
    used=1;
    updated= false;

}

RecMosh& RecMosh::createRecMosh(const int n1, const int l1, const int n2, const int l2, const std::string& path)
{
    RecMosh *r = &maprecmosh.get( n1, l1, n2, l2, path );
    r->use();
    return *r;

}

RecMosh::~RecMosh()
{
    if( updated ) {
        writeToFile();
    }
    cout << n1 << l1 << n2 << l2 << " recmosh removed " << endl;
}

void RecMosh::use()
{
    used++;
//  cout << "USE " << used << endl;
}

void RecMosh::remove()
{
    used--;
    if( used < 1 ) {
        cout << "test" << __FILE__ << __LINE__ << endl;
        // maprecmosh.remove( n1*1e3+ l1*1e2+ n2*10+ l2 );
        //delete this;
    }

}

void RecMosh::writeToFile()
{
    if( updated ) {
        stringstream name;
        name << path << "/recmosh" << n1 << l1 << n2 << l2 << ".dat";
        cout << "[recmosh][RecMosh::writeToFile] write to " << name.str().c_str() << endl;
        ofstream dataFile ( name.str().c_str(), ofstream::out );
        if( !dataFile.is_open() ) {
            cerr << "[recmosh][RecMosh::writeToFile] File could not be opened! Output written to stdout" << endl;
            dataFile.close();

            map<int, double>::iterator it;
            for( it =coefficients.begin(); it!= coefficients.end(); it++ ) {
                cout << (it)->first << " " << (it)->second << endl;
            }
        }
        dataFile << 1000*n1 + 100*l1 + 10*n2 + l2  << endl;

        map<int, double>::iterator it;
        for( it =coefficients.begin(); it!= coefficients.end(); it++ ) {
            dataFile << std::fixed      << std::setw(12) << (it)->first << "   ";
            dataFile << std::scientific << std::setw(12) << (it)->second << endl;
        }

        dataFile.close();
        updated = false;
    }
}

void RecMosh::loadFile(  )
{
    stringstream name;

    name << path << "/recmosh" << n1 << l1 << n2 << l2 << ".dat";
    ifstream dataFile ( name.str().c_str(), ifstream::in );
    if( !dataFile ) {
        cerr << "[recmosh][RecMosh::loadFile] File " << name.str() << " not found! " << endl;
        dataFile.close();
    } else {
        int data;
        dataFile >> data;
        while(!dataFile.eof()) {
            int key;
            double coefficient ;
            dataFile >> key;
            dataFile >> coefficient;
            coefficients[key] = coefficient;
        }
        dataFile.close();
    }
}

double RecMosh::getCoefficient( const int n, const int l, const int N, const int L, const int Lambda )
{
    if( 2*n+l+2*N+L != 2*n1+l1+2*n2+l2) {
        return 0;
    }
    //only coefficients with N L >= nl are computed and stored, so make use of permutation relations
    if( 10*n+l > 10*N+L ) {
        double val = maprecmosh.get( n1, l1, n2, l2, path ).getCoefficient( N, L, n, l, Lambda);
        if( (l1+Lambda)%2 ) return -1*val;
        else return val;
    }

    //
    // Only the files with n1l1 <= n2l2 are calculated and saved,
    // So the n2l2 n1l1 coeff is taken and correct phase factor added
    if( 10*n1+l1 > 10*n2+l2 ) {
        double val = maprecmosh.get( n2, l2, n1, l1, path ).getCoefficient( n, l, N, L, Lambda);
        if( (L+Lambda)%2 ) return -1*val;
        else return val;


    }

    int key = 100000*n + 10000*l + 1000*N + 100*L + Lambda*10;

    map< int, double >::const_iterator it;
    it = coefficients.find(key);
    if( it == coefficients.end()) {
        double result = calculate(n,l,N,L,n1,l1,n2,l2,Lambda);
        //if(abs(result)<1.E-10) result=0.;
        coefficients[key] = result;
        //if( fabs(result) < 1e-10 ) return 0;
        return result;
    } else {
        double result= it->second;
        //if( fabs( result ) < 1e-10 ) return 0;
        return result;
    }

}

/*
 * The moshinskybracket is calculated in a recursive way,
 * See thesis appendix for the expression
 */
double RecMosh::calculate( const int n, const int l, const int N, const int L, 
                            const int n1, const int l1, const int n2, const int l2, const int Lambda) const
{
    cout << "[recmosh][RecMosh::calculate] new recmosh value calculated " << endl;
// 	cout << "<(" << n << l << N << L << ")" << Lambda;
//	cout << "|(0" << l1 << 0 << l2 << ")" << Lambda << "> = " << endl;

    //these are all in principle checked before...
    // tested norms with this commented out and indeed is ok, but no noticeable speed gain, so I leave them in to be sure...
    if( Lambda < fabs(l-L) || Lambda > l+L)
        return 0;
    if( Lambda < fabs(l1-l2) || Lambda > l1+l2 )
        return 0;
    if( 2*n+l + 2*N+L != 2*n1+l1+2*n2+l2) return 0;
//	cerr << " calculating missing coeff. " << endl;
//	cerr << n << l << N << L << "|" << n1 << l1 << n2 << l2 << "; " << Lambda << endl ;
    maprecmosh.get(n1,l1,n2,l2,path).updated= true;
    // n1 not zero -> Eq. A.15 Phd thesis M. Vanhalst
    if( n1 > 0 ) {
// 		cout << "n1>0" << endl;
        double factor = 1./sqrt((n1)*(n1+l1+0.5));
        double sum = 0.;
        for( int na = n-1; na <= n ; na++ ) {
            if( na<0) continue;
            for( int la = l-1; la <= l+1; la++) {
                if( la<0) continue;
                for( int Na = N-1; Na <=N; Na++ ) {
                    if( Na<0) continue;
                    for( int La = L-1; La <= L+1; La++) {
                        if( La<0||Lambda<fabs(la-La)||Lambda>la+La) continue;
                        double me = getMatrixElement(n,l,N,L,na,la,Na,La, Lambda, 1);
                        double moshbr = maprecmosh.get( n1-1, l1, n2, l2, path).getCoefficient( na, la, Na, La, Lambda);
                        // double moshbr = calculate(na,la,Na,La,n1-1,l1,n2,l2,Lambda);
                        sum += me*moshbr;
                    }
                }
            }
        }
        return factor*sum;
    //n1 zero but n2 not zero, see PhD Vanhalst p131 after A.15    
    } else if( n2 > 0 ) {
// 		cout << "n2>0" << endl;
        double factor = 1./sqrt((n2)*(n2+l2+0.5));
        double sum = 0.;
        for( int na = n-1; na <= n; na++ ) {
            if( na<0) continue;
            for( int la = l-1; la <= l+1; la++) {
                if( la<0) continue;
                for( int Na = N-1; Na <=N; Na++ ) {
                    if( Na<0) continue;
                    for( int La = L-1; La <= L+1; La++) {
                        if( La<0||Lambda<fabs(la-La)||Lambda>la+La) continue;
                        double me = getMatrixElement(n,l,N,L,na,la,Na,La, Lambda, 2);
                        double moshbr = maprecmosh.get( n1, l1, n2-1, l2, path).getCoefficient( na, la, Na, La, Lambda);
                        //double moshbr = calculate(na,la,Na,La,n1,l1,n2-1,l2,Lambda);
                        sum += me*moshbr;
                    }
                }
            }
        }
        return factor*sum;
    //n1 and n2 zero -> Eq. A.12 PhD thesis M. Vanhalst    
    } else if( n1 == 0 && n2 == 0) {
        double factor = gsl_sf_fact(l1) * gsl_sf_fact(l2)/ gsl_sf_fact(2*l1) / gsl_sf_fact(2*l2);
        factor *= (2*l+1)*(2*L+1)/pow(2.,l+L);
        factor *= gsl_sf_fact(n+l)/gsl_sf_fact(n) / gsl_sf_fact(2*n+2*l+1);
        factor *= gsl_sf_fact(N+L)/gsl_sf_fact(N)/ gsl_sf_fact(2*N+2*L+1);
        double sign = pow( -1., n+l+L-Lambda);
        double sum = 0;
// 	  	for( int x = 0; x <= 4; x++)
        for( int x = fabs(l-l1); x <= l+l1; x++) {
// 			cout << "x allowed ?" << x << endl;
            if( x < fabs(L - l2) || x > L+l2 ) {
                continue;
            }
            double W = Wc( l, L, l1, l2, Lambda, x);
// 			cout << x << " " << W << endl;
            if( W == 0 ) {
                // cout << "NOOOOOOOOOOO" << endl;
                continue;
            }
// 			cout << "x allowed: " << x << endl;
// 			cout << "-------" << endl;
            double term = 2*x+1;
            term *= A( l1, l, l2, L, x );
            // 		cout << " A" << term << endl;
            term *= W;
            sum += term;
            // 		cout << " ------E " << endl;
        }
// 		cout << "<(" << n << l << N << L << ")" << Lambda;
// 		cout << "|(0" << l1 << 0 << l2 << ")" << Lambda << "> = " << sqrt(factor)*sign*sum << endl;
        if(fabs(sum)<1.E-09) sum=0.;
        return sqrt(factor)*sign*sum;
    }

    cerr << "[recmosh][RecMosh::calculate] IMPOSSIBLE " << endl;
    return 0;
}

double RecMosh::getMatrixElement( const int n, const int l, const int N, const int L, 
                                const int na, const int la, const int Na, const int La, const int Lambda, const int f) const 
{
    if( na == n-1 ) {
        if( la == l) {
            if( Na == N) {
                if( La == L) return 0.5*sqrt(n*(n+l+0.5));
                else return 0;
            } else return 0;
        } else if( la == l+1) {
            if( Na == N-1 ) {
                if ( La == L+1)
                    return pow(-1.,f+1)*sqrt( n*N*(l+1)*(L+1))*pow(-1.,Lambda+L+l)*Wc(l,l+1,L,L+1,1,Lambda);
                else return 0;
            } else if (Na == N) {
                if( La == L-1)
                    return pow(-1,f+1)*sqrt(n*(N+L+0.5)*(l+1)*L)*pow(-1.,Lambda+L+l)*Wc(l,l+1,L,L-1,1,Lambda);
                else return 0;
            } else return 0;
        } else return 0;
    } else if ( na == n) {
        if (la == l) {
            if( Na == N-1) {
                if( La == L) return 0.5*sqrt(N*(N+L+0.5));
                else return 0;
            } else return 0;
        } else if ( la == l-1) {
            if( Na == N-1) {
                if( La == L+1)
                    return pow(-1.,f+1)*sqrt((n+l+0.5)*N*l*(L+1))*pow(-1.,Lambda+L+l)*Wc(l,l-1,L,L+1,1,Lambda);
                else return 0;
            } else if ( Na == N) {
                if( La == L-1)
                    return pow(-1.,f+1)*sqrt((n+l+0.5)*(N+L+0.5)*l*L)*pow(-1,Lambda+L+l)*Wc(l,l-1,L,L-1,1,Lambda);
                else return 0;
            } else return 0;
        } else return 0;
    } else return 0;

}

double RecMosh::Wc( const int a, const int b, const int c, const int d, const int e, const int f)const 
{
    double sign = pow(-1.,a+b+c+d);
    double sixj = gsl_sf_coupling_6j( 2*a, 2*b, 2*e, 2*d, 2*c, 2*f);
    return sign*sixj;
}

double RecMosh::A( const int l1, const int l, const int l2, const int L, const int x )const 
{
    double factor = gsl_sf_fact(l1+l+x+1) ;
    factor*= gsl_sf_fact(l1+l-x);
    factor*= gsl_sf_fact(l1+x-l);
    factor *= gsl_sf_fact(l2+L+x+1)*gsl_sf_fact(l2+L-x)*gsl_sf_fact(l2+x-L);
    factor *= 1./gsl_sf_fact(l+x-l1)/gsl_sf_fact(L+x-l2);
    double sum = 0.;
    for( int q = x; q <= l1+l && q<= l2+L; q++ ) {
// 		cout << "q " << q << endl;
        if( ((l+q-l1)%2) || ((L+q-l2)%2) ) {
            continue;
        }
// 		cout << ((l+q-l1)%2) << " " <<  ((L+q-l2)%2)  <<  q << endl;
        double sign = pow(-1., 0.5*(l+q-l1) );
        double term = gsl_sf_fact(l+q-l1)/gsl_sf_fact((l+q-l1)/2)/gsl_sf_fact((l+l1-q)/2);
        term *= 1./gsl_sf_fact(q-x)/gsl_sf_fact(q+x+1);
        term*= gsl_sf_fact(L+q-l2)/gsl_sf_fact((L+q-l2)/2)/gsl_sf_fact((L+l2-q)/2);
        term *= sign;
// 		cout << " t " << term << endl;
        sum += term;
    }
    return sqrt(factor)*sum;


}

/*
 * Check if appropriate recmosh already exists.
 * If not, create new one
 */
RecMosh& MapRecMosh::get( const int n1, const int l1, const int n2, const int l2, const std::string& path )
{
    int key=  n1*1e3+ l1*1e2+ n2*10+ l2;
    map< int, RecMosh* >::iterator it;
    it= maprecmosh.find( key );
    //it's already in the map
    if( it != maprecmosh.end() ) {
        //it->second->use();
        return *(it->second);
    }
    //not in the map so create a new RecMosh object
    RecMosh* newrecmosh= new RecMosh( n1, l1, n2, l2, path );
    maprecmosh[ key]= newrecmosh;
    return *newrecmosh;
}


void MapRecMosh::remove( int key )
{
    map< int, RecMosh* >::iterator it;
    it= maprecmosh.find( key);
    delete (*it).second;
    maprecmosh.erase( it );
}


MapRecMosh::~MapRecMosh()
{
    cout << "~MapRecMosh "<< maprecmosh.size() << endl;
    map< int, RecMosh* >::iterator it;
    for( it= maprecmosh.begin(); it!= maprecmosh.end(); it++ ) {
        delete it->second;
    }
}

