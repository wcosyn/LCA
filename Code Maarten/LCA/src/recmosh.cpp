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

MapRecMosh RecMosh::maprecmosh = MapRecMosh();

RecMosh::RecMosh(int n1, int l1, int n2, int l2, char* inputPath, char* outputPath)
    : n1(n1),
      l1(l1),
      n2(n2),
      l2(l2),
      coefficients(),
      inputPath(inputPath),
      outputPath(outputPath)

{
    if( n1*10+l1 <= n2*10+l2 )
        loadFile();

    // if n1*10+l1 > n2*10+l2 wordt in getCoefficient de
    // waarde uit n2l2 n1l2 genomen en de juiste voorfactor berekend

// used is 1, so when it is not used any more, it is still saved, so when a next pair or nucleus needs it, it is still loaded, otherwise it again loads the file from input, and updates written out to output are not used anymore.
// Warning, need to be sure that everything gets deleted and written out at the end !!
    used=1;
    updated= false;

}

RecMosh* RecMosh::createRecMosh(int n1, int l1, int n2, int l2, char* inputPath, char* outputPath)
{
    RecMosh* r = maprecmosh.get( n1, l1, n2, l2, inputPath, outputPath );
    r->use();
    return r;

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
        maprecmosh.remove( n1*1e3+ l1*1e2+ n2*10+ l2 );
        //delete this;
    }

}

void RecMosh::writeToFile()
{
    if( updated ) {
        stringstream name;
        name << outputPath << "/recmosh" << n1 << l1 << n2 << l2 << ".dat";
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

    name << inputPath << "/recmosh" << n1 << l1 << n2 << l2 << ".dat";
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

double RecMosh::getCoefficient( int n, int l, int N, int Lambda, int L )
{
    if( 2*n+l+2*N+Lambda != 2*n1+l1+2*n2+l2) {
        return 0;
    }

    if( 10*n+l > 10*N+Lambda ) {
        double val = maprecmosh.get( n1, l1, n2, l2, inputPath, outputPath )->getCoefficient( N, Lambda, n, l, L);
        if( (l1+L)%2 ) return -1*val;
        else return val;
    }

    //
    // Only the files with n1l1 <= n2l2 are calculated and saved,
    // So the n2l2 n1l1 coeff is taken and correct phase factor added
    if( 10*n1+l1 > 10*n2+l2 ) {
        double val = maprecmosh.get( n2, l2, n1, l1, inputPath, outputPath )->getCoefficient( n, l, N, Lambda, L);
        if( (Lambda+L)%2 ) return -1*val;
        else return val;


    }

    int key = 100000*n + 10000*l + 1000*N + 100*Lambda + L*10;

    map< int, double >::iterator it;
    it = coefficients.find(key);
    if( it == coefficients.end()) {
        double result = calculate(n,l,N,Lambda,n1,l1,n2,l2,L);
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
 * The moshinsky braket is calculated in a recursive way,
 * See thesis appendix for the expression
 */
double RecMosh::calculate( int n, int l, int N, int Lambda, int n1, int l1, int n2, int l2, int L)
{
    cout << "[recmosh][RecMosh::calculate] new recmosh value calculated " << endl;
// 	cout << "<(" << n << l << N << Lambda << ")" << L;
//	cout << "|(0" << l1 << 0 << l2 << ")" << L << "> = " << endl;
    if( L < fabs(l-Lambda) || L > l+Lambda)
        return 0;
    if( L < fabs(l1-l2) || L > l1+l2 )
        return 0;
    if( 2*n+l + 2*N+Lambda != 2*n1+l1+2*n2+l2) return 0;
//	cerr << " calculating missing coeff. " << endl;
//	cerr << n << l << N << Lambda << "|" << n1 << l1 << n2 << l2 << "; " << L << endl ;
    updated= true;
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
                    for( int Lambdaa = Lambda-1; Lambdaa <= Lambda+1; Lambdaa++) {
                        if( Lambdaa<0) continue;
                        double me = getMatrixElement(n,l,N,Lambda,na,la,Na,Lambdaa, L, 1);
                        double moshbr = maprecmosh.get( n1-1, l1, n2, l2, inputPath, outputPath)->getCoefficient( na, la, Na, Lambdaa, L);
                        // double moshbr = calculate(na,la,Na,Lambdaa,n1-1,l1,n2,l2,L);
                        sum += me*moshbr;
                    }
                }
            }
        }
        return factor*sum;
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
                    for( int Lambdaa = Lambda-1; Lambdaa <= Lambda+1; Lambdaa++) {
                        if( Lambdaa<0) continue;
                        double me = getMatrixElement(n,l,N,Lambda,na,la,Na,Lambdaa, L, 2);
                        double moshbr = maprecmosh.get( n1, l1, n2-1, l2, inputPath, outputPath)->getCoefficient( na, la, Na, Lambdaa, L);
                        //double moshbr = calculate(na,la,Na,Lambdaa,n1,l1,n2-1,l2,L);
                        sum += me*moshbr;
                    }
                }
            }
        }
        return factor*sum;
    } else if( n1 == 0 && n2 == 0) {
        double factor = gsl_sf_fact(l1) * gsl_sf_fact(l2)/ gsl_sf_fact(2*l1) / gsl_sf_fact(2*l2);
        factor *= (2*l+1)*(2*Lambda+1)/pow(2.,l+Lambda);
        factor *= gsl_sf_fact(n+l)/gsl_sf_fact(n) / gsl_sf_fact(2*n+2*l+1);
        factor *= gsl_sf_fact(N+Lambda)/gsl_sf_fact(N)/ gsl_sf_fact(2*N+2*Lambda+1);
        double sign = pow( -1., n+l+Lambda-L);
        double sum = 0;
// 	  	for( int x = 0; x <= 4; x++)
        for( int x = fabs(l-l1); x <= l+l1; x++) {
// 			cout << "x allowed ?" << x << endl;
            if( x < fabs(Lambda - l2) || x > Lambda+l2 ) {
                continue;
            }
            double W = Wc( l, Lambda, l1, l2, L, x);
// 			cout << x << " " << W << endl;
            if( W == 0 ) {
                // cout << "NOOOOOOOOOOO" << endl;
                continue;
            }
// 			cout << "x allowed: " << x << endl;
// 			cout << "-------" << endl;
            double term = 2*x+1;
            term *= A( l1, l, l2, Lambda, x );
            // 		cout << " A" << term << endl;
            term *= W;
            sum += term;
            // 		cout << " ------E " << endl;
        }
// 		cout << "<(" << n << l << N << Lambda << ")" << L;
// 		cout << "|(0" << l1 << 0 << l2 << ")" << L << "> = " << sqrt(factor)*sign*sum << endl;
        return sqrt(factor)*sign*sum;
    }

    cerr << "[recmosh][RecMosh::calculate] IMPOSSIBLE " << endl;
    return 0;
}

double RecMosh::getMatrixElement( int n, int l, int N, int Lambda, int na, int la, int Na, int Lambdaa, int L, int f)
{
    if( na == n-1 ) {
        if( la == l) {
            if( Na == N) {
                if( Lambdaa == Lambda) return 0.5*sqrt(n*(n+l+0.5));
                else return 0;
            } else return 0;
        } else if( la == l+1) {
            if( Na == N-1 ) {
                if ( Lambdaa == Lambda+1)
                    return pow(-1.,f+1)*sqrt( n*N*(l+1)*(Lambda+1))*pow(-1.,L+Lambda+l)*Wc(l,l+1,Lambda,Lambda+1,1,L);
                else return 0;
            } else if (Na == N) {
                if( Lambdaa == Lambda-1)
                    return pow(-1,f+1)*sqrt(n*(N+Lambda+0.5)*(l+1)*Lambda)*pow(-1.,L+Lambda+l)*Wc(l,l+1,Lambda,Lambda-1,1,L);
                else return 0;
            } else return 0;
        } else return 0;
    } else if ( na == n) {
        if (la == l) {
            if( Na == N-1) {
                if( Lambdaa == Lambda) return 0.5*sqrt(N*(N+Lambda+0.5));
                else return 0;
            } else return 0;
        } else if ( la == l-1) {
            if( Na == N-1) {
                if( Lambdaa == Lambda+1)
                    return pow(-1.,f+1)*sqrt((n+l+0.5)*N*l*(Lambda+1))*pow(-1.,L+Lambda+l)*Wc(l,l-1,Lambda,Lambda+1,1,L);
                else return 0;
            } else if ( Na == N) {
                if( Lambdaa == Lambda-1)
                    return pow(-1.,f+1)*sqrt((n+l+0.5)*(N+Lambda+0.5)*l*Lambda)*pow(-1,L+Lambda+l)*Wc(l,l-1,Lambda,Lambda-1,1,L);
                else return 0;
            } else return 0;
        } else return 0;
    } else return 0;

}

double RecMosh::Wc( int a, int b, int c, int d, int e, int f)
{
    double sign = pow(-1.,a+b+c+d);
    double sixj = gsl_sf_coupling_6j( 2*a, 2*b, 2*e, 2*d, 2*c, 2*f);
    return sign*sixj;
}

double RecMosh::A( int l1, int l, int l2, int Lambda, int x )
{
    double factor = gsl_sf_fact(l1+l+x+1) ;
    factor*= gsl_sf_fact(l1+l-x);
    factor*= gsl_sf_fact(l1+x-l);
    factor *= gsl_sf_fact(l2+Lambda+x+1)*gsl_sf_fact(l2+Lambda-x)*gsl_sf_fact(l2+x-Lambda);
    factor *= 1./gsl_sf_fact(l+x-l1)/gsl_sf_fact(Lambda+x-l2);
    double sum = 0.;
    for( int q = x; q <= l1+l && q<= l2+Lambda; q++ ) {
// 		cout << "q " << q << endl;
        if( ((l+q-l1)%2) || ((Lambda+q-l2)%2) ) {
            continue;
        }
// 		cout << ((l+q-l1)%2) << " " <<  ((Lambda+q-l2)%2)  <<  q << endl;
        double sign = pow(-1., 0.5*(l+q-l1) );
        double term = gsl_sf_fact(l+q-l1)/gsl_sf_fact((l+q-l1)/2)/gsl_sf_fact((l+l1-q)/2);
        term *= 1./gsl_sf_fact(q-x)/gsl_sf_fact(q+x+1);
        term*= gsl_sf_fact(Lambda+q-l2)/gsl_sf_fact((Lambda+q-l2)/2)/gsl_sf_fact((Lambda+l2-q)/2);
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
RecMosh* MapRecMosh::get( int n1, int l1, int n2, int l2, char* inputPath, char* outputPath )
{
    int key=  n1*1e3+ l1*1e2+ n2*10+ l2;
    map< int, RecMosh* >::iterator it;
    it= maprecmosh.find( key );
    if( it != maprecmosh.end() ) {
        //it->second->use();
        return it->second;
    }

    RecMosh* newrecmosh= new RecMosh( n1, l1, n2, l2, inputPath, outputPath );
    maprecmosh[ key]= newrecmosh;
    return newrecmosh;
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

