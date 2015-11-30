#include "threej.h"

threej threej::threejs= threej( 40 );

threej::threej( int lmax ) : lmax( lmax )
{
    cmax= lmax*(274+lmax*(225+lmax*(85+lmax*(15+lmax))))/120 + 1;
    cout << "cmax " << cmax << endl;
    vector_threej= new vector< double >( cmax, 0 );
    bool_threej= new vector< bool >( cmax, false );
    cout << vector_threej->size();
    cout << " " << bool_threej->size();
    cout << endl;
}

threej::~threej()
{
    delete vector_threej;
    delete bool_threej;
    cout << "THREEJ DELETED" << endl;
}

double threej::get( int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3 )
{
    double res;
    get( two_j1, two_j2, two_j3, two_m1, two_m2, two_m3, &res );
    return res;
}

void threej::get( int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3, double* res )
{
    *res= 0;
    if ( triangle_selection_fails(two_j1, two_j2, two_j3)
            || m_selection_fails(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3))
        return;

    double prefactor= 1.;
    double phase= 1.;
    if( (((two_j1+ two_j2+ two_j3)/2)%2) ) {
        phase= -1;
    }

    while( true ) {
        if( two_j1 >= two_j3 ) {
            if( two_j3 >= two_j2 )
                break;
            else {
                std::swap( two_j2, two_j3);
                std::swap( two_m2, two_m3);
                prefactor*= phase;
            }
        } else {
            std::swap( two_j1, two_j3);
            std::swap( two_m1, two_m3);
            prefactor*= phase;
        }
    }
    if( two_m2 < 0 ) {
        two_m1*= -1;
        two_m2*= -1;
        two_m3*= -1;
        prefactor*= phase;
    }
    if( two_m2 == 0 && two_m3 < 0 ) {
        two_m1*= -1;
        two_m2*= -1;
        two_m3*= -1;
        prefactor*= phase;
    }

    int Regge_values[9]= {-two_j1+ two_j2+ two_j3, two_j1-two_j2+ two_j3,
                          two_j1+ two_j2- two_j3, two_j1- two_m1,
                          two_j2 - two_m2, two_j3- two_m3,
                          two_j1+ two_m1, two_j2+ two_m2, two_j3+ two_m3
                         };

    int min= gsl_stats_int_min_index( Regge_values, 1, 9 );
    int max= gsl_stats_int_max_index( Regge_values, 1, 9 );


    double Regge_prefactor= 1;
    // Transposition
    if( !((max - min )%3) ) {
//    cout << "transpose " << phase << endl;
        std::swap( Regge_values[1], Regge_values[3] );
        std::swap( Regge_values[2], Regge_values[6] );
        std::swap( Regge_values[5], Regge_values[7] );
        if( min % 4 == 1 )
            min+= 2;
        else if( min % 4 == 3 )
            min-= 2;
        else if( min % 4 == 2 )
            min= 8-min;
        if( max % 4 == 1 )
            max+= 2;
        else if( max % 4 == 3 )
            max-= 2;
        else if( max % 4 == 2 )
            max= 8-max;
    }

    if( min > 2 ) {
//   cout << "row perm" << endl;
        int row= min/3;
        for( int i= 0; i < 3; i++ )
            std::swap( Regge_values[row*3+i], Regge_values[i] );
        Regge_prefactor*= phase;
        min= min%3;
        max= max%3;
        /*    for( int i= 0; i < 9; i++ )
            {
              if( !(i%3) ) cout << endl;
              cout << Regge_values[i] ;
            }
            cout << endl;
            */
    }

    if( min != 0 ) {
//    cout << "col perm" << endl;
        for( int i= 0; i < 3; i++ )
            std::swap( Regge_values[min+3*i], Regge_values[3*i] );
        Regge_prefactor*= phase;
        if (max == 0 )
            max= min;
        min= 0;
        /*    for( int i= 0; i < 9; i++ )
            {
              if( !(i%3) ) cout << endl;
              cout << Regge_values[i] ;
            }
            cout << endl;
            */
    }

    if( max != 1 ) {
//    cout << "col perm" << endl;
        for( int i= 0; i < 3; i++ )
            std::swap( Regge_values[3*i+1], Regge_values[3*i+2] );
        max= 1;
        Regge_prefactor*= phase;
        /*    for( int i= 0; i < 9; i++ )
            {
              if( !(i%3) ) cout << endl;
              cout << Regge_values[i] ;
            }
            cout << endl;
            */
    }

    if( !( Regge_values[4] < Regge_values[7] ) ) {
        if( ( Regge_values[4] != Regge_values[7] ) ) {
//      cout << "row perm" << endl;
            for( int i= 0; i < 3; i++ )
                std::swap( Regge_values[3+i], Regge_values[6+i] );
            Regge_prefactor*= phase;
        } else {
            if( !( Regge_values[5] < Regge_values[8] ) ) {
//        cout << "row perm" << endl;
                for( int i= 0; i < 3; i++ )
                    std::swap( Regge_values[3+i], Regge_values[6+i] );
                Regge_prefactor*= phase;
            }
        }
        /*    for( int i= 0; i < 9; i++ )
            {
              if( !(i%3) ) cout << endl;
              cout << Regge_values[i] ;
            }
            cout << endl;
            */
    }

    int S= Regge_values[0]/2;
    int L= Regge_values[1]/2;
    int X= Regge_values[3]/2;
    int B= Regge_values[4]/2;
    int T= Regge_values[8]/2;

    if( L < X || X < T || T < B || B < S ) {
        return;
    }

    int c= L*(24+L*(50+L*(35+L*(10+L))))/120 +X*(6+X*(11+X*(6+X)))/24
           +T*(2+T*(3+T))/6 +B*(B+1)/2 + S + 1;

//  cout << "->" << two_j1 << " " << two_j2 << " " << two_j3 << endl;
//  cout << "->" << two_m1 << " " << two_m2 << " " << two_m3 << endl;
//  cout << "->" << "index " << c << endl;

//#pragma omp critical(getthreej)
    {
        while( c >= cmax ) {
//#pragma omp critical(getthreej)
            {
                //    cout << "RESIZE " << c << endl;
                //    cout << two_j1 << " " << two_j2 << " " << two_j3 << endl;
                lmax++;
                cmax= lmax*(274+lmax*(225+lmax*(85+lmax*(15+lmax))))/120 + 1;
                cout << cmax << " " << lmax << endl;
                cout << "max size " << vector_threej->max_size() << " " << bool_threej->max_size() << endl;
                vector_threej->resize( cmax, 0 );
                bool_threej->resize( cmax, false );
            }
        }

        if( (*bool_threej)[c] == true ) {
            *res= (*vector_threej)[c]* prefactor* Regge_prefactor;
            //    double threejval= gsl_sf_coupling_3j( two_j1, two_j2, two_j3, two_m1, two_m2, two_m3)* prefactor;
            //    if( *res != threejval )
            //      cout << "FOUND WRONG " <<  *res << " " << threejval << endl;
        } else {
            double threejval= gsl_sf_coupling_3j( two_j1, two_j2, two_j3, two_m1, two_m2, two_m3);
            //    if( threejval == 0 ) cerr << __FILE__ << __LINE__ << endl;
            *res= threejval* prefactor;
//#pragma omp critical(getthreej)
            {
                (*vector_threej)[c]= Regge_prefactor* threejval;
                (*bool_threej)[c]= true;
            }
        }
    }
}

int threej::triangle_selection_fails(int two_ja, int two_jb, int two_jc)
{
    /*
     * enough to check the triangle condition for one spin vs. the other two
     */
    return ( (two_jb < fabs(two_ja - two_jc)) || (two_jb > two_ja + two_jc) ||
             GSL_IS_ODD(two_ja + two_jb + two_jc) );
}

int threej::m_selection_fails(int two_ja, int two_jb, int two_jc,
                              int two_ma, int two_mb, int two_mc)
{
    return (
               fabs(two_ma) > two_ja
               || fabs(two_mb) > two_jb
               || fabs(two_mc) > two_jc
               || GSL_IS_ODD(two_ja + two_ma)
               || GSL_IS_ODD(two_jb + two_mb)
               || GSL_IS_ODD(two_jc + two_mc)
               || (two_ma + two_mb + two_mc) != 0
           );
}
