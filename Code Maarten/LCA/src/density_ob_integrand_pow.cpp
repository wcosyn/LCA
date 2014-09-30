#include "density_ob_integrand_pow.h"

density_ob_integrand_pow::density_ob_integrand_pow( char* input, char* output, int l1, int l2, int i )
  : input( input ), output( output ), l1( l1 ), l2( l2 ), i( i )
{
  step= 0.1;
  max= 150;
  updated= true;
  loaded= false;
  matrix= new vector< vector< double > > ( max, vector< double >( max, 0 ) );
  loadFile();


}

density_ob_integrand_pow::~density_ob_integrand_pow()
{
  if( updated )
  {
    writeToFile();
  }
  cout << l1 << l2 << i << " density_ob_integrand_pow removed " << endl;
  delete matrix;
}

void density_ob_integrand_pow::calculate()
{
  updated= true;
  gsl_integration_workspace* w= gsl_integration_workspace_alloc( 1e5 );
  gsl_function F;
  F.function= &integrand;
  cout << "Start calculation " << l1 << l2 << i << "..." << endl;
  for( int i1= 0; i1 < max; i1++ )
  {
    for( int i2= 0; i2 < max; i2++ )
    {
      struct params_doip params= { i1*step, i2*step, l1, l2, i };
      F.params= &params;
      double res, abserr;
      int succes = gsl_integration_qagiu( &F, 0, 1e-8, 1e-3, 1e5, w, &res, &abserr );
      if( succes )
      {
        cerr << "integration failed: " << succes << __FILE__ << __LINE__ << endl;
        cerr << res << " " << abserr << endl;
      }
      (*matrix)[i1][i2]= res;
    }
  }
  gsl_integration_workspace_free(w);
  loaded=true;
}

      
double density_ob_integrand_pow::integrand( double r, void* params )
{
  struct params_doip* p = (struct params_doip*) params;
  double k1= (p->k1);
  double k2= (p->k2);
  int l1= (p->l1);
  int l2= (p->l2);
  int i= (p->i);

  gsl_sf_result bessel1;
  int status = gsl_sf_bessel_jl_e(l1, r*k1, &bessel1 );
  if( status )
  {
    if(status == GSL_EUNDRFLW)
    {
      cerr << "Bessel Underflow " << l1 << ", " << r*k1 << endl;
      return 0;
    }
    else 
    {
      cerr << __FILE__ << __LINE__ << " gsl_errno = " << status << endl;
      cerr << l1 << "\t" << r*k1 << endl;
    }
  }

  gsl_sf_result bessel2;
  status = gsl_sf_bessel_jl_e(l2, r*k2, &bessel2 );
  if( status )
  {
    if(status == GSL_EUNDRFLW)
    {
      cerr << "Bessel Underflow " << l2 << ", " << r*k2 << endl;
      return 0;
    }
    else 
    {
      cerr << __FILE__ << __LINE__ << " gsl_errno = " << status << endl;
      cerr << l2 << "\t" << r*k2 << endl;
    }
  }
  gsl_sf_result exp;
  status = gsl_sf_exp_e( -r*r, &exp);
  if( status )
  {
    if(status == GSL_EUNDRFLW)
    {
      return 0;
    }
    else
    {
      cerr << __FILE__ << __LINE__ << " gsl_errno = " << status << endl;
    }
  }
  gsl_sf_result power;
  status= gsl_sf_pow_int_e( r, i, &power );
  if( status )
  {
    cerr << __FILE__ << __LINE__ << " gsl_errno = " << status << endl;
  }
  return bessel1.val* bessel2.val* exp.val* power.val;
}


double density_ob_integrand_pow::get( double k1, double k2 )
{
  if( !loaded )
  {
    calculate();
  }

  int i1= int( k1/step );
  int i2= int( k2/step );
  double delta1= k1-i1*step;
  double delta2= k2-i2*step;

  //
  // Linear Interpolation
  //
  double x1= i1*step;
  double y1= i2*step;
  double z1= (*matrix)[i1][i2] ;
  double x2= x1+step;
  double y2= y1+step;
  double z2= (*matrix)[i1+1][i2+1] ;
  double x3, y3, z3;
  if( delta1 > delta2)
  {
    x3= x2;
    y3= y1;
    z3= (*matrix)[i1+1][i2] ;
  }
  else
  {
    x3= x1;
    y3= y2;
    z3= (*matrix)[i1][i2+1] ;
  }
  double A= y1*(z2-z3)+ y2*(z3-z1)+ y3*(z1-z2);
  double B= z1*(x2-x3)+ z2*(x3-x1)+ z3*(x1-x2);
  double C= x1*(y2-y3)+ x2*(y3-y1)+ x3*(y1-y2);
  double D= -A*x1-B*y1-C*z1;

  double val= -A/C*k1-B/C*k2-D/C;

  return val;
}


void density_ob_integrand_pow::loadFile( )
{
  stringstream name;
  name << input << "/densobpow." << l1 << "." << l2 << "." << i;
  ifstream dataFile ( name.str().c_str(), ifstream::in );
  if( !dataFile )
  {
    cerr << "File " << name.str() << " not found! " << endl;
    dataFile.close();
  }
  else
  {
    cerr << "File " << name.str() << " found! " << endl;
    
    int i=0;
    double data;
    dataFile >> data;
    while ( !dataFile.eof() )
    {
//      cout << i/max << " " <<  i%max << " " << data << endl;
      (*matrix)[i/max][i%max]= data;
      i++;
      dataFile >> data;
    }
    
    /*
    for( int i1= 0; i1< max; i1++ )
    {
      for( int i2= 0; i2< max; i2++ )
      {
        dataFile >> data;
        (*matrix)[i1][i2]= data;
        cout << i1 << i2 << " " << data;
      }
    }
    */
    dataFile.close();
    loaded= true;
  }
}

void density_ob_integrand_pow::writeToFile()
{
  if( updated )
  {
    stringstream name;
    name << output << "/densobpow." << l1 << "." << l2 << "." << i;
    cout << "write to " << name.str().c_str() << endl;
    ofstream dataFile ( name.str().c_str(), ofstream::out );
    if( !dataFile.is_open() )
    {
      cerr << "File could not be opened! Output written to stdout" << endl;
      dataFile.close();
      /*
      vector< vector<double>*>::iterator it1;
      vector< double>::iterator it2;
      for( it1= matrix->begin(); it1 != matrix->end(); it1++ )
      {
        vector< double >* v= *it1;
        for( it2= v->begin(); it2 != v->end(); it2++ )
        {
          cout << setw(14) << setprecision(10) << *it2;
        }
        cout << endl ;
      }
      */

    }
    else
    {
      for( int i= 0; i < matrix->size(); i++ )
      {
        for( int j= 0; j < (*matrix)[i].size(); j++ )
        {
          dataFile << setw(14) << setprecision(7) << (*matrix)[i][j];
        }
        dataFile << endl;
      }

      /*
      vector< vector<double>*>::iterator it1;
      vector< double>::iterator it2;
      for( it1= matrix->begin(); it1 != matrix->end(); it1++ )
      {
        vector< double >* v= *it1;
        for( it2= v->begin(); it2 != v->end(); it2++ )
        {
           dataFile << setw(14) << setprecision(10) << *it2;
//          dataFile << "\t" << *it2;
        }
        dataFile << endl ;
      }
      */
      dataFile.close();
    }
    updated= false;
  }
}

