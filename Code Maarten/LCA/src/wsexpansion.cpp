#include "wsexpansion.h"


WSexpansion::WSexpansion( WSWF* wf, int A, char* outputPath)
		: wf( wf),
		  coeff(),
		  A(A),
                  outputPath( outputPath)
{
	cout << " start expansion ...";
	cout.flush();
	coeff.reserve(100);
	
	two_j = wf->getTwo_j();
	n = wf->getN();
	l = wf->getL();
	two_t= wf->getTwo_t();
	
	for( int n_HO = 0; n_HO <= 20 ; n_HO++ )
	{
		double res;
		integrate( n_HO, l, &res );
		cout << ".";
		cout.flush();
		if( res*res < 1e-5 && n_HO > n) break;

		// We can sum coefficients with same n_HO, l_HO because radial wavefunction only dependent of these 2 quantum numbers.
		//coeff.push_back( sum );
		coeff.push_back( res );

		if( n_HO == 20 ) cerr << " n_HO limit reached !!!!!!" << endl;
	}
		
	double normalisation = normalize();
	while( fabs(normalisation-1.) > 1e-5 )
	{	
		for( unsigned int i = 0; i < coeff.size(); i++)
		{
			coeff[i] *= normalisation;
		}
		normalisation = normalize();
	}
        /*
	double sum = 0;
	for( u_int i = 0; i < coeff.size() ; i++ )
	{
		sum += coeff[i]*coeff[i];
	}
        */
// 	cout << " test " << sum << endl;
		
// 	cout << "!!! " << coeff.size() << endl;
	cout << " expanded und normalized" << endl;
	

}

WSexpansion::~WSexpansion()
{
	
	stringstream name;
	name << outputPath << "/wsexp" << A << n << l << two_j << two_t << ".dat";
 	cout << "write to " << name.str().c_str() << endl;
	ofstream dataFile ( name.str().c_str(), ofstream::out );
	if( !dataFile.is_open() )
	{
		cerr << "File could not be opened! Output written to stdout" << endl;
		dataFile.close();
// 		cout << 1000*n1 + 100*l1 + 10*n2 + l2  << endl;
		for ( u_int i = 0; i < coeff.size(); i++)
		{
			cout << i << "\t" << coeff[i] << endl;
		}
	}
	for ( u_int i = 0; i < coeff.size(); i++)
	{
		dataFile << i << "\t" << coeff[i] << endl;
	}
	
	dataFile.close();
	
}

double WSexpansion::normalize()
{
	
	struct normf_params params = { coeff, l, A};
	
	int status;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc( 10000000 );
	
	gsl_function F;
	F.function = &normf;
	F.params = &params;
	
	double error, result;
	
	status = gsl_integration_qagiu( &F, 0, 1e-3, 1e-5, 10000000, w, &result, &error);
// 	cout << "r " << result << " " << error << endl;
	
	if(status)
	{
		cerr << "failed integration: error status " << status << " in normalization" << endl;
		gsl_integration_workspace_free(w);
		exit(-1);
	}
	gsl_integration_workspace_free(w);
	
	return 1./sqrt(result);
}

double WSexpansion::normf( double x, void * p )
{
	struct normf_params * params = (struct normf_params *)p;
	vector<double> coeff = (params->coeff);
	int l = params->l;
	int A = params->A;
	double sum = 0.;
	for( unsigned int j = 0; j < coeff.size(); j++ )
	{
		sum += coeff[j]*radialwf(j, l, x, A);
	}
	return x*x*sum*sum;
	
}
void WSexpansion::writeToFile( const char* fileName)
{
	ofstream dataFile ( fileName, ofstream::out );
	if( !dataFile.is_open() )
	{
		cout << "datafile not open: " << endl;
	}
	
	for( double x = 0; x < 12; x+=0.01 )
	{
		double res=0;
		for( unsigned int j = 0; j < coeff.size(); j++ )
		{
			res += coeff[j]*radialwf(j, l, x, A);
		}
		dataFile << x << "\t" << res << endl;
// 			cout << i << "\t" << res << endl;
	}	
	dataFile.close();
}

vector<double> WSexpansion::getCoeff()
{
	return coeff;
}
	


void WSexpansion::integrate( int n_HO, int l_HO, double* result)
{
// 	cout << "start integration" << endl;
	
	struct wsexpf_params params = { wf, n_HO, l_HO, A};
	
	int status;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc( 1e6 );
	
	gsl_function F;
	F.function = &wsexpf;
	F.params = &params;
	
	double error;
	
	status = gsl_integration_qagiu( &F, 0.0, 1e-6, 1e-3, 1e5, w, result, &error);

	//status = gsl_integration_qagiu( &f, 0.01, 1e-2, 1e-3, 100000, w, result, &error);
// 	cout << "r expansion " << n_HO << " " << l_HO << " "<<  *result << " " << error << endl;
	
	if(status)
	{
          cerr << "INTEGRATION FAILED" << endl;
	  status = gsl_integration_qagiu( &F, 0.00,  1e-3, 1e-3, 1e5, w, result, &error);
          if(status)
          {
              cerr << n_HO << " " << status << endl;
              gsl_integration_workspace_free(w);
              exit(-1);
            status = gsl_integration_qagiu( &F, 0.01,  1e-4, 1e-3, 1e5, w, result, &error);
            if( status )
            {
              cerr << n_HO << " " << status << endl;
              status = gsl_integration_qagiu( &F, 0.01,  1e-3, 1e-2, 1e5, w, result, &error);
            
            if(status)
            {
              cerr << n_HO << " " << status << endl;
              status = gsl_integration_qag( &F, 0, 20, 1e-3, 1e-2, 1e6, GSL_INTEG_GAUSS61, w, result, &error);
              if(status )
              {
              status = gsl_integration_qag( &F, 0.01, 12, 1e-3, 1e-2, 1e6, GSL_INTEG_GAUSS61, w, result, &error);

              if( status )
              {
              cerr << "failed integration: error status " << status << " in WSexpansion::integrate " << n_HO << endl;
              
              gsl_integration_workspace_free(w);
              exit(-1);
              }
              }
            }
            }
          }
	}
	gsl_integration_workspace_free(w);
	
	
}

double WSexpansion::wsexpf( double x, void * p )
{
	struct wsexpf_params* params = (struct wsexpf_params *)p;
	int n_HO = (params->n_HO);
	int l_HO = (params->l_HO);
	int A = params->A;
	WSWF* wf = (params->wf);
	return x*radialwf(n_HO, l_HO, x, A)* wf->getRadialWF(x);
}

double WSexpansion::radialwf( int n_HO, int l_HO, double x, int A )
{
	return uncorrelatedradialwf( n_HO, l_HO, x, A);
}

