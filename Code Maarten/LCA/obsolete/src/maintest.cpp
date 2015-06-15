#include "nucleuspp.h"
#include "nucleusnp.h"
#include "nucleusnn.h"
#include "wsexpansion.h"
#include "wsnucleuspp.h"
#include "wsnucleusnp.h"
#include <string.h>


void numberOfPairs( char* resultdir, char* inputdir, int A, int Z, int n, int l, int S );

void calcwscom(char* resultdir, char* inputdir, char* name, int A, int Z, int n1, int l1, int two_j1, int n2, int l2, int two_j2, double pmin, double pmax );

void calcwscom(char* resultdir, char* inputdir, char* name, int A, int Z, int n1, int l1, int two_j1, int two_mj1, int n2, int l2, int two_j2, int two_mj2, double pmin, double pmax );

void n1_p( char* inputdir, char* resultdir, char* name, int A, int Z, double pmin, double pmax );

//void n1_p_ws( char* inputdir, char* resultdir, char* name, int A, int Z, double pmin, double pmax );

/*
 * Calculate P_2( P )
 */
void n20_pc( char* inputdir, char* resultdir, char* name, int A, int Z, char* xx );
void n20_pc2( char* inputdir, char* resultdir, char* name, int A, int Z, char* xx );


/*
 * Examples of usage Pair class.
 */
int main( int argc, char *argv[])
{
	gsl_set_error_handler_off();
	cout << "argc " << argc << endl;

	if( argc == 7 )
	{
	  // Calculate P_2( p_c )
	  // inpudir resultdir name A Z  nn/pp/pn
	  n20_pc( argv[1], argv[2], argv[3], atoi( argv[4]), atoi( argv[5]), argv[6]);
	  return 1;
	}

	if( argc == 8 )
	{ 
	  cout << "Calculate one particle momentum distribution in WS and HO base " << endl;
	  cout << "-------------- HO ----------------" << endl;
	  n1_p( argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atof(argv[6]), atof(argv[7]) );
	  cout << "-------------- WS ----------------" << endl;
	  //n1_p_ws( argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atof(argv[6]), atof(argv[7]) );
	  return 1; 
	}

	if( argc == 14  )
	{
	  cout << "Calculates WS two-particle com momentum distribution for combination of 2 shells nlj " << endl;

	  //      calcwscom(char* inputdir, char* resultdir, char* name, int A, int Z, int n1, int l1, int two_j1, int n2, int l2, int two_j2, double pmin, double pmax )
	  calcwscom(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), atoi(argv[10]), atoi(argv[11]), atof(argv[12]), atof(argv[13]) );
	  return 1;

	}
	if( argc == 16  )
	{
	  cout << "Calculates WS two-particle com momentum distribution for combination of 2 shells  nlj m_j " << endl;
	  // calcwscom(char* resultdir, char* inputdir, char* name, int A, int Z, int n1, int l1, int two_j1, int two_mj1, int n2, int l2, int two_j2, int two_mj2, double pmin, double pmax );
	  calcwscom(argv[1], argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8]), atoi(argv[9]), atoi(argv[10]), atoi(argv[11]), atoi(argv[12]), atoi(argv[13]), atof(argv[14]), atof(argv[15]) );
	  return 1;

	}

	if( argc == 5 )
	{
	  cout << "Count pairs with certain quantum numbers " << endl; 
	  cout << "Or count distribution of the pairs for each shell combination" << endl;
	  cout << "inputdir resultdir A Z" << endl;
	  char* inputdir = argv[1];
	  char* resultdir = argv[2];
	  int A = atoi( argv[3] );
	  int Z = atoi( argv[4] );
	  NucleusNP* Xnp = new NucleusNP( inputdir, resultdir, A, Z ); 
	  NucleusNN* Xnn = new NucleusNN( inputdir, resultdir, A, Z ); 
	  NucleusPP* Xpp = new NucleusPP( inputdir, resultdir, A, Z ); 
	  Xnp->printPairsPerShell();
	  cout << "np done ... " << endl;
	  Xnn->printPairsPerShell();
	  cout << "nn done ... " << endl;
	  Xpp->printPairsPerShell();
	  cout << "pp done ... " << endl;
	  delete Xpp;
	  delete Xnp;
	  delete Xnn;

	  return 1;
	  
	}
	/*
	if( argc == 8 ){ numberOfPairs( argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]) ); cout << "test" << endl; return 1; }
	*/
	else 
	{
	  cout << " WRONG NUMBER OF ARGUMENTS" << endl;
	  cerr << " WRONG NUMBER OF ARGUMENTS" << endl;
	  return 0;
	}

}

void n20_pc2( char* inputdir, char* resultdir, char* name, int A, int Z, char* xx )
{
	double dP= 0.02;
	WavefunctionP::mapwfp.setA( A );
	WavefunctionP::mapwfp.setpstep( dP );

	Nucleus* X;
	int S= -1;
	int intmax= int(3/dP);
	if( strcmp( xx, "np")==0 || strcmp( xx, "pn" )==0 )
	{ 
	   X = new NucleusNP( inputdir, resultdir, A, Z ); 
	   cout << "np" << endl;
	   S=0;
	   for( int n= 0; n<5; n++)
	   {
	     X->n20_pcom( name, n, 0, S, 0, intmax );
	   }
	     X->n20_pcom( name, -1, 0, S, 0, intmax );

	   S=1;
	   for( int n= 0; n<5; n++)
	   {
	     X->n20_pcom( name, n, 0, S, 0, intmax );
	   }
	     X->n20_pcom( name, -1, 0, S, 0, intmax );
	}
	else if( strcmp( xx, "pp" )==0 )
	{
	   X = new NucleusPP( inputdir, resultdir, A, Z ); 
	   cout << "pp" << endl;
	   for( int n= 0; n<5; n++)
	   {
	     X->n20_pcom( name, n, 0, S, 0, intmax );
	   }
	     X->n20_pcom( name, -1, 0, S, 0, intmax );
	}
	else if( strcmp( xx, "nn" )==0 )
	{
	   X = new NucleusNN( inputdir, resultdir, A, Z ); 
	   cout << "nn" << endl;
	   for( int n= 0; n<5; n++)
	   {
	     X->n20_pcom( name, n, 0, S, 0, intmax );
	   }
	     X->n20_pcom( name, -1, 0, S, 0, intmax );
	}
	else
	{
	  cerr << "no nn/np/pp specified" << endl; exit( 0 );
	}

	delete X;
}

void n20_pc( char* inputdir, char* resultdir, char* name, int A, int Z, char* xx )
{
	double dP= 0.02;
	WavefunctionP::mapwfp.setA( A );
	WavefunctionP::mapwfp.setpstep( dP );

	Nucleus* X;
	int S= -1;
	int intmax= int(3/dP);
        cout << "intmax " << intmax << endl;
	if( strcmp( xx, "np")==0 || strcmp( xx, "pn" )==0 )
	{ 
	   X = new NucleusNP( inputdir, resultdir, A, Z ); 
	   cout << "np" << endl;
	   S=0;
	   for( int l= 0; l<8; l++)
	   {
	     X->n20_pcom( name, l, S, 0, intmax );
	   }
	   X->n20_pcom( name, -1, S, 0, intmax );

	   S=1;
	   for( int l= 0; l<8; l++)
	   {
	     X->n20_pcom( name, l, S, 0, intmax );
	   }
	   X->n20_pcom( name, -1, S, 0, intmax );
	}
	else if( strcmp( xx, "pp" )==0 )
	{
	   X = new NucleusPP( inputdir, resultdir, A, Z ); 
	   cout << "pp" << endl;
	   for( int l= 0; l<8; l++)
	   {
	     X->n20_pcom( name, l, S, 0, intmax );
	   }
	   X->n20_pcom( name, -1, S, 0, intmax );
	}
	else if( strcmp( xx, "nn" )==0 )
	{
	   X = new NucleusNN( inputdir, resultdir, A, Z ); 
	   cout << "nn" << endl;
	   for( int l= 0; l<8; l++)
	   {
	     X->n20_pcom( name, l, S, 0, intmax );
	   }
	   X->n20_pcom( name, -1, S, 0, intmax );
	}
	else
	{
	  cerr << "no nn/np/pp specified" << endl; exit( 0 );
	}

	delete X;
}


void numberOfPairs( char* inputdir, char* resultdir, int A, int Z, int n, int l, int S )
{
	NucleusNP* X = new NucleusNP( inputdir, resultdir, A, Z ); 
	double pairs= X->getlPairs( n, l, S);
	double totalpairs= X->getlPairs( -1, -1, -1);
	delete X;
//	X.~NucleusPP();
//	WSNucleusPP* WSX =  new WSNucleusPP( inputdir, resultdir, A, Z);
//	WSX->makepairs();
//	double wspairs= WSX->getlPairs( l);
//	double totalwspairs= WSX->getlPairs( -1);
//	delete WSX;
//	WSX.~WSNucleusPP();
	cout << " Pairs of " << A << "(" << Z << ")"  << " with l= " << l << " are " << endl;
	cout << " - HO: " << pairs << " / " << totalpairs << endl;
//	cout << " - WS: " << wspairs << " / " << totalwspairs << endl;
	cout << endl;
	return;
}

void calcwscom(char* inputdir, char* resultdir, char* name, int A, int Z, int n1, int l1, int two_j1, int two_mj1, int n2, int l2, int two_j2, int two_mj2, double pmin, double pmax )
{
	int t1= 1;
	int t2= 1;

	stringstream filename;
	filename << resultdir << "/wscommomgsl_" << n1 << l1 << two_j1 << two_mj1 << "_" << n2 << l2 << two_j2 << two_mj2 << "." << name;
	cout << filename.str().c_str() << endl;
	ofstream file( filename.str().c_str() , ofstream::out );
	file << "# com momentum \t all Pairs \t\t 1S0 pairs " << endl; 
	file << "# P = p1+p2 / sqrt(2) " << endl;

	WSPair* pair= new WSPair( inputdir, inputdir, resultdir, A, n1, l1, two_j1, two_mj1, t1, n2, l2, two_j2, two_mj2, t2 );
	double sum = pair->getSum();
	if( fabs(sum) < 1e-4 ) delete pair;
	else if( sum < 0.9 ) cerr << "CHECK " << __FILE__ << ":" << __LINE__ << ": " << n1 << l1 << two_j1 << two_mj2 << " " << n2 << l2 << two_j2 << two_mj2 << endl;
	else
	{

	double integral= 0;
	double integral0= 0;
	double dP= 0.04;
	int intpmin = pmin/dP;
	int intpmax = pmax/dP;
	WavefunctionP::mapwfp.setA( A );
	WavefunctionP::mapwfp.setpstep( dP );

	for( int intp= intpmin; intp<= intpmax; intp++)
	{
		double P= intp* dP;
		double res0= pair->n20_pcom(true, intp);
		double res= pair->n20_pcom(false, intp);
		file << P*Pair::hbarc;
		file << "\t" << res /Pair::hbarc/Pair::hbarc/Pair::hbarc;
		file << "\t" << res0 /Pair::hbarc/Pair::hbarc/Pair::hbarc;
		file << endl;
		integral+= res* P*P* dP;
		integral0+= res0* P*P* dP;
	}

	cout << "integral over P is " << integral << endl;
	//cerr << "integral over P is " << integral << endl;
	cout << "selected integral over P is " << integral0 << endl;
	//cerr << "selected integral over P is " << integral0 << endl;
	file.close();

	delete pair;
	}

}

void calcwscom(char* inputdir, char* resultdir, char* name, int A, int Z, int n1, int l1, int two_j1, int n2, int l2, int two_j2, double pmin, double pmax )
{
	int t1= 1;
	int t2= 1;
	for( int two_mj1= -two_j1; two_mj1 <= two_j1; two_mj1+=2 )
	{
	for( int two_mj2= -two_j2; two_mj2 <= two_j2; two_mj2+=2 )
	{
		if( t1 == t2 && n1 == n2 && l1 == l2 && two_j1 == two_j2 && two_mj1 <= two_mj2 ) continue;
		calcwscom(  inputdir, resultdir, name, A, Z, n1, l1, two_j1, two_mj1, n2, l2, two_j2, two_mj2, pmin, pmax );
	}
	}

}

void n1_p( char* inputdir, char* resultdir, char* name, int A, int Z, double pmin, double pmax )
{
	double dP= 0.04;
	int intpmin = pmin/dP;
	int intpmax = pmax/dP;
	WavefunctionP::mapwfp.setA( A );
	WavefunctionP::mapwfp.setpstep( dP );

	NucleusNP* X_np = new NucleusNP( inputdir, resultdir, A, Z ); 
	X_np->n1_p( name, intpmin, intpmax );

	delete X_np;
}

/*
void n1_p_ws( char* inputdir, char* resultdir, char* name, int A, int Z, double pmin, double pmax )
{

	WSNucleusNP* X_np = new WSNucleusNP( inputdir, resultdir, A, Z ); 
	X_np->n1_p( name, pmin, pmax );

	delete X_np;
}
*/

