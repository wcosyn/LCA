#include "wswf.h"


WSWF::WSWF(const double ENER0, const double rmatch, int A, int Z, int n, int l, int two_j, int two_t)
// massConstant = 2* (reduced)mass /hbarc / hbarc
    : Rmin(EPS),
      Rmax(R_MAX),
      Rmatchi(rmatch),
      Rmatchf(rmatch),
      Size(NSTEPS+1),
      n(n),
      l(l),
      two_j(two_j),
      two_t(two_t)
{
    cout << "two_j " << two_j << endl;
    // 6 parameters
    double V = 52.06;
    double kappa= 0.639;
    double R = 1.26;
    a= 0.662;
    double lambda = 24.1;
    double R2 = 1.16;

    double mass;

    if( two_t == 1 ) {
        Zc= Z-1;
        mass = 1./( 1./ 938.272 + 1./(A-1)/931.494 );
    } else {
        Zc = 0;
        mass = 1./( 1./ 939.566 + 1./(A-1)/931.494 );

    }
    massConstant = 2.*mass/197.327/197.327;
    R0= R*pow(A, 1./3.);
    Rc= R0;
    RSO= R2*pow(A, 1./3.);

    double factor;
    if( 2*Z == A ) {
        factor = 3;
    } else if ( A-Z > Z ) {
        if( two_t == 1 ) {
            factor = A-2*Z+3;
        } else {
            factor = -( A-2*Z+1)+ 2 ;
        }
    } else if ( A-Z < Z ) {
        if( two_t == 1 ) {
            factor = A-2*Z+1;
        } else {
            factor = -( A-2*Z-1)+ 2 ;
        }
    } else {
        cerr << "ERROR " << __FILE__ << __LINE__ << endl;
        exit(-1);
    }
    V0 = V*( 1+kappa*factor/A);

    double VSa = lambda* V;
    VSO = VSa*0.5/mass/mass*197.327*197.327;

    cout << mass << "\t" << a << "\t" << R0 << "\t" << -1*V0 << "\n" << RSO << "\t" << a << "\t" << VSO << "\n" << Rc << "\t" << Zc << endl;


    RadialWF = new vector<vector<double> >(NSTEPS+1, vector<double>(3,0));

    cout << "Start" << endl;
    double E = ENER0; // Initial E in MeV
    double Emax, Emin, Diff;
    Emax = 1.5*E;
    Emin = E/1.5;
    Rmax = R_MAX;
    Rmin = EPS;
    int count, countmax = 100;

    double DiffMax = diff(Emax);
    cout << "DiffMax " << DiffMax << endl;
    for( count = 0; count <= countmax; count++) {
        E = (Emax+Emin)/2.;
        Diff = diff(E);
        cout <<"E = " << E << ", L-R Log Diver(E) = " << Diff << endl;
// 		if ( isinf(Diff) ) break;
        if( DiffMax*Diff > 0 ) {
            Emax = E;
            DiffMax = Diff;
        } else Emin = E;
        if( fabs(Diff) < ENER_EPS ) break;
    }
    Energy = E;
    if( count < countmax)
        plot();
    cerr << "Final eigenvalue E = " << Energy << endl;
    cout << "Iterations, max = " << count << ", " << countmax << endl;
}

WSWF::WSWF(const char* path )

{
    RadialWF = new vector<vector<double> >(NSTEPS+1, vector<double>(3,0));

    ifstream dataFile ( path  );
    if( !dataFile ) {
        cerr << "File " << path << " could not be opened! " << endl;
        exit(-1);
    }
    dataFile.ignore(std::numeric_limits<int>::max(), '=');
    dataFile >> massConstant;
    dataFile.ignore(std::numeric_limits<int>::max(), '\n');
//		cout << "massConstant " << massConstant << endl;

    dataFile.ignore(std::numeric_limits<int>::max(), '=');
    dataFile >> R0;
    dataFile.ignore(std::numeric_limits<int>::max(), '\n');
// 		cout << "R0 = "  << R0 << endl;

    dataFile.ignore(std::numeric_limits<int>::max(), '=');
    dataFile >> a;
    dataFile.ignore(std::numeric_limits<int>::max(), '\n');
// 		cout << "a = "  << a << endl;


    dataFile.ignore(std::numeric_limits<int>::max(), '=');
    dataFile >> V0;
    dataFile.ignore(std::numeric_limits<int>::max(), '\n');
// 		cout << "V0 = "  << V0 << endl;

    dataFile.ignore(std::numeric_limits<int>::max(), '=');
    dataFile >> RSO;
    dataFile.ignore(std::numeric_limits<int>::max(), '\n');
// 		cout << "RSO = "  << RSO << endl;

    double aSO;
    dataFile.ignore(std::numeric_limits<int>::max(), '=');
    dataFile >> aSO;
    dataFile.ignore(std::numeric_limits<int>::max(), '\n');
// 		cout << "aSO = "  << aSO << endl;

    dataFile.ignore(std::numeric_limits<int>::max(), '=');
    dataFile >> VSO;
    dataFile.ignore(std::numeric_limits<int>::max(), '\n');
// 		cout << "VSO = "  << VSO << endl;
//
    dataFile.ignore(std::numeric_limits<int>::max(), '=');
    dataFile >> Rc;
    dataFile.ignore(std::numeric_limits<int>::max(), '\n');
// 		cout << "Rc = "  << Rc << endl;
//
    dataFile.ignore(std::numeric_limits<int>::max(), '=');
    dataFile >> Zc;
    dataFile.ignore(std::numeric_limits<int>::max(), '\n');
// 		cout << "Zc = "  << Zc << endl;

    dataFile.ignore(std::numeric_limits<int>::max(), '=');
    dataFile >> n;
    dataFile.ignore(std::numeric_limits<int>::max(), '\n');
// 		cout << "n = "  << n << endl;

    dataFile.ignore(std::numeric_limits<int>::max(), '=');
    dataFile >> l;
    dataFile.ignore(std::numeric_limits<int>::max(), '\n');
// 		cout << "l = "  << l << endl;

    dataFile.ignore(std::numeric_limits<int>::max(), '=');
    dataFile >> two_j;
    dataFile.ignore(std::numeric_limits<int>::max(), '\n');
// 		cout << "two_j = "  << two_j << endl;
//
    dataFile.ignore(std::numeric_limits<int>::max(), '=');
    dataFile >> two_t;
    dataFile.ignore(std::numeric_limits<int>::max(), '\n');

    dataFile.ignore(std::numeric_limits<int>::max(), '=');
    dataFile >> Energy;
    dataFile.ignore(std::numeric_limits<int>::max(), '\n');
// 		cout << "Energy = "  << Energy << endl;
//		cout << path << ": " << n << l << two_j << two_t << endl;

    dataFile.ignore(std::numeric_limits<int>::max(), '\n');
    double r, u, phi, R;
    dataFile >> r;
    dataFile >> u;
    dataFile >> phi;
    dataFile >> R;
    Rmin = r;
    (*RadialWF)[0][0] = r;
    (*RadialWF)[0][1] = u;
    (*RadialWF)[0][2] = phi;
    int i = 0;
    while( !dataFile.eof() ) {
        i++;
        dataFile >> r;
        dataFile >> u;
        dataFile >> phi;
        dataFile >> R;
        (*RadialWF)[i][0] = r;
        (*RadialWF)[i][1] = u;
        (*RadialWF)[i][2] = phi;
    }
    Size = i+1;
    Rmax = r;
    dataFile.close();
// 		cout << "WF created" << endl;
}

WSWF::~WSWF()
{
    delete RadialWF;
}

// void WSWF::calculateMomentumDistribution()
// {
// }

/*
 * After first integration found the correct energy level, integration is
 * repeated, normalised and saved to RadialWF
 */
void WSWF::plot( )	// Repeat Integration for Plot
{

    int imatch = (int) ( (Rmatchf - Rmin)*NSTEPS/(Rmax-Rmin) );

    const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, 2);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (1e-6, 0.0);
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (2);
    struct ODE_params params = {Energy, this};
    gsl_odeiv_system sys = {func, jac, 2, &params};

    double h = 1e-5;
    double r = Rmin;
    double y[2] = {pow( Rmin, l+1.), (l+1)*pow(Rmin, l ) };
    (*RadialWF)[0][0]= Rmin;
    (*RadialWF)[0][1]= y[0];
    (*RadialWF)[0][2]= y[1];
    for (int i = 1; i < imatch+1; i++) {
        double ri = Rmin + i * (Rmax-Rmin) / NSTEPS;
        while (r < ri) {
            int status = gsl_odeiv_evolve_apply (e, c, s,&sys,&r, ri, &h,y);
            if (status != GSL_SUCCESS) {
                cerr << "Error in Plot: " << status << endl;
                break;
            }
        }
        (*RadialWF)[i][0]= r;
        (*RadialWF)[i][1]= y[0];
        (*RadialWF)[i][2]= y[1];
    }
    double Y0match = (*RadialWF)[imatch][1];
    double Y1match = (*RadialWF)[imatch][2];


    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);
    s = gsl_odeiv_step_alloc (T, 2);
    c = gsl_odeiv_control_y_new (1e-6, 0.0);
    e = gsl_odeiv_evolve_alloc (2);


    h = -1e-5;
    double kappa = sqrt(-massConstant*Energy);
    gsl_sf_result exp;
    int status = gsl_sf_exp_e(-kappa*Rmax, &exp);
    if (status) {
        cerr << "err: gsl_errno = " << status << endl;
    }
    if(n%2) {
        y[0] = -exp.val;
        y[1] = -kappa*y[0];
    } else {
        y[0] = exp.val;
        y[1] = -kappa*y[0];
    }
    (*RadialWF)[NSTEPS][0]= Rmax;
    (*RadialWF)[NSTEPS][1]= y[0];
    (*RadialWF)[NSTEPS][2]= y[1];
    r = Rmax;
    for (int i = NSTEPS-1; i >= imatch; i--) {
        double ri = Rmin + i * (Rmax-Rmin) / (NSTEPS);

        while (r > ri) {
            int status = gsl_odeiv_evolve_apply (e, c, s,&sys,&r, ri, &h,y);
            if (status != GSL_SUCCESS)
                break;
        }
        (*RadialWF)[i][0]= r;
        (*RadialWF)[i][1]= y[0];
        (*RadialWF)[i][2]= y[1];
    }
    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);

    double normL =  fabs((*RadialWF)[imatch][1]/Y0match);
    double normL2 =  fabs((*RadialWF)[imatch][2]/Y1match);
    for( int i = 0; i < imatch; i++ ) {
        (*RadialWF)[i][1] *= normL;
    }
    for( int i = 0; i < imatch; i++ ) {
        (*RadialWF)[i][2] *= normL2;
    }


    double normalisation = Normalize();
    cout << normalisation << endl;
    while( fabs(normalisation-1.) > 1e-5 ) {
        for( int i = 0; i < NSTEPS; i++) {
            (*RadialWF)[i][1] *= normalisation;
        }
        normalisation = Normalize();
        cout << normalisation << endl;
    }
    cout << "Rmatch " << Rmatchf << endl;
}

int WSWF::writeToFile(const char* output, const char* name )
{
    if( (*RadialWF)[0][1] == 0 )
        return 1;
    stringstream fileWF;
    fileWF << output << "/WF" << name << n <<  l << two_j << two_t << ".dat";
    ofstream dataFile ( fileWF.str().c_str() );
    stringstream fileV;
    fileV << output << "/V" << name << n << l <<  two_j << two_t << ".dat";
    ofstream dataFileV ( fileV.str().c_str() );
    if( !dataFile || !dataFileV) {
        cerr << "File could not be opened! " << endl;
        exit( 1 );
    }
    dataFile.setf(ios::left);
    dataFileV << "# Mass = " << massConstant << " Mev*(hbarc)^2/2 \n";
    dataFileV << "# R = " << R0 << " fm \n# a = " << a << " fm \n# V0 = " << V0 << " MeV " << endl;
    dataFileV<< "# RSO = " << RSO << " fm \n# aSO = " << a << " fm \n# VSO = " << VSO << " MeV fm^2" << endl;
    dataFileV << "# l = " << l << " \n# two_j = " << two_j << " \n# two_t = " << two_t << " \n# Energy = " << Energy << " MeV" << endl;
    dataFileV << "# r \t V_eff [2m/hbar/hbar]" << endl;
    dataFile << "# Mass = " << massConstant << " Mev*(hbarc)^2/2 \n";
    dataFile << "# R = " << R0 << " fm \n# a = " << a << " fm \n# V0 = " << V0 << " MeV " << endl;
    dataFile << "# RSO = " << RSO << " fm \n# aSO = " << a << " fm \n# VSO = " << VSO << " MeV fm^2" << endl;
    dataFile << "# Rc = " << Rc << " fm\n# Zc = " << Zc << endl;
    dataFile << "# n = " << n << " \n# l = " << l << " \n# two_j = " << two_j << " \n# two_t = " << two_t << " \n# Energy = " << Energy << " MeV" << endl;
    dataFile.width(13);
    dataFile << "# r " ;
    dataFile.width(13);
    dataFile << "u ";
    dataFile.width(13);
    dataFile << "Phi";
    dataFile.width(13);
    dataFile << "u/r " << endl;
    double testIntegral = 0;
    for( int i = 0; i < NSTEPS; i++) {
        dataFile.width(13);
        dataFile << (*RadialWF)[i][0];
        dataFile.width(13);
        dataFile <<  (*RadialWF)[i][1];
        dataFile.width(13);
        dataFile <<  (*RadialWF)[i][2];
        dataFile.width(13);
        dataFile << (*RadialWF)[i][1]/(*RadialWF)[i][0] << endl;
        dataFileV << (*RadialWF)[i][0] << " \t" << -k2((*RadialWF)[i][0], 0)/massConstant << endl;
        testIntegral += (*RadialWF)[i][1]*(*RadialWF)[i][1]*(Rmax-Rmin)/NSTEPS;
    }
    dataFile.close();
    dataFileV.close();
    cout << "Integral is approx " << testIntegral << endl;
    cout << "WF in " << fileWF.str().c_str() << ", V in " << fileV.str().c_str() << endl;
    return 0;
}

void WSWF::write_momentum_density( const char* outputPath, const char* name )
{
    stringstream fileWF;
    fileWF << outputPath << "/MomDens" << name << n <<  l << two_j << two_t << ".dat";
    ofstream dataFile ( fileWF.str().c_str() );
    if( !dataFile) {
        cerr << "File could not be opened! " << endl;
        exit( 1 );
    }
    dataFile.setf(ios::left);
    dataFile << "# Mass = " << massConstant << " Mev*(hbarc)^2/2 \n";
    dataFile << "# R = " << R0 << " fm \n# a = " << a << " fm \n# V0 = " << V0 << " MeV " << endl;
    dataFile << "# RSO = " << RSO << " fm \n# aSO = " << a << " fm \n# VSO = " << VSO << " MeV fm^2" << endl;
    dataFile << "# Rc = " << Rc << " fm\n# Zc = " << Zc << endl;
    dataFile << "# n = " << n << " \n# l = " << l << " \n# two_j = " << two_j << " \n# two_t = " << two_t << " \n# Energy = " << Energy << " MeV" << endl;
    dataFile.width(13);
    dataFile << "# k " ;
    dataFile.width(13);
    dataFile << "n_2 " << endl;
    double testIntegral = 0;
    double stepsize = 0.05;
    double nsteps = (3)/stepsize;

    for( int i = 0; i < nsteps; i++) {
        double k = i* stepsize;
        double result = calculate_momentum_density( k);
        dataFile.width(13);
        dataFile << k;
        dataFile.width(13);
        dataFile << result;
        dataFile << endl;
        testIntegral += result*k*k* stepsize;
    }
    dataFile.close();
    cout << "Momentum Density Integral is approx " << testIntegral << endl;
    cout << "Mom Dist in " << fileWF.str().c_str() << endl;
}


double  WSWF::calculate_momentum_density( double k )
{
    //= Integrate( j_l R r^2 /sqrt(0.5 pi) )^2;
    double res_fourier = fourier( k );
    return res_fourier* res_fourier;
}

double WSWF::fourier( double k)
{
    int status = 0;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);

    double result, error;

    struct fourier_integrand_params params = {this, k };
    gsl_function F;
    F.function = &fourier_integrand;
    F.params= &params;
    status = gsl_integration_qagiu (&F,0 , 1e-3, 1e-3, 100000, w,&result, &error);
    //status = gsl_integration_qag(&F,0, Rmax, 1e-2, 1e-2, 10000,GSL_INTEG_GAUSS61, w,&result, &error);

    if(status) {

        status = gsl_integration_qagiu (&F,0 , 1e-3, 1e-3, 100000, w,&result, &error);
        if(status) {
            status = gsl_integration_qag(&F,0, Rmax, 1e-2, 1e-2, 10000,GSL_INTEG_GAUSS61, w,&result, &error);
            if( status ) {
                cerr << "fourier integrand error status " << status << " " << k  << endl;

                gsl_integration_workspace_free (w);
                exit(-1);
            }
        }
    }
    //cout << "Integral result: " << result << " " << error << endl;
    gsl_integration_workspace_free (w);
    return result;


}

double WSWF::fourier_integrand( double r, void* params )
{
    struct fourier_integrand_params* p = (struct fourier_integrand_params*) params;
    WSWF* wswf= (p->wswf);
    double mom = p->mom;
    int l= wswf->getL();

    double radial = wswf->getRadialWF( r );
    gsl_sf_result bessel;
    int status = gsl_sf_bessel_jl_e(l, r*mom, &bessel );
    if( status ) {
        if(status == GSL_EUNDRFLW) {
            //cerr << "Bessel Underflow" << endl;
            return 0;
        } else cerr << "failed bessel, gsl_errno = " << status << endl;
    }

//        cerr << r << "\t" << radial * bessel.val* r << endl;
    return bessel.val* radial* r/ sqrt(0.5*M_PI);
}


double WSWF::Normalize(  )
{
    int status = 0;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (100000);

    double result, error;

    gsl_function F;
    F.function = &WSIntegral;
    F.params = this;
    status = gsl_integration_qagiu (&F,0, 1e-3, 1e-3, 100000, w,&result, &error);
// 	status = gsl_integration_qag(&F,0, Rmax, 1e-4, 1e-4, 10000,GSL_INTEG_GAUSS61, w,&result, &error);
    gsl_integration_workspace_free (w);

    if(status) {
        cerr << "error status " << status << endl;
        return 1;
    }
    cout << "Integral result: " << result << " " << error << endl;
    return 1./sqrt(result);

}

double WSWF::diff( double energy )
{
    Rmatchf = Rmatchi;
    double left[2];
    double right[2];
    int status = integrate(left, Rmin, Rmatchf, energy);
    status += integrate(right,Rmax, Rmatchf, energy );
    double l = left[1]/left[0];
    double r = right[1]/right[0];
    double diff = (l-r)/(l+r);
    while ( status || isinf(diff) ) {
        cout << status << " " << diff << endl;
        Rmatchf *= 0.95;
        status = integrate(left, Rmin, Rmatchf, energy);
        status += integrate(right,Rmax, Rmatchf, energy );
        l = left[1]/left[0];
        r = right[1]/right[0];
        diff = (l-r)/(l+r);

    }

    return diff;
}


int WSWF::integrate( double solution[], double rstart, double rend, double energy)
{
    const gsl_odeiv_step_type * T
        = gsl_odeiv_step_rk8pd;

    gsl_odeiv_step * s
        = gsl_odeiv_step_alloc (T, 2);
    gsl_odeiv_control * c
        = gsl_odeiv_control_y_new (1e-6, 0.0);
    gsl_odeiv_evolve * e
        = gsl_odeiv_evolve_alloc (2);

    struct ODE_params params = {energy, this};
    gsl_odeiv_system sys = {func, jac, 2, &params};

    double y[2];
    double h;
    if( rstart < rend ) {
        h = 1e-5;
        y[0] = pow( rstart, l+1.);
        y[1] = (l+1.) * pow(rstart, l );
    } else {
        h = -1e-5;
        double kappa = sqrt(-massConstant*energy);
        gsl_sf_result exp;
        int status = gsl_sf_exp_e(-kappa*rend, &exp);
        if (status) {
            if(status == GSL_EUNDRFLW) {
                exp.val = 0;
            } else cerr << "failed, gsl_errno = " << status << endl;
        }
        if(n%2) {
            y[0] = -exp.val;
            y[1] = -kappa*y[0];
        } else {
            y[0] = exp.val;
            y[1] = -kappa*y[0];
        }

    }

    int status;
    double r = rstart;
    for (int i = 1; i <= NSTEPS/10; i++) {
        double ri = rstart + i * (rend - rstart) / (NSTEPS/10);

        int sign = (h > 0) - (h < 0);
        while (sign*r < sign*ri) {
            status = gsl_odeiv_evolve_apply (e, c, s,&sys,&r, ri, &h,y);
            if (status != GSL_SUCCESS) {
                cerr << "gsl_odeiv_evolve_apply err: " << status << endl;
                return status;
            }
        }
// 			cout << r << " " << y[0] << " " <<  y[1] << endl;

// 		if( !isnan(y[1]) && rstart<rend)
// 			cout << i << " " << y[0] << " " << y[1] << endl;
    }

    solution[0] = y[0];
    solution[1] = y[1];

    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);
    return status;
}

// gives R * R * r * r = u * u
double WSWF::WSIntegral( double x, void *p)
{
    WSWF* wswf = (WSWF *) p;
    double wf = wswf->getRadialWF(x);
    return wf*wf;
}

//
// Return u_nlj, not R_nlj !!! R = u/r
//
double WSWF::getRadialWF( double r)
{
    if( r >= Rmax ) {
        return 0;
    }
    if( r < Rmin ) {
        return pow( r, getL()+1 ) * (*RadialWF)[0][1] / pow( Rmin, getL()+1 ) ;
        /*
        if( (*RadialWF)[0][1] < 0.0001 )
        		return 0;
        	return (*RadialWF)[0][1];
                      */
    }
    int i = (int) ( (r-Rmin)*Size/(Rmax-Rmin) );
    while( r > (*RadialWF)[i][0] )
        i++;
    while( r < (*RadialWF)[i-1][0] )
        i--;
    if( i >= Size ) cerr << " i >= Size" << endl;
    double diffy, diffx, diffr, result;
    diffx = (*RadialWF)[i-1][0]-(*RadialWF)[i][0];
// 	cout << r<< " " << (*RadialWF)[i-1][0] << " " << (*RadialWF)[i][0] << " " <<  i << endl;
    if( (*RadialWF)[i-1][1] > (*RadialWF)[i][1] ) {
        diffy = (*RadialWF)[i-1][1] - (*RadialWF)[i][1];
        diffr = (*RadialWF)[i][0]-r;
        result = (*RadialWF)[i][1] + diffy*diffr/diffx;
    } else {
        diffy = (*RadialWF)[i][1] - (*RadialWF)[i-1][1];
        diffr = r-(*RadialWF)[i-1][0];
        result = (*RadialWF)[i-1][1] + diffy*diffr/diffx;
    }
    if( isnan( result ) ) {
        cerr << "result is NaN " << endl;
        cerr << r << " " << i << " " << diffx << " " << diffy << "  " << diffr << endl;
    }
    if( isinf( result ) ) {
        cerr << "result is inf " << endl;
        cerr <<  " " << i << " " << diffx << " " << diffy << "  " << diffr << endl;
        cerr << (*RadialWF)[i][0] << " " << r << " " << (*RadialWF)[i-1][0] << endl;
        cerr << (*RadialWF)[i][1] << " " << result << " " << (*RadialWF)[i-1][1] << endl;
        cerr << "Rmin " << Rmin << endl;
    }
    return result;
}

int WSWF::func (double r, const double y[], double f[], void *params)
{
    struct ODE_params* p = (struct ODE_params *) params;
    double energy = p->energy;
    WSWF* wswf = p->wswf;
    f[0] = y[1];
    f[1] = -wswf->k2(r, energy)*y[0];
    if( isnan(f[0]) || isnan(f[1] ) ) {
        cerr << "NaN error in func " << y[0] << " " << y[1] << ": " << f[0] << " " << f[1] <<endl;
        return GSL_FAILURE;
    }
    return GSL_SUCCESS;
}

int WSWF::jac (double r, const double y[], double *dfdy, double dfdr[], void *params)
{
    struct ODE_params* p = (struct ODE_params *) params;
    double energy = p->energy;
    WSWF* wswf = p->wswf;
    gsl_matrix_view dfdy_mat
        = gsl_matrix_view_array (dfdy, 2, 2);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, 1.0);
    gsl_matrix_set (m, 1, 0, -wswf->k2(r, energy));
    gsl_matrix_set (m, 1, 1, 0.0);
    dfdr[0] = 0.0;
    dfdr[1] = -wswf->dkdr(r, energy)*y[0];
    if( isnan(dfdr[0]) || isnan(dfdr[1] ) ) {
        cerr << "NaN error is jac" << endl;
        return GSL_FAILURE;
    }
    return GSL_SUCCESS;
}

double WSWF::k2( double r, double energy)
{
    gsl_sf_result exp0;
    int status = gsl_sf_exp_e((r-R0)/a, &exp0);
    if (status) {
        if(status == GSL_EUNDRFLW) {
            exp0.val = 0;
        } else cerr << "failed, gsl_errno = " << status << endl;
    }
    gsl_sf_result expSO;
    status = gsl_sf_exp_e((r-RSO)/a, &expSO);
    if (status) {
        if(status == GSL_EUNDRFLW) {
            expSO.val = 0;
        } else cerr << "failed, gsl_errno = " << status << endl;
    }
    double V = -V0/(1.+exp0.val);
    double Vbar = l*(l+1.) / r / r / massConstant;
    double ls = two_j/2.*(two_j/2.+1.) - l*(l+1.) - 3./4.;
// 	double Vso = -1./r/massConstant/197/197 * VSO/(1.+expSO.val)/(1.+expSO.val) * expSO.val / aSO * ls/2.;
    double Vso = -VSO/r/(1.+expSO.val)/(1.+expSO.val) * expSO.val / a * ls/2.;
    double Vc = getVc( r );
    double Veff = V + Vbar + Vso + Vc;
    // cerr << r << "\t"  << Veff << "\t" << V << "\t" << Vbar << "\t" << Vso << "\t" << Vc << endl;
    double result = massConstant * (energy-Veff);
    return result;

}


double WSWF::dkdr( double r, double energy)
{
    gsl_sf_result exp0;
    int status = gsl_sf_exp_e((r-R0)/a, &exp0);
    if (status) {
        if(status == GSL_EUNDRFLW) {
            exp0.val = 0.;
        } else cerr << "failed, gsl_errno = " << status << endl;
    }
    gsl_sf_result expSO;
    status = gsl_sf_exp_e((r-RSO)/a, &expSO);
    if (status) {
        if(status == GSL_EUNDRFLW) {
            expSO.val = 0.;
        } else cerr << "failed, gsl_errno = " << status << endl;
    }
    double dVdr = V0/ ( 1.+exp0.val) / ( 1. + exp0.val ) * exp0.val / a;
    double dVbardr = -2.*l*(l+1.)/r/r/r/massConstant;
    double VSO1 = 1./r/r*expSO.val / (1.+expSO.val)/(1.+expSO.val);
    double VSO2 = 2./r/(1.+expSO.val)/(1.+expSO.val)/(1.+expSO.val) * expSO.val*expSO.val / a;
    double VSO3 = -1./r/(1.+expSO.val)/(1.+expSO.val) * expSO.val / a;
    double ls = two_j/2.*(two_j/2.+1.) - l*(l+1.) - 3./4.;
    double dVSOdr = ls/2. * VSO / a * ( VSO1 + VSO2 + VSO3);
    double dVcdr = getdVcdr( r );
    double result = -1.* massConstant*( dVdr + dVbardr + dVSOdr+ dVcdr );
    return result;
}

double WSWF::getVc( double r )
{
    double factor = 1.* Zc /137. * 197.;
    if ( r > Rc ) {
        return factor/r;
    } else {
        return factor*( 3*Rc*Rc - r*r) / (2*Rc*Rc*Rc);
    }
}

double WSWF::getdVcdr( double r )
{
    double factor = -1.* Zc/137. * 197.;
    if ( r > Rc ) {
        return factor/r/r;
    } else {
        return factor*r/ (Rc*Rc*Rc);
    }
}

double WSWF::getLSCoupling( int two_ml, int two_ms, int two_mj )
{
    if (two_ml + two_ms != two_mj  ) {
        return 0;
    }
    return pow( -1., (2*l-1+two_mj)/2 ) * sqrt(two_j+1) * gsl_sf_coupling_3j( 2*l, 1, two_j, two_ml, two_ms, -two_ml);

}
