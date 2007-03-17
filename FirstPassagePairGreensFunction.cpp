//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include <iostream>
#include <exception>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
//#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_roots.h>

#include "factorial.hpp"
#include "bessel.hpp"

#include "HalfOrderBesselGenerator.hpp"

#include "FirstPassagePairGreensFunction.hpp"



FirstPassagePairGreensFunction::
FirstPassagePairGreensFunction( const Real D, 
				const Real kf, 
				const Real Sigma )
    :
    PairGreensFunction( D, kf, Sigma ),
    h( getkf() / ( 4.0 * M_PI * getSigma() * getSigma() * getD() ) ),
    hsigma_p_1( 1.0 + h * getSigma() ),
    a( INFINITY )
{
    this->alphaTable.reserve( 32 );
    this->expTable.reserve( 32 );
    this->psurvTable.reserve( 32 );
}

FirstPassagePairGreensFunction::~FirstPassagePairGreensFunction()
{
    ; // do nothing
}

void FirstPassagePairGreensFunction::seta( Real a )
{
    assert( a > this->getSigma() );

    this->a = a;

    this->alphaTable.clear();
    this->expTable.clear();
    this->psurvTable.clear();
}



const Real 
FirstPassagePairGreensFunction::f_alpha0( const Real alpha ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real alpha_a_m_sigma( alpha * ( a - sigma ) );
    const Real hsigma_p_1( this->hsigma_p_1 );

    Real sin_alpha_a_m_sigma;
    Real cos_alpha_a_m_sigma;
    sincos( alpha_a_m_sigma, &sin_alpha_a_m_sigma, &cos_alpha_a_m_sigma );

    const Real term1( alpha * sigma * cos_alpha_a_m_sigma );
    const Real term2( hsigma_p_1 * sin_alpha_a_m_sigma );

    const Real result( term1 + term2 );

    return result;
}

const Real 
FirstPassagePairGreensFunction::f_alpha0_aux( const Real alpha ) const

{
    const Real a( geta() );
    const Real sigma( getSigma() );

    const Real term1( ( a - sigma ) * alpha );

    const Real angle( this->hsigma_p_1 / ( sigma * alpha ) );
    const Real term2( std::atan( angle ) );

    const Real result( term1 - term2 );

    return result;
}



const Real 
FirstPassagePairGreensFunction::
f_alpha0_aux_F( const Real alpha,
		const f_alpha0_aux_params* const params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real value( params->value );

    return gf->f_alpha0_aux( alpha ) - value;
}


const Real 
FirstPassagePairGreensFunction::alpha0_i( const Int i ) const
{
    assert( i >= 0 );

    const Real sigma( this->getSigma() );

    const Real target( i * M_PI + M_PI_2 );
    f_alpha0_aux_params params = { this, target };


    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &f_alpha0_aux_F ),
	    &params 
	};


    // We know the range of the solution from - Pi/2 <= atan <= Pi.
    const Real rangeFactor( M_PI / ( a - sigma ) );
    Real low( i * rangeFactor );
    Real high( (i+1) * rangeFactor );

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, low, high );

    const unsigned int maxIter( 100 );

    unsigned int j( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );

        low = gsl_root_fsolver_x_lower( solver );
        high = gsl_root_fsolver_x_upper( solver );
	int status( gsl_root_test_interval( low, high, 0.0, 1e-12 ) );
        //	printf("%g %g\n", low, high );


	if( status == GSL_CONTINUE )
	{
	    if( j >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "alpha0_i: failed to converge." 
			  << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++j;
    }

    // printf("%d\n",j);

    const Real alpha( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );
  
    return alpha;
}



const Real 
FirstPassagePairGreensFunction::p_survival_i( const Real alpha,
					      const Real r0 ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );


    Real num1;
    {
	const Real angle_a( alpha * ( a - sigma ) );
	Real sin_a;
	Real cos_a;
	sincos( angle_a, &sin_a, &cos_a );
	num1 = alpha * sigmasq * h - 
	    alpha * ( a - sigma + a * h * sigma ) * cos_a +
	    ( hsigma_p_1 + a * sigma * alphasq ) * sin_a ;
    }

    const Real num2( num_r0( alpha, r0 ) );

    const Real den( r0 * alpha * 
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( num1 * num2 / den );

    return result;
}


const Real 
FirstPassagePairGreensFunction::p_leavea_i( const Real alpha,
					    const Real r0 ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    Real num1;
    {
	const Real angle_a( alpha * ( a - sigma ) );
	Real sin_a;
	Real cos_a;
	sincos( angle_a, &sin_a, &cos_a );

	num1 = - alpha * ( a - sigma + a * h * sigma ) * cos_a
	    + ( hsigma_p_1 + a * sigma * alphasq ) * sin_a;
    }
    
    Real num2;
    {
	const Real angle_r0( alpha * ( r0 - sigma ) );
	Real sin_r0;
	Real cos_r0;
	sincos( angle_r0, &sin_r0, &cos_r0 );
	
	num2 = alpha * sigma * cos_r0 + hsigma_p_1 * sin_r0;
    }
    
    
    const Real den( r0 * alpha *
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( num1 * num2 / den );

    return result;
}


const Real 
FirstPassagePairGreensFunction::p_leaves_i( const Real alpha,
					    const Real r0 ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );


    Real num;
    {
	const Real angle_r0( alpha * ( r0 - sigma ) );
	Real sin_r0;
	Real cos_r0;
	sincos( angle_r0, &sin_r0, &cos_r0 );
	
	num = h * sigmasq * ( alpha * sigma * cos_r0 + hsigma_p_1 * sin_r0 );
    }
		      
    const Real den( r0 *
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( num / den );
	
    return result;
}


const Real
FirstPassagePairGreensFunction::num_r0( const Real alpha,
					const Real r0 ) const
{
    const Real sigma( getSigma() );
    const Real angle_r0( alpha * ( r0 - sigma ) );
    Real sin_r0;
    Real cos_r0;
    sincos( angle_r0, &sin_r0, &cos_r0 );

    const Real hsigma_p_1( this->hsigma_p_1 );
    const Real result( alpha * sigma * cos_r0 + hsigma_p_1 * sin_r0 );

    return result;
}

const Real
FirstPassagePairGreensFunction::p_int_r_i( const Real r,
					   const Real alpha,
					   const Real r0,
					   const Real num_r0 ) const
{
    // NOTE: identical to p_survival_i with a -> r.

    const Real sigma( getSigma() );

    const Real angle_r( alpha * ( r - sigma ) );
    Real sin_r;
    Real cos_r;
    sincos( angle_r, &sin_r, &cos_r );  // do sincos here; latency. 

    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real num1( alpha * sigmasq * h - 
		     alpha * ( r - sigma + r * h * sigma ) * cos_r +
		     ( hsigma_p_1 + r * sigma * alphasq ) * sin_r );


//    const Real num2( num_r0( alpha, r0 ) );
    const Real num2( num_r0 );

    const Real den( r0 * alpha * 
		    ( ( r - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( r + r * h * sigma - h * sigmasq ) ) );

    const Real result( num1 * num2 / den );

    return result;
}


void 
FirstPassagePairGreensFunction::updateAlphaTable0( const Real t ) const
{
    const Real a( geta() );

    RealVector& alphaTable_0( this->getAlphaTable( 0 ) );
    alphaTable_0.clear();

    const Real alpha0_0( this->alpha0_i( 0 ) );
    alphaTable_0.push_back( alpha0_0 );

    const Real Dt( this->getD() * t );

    const Real alpha_cutoff( sqrt( ( - log( 1e-10 ) / Dt )
				   + alpha0_0 * alpha0_0 ) );


//    printf("%g %g\n", alphaTable_0[0], alpha_cutoff );


    const Int maxIter( 2000 );

    Int i( 1 );
    while( true )
    {
	const Real alpha_i( this->alpha0_i( i ) );
	alphaTable_0.push_back( alpha_i );

	if( alpha_i > alpha_cutoff )
	{
//	    printf("%d %g %g\n", 
//		   i, alpha_i, std::exp( - getD() * t * alpha_i * alpha_i ) 
//		   / alpha_i );
	    break;
	}

	if( i >= maxIter )
	{
	    std::cerr << "alphaTable_0: max iteration reached." << std::endl;
	    throw std::exception();
	}

	++i;
    }
}


void 
FirstPassagePairGreensFunction::updatePsurvTable( const Real r0 ) const
{
    const Real a( geta() );

    const RealVector& alphaTable_0( this->alphaTable[0] );
    RealVector& psurvTable( this->psurvTable );
    psurvTable.clear();
    psurvTable.reserve( alphaTable_0.size() );

    for( RealVector::const_iterator i( alphaTable_0.begin() );
	 i != alphaTable_0.end(); ++i )
    {
	const Real alpha( *i );
	psurvTable.push_back( p_survival_i( alpha, r0 ) );
    }
}


void 
FirstPassagePairGreensFunction::updateNum_r0Table( RealVector& num_r0Table,
						   const Real r0 ) const
{
    const Real a( geta() );

    const RealVector& alphaTable_0( this->alphaTable[0] );
    const RealVector& expTable( this->expTable );

    num_r0Table.clear();
    num_r0Table.reserve( alphaTable_0.size() );

    assert( alphaTable_0.size() >= expTable.size() );

    for( unsigned int j( 0 ); j < expTable.size(); ++j )
    {
	const Real alpha( alphaTable_0[j] );
	num_r0Table.push_back( num_r0( alpha, r0 ) );
    }
}


void 
FirstPassagePairGreensFunction::updateExpTable( const Real t ) const
{
    const RealVector& alphaTable_0( this->getAlphaTable( 0 ) );

    RealVector& expTable( this->expTable );
    expTable.clear();

    const Real mDt( - getD() * t );
    const Real alpha0( alphaTable_0[0] );
    const Real value0( 2.0 * exp( mDt * alpha0 * alpha0 ) / alpha0 );
    expTable.push_back( value0 );

    const Real factor( 1.0 / value0 );

    for( RealVector::size_type j( 1 ); j < alphaTable_0.size(); ++j )
    {
	const Real alpha( alphaTable_0[j] );
	const Real value( 2.0 * std::exp( mDt * alpha * alpha ) / alpha );
	expTable.push_back( value );

	if( value * factor < this->CUTOFF )
	{
	    // printf("cutoff at %d; %g\n", j, value );
	    break;
	}
    }

}


const Real 
FirstPassagePairGreensFunction::p_survival( const Real t,
					    const Real r0 ) const
{
    Real p( 0.0 );

    const RealVector& alphaTable_0( this->getAlphaTable( 0 ) );
    const RealVector& expTable( this->expTable );
    const RealVector& psurvTable( this->psurvTable );

    assert( alphaTable_0.size() >= expTable.size() );

    for( RealVector::size_type j( 0 ); j < expTable.size(); ++j )
    {
	const Real value( psurvTable[j] );
	p += value * expTable[j];
    }

    return p;
}

const Real 
FirstPassagePairGreensFunction::p_leaves( const Real t,
					  const Real r0 ) const
{
    Real p( 0.0 );

    const RealVector& alphaTable_0( this->getAlphaTable( 0 ) );
    const RealVector& expTable( this->expTable );

    this->updateExpTable( t );

    for( unsigned int j( 0 ); j < expTable.size(); ++j )
    {
	const Real value( p_leaves_i( alphaTable_0[j], r0 ) );

	p += value * expTable[j];
    }

    return p;
}

const Real 
FirstPassagePairGreensFunction::p_leavea( const Real t,
					  const Real r0 ) const
{
    Real p( 0.0 );

    const RealVector& alphaTable_0( this->getAlphaTable( 0 ) );
    const RealVector& expTable( this->expTable );

    this->updateExpTable( t );

    for( unsigned int j( 0 ); j < expTable.size(); ++j )
    {
	const Real value( p_leavea_i( alphaTable_0[j], r0 ) );

	p += value * expTable[j];
    }

    return p;
}



const Real 
FirstPassagePairGreensFunction::p_int_r( const Real r,
					 const Real t,
					 const Real r0,
					 const RealVector& num_r0Table ) const
{
    Real p( 0.0 );

    const RealVector& alphaTable_0( this->getAlphaTable( 0 ) );
    const RealVector& expTable( this->expTable );

    assert( alphaTable_0.size() >= expTable.size() );

    for( RealVector::size_type j( 0 ); j < expTable.size(); ++j )
    {
	const Real alpha( alphaTable_0[j] );
	const Real value( p_int_r_i( r, alpha, r0, num_r0Table[j] ) );
	p += value * expTable[j];
    }

    return p;
}





const Real
FirstPassagePairGreensFunction::p_survival_F( const Real t,
					      const p_survival_params* params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real r0( params->r0 );
    const Real rnd( params->rnd );

    gf->updateExpTable( t );

    return gf->p_survival( t, r0 ) - rnd;
}


const Real
FirstPassagePairGreensFunction::p_int_r_F( const Real r,
					   const p_int_r_params* params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real t( params->t );
    const Real r0( params->r0 );
    const Real psurv( params->psurv );
    const RealVector& num_r0Table( params->num_r0Table );
    const Real rnd( params->rnd );

    return ( gf->p_int_r( r, t, r0, num_r0Table ) / psurv ) - rnd;
}


const Real FirstPassagePairGreensFunction::drawTime( const Real rnd, 
						     const Real r0 ) const
{
    assert( rnd <= 1.0 && rnd >= 0.0 );
    assert( r0 > 0.0 );

    p_survival_params params = { this, r0, rnd };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_survival_F ),
	    &params 
	};

    Real low( 1e-5 );
    Real high( 10.0 );

    this->updateAlphaTable0( low );
    this->updatePsurvTable( r0 );

    //FIXME: adjust high here.

    // adjust low to make sure tha f( low ) and f( high ) straddle.
    const Real highvalue( GSL_FN_EVAL( &F, high ) );
    while( GSL_FN_EVAL( &F, low ) * highvalue >= 0.0 )
    {
	low *= .1;
	//printf("drawTime: adjusting low: %g\n",low);

	if( fabs( low ) <= 1e-50 )
	{
	    std::cerr << "Couldn't adjust low. (" << low <<
		")" << std::endl;
	    throw std::exception();
	    
	}
	this->updateAlphaTable0( low );
	this->updatePsurvTable( r0 );
    }

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, low, high );

    const unsigned int maxIter( 100 );

    unsigned int i( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );
	low = gsl_root_fsolver_x_lower( solver );
	high = gsl_root_fsolver_x_upper( solver );
	int status( gsl_root_test_interval( low, high, .0, this->CUTOFF ) );

	if( status == GSL_CONTINUE )
	{
	    if( i >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "drawTime: failed to converge." << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++i;
    }
  
    // printf("%d\n", i );


    Real t( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );

#ifndef NDEBUG
    Real p( p_survival( t, r0 ) );
    Real leaves( p_leaves( t, r0 ) );
    Real leavea( p_leavea( t, r0 ) );
    const Real error( p - leaves - leavea );
    if( fabs( error ) >= this->CUTOFF )
    {
	printf("ERROR >= CUTOFF: %g; %g %g %g\n", error, p, leaves, leavea );
    }
#endif

    return t;
}

const bool FirstPassagePairGreensFunction::drawEventType( const Real rnd, 
							  const Real r0,
							  const Real t ) const
{
    // NOTE: The huge assumption for this method to operate correctly is that
    // drawTime() was called immediately before invokation of this with
    // the same parameter r0.

    Real p( p_survival( t, r0 ) );
    Real leaves( p_leaves( t, r0 ) );

    Real value( leaves / p );

//    Real leavea( p - leaves );
//    printf("%g %g %g %g\n", value, p, leaves, leavea);

    if( rnd < value )  
    {
	return true;   // leaves
    }
    else 
    {
	return false;  // leavea
    }
}


const Real FirstPassagePairGreensFunction::drawR( const Real rnd, 
						  const Real r0, 
						  const Real t ) const
{
    assert( rnd <= 1.0 && rnd >= 0.0 );
    assert( r0 > 0.0 );

    this->updateAlphaTable0( t );
    this->updateExpTable( t );
    this->updatePsurvTable( r0 );

    const Real psurv( p_survival( t, r0 ) );

    RealVector num_r0Table;
    updateNum_r0Table( num_r0Table, r0 );

    p_int_r_params params = { this, t, r0, psurv, num_r0Table, rnd };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_int_r_F ),
	    &params 
	};

    Real low( this->getSigma() );
    Real high( this->geta() );

    //const Real lowvalue( GSL_FN_EVAL( &F, low*1.1 ) );
    //    const Real highvalue( GSL_FN_EVAL( &F, high*0.9 ) );
    //printf("%g %g %g\n", lowvalue, highvalue, psurv );

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, low, high );

    const unsigned int maxIter( 100 );

    unsigned int i( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );
	low = gsl_root_fsolver_x_lower( solver );
	high = gsl_root_fsolver_x_upper( solver );
	int status( gsl_root_test_interval( low, high, .0, this->CUTOFF ) );

	if( status == GSL_CONTINUE )
	{
	    if( i >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "drawR: failed to converge." << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++i;
    }
  
    //printf("%d\n", i );


    Real r( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );

    return r;
}


const Real FirstPassagePairGreensFunction::f_alpha( const Real alpha,
						    const Int n ) const
{
    const Real a( this->geta() );
    const Real aAlpha( a * alpha );
    const Real sigmaAlpha( getSigma() * alpha );
    const Real hSigma( geth() * getSigma() );
    const Real realn( static_cast<Real>( n ) );

    const Real hSigma_m_n( hSigma - realn );

#if 1
    // Numerical recipes
    Real tmp, jas1, yas1, jas2, yas2, jaa, yaa;
    bessjy( sigmaAlpha, realn + 0.5, &jas1, &yas1, &tmp, &tmp );
    bessjy( sigmaAlpha, realn + 1.5, &jas2, &yas2, &tmp, &tmp );
    bessjy( aAlpha, realn + 0.5, &jaa, &yaa, &tmp, &tmp );
#else
    // GSL
    const Real factors( sqrt( sigmaAlpha * M_2_PI ) );
    const Real factora( sqrt( aAlpha * M_2_PI ) );
    const Real jas1( gsl_sf_bessel_jl( n, sigmaAlpha ) * factors );
    const Real yas1( gsl_sf_bessel_yl( n, sigmaAlpha ) * factors );
    const Real jas2( gsl_sf_bessel_jl( n+1, sigmaAlpha ) * factors );
    const Real yas2( gsl_sf_bessel_yl( n+1, sigmaAlpha ) * factors );
    const Real jaa( gsl_sf_bessel_jl( n, aAlpha ) * factora );
    const Real yaa( gsl_sf_bessel_yl( n, aAlpha ) * factora );
#endif

//    printf("b %g %g %g %g %g %g\n", jas1, jas2, yas1, yas2, yaa, jaa );
    const Real term1( ( hSigma_m_n * jas1 + sigmaAlpha * jas2 ) * yaa );
    const Real term2( ( hSigma_m_n * yas1 + sigmaAlpha * yas2 ) * jaa );

//    printf("s %g %g %g %g\n", hSigma_m_n * jas1 * yaa, sigmaAlpha * jas2 * yaa,
//	   hSigma_m_n * yas1 * jaa, sigmaAlpha * yas2 * jaa);

//    printf("t %g %g %g %g\n", alpha, term1, term2, term1-term2 );// cos(f_alpha_aux( alpha,n )) );
    const Real result( term1 - term2 );
    
    return result;
}

inline const Real G( const unsigned int n, const unsigned int k )
{
    //    std::cerr << n << ' ' << k << std::endl;
    //    return gsl_sf_fact( n + k ) / ( gsl_sf_fact( k ) * 
    //                                    gsl_sf_fact( n - k ) );    
    return factorial( n + k ) * ( factorial_r( k ) * factorial_r( n - k ) );
}


const Real FirstPassagePairGreensFunction::P( const Int n,
					       const Real x )
{
    Real result( 0.0 );

    Real sx2( 1.0 );
    Int term1( 1 );

    const Real x2sq_r( 1.0 / gsl_pow_2( x + x ) );
    const unsigned int maxm( n / 2 );
    for( unsigned int m( 0 ); m <= maxm; ++m )
    {
	const Real value( term1 * sx2 * G( n, 2 * m ) );
	result += value;

	term1 = - term1;
	sx2 *= x2sq_r;
    }

    return result;
}

const boost::tuple<Real,Real>
FirstPassagePairGreensFunction::P2( const Int n, const Real x )
{
    Real result( 0.0 );
    Real resultp( 0.0 );

    Real sx2( 1.0 );
    Int term1( 1 );

    const Real x2sq_r( 1.0 / gsl_pow_2( x + x ) );
    const unsigned int np1( n + 1 );
    const unsigned int maxm( n / 2 );
    for( unsigned int m( 0 ); m <= maxm; ++m )
    {
	const Real sx2p( term1 * sx2 );
	const unsigned int m2( 2 * m );
	const Real value( sx2p * G( n, m2 ) );
	result += value;

	const Real valuep( sx2p * G( np1, m2 ) );
	resultp += valuep;

	term1 = - term1;
	sx2 *= x2sq_r;
    }

    if( n % 2 )
    {
	resultp += term1 * sx2 * G( np1, np1 );
    }


    return boost::make_tuple( result, resultp );
}


const Real FirstPassagePairGreensFunction::Q( const Int n,
					      const Real x )
{
    Real result( 0.0 );

    Real sx2( 1.0 / ( x + x ) );
    Int term1( 1 );

    const Real x2sq( sx2 * sx2 );
    const unsigned int maxm( (n+1)/2 ); // sum_(0)^((n-1)/2)
    for( unsigned int m( 0 ); m < maxm; ++m )
    {
	const Real value( term1 * sx2 * G( n, 2 * m + 1 ) );
	result += value;

	term1 = - term1;  // (-1)^m
	sx2 *= x2sq;
    }

    return result;
}

const boost::tuple<Real,Real>
FirstPassagePairGreensFunction::Q2( const Int n, const Real x )
{
    Real result( 0.0 );
    Real resultp( 0.0 );

    Real sx2( 1.0 / ( x + x ) );
    Int term1( 1 );  // (-1)^m

    const Real x2sq( sx2 * sx2 );
    const unsigned int np1( n + 1 );
    const unsigned int maxm( (n+1)/2 ); // sum_(0)^((n-1)/2)
    for( unsigned int m( 0 ); m < maxm; ++m )
    {
	const Real sx2p( term1 * sx2 );
	const unsigned int m2p1( 2 * m + 1 );
	const Real value( sx2p * G( n, m2p1 ) );
	result += value;

	const Real valuep( sx2p * G( np1, m2p1 ) );
	resultp += valuep;

	term1 = - term1; // (-1)^m
	sx2 *= x2sq;
    } 


    if( !( n % 2 ) )
    {
	resultp += term1 * sx2 * G( np1, np1 );
    }


    return boost::make_tuple( result, resultp );
}


const Real 
FirstPassagePairGreensFunction::f_alpha_aux( const Real alpha, 
					     const Int n ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );

    const Real aAlpha( a * alpha );
    const Real sigmaAlpha( sigma * alpha );

    const Real realn( static_cast<Real>( n ) );
    const Real n_m_hSigma( n - h * sigma );

    /*(a - s) u - 
      ArcTan[( P[n, a u] ((-n + h s) P[n, s u] + s u Q[1 + n, s u]) -
               Q[n, a u] (s u P[1 + n, s u] + (n - h s) Q[n, s u]) )/
             ( Q[n, a u] ((-n + h s) P[n, s u] + s u Q[1 + n, s u]) + 
               P[n, a u] (s u P[1 + n, s u] + (n - h s) Q[n, s u]) )]
    */

    const Real Pa( P( n, aAlpha ) );
    const Real Qa( Q( n, aAlpha ) );

    Real Ps;
    Real Psp;
    boost::tie( Ps, Psp ) = P2( n, sigmaAlpha );

    Real Qs;
    Real Qsp;
    boost::tie( Qs, Qsp ) = Q2( n, sigmaAlpha );

    const Real n_m_hSigmaPs( n_m_hSigma * Ps );
    const Real n_m_hSigmaQs( n_m_hSigma * Qs );
    const Real sigmaAlphaPsp( sigmaAlpha * Psp );
    const Real sigmaAlphaQsp( sigmaAlpha * Qsp );

    const Real t1( Pa * ( sigmaAlphaQsp - n_m_hSigmaPs ) ); 
    const Real t2( Qa * ( n_m_hSigmaQs  + sigmaAlphaPsp ) );
    const Real t3( Qa * sigmaAlphaQsp + Pa * sigmaAlphaPsp );
    const Real t4( Pa * n_m_hSigmaQs  - Qa * n_m_hSigmaPs );

    const Real angle( (t1 - t2) / (t3 + t4) );

    const Real term1( ( a - sigma ) * alpha );
    const Real term2( std::atan( angle ) );

    const Real result( term1 - term2 );

    return result;
}


const Real 
FirstPassagePairGreensFunction::
f_alpha_aux_F( const Real alpha,
	       const f_alpha_aux_params* const params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Int n( params->n );
    const Real value( params->value );


    return gf->f_alpha_aux( alpha, n ) - value;
}


const Real 
FirstPassagePairGreensFunction::alpha_i( const Int i, const Int n, 
					 gsl_root_fsolver* const solver ) const
{
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    const Real target( M_PI * i + M_PI_2 );

    const Real factor( 1.0 / ( a - sigma ) );
    Real low( (target - M_PI_2) * factor );
    Real high( (target + M_PI_2) * factor );

    f_alpha_aux_params params = { this, n, target };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>
	    ( &FirstPassagePairGreensFunction::f_alpha_aux_F ),
	    &params 
	};

    gsl_root_fsolver_set( solver, &F, low, high );

    const unsigned int maxIter( 100 );
    unsigned int k( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );
	
	low = gsl_root_fsolver_x_lower( solver );
	high = gsl_root_fsolver_x_upper( solver );
	int status( gsl_root_test_interval( low, high, 1e-15, 1e-15 ) );
	
	if( status == GSL_CONTINUE )
	{
	    if( k >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "alpha_i: failed to converge." 
			  << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}
	
	++k;
    }
    
    const Real alpha( gsl_root_fsolver_root( solver ) );

    return alpha;
}


void FirstPassagePairGreensFunction::updateAlphaTable( const Int n,
						       const Real t ) const
{
    assert( n >= 0 );

    if( n == 0 )
    {
	this->updateAlphaTable0( t );
	return;
    }

    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    Real target( M_PI_2 );

    // We know the range of the solution from - Pi/2 <= atan <= Pi/2.
    const Real alphaMid( target / ( a - sigma ) );
    const Real alphaHalfRange( M_PI_2 / ( a - sigma ) );
    Real low( alphaMid - alphaHalfRange * .9 ); // avoid zero.
    Real high( alphaMid + alphaHalfRange );


    // Here we find the interval where the first positive root is in.
    // We find the first pair of alpha
    // ( Pi * offset + Pi/2 ) +- Pi/2 / ( a - sigma )
    // where the values of f_alpha() straddle.
    // The assumption is the interval between roots is not much
    // smaller than Pi / ( a - sigma ).

    Int offset( 0 );
    while( true ) // this can be much faster if better initial guess is given.
    {
	Real lowvalue( f_alpha(low,n) );
	Real highvalue( f_alpha(high,n) );

	if( lowvalue * highvalue < 0 ) // low and high straddle?
	{
	    break;
	}

//	printf("lh: %g %g %g %g\n", low, high, lowvalue, highvalue );
	++offset;
	target = M_PI * offset + M_PI_2;
	low = (target - M_PI_2) / (a-sigma);
	high = (target + M_PI_2) / (a-sigma);
    }


    RealVector& alphaTable_n( this->getAlphaTable( n ) );
    alphaTable_n.clear();

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    const Real alphan_0( alpha_i( offset, n, solver ) );
    alphaTable_n.push_back( alphan_0 );

    const Real Dt( this->getD() * t );
    const Real alpha_cutoff( sqrt( ( - log( 1e-10 ) / Dt )
				   + alphan_0 * alphan_0 ) ); 
    
    const unsigned int MAXI( offset + 500 );
    for( unsigned int i( offset+1 ); i <= MAXI; ++i )
    {
	const Real alpha_i( this->alpha_i( i, n, solver ) );
	alphaTable_n.push_back( alpha_i );

	// cutoff
	if( alpha_i > alpha_cutoff )
	{
//	    printf("alpha %d %g %g\n", 
//		   i-offset, alpha_i, std::exp( - Dt * alpha_i * alpha_i ) );
	    break;
	}
    }

    gsl_root_fsolver_free( solver );
}


const Real FirstPassagePairGreensFunction::p_n_alpha( const Real alpha,
						      const Int n,
						      const Real r,
						      const Real r0, 
						      const Real t ) const
{
    const Real Dt( this->getD() * t );
    const Real sigma( this->getSigma() );
    const Real h( this->geth() );

    const Real alphasq( alpha * alpha );

    const Real aAlpha( a * alpha );
    const Real sigmaAlpha( sigma * alpha );
    const Real hSigma( geth() * getSigma() );
    const Real realn( static_cast<Real>( n ) );
    const Real hSigma_m_n( hSigma - realn );

    const Real term1( alphasq * exp( - Dt * alphasq ) );
    const Real np( realn + 0.5 );

    Real jas1, yas1, jas2, yas2, jaa, yaa, jar, yar, jar0, yar0, _;
    bessjy( sigmaAlpha, np,       &jas1, &yas1, &_, &_ );
    bessjy( sigmaAlpha, np + 1.0, &jas2, &yas2, &_, &_ );
    bessjy( aAlpha,     np,       &jaa,  &yaa,  &_, &_ );
    bessjy( r * alpha,  np,       &jar,  &yar,  &_, &_ );
    bessjy( r0 * alpha, np,       &jar0, &yar0, &_, &_ );

    const Real f_1( hSigma_m_n * jas1 + sigmaAlpha * jas2 );
    const Real f_2( hSigma_m_n * yas1 + sigmaAlpha * yas2 );

    const Real num( ( f_1 * yar - f_2 * jar ) * 
		    ( f_1 * yar0 - f_2 * jar0 ) );

    const Real E1( realn + realn * realn - 
		   sigma * ( h + h * h * sigma + sigma * alphasq ) );

    const Real E2( ( f_1 * f_1 + f_2 * f_2 ) / 
		   ( jaa * jaa + yaa * yaa ) );

//    const Real E2( ( hSigma_m_n * jas1 + sigmaAlpha * jas2 ) / 
//		   jaa * jaa );

    const Real den( E1 + E2 );

    const Real result( term1 * num / den );

    return result;
}


const Real 
FirstPassagePairGreensFunction::p_n( const Int n,
				     const Real r,
				     const Real r0, 
				     const Real t ) const
{
    Real p( 0.0 );

    updateAlphaTable( n, t );

    const RealVector& alphaTable_n( this->getAlphaTable( n ) );

    for( unsigned int i( 0 ); i < alphaTable_n.size(); ++i )
    {
	const Real alpha( alphaTable_n[i] );
	const Real value( p_n_alpha( alpha, n, r, r0, t ) );
	p += value;

	//printf("p_n %d %d %g %g\n", 
	// n, i, value, p_n_alpha( alpha, n, r, r0, t ) );

	if( fabs( value ) < fabs( p * this->CUTOFF ) )
	{
	    break;
	}


    }

    const Real factor( ( 1 + 2 * n ) * M_PI / ( 8.0 * sqrt( r * r0 ) ) );

    p *= factor;


    return p;
}


void
FirstPassagePairGreensFunction::makep_nTable( const Real r, 
					      const Real r0, 
					      const Real t,
					      RealVector& p_nTable ) const
{
    const unsigned NMAX( 1000 );
    const Real truncationTolerance( 1e-8 );

    p_nTable.clear();

    const Real p_0( this->p_n( 0, r, r0, t ) );
    p_nTable.push_back( p_0 );

    const Real threshold( truncationTolerance * p_0 );

    Real p_n_prev( p_0 );
    unsigned int n( 1 );
    while( true )
    {
	Real p_n( this->p_n( n, r, r0, t ) );

	if( p_n < 0.0 )
	{
#ifndef NDEBUG
	    printf("makep_nTable: p_n < 0 %g \n", p_n );
#endif // NDEBUG
	    p_n = 0.0;
	}

	p_nTable.push_back( p_n );

	
	// truncate when converged enough.
	if( fabs( p_n ) < threshold &&
	    fabs( p_n ) < fabs( p_n_prev ) )
	{
	    break;
	}
	
	if( n >= NMAX )
	{
	    std::cerr << "p_n didn't converge." << std::endl;
	    break;
	}
	
	p_n_prev = p_n;
	++n;
    }

}

const Real 
FirstPassagePairGreensFunction::p_theta( const Real theta,
					 const Real r, 
					 const Real r0, 
					 const Real t, 
					 const RealVector& p_nTable ) const
{
    Real p( 0.0 );

    const unsigned int tableSize( p_nTable.size() );

    RealVector LgndTable( tableSize );

    gsl_sf_legendre_Pl_array( tableSize - 1, cos( theta ), 
			      &LgndTable[0] );

    for( RealVector::size_type n( 0 ); n < tableSize; ++n )
    {
	p += p_nTable[n] * LgndTable[n];
//	printf("%d %g %g %g\n",n,p_nTable[n],LgndTable[n],p);
    }

    return p;
}


const Real 
FirstPassagePairGreensFunction::drawTheta( const Real rnd,
					   const Real r, 
					   const Real r0, 
					   const Real t ) const
{
    RealVector p_nTable;
    makep_nTable( r, r0, t, p_nTable );


    const unsigned int tableSize( 100 );
    const Real thetaStep( M_PI / tableSize );

    RealVector pTable( tableSize );
    
    // pTable[0] = 0.0;
    Real p_prev( 0.0 );

    unsigned int i( 1 );
    while( true )
    {
	const Real theta( thetaStep * i );

	Real p( this->p_theta( theta, r, r0, t, p_nTable ) );
	if( p < 0.0 )
	{
	    printf("drawTheta: p<0 %g\n", p );
	    p = 0.0;
	}


	const Real value( ( p_prev + p ) * 0.5 );
	pTable[i] = pTable[i-1] + value;

//	printf("p %g %g %g\n", theta, pTable[i], p );

	if( value < pTable[i] * std::numeric_limits<Real>::epsilon() ||
	    i >= tableSize - 1 )
	{
	    break;   // pTable is valid in [0,i].
	}

	p_prev = p;
	++i;
    }

//    printf("p_int %g %g\n", pTable[i],
//	   pTable[i] * 4.0 * M_PI * r * r / thetaStep );


    Real theta;

    const Real targetPoint( rnd * pTable[i] );
    const size_t lowerBound( gsl_interp_bsearch( &pTable[0], targetPoint, 
						 0, i ) );
    const Real low( lowerBound * thetaStep );
    
    if( pTable[lowerBound+1] - pTable[lowerBound] != 0.0 )
    {
	theta = low + thetaStep * ( targetPoint - pTable[lowerBound] ) / 
	    ( pTable[lowerBound+1] - pTable[lowerBound] );
    }
    else
    {
	// this can happen when rnd is equal to or is too close to 1.0.
	theta = low;
    }


    return theta;
}
    
