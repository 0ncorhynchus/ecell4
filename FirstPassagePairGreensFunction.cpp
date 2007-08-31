//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include <iostream>
#include <stdexcept>
#include <vector>
#include <sstream>

#include <boost/bind.hpp>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
//#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sum.h>

#include "factorial.hpp"
#include "bessel.hpp"

//#include "HalfOrderBesselGenerator.hpp"

#include "FirstPassagePairGreensFunction.hpp"



FirstPassagePairGreensFunction::
FirstPassagePairGreensFunction( const Real D, 
				const Real kf, 
				const Real Sigma )
    :
    PairGreensFunction( D, kf, Sigma ),
//    alphaTable( this->MAX_ORDER+1 ),
    h( kf / ( 4.0 * M_PI * Sigma * Sigma * D ) ),
    hsigma_p_1( 1.0 + h * Sigma ),
    a( INFINITY )
{
    ; // do nothing
}

FirstPassagePairGreensFunction::~FirstPassagePairGreensFunction()
{
    ; // do nothing
}

void FirstPassagePairGreensFunction::seta( const Real a )
{
    THROW_UNLESS( std::invalid_argument, a >= this->getSigma() );

    this->a = a;

    clearAlphaTable();
}

//
// Alpha-related methods
//

void FirstPassagePairGreensFunction::clearAlphaTable() const
{
    std::for_each( this->alphaTable.begin(), this->alphaTable.end(),
		   boost::mem_fn( &RealVector::clear ) );
    std::fill( this->alphaOffsetTable.begin(), this->alphaOffsetTable.end(),
	       -1 );

}


const Real 
FirstPassagePairGreensFunction::f_alpha0( const Real alpha ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );

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
    const Real a( this->geta() );
    const Real sigma( this->getSigma() );

    const Real term1( ( a - sigma ) * alpha );

    const Real angle( this->hsigma_p_1 / ( sigma * alpha ) );
    const Real term2( std::atan( angle ) );

    const Real result( term1 - term2 );

    return result;
}




const Real 
FirstPassagePairGreensFunction::
f_alpha0_aux_F( const Real alpha, const f_alpha0_aux_params* const params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real value( params->value );

    return gf->f_alpha0_aux( alpha ) - value;
//    return gf->f_alpha0( alpha );
}


const Real 
FirstPassagePairGreensFunction::alpha0_i( const Integer i ) const
{
    THROW_UNLESS( std::out_of_range, i >= 0 );

    const Real a( this->geta() );
    const Real sigma( this->getSigma() );


    const Real target( i * M_PI + M_PI_2 );
    f_alpha0_aux_params params = { this, target };


    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &f_alpha0_aux_F ),
	    &params 
	};


    // We know the range of the solution from - Pi/2 <= atan <= Pi/2.
    const Real interval( M_PI / ( a - sigma ) );
    Real low( i * interval );// + std::numeric_limits<Real>::epsilon() );
    Real high( (i+1) * interval );

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
	const int status( gsl_root_test_interval( low, high, 0.0, 1e-15 ) );

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

    const Real alpha( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );
  
//    printf("alphavalue %g, alpha %g\n",GSL_FN_EVAL( &F, alpha ), alpha);

    return alpha;
}


void
FirstPassagePairGreensFunction::updateAlphaTable0( const Real t ) const
{
    RealVector& alphaTable_0( this->getAlphaTable( 0 ) );
    alphaTable_0.clear();
    alphaTable_0.reserve( MAX_ALPHA_SEQ );

    const Real alpha0_0( this->alpha0_i( 0 ) );
    alphaTable_0.push_back( alpha0_0 );

    const Real Dt( this->getD() * t );

    const Real alpha_cutoff( sqrt( ( - log( ALPHA_CUTOFF ) / Dt )
				   + alpha0_0 * alpha0_0 ) );


//    printf("%g %g\n", alpha0_0, alpha_cutoff );

    unsigned int i( 1 );
    while( true )
    {
	const Real alpha0_i( this->alpha0_i( i ) );
	alphaTable_0.push_back( alpha0_i );

	if( alpha0_i > alpha_cutoff )
	{
	    break;
	}

	++i;

	if( i >= MAX_ALPHA_SEQ )
	{
	    break;
	}
    }
}


const Real FirstPassagePairGreensFunction::f_alpha( const Real alpha,
						    const Integer n ) const
{
    const Real a( this->geta() );
    const Real aAlpha( a * alpha );
    const Real sigmaAlpha( getSigma() * alpha );
    const Real hSigma( geth() * getSigma() );
    const Real realn( static_cast<Real>( n ) );

    const Real hSigma_m_n( hSigma - realn );
#if 0
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


const Real FirstPassagePairGreensFunction::P( const Integer n,
					       const Real x )
{
    Real result( 0.0 );

    Real sx2( 1.0 );
    Integer term1( 1 );

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
FirstPassagePairGreensFunction::P2( const Integer n, const Real x )
{
    Real result( 0.0 );
    Real resultp( 0.0 );

    Real sx2( 1.0 );
    Integer term1( 1 );

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


const Real FirstPassagePairGreensFunction::Q( const Integer n,
					      const Real x )
{
    Real result( 0.0 );

    Real sx2( 1.0 / ( x + x ) );
    Integer term1( 1 );

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
FirstPassagePairGreensFunction::Q2( const Integer n, const Real x )
{
    Real result( 0.0 );
    Real resultp( 0.0 );

    Real sx2( 1.0 / ( x + x ) );
    Integer term1( 1 );  // (-1)^m

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
					     const Integer n ) const
{
    if( alpha == 0.0 )
    {
	return -1.0;
    }

    const Real a( geta() );
    const Real sigma( getSigma() );

    const Real aAlpha( a * alpha );
    const Real sigmaAlpha( sigma * alpha );

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

    //printf("aux %g %g %g %g\n",alpha,term1, term2, result );

    return result;
}


const Real 
FirstPassagePairGreensFunction::
f_alpha_aux_F( const Real alpha,
	       const f_alpha_aux_params* const params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Integer n( params->n );
    const Real value( params->value );


    return gf->f_alpha_aux( alpha, n ) - value;
//    return gf->f_alpha( alpha, n );
}


const Real 
FirstPassagePairGreensFunction::alpha_i( const Integer i, const Integer n, 
					 gsl_root_fsolver* const solver ) const
{
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    const Real target( M_PI * i + M_PI_2 );

    const Real factor( 1.0 / ( a - sigma ) );
    Real low( ( target - M_PI_2 ) * factor );
    Real high( ( target + M_PI_2 ) * factor );

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
	const int status( gsl_root_test_interval( low, high, 1e-6, 1e-21 ) );
	
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


const unsigned int
FirstPassagePairGreensFunction::alphaOffset( const unsigned int n ) const
{
    if( this->alphaOffsetTable[n] >= 0 )
    {
	return this->alphaOffsetTable[n];
    }

    unsigned int offset( 0 );

    if( n == 0 )
    {
	return offset;
    }

    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    Real target( M_PI_2 );
    const Real factor( 1.0 / (a-sigma) );

    // We know the range of the solution from - Pi/2 <= atan <= Pi/2.
    const Real alphaMid( target * factor );
    const Real alphaHalfRange( M_PI_2 * factor );
    Real low( alphaMid - alphaHalfRange * (1.0 - 1e-3) ); // avoid zero.
    Real high( alphaMid + alphaHalfRange );


    // Here we find the interval where the first positive root is in.
    // We find the first pair of alpha
    // ( Pi * offset + Pi/2 ) +- Pi/2 / ( a - sigma )
    // where the values of f_alpha() straddle.
    // The assumption is the interval between roots is not much
    // smaller than Pi / ( a - sigma ).


    Real lowvalue( f_alpha(low,n) );
    Real highvalue( f_alpha(high,n) );

    while( true ) // this can be much faster if better initial guess is given.
    {

	if( lowvalue * highvalue < 0 ) // low and high straddle?
	{
	    break;
	}

//	printf("lh: %d %d %g %g %g %g\n", 
//	       n, offset, low, high, f_alpha( low, n ), f_alpha( high, n ) );
	++offset;
	target = M_PI * offset + M_PI_2;
	low = ( target - M_PI_2 ) * factor;
	high = ( target + M_PI_2 ) * factor;

	lowvalue = highvalue;
	highvalue = f_alpha( high, n );
    }

    this->alphaOffsetTable[n] = offset;

    return offset;
}


void
FirstPassagePairGreensFunction::updateAlphaTable( const unsigned int n,
						  const Real t ) const
{
    THROW_UNLESS( std::range_error, n >= 0 && n <= this->MAX_ORDER );

    if( n == 0 )
    {
	this->updateAlphaTable0( t );
	return;
    }

    const unsigned int offset( alphaOffset( n ) );

    RealVector& alphaTable_n( this->getAlphaTable( n ) );
    alphaTable_n.clear();
    alphaTable_n.reserve( MAX_ALPHA_SEQ );

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    const Real alphan_0( alpha_i( offset, n, solver ) );
    const Real alphan_0_sq( alphan_0 * alphan_0 );

    alphaTable_n.push_back( alphan_0 );

    const Real Dt( this->getD() * t );

    const Real threshold( this->ALPHA_CUTOFF * 
			  alphan_0_sq * exp( - Dt * alphan_0_sq ) );
   
    const unsigned int end( offset + MAX_ALPHA_SEQ );
    unsigned int i( offset + 1 );
    while( true )
    {
	const Real alpha_i( this->alpha_i( i, n, solver ) );

//	printf("alpha %d %d %g %g %g\n", n, i, alpha_i, f_alpha(alpha_i+1,n),
//	       f_alpha(alpha_i-1,n));

	alphaTable_n.push_back( alpha_i );

	// cutoff
	const Real alpha_i_sq( alpha_i * alpha_i );
	if( alpha_i_sq * exp( - Dt * alpha_i_sq )  < threshold )
	{
	    break;
	}

	if( i >= end )
	{
	    std::cerr << "alphaTable (" << n << 
		"): didn't converge. t = " << t  << ", " 
		      << dump() << std::endl;
	    break;
	}

	++i;
    }

    gsl_root_fsolver_free( solver );
}




const Real 
FirstPassagePairGreensFunction::p_0_i( const Real alpha, 
				       const Real r,
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
	const Real angle_r( alpha * ( r - sigma ) );
	Real sin_r;
	Real cos_r;
	sincos( angle_r, &sin_r, &cos_r );
	num1 = alpha * sigma * cos_r + hsigma_p_1 * sin_r ;
    }

    const Real num2( num_r0( alpha, r0 ) );

    const Real den( 2 * M_PI * r * r0 * 
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( num1 * num2 / den );

    return result;
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

	/*
	num1 = alpha * sigmasq * h - 
	    alpha * ( a - sigma + a * h * sigma ) * cos_a +
	    ( hsigma_p_1 + a * sigma * alphasq ) * sin_a ;
	*/

	num1 = alpha * sigmasq * h - 
	    alpha * ( a + a * h * sigma ) * cos_a +
	    ( a * sigma * alphasq ) * sin_a ;
    }

    const Real num2( num_r0( alpha, r0 ) );

    const Real den( r0 * alphasq * 
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( 2 * num1 * num2 / den );

    return result;
}


const Real 
FirstPassagePairGreensFunction::dp_survival_i( const Real alpha,
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

	num1 = alpha * ( - a * hsigma_p_1 * cos_a +
			 sigma * ( h * sigma + a * alpha * sin_a ) );
    }

    const Real num2( num_r0( alpha, r0 ) );

    const Real den( r0 * 
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( - 2 * getD() * num1 * num2 / den );

    return result;
}


const Real 
FirstPassagePairGreensFunction::leavea_i( const Real alpha,
					  const Real r0 ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real D( getD() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    Real num1;
    {
	const Real angle_a( alpha * ( a - sigma ) );
	Real sin_a;
	Real cos_a;
	sincos( angle_a, &sin_a, &cos_a );

	num1 = alpha * ( hsigma_p_1 * cos_a - sigma * alpha * sin_a );
    }

    const Real num2( num_r0( alpha, r0 ) );
    
    const Real den( 2 * a * M_PI * r0 *
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( D * num1 * num2 / den );

    return result;
}


const Real 
FirstPassagePairGreensFunction::leaves_i( const Real alpha,
					  const Real r0 ) const
{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real D( getD() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real num( h * alpha * num_r0( alpha, r0 ) );
		      
    const Real den( 2 * M_PI * r0 *
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( - D * num / den );
	
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

//	num1 = - ( a - sigma + a * h * sigma ) * alpha * cos_a + 
//	    ( hsigma_p_1 + a * sigma * alphasq ) * sin_a;
	num1 = - ( a * h * sigma + a ) * cos_a + a * sigma * alpha * sin_a;
    }

    const Real num2( num_r0( alpha, r0 ) );
    
    const Real den( r0 * alpha *
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( 2 * num1 * num2 / den );

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

    const Real num( h * sigmasq * num_r0( alpha, r0 ) );
		      
    const Real den( r0 * alpha *
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( 2 * num / den );
	
    return result;
}


const Real FirstPassagePairGreensFunction::num_r0( const Real alpha,
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
    const Real sigma( getSigma() );

    const Real angle_r( alpha * ( r - sigma ) );
    Real sin_r;
    Real cos_r;
    sincos( angle_r, &sin_r, &cos_r );  // do sincos here; latency. 

    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real hsigma( h * sigma );

    const Real num1( alpha * ( hsigma * sigma - hsigma * r * cos_r +
			       ( sigma - r ) * cos_r ) +
		     ( hsigma_p_1 + r * sigma * alphasq ) * sin_r );

    const Real num2( num_r0 );

    const Real den( r0 * alphasq * 
		    ( ( a - sigma ) * sigmasq * alphasq +
		      hsigma_p_1 * ( a + a * h * sigma - h * sigmasq ) ) );

    const Real result( 2 * num1 * num2 / den );


    return result;
}


void 
FirstPassagePairGreensFunction::
createPsurvTable( RealVector& psurvTable, const Real r0 ) const
{
    const RealVector& alphaTable_0( this->getAlphaTable( 0 ) );

    psurvTable.clear();
    psurvTable.reserve( alphaTable_0.size() );

    std::transform( alphaTable_0.begin(), alphaTable_0.end(),
		    std::back_inserter( psurvTable ),
		    boost::bind( &FirstPassagePairGreensFunction::p_survival_i,
				 this, _1, r0 ) );
}


void 
FirstPassagePairGreensFunction::createNum_r0Table( RealVector& num_r0Table,
						   const Real r0 ) const
{
    const RealVector& alphaTable_0( this->alphaTable[0] );

    num_r0Table.clear();
    num_r0Table.reserve( alphaTable_0.size() );

    std::transform( alphaTable_0.begin(), alphaTable_0.end(),
		    std::back_inserter( num_r0Table ),
		    boost::bind( &FirstPassagePairGreensFunction::num_r0,
				 this, _1, r0 ) );
}

const Real 
FirstPassagePairGreensFunction::p_0_i_exp( const unsigned int i,
					   const Real t,
					   const Real r,
					   const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * p_0_i( alpha, r, r0 );
}


const Real 
FirstPassagePairGreensFunction::p_survival_i_exp( const unsigned int i,
						  const Real t,
						  const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * 
	p_survival_i( alpha, r0 );
}

const Real 
FirstPassagePairGreensFunction::
p_survival_i_exp_table( const unsigned int i,
			const Real t,
			const Real r0,
			const RealVector& psurvTable ) const
{

    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * psurvTable[i];
}

const Real 
FirstPassagePairGreensFunction::dp_survival_i_exp( const unsigned int i,
						   const Real t,
						   const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * 
	dp_survival_i( alpha, r0 );
}

const Real 
FirstPassagePairGreensFunction::leavea_i_exp( const unsigned int i,
					      const Real t,
					      const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * leavea_i( alpha, r0 );
}

const Real 
FirstPassagePairGreensFunction::leaves_i_exp( const unsigned int i,
					      const Real t,
					      const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * leaves_i( alpha, r0 );
}

const Real 
FirstPassagePairGreensFunction::p_leavea_i_exp( const unsigned int i,
						const Real t,
						const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * p_leavea_i( alpha, r0 );
}

const Real 
FirstPassagePairGreensFunction::p_leaves_i_exp( const unsigned int i,
						const Real t,
						const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * p_leaves_i( alpha, r0 );
}

const Real 
FirstPassagePairGreensFunction::
p_int_r_i_exp( const unsigned int i,
	       const Real t,
	       const Real r,
	       const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * 
	p_int_r_i( r, alpha, r0, num_r0( alpha, r0 ) );
}

const Real 
FirstPassagePairGreensFunction::
p_int_r_i_exp_table( const unsigned int i,
		     const Real t,
		     const Real r,
		     const Real r0,
		     const RealVector& num_r0Table ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * 
	p_int_r_i( r, alpha, r0, num_r0Table[i] );
}

const Real 
FirstPassagePairGreensFunction::
funcSum( boost::function<const Real( const unsigned int i )> f,
	 const size_t max_i,
	 const Real tolerance )
{
    Real sum( 0.0 );

    RealVector pTable;
    const Real p_0( f( 0 ) );
    if ( p_0 == 0.0 )
    {
	return 0.0;
    }

    const Real threshold( fabs( p_0 * tolerance ) );
    pTable.push_back( p_0 );

    bool extrapolationNeeded( true );

    RealVector::size_type i( 1 ); 
    while( i <= max_i )
    {
	const Real p_i( f( i ) );
	pTable.push_back( p_i );
	
	if( threshold >= fabs( p_i ) ) // '=' is important when p0 is so small.
	{
	    extrapolationNeeded = false;
	    break;
	}
	
	++i;
    }

    if( ! extrapolationNeeded )
    {
	sum = std::accumulate( pTable.begin(), pTable.end(), 0.0 );
    }
    else
    {
	// std::cerr << "Using series acceleration." << std::endl;
	Real error;

	gsl_sum_levin_u_workspace* 
	    workspace( gsl_sum_levin_u_alloc( i ) );
	gsl_sum_levin_u_accel( &pTable[0], pTable.size(), workspace, 
        &sum, &error );
/*	gsl_sum_levin_utrunc_workspace* 
	    workspace( gsl_sum_levin_utrunc_alloc( i ) );
	gsl_sum_levin_utrunc_accel( &pTable[0], pTable.size(), workspace, 
        &sum, &error );*/
	if( fabs( error ) >= fabs( sum * tolerance ) )
	{
	    std::cerr << "Series acceleration error; "
		      << fabs( error ) << " (rel error: " 
		      << fabs( error / sum ) << "), terms_used = " 
		      << workspace->terms_used << " (" 
		      << pTable.size() << " given)." << std::endl;
	}

	gsl_sum_levin_u_free( workspace );
//	gsl_sum_levin_utrunc_free( workspace );
    }

    return sum;
}

const Real 
FirstPassagePairGreensFunction::p_0( const Real t,
				     const Real r,
				     const Real r0 ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_0_i_exp,
					this,
					_1, t, r, r0 ),
			   this->MAX_ALPHA_SEQ ) );
    return p;
}


const Real 
FirstPassagePairGreensFunction::p_survival( const Real t,
					    const Real r0 ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_survival_i_exp, 
					this,
					_1, t, r0 ), 
			   this->MAX_ALPHA_SEQ ) );
    return p;
}

const Real 
FirstPassagePairGreensFunction::
p_survival_table( const Real t,
		  const Real r0,
		  const RealVector& psurvTable ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_survival_i_exp_table, 
					this,
					_1, t, r0, psurvTable ),
			   psurvTable.size()-1 ) );

    return p;
}

const Real 
FirstPassagePairGreensFunction::dp_survival( const Real t,
					     const Real r0 ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					dp_survival_i_exp, 
					this,
					_1, t, r0 ),
			   this->MAX_ALPHA_SEQ ) );
    return p;
}


const Real 
FirstPassagePairGreensFunction::leaves( const Real t,
					const Real r0 ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					leaves_i_exp,
					this,
					_1, t, r0 ),
			   this->MAX_ALPHA_SEQ ) );
    return p;
}

const Real 
FirstPassagePairGreensFunction::leavea( const Real t,
					const Real r0 ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					leavea_i_exp,
					this,
					_1, t, r0 ),
			   this->MAX_ALPHA_SEQ ) );
    return p;
}


const Real 
FirstPassagePairGreensFunction::p_leaves( const Real t,
					  const Real r0 ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_leaves_i_exp,
					this,
					_1, t, r0 ),
			   this->MAX_ALPHA_SEQ ) );
    return p;
}


const Real 
FirstPassagePairGreensFunction::p_leavea( const Real t,
					  const Real r0 ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_leavea_i_exp,
					this,
					_1, t, r0 ),
			   this->MAX_ALPHA_SEQ ) );
    return p;
}

const Real 
FirstPassagePairGreensFunction::p_int_r( const Real r,
					 const Real t,
					 const Real r0 ) const

{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_int_r_i_exp,
					this,
					_1, t, r, r0 ),
			   this->MAX_ALPHA_SEQ ) );
    return p;
}

const Real 
FirstPassagePairGreensFunction::
p_int_r_table( const Real r,
	       const Real t,
	       const Real r0,
	       const RealVector& num_r0Table ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_int_r_i_exp_table,
					this,
					_1, t, r, r0, num_r0Table ),
			   this->MAX_ALPHA_SEQ ) );
    return p;
}


const Real
FirstPassagePairGreensFunction::p_survival_F( const Real t,
					      const p_survival_params* params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real r0( params->r0 );
    const RealVector& psurvTable( params->psurvTable );
    const Real rnd( params->rnd );

    return rnd - gf->p_survival_table( t, r0, psurvTable );
}


const Real
FirstPassagePairGreensFunction::p_int_r_F( const Real r,
					   const p_int_r_params* params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real t( params->t );
    const Real r0( params->r0 );
    const RealVector& num_r0Table( params->num_r0Table );
    const Real rnd( params->rnd );

    return gf->p_int_r_table( r, t, r0, num_r0Table ) - rnd;
}


const Real FirstPassagePairGreensFunction::drawTime( const Real rnd, 
						     const Real r0 ) const
{
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    THROW_UNLESS( std::invalid_argument, rnd <= 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= sigma && r0 <= a );

    if( r0 == a || a == sigma )
    {
	return 0.0;
    }

    Real low( 1e-5 );
    Real high( 1.0 );

    this->updateAlphaTable0( high );

    RealVector psurvTable;
    this->createPsurvTable( psurvTable, r0 );

    p_survival_params params = { this, r0, psurvTable, rnd };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_survival_F ),
	    &params 
	};

    // adjust high and low to make sure that f( low ) and f( high ) straddle.
    while( GSL_FN_EVAL( &F, high ) < 0.0 )
    {
	high *= 10;
	printf( "drawTime: adjusting high: %g\n", high );
	if( fabs( high ) >= 1e10 )
	{
	    std::cerr << "Couldn't adjust high. F(" << high <<
		") = " << GSL_FN_EVAL( &F, high ) << "; r0 = " << r0 << 
		", " << dump() << std::endl;
	    throw std::exception();
	}
//	this->updateAlphaTable0( high );
    }

    this->updateAlphaTable0( low );
    this->createPsurvTable( psurvTable, r0 );

    while( GSL_FN_EVAL( &F, low ) > 0.0 )
    {
	low *= .1;
	printf( "drawTime: adjusting low: %g\n",low );

	if( fabs( low ) <= this->MIN_T )
	{
	    std::cerr << "Couldn't adjust low.  Returning MIN_T (= "
		      << this->MIN_T << "); F(" << low <<
		") = " << GSL_FN_EVAL( &F, low ) << "; r0 = " << r0 << ", "
		      << dump() << std::endl;
	    return this->MIN_T;
	}
	this->updateAlphaTable0( low );
	this->createPsurvTable( psurvTable, r0 );
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

	const int status( gsl_root_test_interval( low, high, this->MIN_T, 
						  this->TOLERANCE ) );

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

    return t;
}

const EventType 
FirstPassagePairGreensFunction::drawEventType( const Real rnd, 
					       const Real r0,
					       const Real t ) const
{
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    THROW_UNLESS( std::invalid_argument, rnd <= 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 > sigma && r0 < a );
    THROW_UNLESS( std::invalid_argument, t > 0.0 );

    // leaves() and leavea() call updateAlphaTable0(). redundant?
    const Real reaction( leaves( t, r0 ) * 
			 4.0 * M_PI * sigma * sigma );
    const Real escape( leavea( t, r0 ) *
		       4.0 * M_PI * a * a );

    const Real value( reaction / ( reaction + escape ) );

    if( rnd < value )  
    {
	return REACTION;   // leaves
    }
    else 
    {
	return ESCAPE;     // leavea
    }
}


const Real FirstPassagePairGreensFunction::drawR( const Real rnd, 
						  const Real r0, 
						  const Real t ) const
{
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    THROW_UNLESS( std::invalid_argument, rnd <= 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= sigma && r0 < a );

    if( t == 0.0 )
    {
	return r0;
    }

    const Real psurv( p_survival( t, r0 ) );

    RealVector num_r0Table;
    createNum_r0Table( num_r0Table, r0 );

    p_int_r_params params = { this, t, r0, num_r0Table, rnd * psurv };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_int_r_F ),
	    &params 
	};

    Real low( sigma );
    Real high( a );

    //const Real lowvalue( GSL_FN_EVAL( &F, low  ) );
    const Real highvalue( GSL_FN_EVAL( &F, high ) );

    if( highvalue < 0.0 )
    {
	//printf( "drawR: highvalue < 0.0 (%g). returning a.\n", highvalue );
	return a;
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
	const int status( gsl_root_test_interval( low, high, 1e-15, 
						  this->TOLERANCE ) );

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


    const Real r( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );

    return r;
}




const Real FirstPassagePairGreensFunction::p_n_alpha( const unsigned int i,
						      const unsigned int n,
						      const Real r,
						      const Real r0, 
						      const Real t ) const
{
    const Real Dt( this->getD() * t );
    const Real sigma( this->getSigma() );
    const Real h( this->geth() );

    const Real alpha( this->getAlpha( n, i ) );
    const Real alphasq( alpha * alpha );

    const Real aAlpha( a * alpha );
    const Real sigmaAlpha( sigma * alpha );
    const Real hSigma( geth() * getSigma() );
    const Real realn( static_cast<Real>( n ) );
    const Real hSigma_m_n( hSigma - realn );

    const Real term1( alphasq * exp( - Dt * alphasq ) );

#if 0
    const Real np( realn + 0.5 );
    Real jas1, yas1, jas2, yas2, jaa, yaa, jar, yar, jar0, yar0, _;
    bessjy( sigmaAlpha, np,       &jas1, &yas1, &_, &_ );
    bessjy( sigmaAlpha, np + 1.0, &jas2, &yas2, &_, &_ );
    bessjy( aAlpha,     np,       &jaa,  &yaa,  &_, &_ );
    bessjy( r * alpha,  np,       &jar,  &yar,  &_, &_ );
    bessjy( r0 * alpha, np,       &jar0, &yar0, &_, &_ );
#else
    // GSL
    const Real factors( sqrt( sigmaAlpha * M_2_PI ) );
    const Real factora( sqrt( aAlpha * M_2_PI ) );
    const Real factorr( sqrt( r * alpha * M_2_PI ) );
    const Real factorr0( sqrt( r0 * alpha * M_2_PI ) );
    const Real jas1( gsl_sf_bessel_jl( n,   sigmaAlpha ) * factors );
    const Real yas1( gsl_sf_bessel_yl( n,   sigmaAlpha ) * factors );
    const Real jas2( gsl_sf_bessel_jl( n+1, sigmaAlpha ) * factors );
    const Real yas2( gsl_sf_bessel_yl( n+1, sigmaAlpha ) * factors );
    const Real jaa(  gsl_sf_bessel_jl( n,   aAlpha )     * factora );
    const Real yaa(  gsl_sf_bessel_yl( n,   aAlpha )     * factora );
    const Real jar(  gsl_sf_bessel_jl( n,   r * alpha )  * factorr );
    const Real yar(  gsl_sf_bessel_yl( n,   r * alpha )  * factorr );
    const Real jar0( gsl_sf_bessel_jl( n,   r0 * alpha ) * factorr0 );
    const Real yar0( gsl_sf_bessel_yl( n,   r0 * alpha ) * factorr0 );
#endif

    const Real J( hSigma_m_n * jas1 + sigmaAlpha * jas2 );
    const Real Y( hSigma_m_n * yas1 + sigmaAlpha * yas2 );
    const Real falpha_r( - J * yar + Y * jar );
    const Real falpha_r0( - J * yar0 + Y * jar0 );

    const Real num( falpha_r * falpha_r0 );

    const Real E1( realn + realn * realn - 
		   sigma * ( h + h * h * sigma + sigma * alphasq ) );

    const Real E2( ( J * J + Y * Y ) / ( jaa * jaa + yaa * yaa ) );

    // Following is simpler, but jaa sometimes gets irregularly huge,
    // perhaps because of some error in alpha, so E2 above is safer for now.
    // const Real E2a( J * J / ( jaa * jaa ) );

    const Real den( E1 + E2 );

    const Real result( term1 * num / den );

    return result;
}



const Real 
FirstPassagePairGreensFunction::p_n( const Integer n,
				     const Real r,
				     const Real r0, 
				     const Real t ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_n_alpha,
					this,
					_1, n, r, r0, t ),
			   this->MAX_ALPHA_SEQ ) );

    const Real factor( ( 1 + 2 * n ) * M_PI / ( 8.0 * sqrt( r * r0 ) ) );

    return p * factor;
}

void
FirstPassagePairGreensFunction::makep_nTable( RealVector& p_nTable,
					      const Real r, 
					      const Real r0, 
					      const Real t ) const
{
    p_nTable.clear();

    const Real p_0( this->p_n( 0, r, r0, t ) );
    p_nTable.push_back( p_0 );

    const Real threshold( fabs( this->TOLERANCE * p_0 ) );

    Real p_n_prev_abs( fabs( p_0 ) );
    unsigned int n( 1 );
    while( true )
    {
	Real p_n( this->p_n( n, r, r0, t ) );

	if( ! std::isnormal( p_n ) )
	{
	    std::cerr << "makep_nTable: invalid value; " <<
		p_n << "( n= " << n << ")." << std::endl;
//	    p_n = 0.0;
	    break;
	}
//	printf("%d p_n %g\n", n, p_n );

	p_nTable.push_back( p_n );

	const Real p_n_abs( fabs( p_n ) );
	// truncate when converged enough.
	if( p_n_abs < threshold &&
	    p_n_abs < p_n_prev_abs )
	{
	    break;
	}
	
	++n;

	if( n >= this->MAX_ORDER )
	{
	    // std::cerr << "p_n didn't converge." << std::endl;
	    break;
	}
	
	p_n_prev_abs = p_n_abs;
    }

}


const Real 
FirstPassagePairGreensFunction::dp_n_alpha_at_a( const unsigned int i,
						 const unsigned int n,
						 const Real r0, 
						 const Real t ) const
{
    const Real mDt( - this->getD() * t );
    const Real sigma( this->getSigma() );
    const Real h( this->geth() );
    const Real alpha( this->getAlpha( n, i ) );

    const Real alphasq( alpha * alpha );

    const Real aAlpha( a * alpha );
    const Real sigmaAlpha( sigma * alpha );
    const Real hSigma( geth() * getSigma() );
    const Real realn( static_cast<Real>( n ) );
    const Real hSigma_m_n( hSigma - realn );

    const Real term1( alphasq * exp( mDt * alphasq ) );

#if 0  
    const Real np( realn + 0.5 );
    Real jas1, yas1, jas2, yas2, jaa1, yaa1, jaa2, yaa2, jar0, yar0, _;
    bessjy( sigmaAlpha, np,       &jas1, &yas1, &_, &_ );
    bessjy( sigmaAlpha, np + 1.0, &jas2, &yas2, &_, &_ );
    bessjy( aAlpha,     np,       &jaa1, &yaa1, &_, &_ );
    bessjy( r0 * alpha, np,       &jar0, &yar0, &_, &_ );
#else
    // GSL
    const Real factors( sqrt( sigmaAlpha * M_2_PI ) );
    const Real factora( sqrt( aAlpha * M_2_PI ) );
    const Real factorr0( sqrt( r0 * alpha * M_2_PI ) );
    const Real jas1( gsl_sf_bessel_jl( n,   sigmaAlpha ) * factors );
    const Real yas1( gsl_sf_bessel_yl( n,   sigmaAlpha ) * factors );
    const Real jas2( gsl_sf_bessel_jl( n+1, sigmaAlpha ) * factors );
    const Real yas2( gsl_sf_bessel_yl( n+1, sigmaAlpha ) * factors );
    const Real jaa1( gsl_sf_bessel_jl( n,   aAlpha )     * factora );
    const Real yaa1( gsl_sf_bessel_yl( n,   aAlpha )     * factora );
    const Real jar0( gsl_sf_bessel_jl( n,   r0 * alpha ) * factorr0 );
    const Real yar0( gsl_sf_bessel_yl( n,   r0 * alpha ) * factorr0 );
#endif

    const Real J( hSigma_m_n * jas1 + sigmaAlpha * jas2 );
    const Real Y( hSigma_m_n * yas1 + sigmaAlpha * yas2 );

    const Real dfalpha_r( - 2.0 * J / ( a * M_PI * jaa1 ) );
    const Real falpha_r0( - J * yar0 + Y * jar0 );

    const Real num( dfalpha_r * falpha_r0 );

    const Real E1( realn + realn * realn - 
		   sigma * ( h + h * h * sigma + sigma * alphasq ) );

    const Real E2( ( J * J + Y * Y ) / 
		   ( jaa1 * jaa1 + yaa1 * yaa1 ) );

//    const Real E2( ( J * J ) / 
//		   ( jaa1 * jaa1 ) );

    const Real den( E1 + E2 );

//    printf("f n1 n2 d %g %g %g %g\n",falpha_r0, num1, num2, den );

    const Real result( term1 * num / den );

    return result;
}

const Real 
FirstPassagePairGreensFunction::dp_n_at_a( const Integer n,
					   const Real r0, 
					   const Real t ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					dp_n_alpha_at_a,
					this,
					_1, n, r0, t ),
			   this->MAX_ALPHA_SEQ ) );

    const Real factor( getD() * ( 1 + 2 * n ) * M_PI / 
		       ( 8.0 * sqrt( geta() * r0 ) ) );


    return p * factor;
}


void
FirstPassagePairGreensFunction::makedp_n_at_aTable( RealVector& p_nTable,
						    const Real r0, 
						    const Real t ) const
{
    p_nTable.clear();

    const Real p_0( this->dp_n_at_a( 0, r0, t ) );
    p_nTable.push_back( p_0 );

    const Real threshold( fabs( this->TOLERANCE * p_0 ) );

    Real p_n_prev_abs( fabs( p_0 ) );
    unsigned int n( 1 );
    while( true )
    {
	Real p_n( this->dp_n_at_a( n, r0, t ) );

	if( ! std::isnormal( p_n ) )
	{
	    std::cerr << "makedp_n_at_aTable: invalid value; " <<
		p_n << "( n= " << n << ")." << std::endl;
//	    p_n = 0.0;
	    break;
	}
//	printf("p_n %g\n",p_n );

	p_nTable.push_back( p_n );

	const Real p_n_abs( fabs( p_n ) );
	// truncate when converged enough.
	if( p_n_abs < threshold &&
	    p_n_abs < p_n_prev_abs )
	{
	    break;
	}
	
	++n;

	if( n >= this->MAX_ORDER )
	{
	    std::cerr << "dp_n_at_a didn't converge." << std::endl;
	    break;
	}
	
	p_n_prev_abs = p_n_abs;
    }

}

const Real 
FirstPassagePairGreensFunction::p_theta( const Real theta,
					 const Real r, 
					 const Real r0, 
					 const Real t ) const 
{
    {
	const Real sigma( this->getSigma() );
	const Real a( this->geta() );
	
	THROW_UNLESS( std::invalid_argument, theta >= 0.0 && theta <= M_PI );
        // r \in ( sigma, a );  not defined at r == sigma and r == a.
	THROW_UNLESS( std::invalid_argument, r > sigma && r < a );
	THROW_UNLESS( std::invalid_argument, r0 > sigma && r0 < a );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    }

    if( t == 0.0 )
    {
	return 0.0;
    }

    RealVector p_nTable;

    makep_nTable( p_nTable, r, r0, t );

    const Real p( p_theta_table( theta, r, r0, t, p_nTable ) );

    return p;
}


const Real 
FirstPassagePairGreensFunction::dp_theta( const Real theta,
                                          const Real r, 
                                          const Real r0, 
                                          const Real t ) const 
{
    {
	const Real sigma( this->getSigma() );
	const Real a( this->geta() );
	
	THROW_UNLESS( std::invalid_argument, theta >= 0.0 && theta <= M_PI );

        // r \in [ sigma, a ]  ;  unlike p_theta,
        // defined at r == sigma and r == a.
	THROW_UNLESS( std::invalid_argument, r >= sigma && r <= a );
	THROW_UNLESS( std::invalid_argument, r0 > sigma && r0 < a );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    }

    if( t == 0.0 )
    {
	return 0.0;
    }

    RealVector p_nTable;

    makedp_n_at_aTable( p_nTable, r0, t );

    const Real p( p_theta_table( theta, r, r0, t, p_nTable ) );

    return p;
}


const Real 
FirstPassagePairGreensFunction::
p_theta_i( const unsigned int n,
	   const RealVector& p_nTable, const RealVector& lgndTable ) const
{
    return p_nTable[n] * lgndTable[n];
}

const Real
FirstPassagePairGreensFunction::
p_theta_table( const Real theta,
	       const Real r, 
	       const Real r0, 
	       const Real t, 
	       const RealVector& p_nTable ) const
{
    const unsigned int tableSize( p_nTable.size() );

    RealVector lgndTable( tableSize );

    Real sin_theta;
    Real cos_theta;
    sincos( theta, &sin_theta, &cos_theta );
    gsl_sf_legendre_Pl_array( tableSize-1, cos_theta, &lgndTable[0] );

    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					p_theta_i,
					this,
					_1, p_nTable, lgndTable ),
			   tableSize-1 ) );

    return p * sin_theta;
}


void
FirstPassagePairGreensFunction::
make_p_thetaTable( RealVector& pTable,
		   const Real r, 
		   const Real r0, 
		   const Real t,
		   const unsigned int n,
		   const RealVector& p_nTable ) const
{
    const Real thetaStep( M_PI / n );

    pTable.push_back( 0.0 );

    Real p_prev( 0.0 );
    unsigned int i( 1 );
    while( true )
    {
	const Real theta( thetaStep * i );

	Real p( this->p_theta_table( theta, r, r0, t, p_nTable ) );

	if( p < 0.0 )
	{
	    printf("drawTheta: p<0 %g\n", p );
	    p = 0.0;
	}

	const Real value( ( p_prev + p ) * 0.5 );
	pTable.push_back( *( pTable.end() - 1 ) + value );

//	printf("p %g %g %g\n", theta, pTable[i], p );

	if( //value < pTable[i] * std::numeric_limits<Real>::epsilon() ||
	    i >= n - 1 )
	{
	    break;   // pTable is valid in [0,i].
	}

	p_prev = p;
	++i;
    }

}


const Real 
FirstPassagePairGreensFunction::ip_theta( const Real theta,
					  const Real r, 
					  const Real r0, 
					  const Real t ) const
{
    {
	const Real sigma( this->getSigma() );
	const Real a( this->geta() );
	
	THROW_UNLESS( std::invalid_argument, theta >= 0.0 && theta <= M_PI );
        // r \in ( sigma, a )
	THROW_UNLESS( std::invalid_argument, r > sigma && r < a );
	THROW_UNLESS( std::invalid_argument, r0 > sigma && r0 < a );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    }

    if( t == 0.0 || theta == 0.0 )
    {
	return 0.0;
    }

    RealVector p_nTable;

    makep_nTable( p_nTable, r, r0, t );

    const Real p( ip_theta_table( theta, r, r0, t, p_nTable ) );

    return p;
}


const Real 
FirstPassagePairGreensFunction::idp_theta( const Real theta,
                                           const Real r, 
                                           const Real r0, 
                                           const Real t ) const
{
    {
	const Real sigma( this->getSigma() );
	const Real a( this->geta() );
	
	THROW_UNLESS( std::invalid_argument, theta >= 0.0 && theta <= M_PI );
        // r \in [ sigma, a ]
	THROW_UNLESS( std::invalid_argument, r >= sigma && r <= a );
	THROW_UNLESS( std::invalid_argument, r0 > sigma && r0 < a );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    }

    if( t == 0.0 || theta == 0.0 )
    {
	return 0.0;
    }

    RealVector p_nTable;

    makedp_n_at_aTable( p_nTable, r0, t );

    const Real p( ip_theta_table( theta, r, r0, t, p_nTable ) );

    return p;
}

const Real 
FirstPassagePairGreensFunction::
ip_theta_i( const unsigned int n,
	    const RealVector& p_nTable, 
	    const RealVector& lgndTable1 ) const
{
    // lgndTable1 is offset by 1; lgndTable1[0] is for n=-1.

    const Real lgnd_n_m1( lgndTable1[n] );   // n-1
    const Real lgnd_n_p1( lgndTable1[n+2] ); // n+1
    
    return p_nTable[n] * ( lgnd_n_m1 - lgnd_n_p1 ) / ( 1.0 + 2.0 * n );
}


const Real 
FirstPassagePairGreensFunction::
ip_theta_table( const Real theta,
		const Real r, 
		const Real r0, 
		const Real t,	 
		const RealVector& p_nTable ) const
{
    const unsigned int tableSize( p_nTable.size() );

    RealVector pTable;
    pTable.reserve( tableSize );

    const Real cos_theta( cos( theta ) );

    // LgndTable is offset by 1 to incorporate the n=-1 case.
    // For ex: LgndTable[0] is for n=-1, lgndTable[1] is for n=0 ...

    RealVector lgndTable1( tableSize + 2 );
    lgndTable1[0] = 1.0;  // n = -1
    gsl_sf_legendre_Pl_array( tableSize, cos_theta, &lgndTable1[1] );


    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction::
					ip_theta_i,
					this,
					_1, p_nTable, lgndTable1 ),
			   tableSize-1 ) );

    return p;
}

const Real
FirstPassagePairGreensFunction::ip_theta_F( const Real theta,
					    const ip_theta_params* params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real r( params->r );
    const Real r0( params->r0 );
    const Real t( params->t );
    const RealVector& p_nTable( params->p_nTable );
    const Real value( params->value );

    return gf->ip_theta_table( theta, r, r0, t, p_nTable ) - value;
}



const Real 
FirstPassagePairGreensFunction::drawTheta( const Real rnd,
					   const Real r, 
					   const Real r0, 
					   const Real t ) const
{
    Real theta;

    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    // input parameter range checks.
    THROW_UNLESS( std::invalid_argument, rnd <= 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= sigma && r0 < a );
    THROW_UNLESS( std::invalid_argument, r >= sigma && r <= a );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    // t == 0 means no move.
    if( t == 0.0 )
    {
	return 0.0;
    }

    RealVector p_nTable;

    if( r == geta() || r < 0.0 )
    {
	puts("dp");
	makedp_n_at_aTable( p_nTable, r0, t );
    }
    else
    {
	makep_nTable( p_nTable, r, r0, t );
    }

#if 1
    // root finding with the integrand form.

    const Real ip_theta_pi( ip_theta_table( M_PI, r, r0, t, p_nTable ) );

    ip_theta_params params = { this, r, r0, t, p_nTable, rnd * ip_theta_pi };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &ip_theta_F ),
	    &params 
	};

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, 0.0, M_PI );

    const unsigned int maxIter( 100 );

    unsigned int i( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );
	const Real low( gsl_root_fsolver_x_lower( solver ) );
	const Real high( gsl_root_fsolver_x_upper( solver ) );
	const int status( gsl_root_test_interval( low, high, 1e-15, 
						  this->TOLERANCE ) );

	if( status == GSL_CONTINUE )
	{
	    if( i >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "drawTheta: failed to converge." << std::endl;
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

    theta = gsl_root_fsolver_root( solver );
    gsl_root_fsolver_free( solver );
    

#else
    // bisection after numerical integration.

    const RealVector::size_type thetaResolution( 200 );
    const Real thetaStep( M_PI / thetaResolution );

    RealVector pTable;
    pTable.reserve( thetaResolution );

    make_p_thetaTable( pTable, r, r0, t, thetaResolution, p_nTable );
    const RealVector::size_type tableLast( pTable.size()-1 );

    // debug
    /*
    const Real psurv( p_survival( t, r0 ) );
    const Real p0r( p_0( t, r, r0 ) * 4.0 * M_PI * r * r );
    const Real ip( ip_theta_table( M_PI, r, r0, t, p_nTable ) );
    printf("theta %g %g %g\n", p0r, pTable[tableLast] * 
           2 * M_PI * r * r * thetaStep, ip * 2 * M_PI * r * r );
    */

    const Real targetPoint( rnd * pTable[tableLast] );
    const size_t lowerBound( gsl_interp_bsearch( &pTable[0], targetPoint, 
						 0, tableLast ) );
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
#endif


    return theta;
}


//
// debug
//

const std::string FirstPassagePairGreensFunction::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << ", sigma = " << this->getSigma() <<
	", a = " << this->geta() <<
	", kf = " << this->getkf() <<
	", h = " << this->geth() << std::endl;
    return ss.str();
}    
