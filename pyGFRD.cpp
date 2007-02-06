#include <exception>

#include <gsl/gsl_math.h>

#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <numpy/arrayobject.h>

#include "DynamicPriorityQueue.hpp"

class PyEvent
{
  
public:
  
  PyEvent( const double time )
    :
    time( time )
  {
    ; // do nothing
  }
  
  void setTime( const double time )
  {
    this->time = time;
  }
  
  const double getTime() const
  {
    return this->time;
  }
  
  
  
  const bool operator< ( const PyEvent& rhs ) const
  {
    if( getTime() < rhs.getTime() )
      {
	return true;
      }
    else
      {
	return false;
      }
  }
  
  const bool operator!= ( const PyEvent& rhs ) const
  {
    if( getTime() == rhs.getTime() )
      {
	return false;
      }
    else
      {
	return true;
      }
  }
  
  
  // dummy, because DynamicPriorityQueue requires this. better without.
  PyEvent()
  {
    ; // do nothing
  }
  
  
private:
  
  double             time;
};


typedef DynamicPriorityQueue<PyEvent> PyEventDynamicPriorityQueue;


#include "PairGreensFunction.hpp"
#include "PlainPairGreensFunction.hpp"

#include "FirstPassageGreensFunction.hpp"

using namespace boost::python;

const double distanceSq( const double* const p1, const double* const p2 )
{
    return gsl_pow_2( p1[0] - p2[0] ) 
	+ gsl_pow_2( p1[1] - p2[1] ) 
	+ gsl_pow_2( p1[2] - p2[2] );
}

const double distanceSq_( PyArrayObject* o1, PyArrayObject* o2 )
{
    // check type, dimension, size.
    if( o1->descr->type_num != PyArray_DOUBLE || 
	o2->descr->type_num != PyArray_DOUBLE ||
	o1->nd != 1 || o2->nd != 1 ||
	o1->dimensions[0] != 3 || o2->dimensions[0] != 3 )
    {
	throw std::exception();
    }

    const double* const 
	p1( reinterpret_cast<const double* const>( o1->data ) );
    const double* const 
	p2( reinterpret_cast<const double* const>( o2->data ) );
  
    return distanceSq( p1, p2 );
}


void* extract_pyarray(PyObject* x)
{
    return PyObject_TypeCheck(x, &PyArray_Type) ? x : 0;
}


BOOST_PYTHON_MODULE( _gfrd )
{
    import_array();

    converter::registry::insert( &extract_pyarray, type_id<PyArrayObject>());


    class_<PlainPairGreensFunction>( "PlainPairGreensFunction",
				     init<const Real, 
				     const Real, 
				     const Real>() )
	.def( "getD", &PlainPairGreensFunction::getD )
	.def( "getkf", &PlainPairGreensFunction::getkf )
	.def( "getSigma", &PlainPairGreensFunction::getSigma )
	.def( "drawTime", &PlainPairGreensFunction::drawTime )
	.def( "drawR", &PlainPairGreensFunction::drawR )
	.def( "drawTheta", &PlainPairGreensFunction::drawTheta )
	;

    class_<FirstPassageGreensFunction>( "FirstPassageGreensFunction",
					init<const Real>() )
	.def( "getD", &FirstPassageGreensFunction::getD )
	.def( "drawTime", &FirstPassageGreensFunction::drawTime )
	.def( "drawR", &FirstPassageGreensFunction::drawR )
	.def( "p_survival", &FirstPassageGreensFunction::p_survival )
	.def( "p_r_int", &FirstPassageGreensFunction::p_r_int )
	.def( "p_r_fourier", &FirstPassageGreensFunction::p_r_fourier )
	;

    def( "distanceSq", &distanceSq_ );

}
