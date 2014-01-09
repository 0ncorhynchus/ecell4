#define BOOST_TEST_MODULE "EventScheduler_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../EventScheduler.hpp"

using namespace ecell4;
using namespace ecell4::lattice;

BOOST_AUTO_TEST_CASE(EventScheduler_test_constructor)
{
    EventScheduler scheduler();
}
