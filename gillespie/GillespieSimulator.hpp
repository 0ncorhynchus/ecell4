#ifndef __GILLESPIESIMULATOR_HPP
#define __GILLESPIESIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Simulator.hpp>

#include "GillespieWorld.hpp"

namespace ecell4 
{

namespace gillespie 
{

class GillespieSimulator 
	: 
		public Simulator 
{
public:
	GillespieSimulator(
		boost::shared_ptr<NetworkModel> model, boost::shared_ptr<GillespieWorld> world,
		RandomNumberGenerator &rng)
		: model_(model), world_(world), rng_(rng)
	{
		this->num_steps_ = 0;
		this->initialize();	// calucate the time the first reaction occurs.
	}
		
	Integer num_steps(void) const;
	void step(void) ;
	bool step(Real const & upto);

	Real t(void) const;
	void set_t(Real const &t);
	Real dt(void) const;

	void initialize(void);	// re-calcurate the next reaction.
	RandomNumberGenerator &rng(void);
					
protected:
	boost::shared_ptr<NetworkModel> model_;
	boost::shared_ptr<GillespieWorld> world_;
	
	Integer num_steps_;
	RandomNumberGenerator &rng_;

	Real dt_;	
	int next_reaction_num_; 	// the index of the next reaction.
	void calc_next_reaction_(void);
};

}

}	// ecell4

#endif //__GILLESPIESIMULATOR_HPP
