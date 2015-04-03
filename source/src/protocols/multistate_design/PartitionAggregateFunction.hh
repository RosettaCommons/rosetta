// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PartitionAggregateFunction.hh
/// @brief
/// @author Colin A. Smith

#ifndef INCLUDED_protocols_multistate_design_PartitionAggregateFunction_hh
#define INCLUDED_protocols_multistate_design_PartitionAggregateFunction_hh


#include <protocols/multistate_design/MultiStateFitnessFunction.fwd.hh>

#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <protocols/multistate_design/MultiStateAggregateFunction.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace multistate_design {

class PartitionAggregateFunction : public MultiStateAggregateFunction {

public:
	typedef utility::pointer::shared_ptr< PartitionAggregateFunction > OP;
	typedef utility::pointer::shared_ptr< PartitionAggregateFunction const > COP;

	PartitionAggregateFunction();
	PartitionAggregateFunction(core::Real temp, core::Real anchor_offset, bool const compare_to_ground_state=false);
	virtual ~PartitionAggregateFunction();

	/*
	PartitionAggregateFunction() : MultiStateAggregateFunction<T>(), temp_(1), anchor_offset_(0), compare_all_to_ground_state_( false ) {}
	PartitionAggregateFunction(core::Real const temp, core::Real const anchor_offset, bool const compare_to_ground_state = false ) :
		MultiStateAggregateFunction<T>(), temp_(temp), anchor_offset_(anchor_offset), compare_all_to_ground_state_( compare_to_ground_state ) {}
*/

	virtual core::Real temp() const;
	virtual void set_temp( core::Real temp );

	virtual core::Real anchor_offset() const;
	virtual void set_anchor_offset( core::Real offset );

	virtual
	core::Real
	evaluate(
		utility::vector1<core::Real> const & single_state_fitnesses,
		MultiStateFitnessFunction & fitness_function
	) const;

private:
	core::Real temp_;
	core::Real anchor_offset_;

	// if states are based on different starting structures, their ground states are different.
	// in this case, the energies are taken as differences from the ground state.
	bool const compare_all_to_ground_state_; // const to make sure that it doesn't change during a run
};

	/*

/// SJF see https://wiki.rosettacommons.org/index.php/IDDocumentation#ProteinInterfaceMS
/// for explanations on the parameters
template <typename T>
core::Real
PartitionAggregateFunction<T>::evaluate(
	utility::vector1<core::Real> const & single_state_fitnesses,
	MultiStateFitnessFunction<T> & fitness_function
) const
{
	using namespace ObjexxFCL::format;

	utility::vector1<SingleStateCOP> single_states(fitness_function.const_states());
	runtime_assert(single_state_fitnesses.size() == single_states.size());

/// SJF the temperature plays the role of scaling the relative contributions of states' energies to
/// fitness. At low temperatures, even mild energy differences between positive and negative states will
/// translate to large fitness gains.

	core::Real const inv_temp( -1.0 / temp_ );
	core::Real numer(0.), denom(0.);

	// 'normalize' is just a constant value to subtract from energies of large magnitude, in order to take exponents of
	// smaller numbers (exact value used should not be important for the calculations)
	core::Real const normalize( compare_all_to_ground_state_ ? 0 : single_state_fitnesses[ 1 ] );
	for (core::Size i = 1; i <= single_state_fitnesses.size(); ++i) {
/// SJF ground_state_offset changes the energy of each design as used by the fitness function by the value of the best-score energy of the single state. This means that designs will be judged by the differences in energy that they imply relative to the starting 'best-score' design. This is useful if you have multiple (e.g., homologous) pdb files of the starting structures, each with a different starting energy. ground_state_offset overrides normalize and has a different meaning from it.
		core::Real const ground_state_offset( compare_all_to_ground_state_ ? single_states[ i ]->best_score() : 0 );
		TR(basic::t_trace) << "State fitness " << F(8,2,single_state_fitnesses[i]);
		TR(basic::t_trace) << "State ground-state offset " << F(8,2,ground_state_offset);
		core::Real const exp_term( std::exp( ( single_state_fitnesses[i] - ground_state_offset - normalize ) * inv_temp ) );
		TR(basic::t_trace) << " exp. term " << F(6,2,exp_term);
		denom += exp_term;
		if ( single_states[i]->is_positive_state() ) {
			TR(basic::t_trace) << " (POSITIVE STATE)";
			numer += exp_term;
/// SJF The affinity anchor is not strictly part of the multi-state design framework. The idea
/// is to have some 'memory' of the energy of the best sequence for the target state. If you make large gains
/// in
/// specificity but the target state's energy is losing in a big way, the anchor_term will make the
/// fitness suffer. Intuitively, anchor_offset_ is a way to state what is an acceptable hit in stability,
/// beyond which fitness will suffer.
/// This is important of course in specificity designs, where large gains in specificity can come at
/// a cost to stability. In other scenarios, such as binding, an appropriate anchor_offset_ value might be
/// approximately the difference between the positive and negative states. This way, the energy of the
/// negative state will make decisive contributions to the fitness.
			// also add 'affinity anchor(s)' to denominator for each positive state
			core::Real const anchor_term( std::exp( ( single_states[i]->best_score() - ground_state_offset + anchor_offset_ - normalize ) * inv_temp ) );
			TR(basic::t_trace) << " anchor exp. term " << F(6,2,anchor_term);
			denom += anchor_term;
		}//is positive state
		TR(basic::t_trace) << std::endl;
	}//for i
	TR(basic::t_trace) << "numer " << F(5,2,numer) << " / denom " << F(5,2,denom) << std::endl;
	// flip sign on the Boltzmann probability for proper ranking elsewhere (more negative is better)
	core::Real const prob_target( numer/denom );
	TR(basic::t_debug) << "Boltzmann prob. for target state(s) vs. competitor(s): "
	                        << ObjexxFCL::format::F(5,2,prob_target) << std::endl;

	return -1.*prob_target; // in genetic algorithm, better fitnesses are more negative
}

template <typename T>
basic::Tracer
PartitionAggregateFunction<T>::TR("protocols.multistate_design.PartitionAggregateFunction");
*/

} // namespace multistate_design
} // namespace protocols

#endif
