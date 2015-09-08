// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/IntermolEnergy.cc
/// @brief  Cost of bringing two chains together.
/// @author Rhiju Das


// Unit headers
#include <core/scoring/methods/IntermolEnergy.hh>
#include <core/scoring/methods/IntermolEnergyCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/scoring/EnergyMap.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

// C++
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "core.scoring.methods.IntermolEnergy" );

/////////////////////////////////////////////////////////////////////////////////////
//
// Created in attempt to fit 'Turner rules' for METHODS.
//  assumes 1 M standard state -- later will allow
//  specification of strand concentration(s) from
//  command line.
//
/////////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the IntermolEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
IntermolEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new IntermolEnergy );
}

ScoreTypes
IntermolEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( intermol );
	return sts;
}


/// c-tor
IntermolEnergy::IntermolEnergy() :
	parent( methods::EnergyMethodCreatorOP( new IntermolEnergyCreator ) ),
	penalty_at_1M_( 2.30 ), // calibrated outside.
	log_conc_( std::log( basic::options::option[ basic::options::OptionKeys::score::conc ]() ) ) // in kT
{}

/// clone
methods::EnergyMethodOP
IntermolEnergy::clone() const
{
	return methods::EnergyMethodOP( new IntermolEnergy );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
void
IntermolEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {

	using namespace core::pose::full_model_info;
	make_sure_full_model_info_is_setup( pose );
	Size const num_chains_frozen = get_num_chains_frozen( pose );
	totals[ intermol ] = num_chains_frozen * ( penalty_at_1M_ - log_conc_ );
} // finalize_total_energy


///////////////////////////////////////////////////////////////////////////////
// quick graph traversal lets us figure out number of connected components.
// this may repeat work in LoopGraph and LoopClose Energy, in which case unify!
Size
IntermolEnergy::get_num_chains_frozen( pose::Pose const & pose ) const {
	using namespace core::pose::full_model_info;
	Size const num_chains = const_full_model_info( pose ).cutpoint_open_in_full_model().size() + 1;

	utility::vector1< std::pair< Size, Size > > chain_connections = get_chain_connections( pose );
	utility::vector1< Size > connection_domains = get_connection_domains( chain_connections, num_chains  );

	// there used to be an alternative way to calculate # chains frozen, but I discovered
	// a bug in it & deprecated it -- rhiju.
	// runtime_assert( num_subgraphs == max( connection_domains ) );

	return ( num_chains - max( connection_domains ) );
}

///////////////////////////////////////////////////////////////////////////////
void
IntermolEnergy::eval_atom_derivative(
	id::AtomID const &,
	pose::Pose const &,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const &,
	Vector &,
	Vector &
) const
{
	// no op.
} // eval atom derivative

core::Size
IntermolEnergy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core
