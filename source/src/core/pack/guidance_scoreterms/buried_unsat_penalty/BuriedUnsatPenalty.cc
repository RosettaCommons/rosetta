// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/guidance_scoreterms/buried_unsat_penalty/BuriedUnsatPenalty.cc
/// @brief An EnergyMethod that gives a penalty for buried unsatisfied hydrogen bond donors and acceptors.
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

// Unit headers
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/BuriedUnsatPenalty.hh>
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/BuriedUnsatPenaltyCreator.hh>
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraph.hh>
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphContainer.hh>
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphOptions.hh>

// Package headers
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/chemical/AtomType.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/numbers.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/memory.hh>
#include <utility/graph/Graph.hh>
#include <utility/graph/graph_util.hh>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace buried_unsat_penalty {

static basic::Tracer TR("core.pack.guidance_scoreterms.buried_unsat_penalty.BuriedUnsatPenalty");

/// @brief This must return a fresh instance of the BuriedUnsatPenalty class, never an instance already in use.
///
core::scoring::methods::EnergyMethodOP
BuriedUnsatPenaltyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const &options ) const
{
	return core::scoring::methods::EnergyMethodOP( utility::pointer::make_shared< BuriedUnsatPenalty >( options ) );
}

/// @brief Defines the score types that this energy method calculates.
///
core::scoring::ScoreTypes
BuriedUnsatPenaltyCreator::score_types_for_method() const
{
	core::scoring::ScoreTypes sts;
	sts.push_back( core::scoring::buried_unsatisfied_penalty );
	return sts;
}

/// @brief Options constructor.
///
BuriedUnsatPenalty::BuriedUnsatPenalty ( core::scoring::methods::EnergyMethodOptions const &options ) :
	parent1( core::scoring::methods::EnergyMethodCreatorOP( utility::pointer::make_shared< BuriedUnsatPenaltyCreator >() ) ),
	parent2( ),
	disabled_(false),
	unsat_graph_(nullptr),
	graph_options_( utility::pointer::make_shared< graph::BuriedUnsatPenaltyGraphOptions > (
	options.buried_unsatisfied_penalty_cone_angle_exponent(),
	options.buried_unsatisfied_penalty_cone_angle_shift_factor(),
	options.buried_unsatisfied_penalty_cone_dist_exponent(),
	options.buried_unsatisfied_penalty_cone_dist_midpoint(),
	options.buried_unsatisfied_penalty_burial_threshold(),
	options.buried_unsatisfied_penalty_hbond_energy_threshold()
	)
	),
	hbond_options_( utility::pointer::make_shared< core::scoring::hbonds::HBondOptions >( options.hbond_options() ) ),
	symmetric_(false),
	nres_(0),
	num_symmetric_copies_(0)
{}

/// @brief Default destructor.
///
BuriedUnsatPenalty::~BuriedUnsatPenalty() {}

/// @brief Clone: create a copy of this object, and return an owning pointer
/// to the copy.
core::scoring::methods::EnergyMethodOP BuriedUnsatPenalty::clone() const {
	return core::scoring::methods::EnergyMethodOP( utility::pointer::make_shared< BuriedUnsatPenalty >(*this) );
}

/// @brief BuriedUnsatPenalty is context-independent and thus indicates that no context graphs need to be maintained by
/// class Energies.
void BuriedUnsatPenalty::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const
{
	//Do nothing.
	return;
}

/// @brief BuriedUnsatPenalty is version 1.0 right now.
///
core::Size BuriedUnsatPenalty::version() const
{
	return 1; // Initial versioning
}

/// @brief Actually calculate the total energy.
/// @details Called by the scoring machinery.  The update_residue_neighbors() function of the pose
/// must be called first.
void BuriedUnsatPenalty::finalize_total_energy(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & /*sfxn*/,
	core::scoring::EnergyMap & totals
) const {
	if ( disabled_ ) return; //Do nothing during minimization trajectories.
	totals[ core::scoring::buried_unsatisfied_penalty ] += calculate_penalty_once_from_scratch( pose );
}

/// @brief Calculate the total energy given a vector of const owning pointers to residues.
/// @details Called directly by the ResidueArrayAnnealingEvaluator during packer runs.  Requires
/// that set_up_residuearrayannealablenergy_for_packing() be called first.
/// @note Outside of the context of the packer, this doesn't behave as expected (and finalize_total_energy() should
/// be called instead).
core::Real
BuriedUnsatPenalty::calculate_energy(
	utility::vector1< core::conformation::ResidueCOP > const &resvect,
	core::Size const //substitution_position
) const {

	if ( disabled_ ) return 0.0; //Do nothing during minimization trajectories.

	core::Size unsat_acceptor_count_lastconsidered = 0;
	core::Size unsat_donor_count_lastconsidered = 0;
	core::Size unsat_acceptor_and_donor_count_lastconsidered = 0;
	core::Size oversat_acceptor_count_lastconsidered = 0;
	core::Size oversat_donor_count_lastconsidered = 0;
	core::Size oversat_acceptor_and_donor_count_lastconsidered = 0;

	debug_assert( resvect.size() >= nres_ );
	graph::BuriedUnsatPenaltyGraphOP current_state_unsat_graph( utility::pointer::make_shared< graph::BuriedUnsatPenaltyGraph >( nres_, graph_options_, hbond_options_ ) );
	current_state_unsat_graph->set_always_rotamer_one(true);
	for ( core::Size i(1), imax( current_state_unsat_graph->num_nodes()); i<=imax; ++i ) {
		current_state_unsat_graph->copy_node_and_connected_edges( i, *unsat_graph_, unsat_graph_->get_node_index( resvect[i] ) );
	}

	current_state_unsat_graph->compute_unsats_all_nodes();
	for ( core::Size i(1), imax( current_state_unsat_graph->num_nodes() ); i<=imax; ++i ) {
		graph::BuriedUnsatPenaltyNode const &curnode( *(static_cast< graph::BuriedUnsatPenaltyNode const * >( current_state_unsat_graph->get_node(i) )) );
		curnode.increment_counts( unsat_acceptor_count_lastconsidered, unsat_donor_count_lastconsidered, unsat_acceptor_and_donor_count_lastconsidered, oversat_acceptor_count_lastconsidered, oversat_donor_count_lastconsidered, oversat_acceptor_and_donor_count_lastconsidered);
	}

	if ( symmetric_ ) {
		unsat_acceptor_count_lastconsidered *= num_symmetric_copies_;
		unsat_donor_count_lastconsidered *= num_symmetric_copies_;
		unsat_acceptor_and_donor_count_lastconsidered *= num_symmetric_copies_;
		oversat_acceptor_count_lastconsidered *= num_symmetric_copies_;
		oversat_donor_count_lastconsidered *= num_symmetric_copies_;
		oversat_acceptor_and_donor_count_lastconsidered *= num_symmetric_copies_;
	}

	return compute_penalty( unsat_acceptor_count_lastconsidered, unsat_donor_count_lastconsidered, unsat_acceptor_and_donor_count_lastconsidered, oversat_acceptor_count_lastconsidered, oversat_donor_count_lastconsidered, oversat_acceptor_and_donor_count_lastconsidered );
}

/// @brief What to do when a substitution that was considered is accepted.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
BuriedUnsatPenalty::commit_considered_substitution() {
	//GNDN
}

/// @brief Get a summary of all loaded data.
///
void BuriedUnsatPenalty::report() const {
	if ( !TR.Debug.visible() ) return; //Do nothing if I don't have a tracer.

	TR.Debug << std::endl << "The BuriedUnsatPenalty object has no loaded data (aside from graphs that are generated and updated on-the-fly)." << std::endl;
}

/// @brief Cache data from the pose in this EnergyMethod in anticipation of scoring.
///
void
BuriedUnsatPenalty::set_up_residuearrayannealableenergy_for_packing (
	core::pose::Pose & pose,
	core::pack::rotamer_set::RotamerSets const &rotamersets,
	core::scoring::ScoreFunction const & /*sfxn*/
) {
	TR << "Setting up BuriedUnsatPenalty energy method for packing." << std::endl;

	// Store information for symmetry:
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		symmetric_ = true;
		core::conformation::symmetry::SymmetricConformationCOP symmconf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( pose.conformation_ptr() ) );
		debug_assert( symmconf != nullptr );
		core::conformation::symmetry::SymmetryInfoCOP symminfo( symmconf->Symmetry_Info() );
		debug_assert( symminfo != nullptr );
		nres_ = symminfo->num_independent_residues();
		num_symmetric_copies_ = symminfo->num_bb_clones() + 1;
	} else {
		symmetric_ = false;
		nres_ = pose.total_residue();
		num_symmetric_copies_ = 1;
	}

	if ( !pose.energies().data().has( core::scoring::EnergiesCacheableDataType::BURIED_UNSAT_HBOND_GRAPH ) ) {
		TR << "Creating new BuriedUnsatPenaltyGraph and caching it in the pose." << std::endl;
		graph::BuriedUnsatPenaltyGraphOP unsat_graph( utility::pointer::make_shared< graph::BuriedUnsatPenaltyGraph >( graph_options_, hbond_options_ ) );
		unsat_graph->initialize_graph_for_packing( pose, rotamersets );
		graph::BuriedUnsatPenaltyGraphContainerOP unsat_graph_container( utility::pointer::make_shared< graph::BuriedUnsatPenaltyGraphContainer >( unsat_graph ) );
		pose.energies().data().set( core::scoring::EnergiesCacheableDataType::BURIED_UNSAT_HBOND_GRAPH, unsat_graph_container );
		unsat_graph_ = unsat_graph;
	} else {
		TR << "Fetching BuriedUnsatPenaltyGraph from cache in the pose." << std::endl;
		unsat_graph_ = (utility::pointer::static_pointer_cast< graph::BuriedUnsatPenaltyGraphContainer const >( pose.energies().data().get_ptr( core::scoring::EnergiesCacheableDataType::BURIED_UNSAT_HBOND_GRAPH ) ))->graph();
	}

	debug_assert( unsat_graph_ != nullptr );
	TR << "Successfully set up BuriedUnsatPenalty energy method for packing." << std::endl;
	TR.flush();
}

/// @brief Disable this scoreterm during minimization trajectory.
void
BuriedUnsatPenalty::setup_for_minimizing(
	core::pose::Pose &,//pose,
	core::scoring::ScoreFunction const & ,//sfxn,
	core::kinematics::MinimizerMapBase const & //minmap
) const {
	TR << "Disabling BuriedUnsatPenalty scoreterm during minimization trajectory." << std::endl;
	disabled_ = true;
}

/// @brief Re-enable this scoreterm after a minimization trajectory.
void
BuriedUnsatPenalty::finalize_after_minimizing(
	core::pose::Pose &// pose
) const {
	TR << "Re-enabling BuriedUnsatPenalty scoreterm after minimization trajectory." << std::endl;
	disabled_ = false;
}

//////////////////////////////////Public member functions specific to this energy method//////////////////////////////////////

/// @brief Given the counts of various unsaturateds, return a penalty value.
core::Real
BuriedUnsatPenalty::compute_penalty(
	core::Size const unsat_acceptor_count,
	core::Size const unsat_donor_count,
	core::Size const unsat_acceptor_and_donor_count,
	core::Size const oversat_acceptor_count,
	core::Size const oversat_donor_count,
	core::Size const oversat_acceptor_and_donor_count
) const {
	//The following could be computed in different, more complex ways in the future:
	return static_cast< core::Real >( std::pow( 5 * ( unsat_acceptor_count + unsat_donor_count + unsat_acceptor_and_donor_count + oversat_acceptor_count + oversat_donor_count + oversat_acceptor_and_donor_count ), 2 ) );
}


/// @brief Given a pose, calculate the penalty energy.
core::Real
BuriedUnsatPenalty::calculate_penalty_once_from_scratch(
	core::pose::Pose const &pose
) const {

	graph::BuriedUnsatPenaltyGraphOP curstate_graph( utility::pointer::make_shared< graph::BuriedUnsatPenaltyGraph >(graph_options_, hbond_options_ ) );
	curstate_graph->initialize_graph_for_scoring( pose );

	curstate_graph->compute_unsats_all_nodes();

	core::Size unsat_acceptor_count(0);
	core::Size unsat_donor_count(0);
	core::Size unsat_acceptor_and_donor_count(0); //The above two don't count hydroxyls, which must EITHER donate OR accept (though both is fine, too).
	core::Size oversat_acceptor_count(0);
	core::Size oversat_donor_count(0);
	core::Size oversat_acceptor_and_donor_count(0); //Again, hydroxyls.

#ifndef NDEBUG
	//Check for debug mode.
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformationCOP symmconf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( pose.conformation_ptr() ) );
		debug_assert( symmconf != nullptr );
		core::conformation::symmetry::SymmetryInfoCOP symminfo( symmconf->Symmetry_Info() );
		debug_assert( symminfo != nullptr );
		debug_assert( curstate_graph->num_nodes() == symminfo->num_independent_residues() );
	} else {
		debug_assert( curstate_graph->num_nodes() == pose.total_residue() );
	}
#endif

	for ( core::Size i(1), imax(curstate_graph->num_nodes()); i<=imax; ++i ) {
		graph::BuriedUnsatPenaltyNode const & curnode( *( static_cast< graph::BuriedUnsatPenaltyNode const * >( curstate_graph->get_node( i ) ) ) );
		curnode.increment_counts(unsat_acceptor_count, unsat_donor_count, unsat_acceptor_and_donor_count, oversat_acceptor_count, oversat_donor_count, oversat_acceptor_and_donor_count );
	}

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetricConformationCOP symmconf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( pose.conformation_ptr() ) );
		debug_assert( symmconf != nullptr );
		core::conformation::symmetry::SymmetryInfoCOP symminfo( symmconf->Symmetry_Info() );
		debug_assert( symminfo != nullptr );
		core::Size const symm_multiplier( symminfo->num_bb_clones() + 1 );
		unsat_acceptor_count *= symm_multiplier;
		unsat_donor_count *= symm_multiplier;
		unsat_acceptor_and_donor_count *= symm_multiplier;
		oversat_acceptor_count *= symm_multiplier;
		oversat_donor_count *= symm_multiplier;
		oversat_acceptor_and_donor_count *= symm_multiplier;
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Computing penalty once from scratch.  Found unsat_acceptor=" << unsat_acceptor_count << " unsat_donor=" << unsat_donor_count << " unsat_acceptor_and_donor=" << unsat_acceptor_and_donor_count << " oversat_acceptor=" << oversat_acceptor_count << " oversat_donor=" << oversat_donor_count << " oversat_acceptor_and_donor=" << oversat_acceptor_and_donor_count << std::endl;
	}

	return compute_penalty( unsat_acceptor_count, unsat_donor_count, unsat_acceptor_and_donor_count, oversat_acceptor_count, oversat_donor_count, oversat_acceptor_and_donor_count );
}

/// @brief Provide Pymol commands to colour the pose grey, non-buried donor and acceptor groups cyan, and buried acceptor
/// and donor groups orange.  Useful for debugging degree of burial.
/// @details To use, pass in a pose.  If this graph contains residues corresponding to those in the pose, commands for colouring
/// them will be written out.  Calls BuriedUnsatPenaltyGraph::provide_pymol_commands_to_show_groups().  Must be called only after
/// set_up_residuearrayannealableenergy_for_packing.
void
BuriedUnsatPenalty::provide_pymol_commands_to_show_groups(
	std::ostream &out,
	core::pose::Pose const &pose
) const {
	runtime_assert( unsat_graph_ != nullptr );
	unsat_graph_->provide_pymol_commands_to_show_groups(out, pose);
}

//////////////////////////////////PRIVATE FUNCTIONS//////////////////////////////////////


} // buried_unsat_penalty
} // guidance_scoreterms
} // pack
} // core
