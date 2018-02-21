// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbnet_energy/HBNetEnergy.cc
/// @brief An EnergyMethod that gives a bonus for hydrogen bond networks, which ramps nonlinearly with the size of the
/// networks.
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.  It has also been modified to permit sub-regions of a pose to be scored.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

// Unit headers
#include <core/scoring/hbnet_energy/HBNetEnergy.hh>
#include <core/scoring/hbnet_energy/HBNetEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/chemical/AtomType.hh>

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
#include <utility/graph/Graph.hh>
#include <utility/graph/graph_util.hh>

namespace core {
namespace scoring {
namespace hbnet_energy {

// The square of the maximum hydrogen bond length (2.5 A).
#define HBNETENERGY_HBOND_DIST_CUTOFF_SQ 6.25

static basic::Tracer TR("core.scoring.hbnet_energy.HBNetEnergy");

/// @brief This must return a fresh instance of the HBNetEnergy class, never an instance already in use.
///
core::scoring::methods::EnergyMethodOP
HBNetEnergyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const &options ) const
{
	return core::scoring::methods::EnergyMethodOP( new HBNetEnergy( options ) );
}

/// @brief Defines the score types that this energy method calculates.
///
ScoreTypes
HBNetEnergyCreator::score_types_for_method() const
{
	ScoreTypes sts;
	sts.push_back( hbnet );
	return sts;
}

/// @brief Options constructor.
///
HBNetEnergy::HBNetEnergy ( core::scoring::methods::EnergyMethodOptions const &options ) :
	parent1( core::scoring::methods::EnergyMethodCreatorOP( new HBNetEnergyCreator ) ),
	parent2( ),
	neighbour_graph_(),
	hbonds_graph_(),
	hbonds_graph_last_accepted_(),
	bb_hbonds_graph_(),
	symm_info_(),
	ramping_type_(HBNetEnergyRampQuadratic),
	max_network_size_(options.hbnet_max_network_size())
{
	set_hbnet_energy_ramping( options.hbnet_bonus_function_ramping() );
}

/// @brief Default destructor.
///
HBNetEnergy::~HBNetEnergy() = default;

/// @brief Clone: create a copy of this object, and return an owning pointer
/// to the copy.
core::scoring::methods::EnergyMethodOP HBNetEnergy::clone() const {
	return core::scoring::methods::EnergyMethodOP( new HBNetEnergy(*this) );
}

/// @brief HBNetEnergy is context-independent and thus indicates that no context graphs need to be maintained by
/// class Energies.
void HBNetEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const
{
	//Do nothing.
	return;
}

/// @brief HBNetEnergy is version 1.0 right now.
///
core::Size HBNetEnergy::version() const
{
	return 1; // Initial versioning
}

/// @brief Actually calculate the total energy
/// @details Called by the scoring machinery.  The update_residue_neighbors() function of the pose
/// must be called first.
void HBNetEnergy::finalize_total_energy(
	core::pose::Pose & pose,
	ScoreFunction const & /*sfxn*/,
	EnergyMap & totals
) const {
	initialize_graphs(pose);

	core::Size nres(pose.total_residue());

	utility::vector1< core::conformation::ResidueCOP > reslist( nres );

	for ( core::Size ir(1); ir<=nres; ++ir ) {
		reslist[ir] = pose.conformation().residue_cop(ir);
	}

	totals[ hbnet ] += calculate_energy(reslist);
}

/// @brief Calculate the total energy given a vector of const owning pointers to residues.
/// @details Called directly by the ResidueArrayAnnealingEvaluator during packer runs.  Requires
/// that set_up_residuearrayannealablenergy_for_packing() be called first.
core::Real
HBNetEnergy::calculate_energy(
	utility::vector1< core::conformation::ResidueCOP > const & resvect,
	core::Size const substitution_position /*= 0*/
) const {

	if ( substitution_position == 0 ) { //If no substitution position is given, iterate over all positions

		//Clear the edges in the hbonds graph:
		hbonds_graph_.drop_all_edges();

		//Iterate over all edges in the neighbour graph.
		for (  utility::graph::EdgeListConstIterator edgeiter( neighbour_graph_.edge_list_begin() ); edgeiter != neighbour_graph_.edge_list_end(); ++edgeiter ) {
			core::Size const res1( (*edgeiter)->get_first_node_ind() ), res2( (*edgeiter)->get_second_node_ind() );
			if ( has_hbond( *(resvect[res1]), *(resvect[res2]), res1, res2 ) || has_hbond( *(resvect[res2]), *(resvect[res1]), res2, res1 ) ) {
				if ( symm_info_ != nullptr ) { symmetrize_hbonds_graph( res1, res2, resvect ); }
				else { hbonds_graph_.add_edge(res1, res2); }
			}
		}
		hbonds_graph_last_accepted_ = hbonds_graph_;
	} else { //If a substitution position is given, only clear edges to that node and repopulate those edges.
		hbonds_graph_ = hbonds_graph_last_accepted_;

		//Clear the edges to node "substitution_position":
		hbonds_graph_.drop_all_edges_for_node( substitution_position );
		if ( symm_info_ != nullptr ) { drop_all_edges_for_symmetric_nodes( substitution_position ); }

		//Iterate through all edges connected to the substitution position node:
		for ( utility::graph::EdgeListConstIterator edgeiter( neighbour_graph_.get_node( substitution_position )->edge_list_begin() ); edgeiter != neighbour_graph_.get_node( substitution_position)->edge_list_end(); ++edgeiter ) {
			core::Size const res2( (*edgeiter)->get_first_node_ind() != substitution_position ? (*edgeiter)->get_first_node_ind() : (*edgeiter)->get_second_node_ind() );
			if ( has_hbond( *(resvect[substitution_position]), *(resvect[res2]), substitution_position, res2 ) || has_hbond( *(resvect[res2]), *(resvect[substitution_position]), res2, substitution_position ) ) {
				if ( symm_info_ != nullptr ) { symmetrize_hbonds_graph( substitution_position, res2, resvect ); }
				else { hbonds_graph_.add_edge(substitution_position, res2); }
			}
		}
	}

	//Count the sizes of all islands ("connected components") in the hbonds graph:
	utility::vector1< std::pair< core::Size, core::Size > > const connected_components( utility::graph::find_connected_components( hbonds_graph_ ) ); //First entry is size of island, second entry is a representative vertex in the island.

	//Score each island nonlinearly.  Accumulate a total score.
	core::Real accumulator(0.0);
	for ( core::Size i(1), imax(connected_components.size()); i<=imax; ++i ) { //Iterate over all islands.
		if ( connected_components[i].second > 1 ) accumulator -= bonus_function( (connected_components[i]).second );
	}

	return accumulator;
}

/// @brief What to do when a substitution that was considered is accepted.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
HBNetEnergy::commit_considered_substitution() {
	hbonds_graph_last_accepted_ = hbonds_graph_;
}

/// @brief Get a summary of all loaded data.
///
void HBNetEnergy::report() const {
	if ( !TR.Debug.visible() ) return; //Do nothing if I don't have a tracer.

	TR.Debug << std::endl << "The HBNetEnergy object has no loaded data (aside from graphs that are generated and updated on-the-fly)." << std::endl;
}

/// @brief Cache data from the pose in this EnergyMethod in anticipation of scoring.
///
void
HBNetEnergy::set_up_residuearrayannealableenergy_for_packing (
	core::pose::Pose const & pose,
	core::pack::rotamer_set::RotamerSets const &/*rotamersets*/,
	core::scoring::ScoreFunction const & /*sfxn*/
) {
	initialize_graphs(pose);
}

/// @brief Set the way that HBNetEnergy scales with network size,
/// by string.
void
HBNetEnergy::set_hbnet_energy_ramping (
	std::string const &ramping_string
) {
	HBNetEnergyRamping const ramping_enum( hbnet_energy_ramping_enum_from_string( ramping_string ) );
	runtime_assert_string_msg( ramping_enum != HBNetEnergyRampINVALID, "Error in HBNetEnergy::set_hbnet_energy_ramping(): \"" + ramping_string + "\" is not a valid ramping type." );
	set_hbnet_energy_ramping( ramping_enum );
}

/// @brief Set the way that HBNetEnergy scales with network size,
/// by enum.
void
HBNetEnergy::set_hbnet_energy_ramping (
	HBNetEnergyRamping const ramping_enum
) {
	runtime_assert_string_msg( ramping_enum != HBNetEnergyRampINVALID, "Error in HBNetEnergy::set_hbnet_energy_ramping(): an invalid ramping type was passed to this function." );
	ramping_type_ = ramping_enum;
}

/// @brief Given a string for an HBNetEnergyRamping type, return the corresponding enum.
/// @details Returns HBNetEnergyRampINVALID if the string isn't recognized.
HBNetEnergyRamping
HBNetEnergy::hbnet_energy_ramping_enum_from_string(
	std::string const &ramping_string
) const {
	for ( core::Size i(1); i<static_cast<core::Size>( HBNetEnergyRamp_end_of_list ); ++i ) {
		if ( !ramping_string.compare( hbnet_energy_ramping_string_from_enum( static_cast< HBNetEnergyRamping >(i) ) ) ) {
			return static_cast< HBNetEnergyRamping >(i);
		}
	}
	return HBNetEnergyRampINVALID;
}

/// @brief Given an enum for an HBNetEnergyRamping type, return the corresponding string.
/// @details Returns "INVALID" if the string isn't recognized.
std::string
HBNetEnergy::hbnet_energy_ramping_string_from_enum(
	HBNetEnergyRamping const ramping_enum
) const {
	switch(ramping_enum) {
	case HBNetEnergyRampQuadratic :
		return "quadratic";
	case HBNetEnergyRampLinear :
		return "linear";
	case HBNetEnergyRampLogarithmic :
		return "logarithmic";
	case HBNetEnergyRampSquareRoot :
		return "squareroot";
	default :
		return "INVALID";
	}
}

/// @brief Set the maximum network size, beyond which there is no bonus for making a network bigger.
/// @details A value of "0" (the default) means no max.
void HBNetEnergy::max_network_size( core::Size const setting ) { max_network_size_ = setting; }

//////////////////////////////////PRIVATE FUNCTIONS//////////////////////////////////////

/// @brief Initializes the neighbor_graph_ and hbonds_graph_ objects, given a pose.
///
void
HBNetEnergy::initialize_graphs(
	core::pose::Pose const &pose
) const {
	//The neighbour graph is initialized from the pose.  Note that this assumes that the pose
	//has an up-to-date neighbor graph.
	neighbour_graph_ = pose.energies().tenA_neighbor_graph(); //This is mutable.

	//The hydrogen bonds graph starts out with a node for every residue, but no hydrogen bonds.
	hbonds_graph_.drop_all_edges(); //This is also mutable.
	hbonds_graph_.set_num_nodes( pose.total_residue() );

	initialize_bb_hbonds_graph( pose );

	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::conformation::symmetry::SymmetricConformationCOP symmconf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( pose.conformation_ptr() ) );
		debug_assert( symmconf != nullptr );
		symm_info_ = symmconf->Symmetry_Info();
	} else {
		symm_info_ = nullptr;
	}

}

/// @brief Set up the bb_hbonds_graph_ object.
/// @details If X and Y share an edge, then X is bb-bb hydrogen bonded to Y, or to a residue that is covalently bonded to Y,
/// or a residue covalently bonded to X is bb-bb hydrogen bonded to Y.
void
HBNetEnergy::initialize_bb_hbonds_graph(
	core::pose::Pose const &pose
) const {
	core::Size const nres(pose.total_residue());
	bb_hbonds_graph_.drop_all_edges();
	bb_hbonds_graph_.set_num_nodes(nres);

	//Iterate over all edges in the neighbour graph.
	for (  utility::graph::EdgeListConstIterator edgeiter( neighbour_graph_.edge_list_begin() ); edgeiter != neighbour_graph_.edge_list_end(); ++edgeiter ) {
		core::Size const res1( (*edgeiter)->get_first_node_ind() ), res2( (*edgeiter)->get_second_node_ind() );
		if ( has_bb_hbond_connection( pose, pose.residue(res1), pose.residue(res2) ) || has_bb_hbond_connection( pose, pose.residue(res2), pose.residue(res1) ) ) {
			bb_hbonds_graph_.add_edge(res1, res2);
		}
	}

}


/// @brief Determines whether a hydrogen bond exists between res1 and res2.
/// @details Ignores backbone-backbone hydrogen bonds.  Also, this is directional: it considers donors in res1 to acceptors in res2.
/// Call twice for bidirectional hbonding.
bool
HBNetEnergy::has_hbond(
	core::conformation::Residue const &res1,
	core::conformation::Residue const &res2,
	core::Size const res1_index,
	core::Size const res2_index
) const {
	for ( core::Size ia(1), iamax(res1.natoms()); ia<=iamax; ++ia ) {
		if ( !res1.atom_type(ia).is_polar_hydrogen() || !res1.atom_type( res1.atom_base(ia) ).is_donor() ) continue;
		bool const one_is_bb( res1.atom_base(ia) < res1.first_sidechain_atom() );
		if ( one_is_bb && bb_hbonds_graph_.find_edge(res1_index, res2_index) != nullptr ) continue; //If either of the atoms is backbone and this is a closely backbone-connected pair, don't allow an hbond edge
		for ( core::Size ja(1), jamax(res2.natoms()); ja<=jamax; ++ja ) {
			bool const two_is_bb( ja < res2.first_sidechain_atom() );
			if ( one_is_bb && two_is_bb ) continue; //Ignore backbone-backbone.
			if ( two_is_bb && bb_hbonds_graph_.find_edge(res1_index, res2_index) != nullptr ) continue; //If either of the atoms is backbone and this is a closely backbone-connected pair, don't allow an hbond edge

			if ( !res2.atom_type(ja).is_acceptor() ) continue;
			if ( res1.xyz(ia).distance_squared( res2.xyz(ja) ) < HBNETENERGY_HBOND_DIST_CUTOFF_SQ ) return true;
		}
	}
	return false;
}

/// @brief Determines whether a pair of residues has a close backbone hydrogen bond connection.
/// @details Returns true if (a) the residues have a bb-bb hbond, (b) res1 has a bb-bb hbond to a residue that is covalently bonded
/// to res2, (c) res2 has a bb-bb hbond to a residue that is covalently bonded to res1, or (d) if res1 and res2 are directly covalently bonded.
/// @note This only considers res1 (and its neighbors) as donor and res2 (and its neighbors) as acceptor.  Call twice to check two-way hydrogen bonds.
bool
HBNetEnergy::has_bb_hbond_connection(
	core::pose::Pose const &pose,
	core::conformation::Residue const &res1,
	core::conformation::Residue const &res2
) const {

	//Case 1a: direct covalent bond
	if ( res1.has_lower_connect() && res1.connected_residue_at_lower() == res2.seqpos() ) return true;
	if ( res1.has_upper_connect() && res1.connected_residue_at_upper() == res2.seqpos() ) return true;

	//Case 1b: direct hydrogen bond connection
	if ( has_bb_hbond(res1, res2) ) return true;

	// Case 2a: res1 i-1 connected to res2.
	if ( res1.has_lower_connect() ) {
		core::Size const res1_lower( res1.connected_residue_at_lower() );
		if ( res1_lower != 0 && has_bb_hbond(pose.residue(res1_lower), res2) ) return true;
	}

	// Case 2b: res1 i+1 connected to res2.
	if ( res1.has_upper_connect() ) {
		core::Size const res1_upper( res1.connected_residue_at_upper() );
		if ( res1_upper != 0 && has_bb_hbond(pose.residue(res1_upper), res2) ) return true;
	}

	// Case 3a: res2 i-1 connected to res1.
	if ( res2.has_lower_connect() ) {
		core::Size const res2_lower( res2.connected_residue_at_lower() );
		if ( res2_lower != 0 && has_bb_hbond( res1, pose.residue(res2_lower) ) ) return true;
	}

	// Case 3a: res2 1+1 connected to res1.
	if ( res2.has_upper_connect() ) {
		core::Size const res2_upper( res2.connected_residue_at_upper() );
		if ( res2_upper != 0 && has_bb_hbond( res1, pose.residue(res2_upper) ) ) return true;
	}

	return false;
}


/// @brief Determines whether a pair of residues has a direct backbone hydrogen bond connection.
/// @details Returns true if the residues have a bb-bb hbond.
/// @note This only considers res1 as donor and res2 as acceptor.  Call twice to check two-way hydrogen bonds.
bool
HBNetEnergy::has_bb_hbond(
	core::conformation::Residue const &res1,
	core::conformation::Residue const &res2
) const {
	for ( core::Size ia(1), iamax(res1.natoms()); ia<=iamax; ++ia ) {
		if ( !res1.atom_type(ia).is_polar_hydrogen() || !res1.atom_type( res1.atom_base(ia) ).is_donor() ) continue;
		if ( res1.atom_base(ia) >= res1.first_sidechain_atom() ) continue; //The first is not a backbone hbond donor.
		for ( core::Size ja(1), jamax(res2.first_sidechain_atom()); ja<jamax; ++ja ) {
			if ( !res2.atom_type(ja).is_acceptor() ) continue;
			if ( res1.xyz(ia).distance_squared( res2.xyz(ja) ) < HBNETENERGY_HBOND_DIST_CUTOFF_SQ ) return true;
		}
	}
	return false;
}


/// @brief Given the size of a connected component in a graph, return a value to pass to the bonus function accumulator.
/// @details The value returned by this function is more POSITIVE for bigger bonuses; it should be SUBTRACTED from the accumulator
/// to yield a score that gets more negative with greater favourability.
core::Real
HBNetEnergy::bonus_function(
	core::Size count_in
) const {
	if ( max_network_size_ > 0  && count_in > max_network_size_ ) count_in = max_network_size_;

	switch( ramping_type_ ) {
	case HBNetEnergyRampQuadratic :
		return static_cast<core::Real>( count_in*count_in );
	case HBNetEnergyRampLinear :
		return static_cast< core::Real >( count_in );
	case HBNetEnergyRampLogarithmic :
		return std::log( static_cast<core::Real>( count_in + 1 ) ); /*Plus one to get around ln(0) being undefined.*/
	case HBNetEnergyRampSquareRoot :
		return std::sqrt( static_cast<core::Real>( count_in ) );
	default :
		utility_exit_with_message("Error in HBNetEnergy::bonus_function(): The energy ramping was incorrectly set!");
	}
	return 0; //To keep compiler happy.
}

/// @brief Given the index of an asymmetric node in the hbonds_graph_ object, drop all of the edges for all corresponding symmetric nodes.
/// @details Assumes symm_info_ points to something.  Check symm_info != nullptr before calling this function.
void
HBNetEnergy::drop_all_edges_for_symmetric_nodes(
	core::Size const asymm_node
) const {
	debug_assert( hbonds_graph_.num_nodes() >= asymm_node && asymm_node > 0 );
	debug_assert( symm_info_ != nullptr );
	debug_assert( symm_info_->bb_is_independent( asymm_node ) );

	for ( core::Size i(1), imax(symm_info_->num_total_residues()); i<=imax; ++i ) { //Loop through all residues
		if ( symm_info_->bb_follows( i ) == asymm_node ) {
			hbonds_graph_.drop_all_edges_for_node( i );
		}
	}
}

/// @brief Given residue indices node1 and node2 in the asymmetric unit, add an edge between the corresponding nodes in all symmetric units.
/// @details Assumes symm_info_ points to something.  Check symm_info != nullptr before calling this function.
void
HBNetEnergy::symmetrize_hbonds_graph(
	core::Size const node1,
	core::Size const node2,
	utility::vector1< core::conformation::ResidueCOP > const & resvect
) const {
	core::Size const nres(resvect.size());
	debug_assert(symm_info_ != nullptr);
	debug_assert( nres >= node1 && node1 > 0 );
	debug_assert( nres >= node2 && node2 > 0 );

	//The equivalent residues in the asymmetric unit:
	core::Size const masternode1( symm_info_->bb_is_independent(node1) ? node1 : symm_info_->bb_follows(node1) );
	core::Size const masternode2( symm_info_->bb_is_independent(node2) ? node2 : symm_info_->bb_follows(node2) );

	core::conformation::ResidueCOP node1res( resvect[node1] ), node2res( resvect[node2] );
	numeric::xyzVector< core::Real > xform( node2res->xyz(1) - node1res->xyz(1) );
	core::Real const lensq( xform.length_squared() );

	//Get all the copies of node1 and node2:
	utility::vector1< core::Size > node1_copies(symm_info_->bb_clones(masternode1)), node2_copies(symm_info_->bb_clones(masternode2));
	node1_copies.push_back(masternode1);
	node2_copies.push_back(masternode2);

	for ( core::Size i(1), imax(node1_copies.size()); i<=imax; ++i ) {
		for ( core::Size j(1), jmax(node2_copies.size()); j<=jmax; ++j ) {
			numeric::xyzVector< core::Real > xform2( resvect[node2_copies[j]]->xyz(1) - resvect[node1_copies[i]]->xyz(1) );
			if ( std::abs( xform2.length_squared() - lensq ) < 1e-4 ) {
				if ( !hbonds_graph_.find_edge(node1_copies[i], node2_copies[j]) ) {
					hbonds_graph_.add_edge( node1_copies[i], node2_copies[j] );
				}
			}
		}
	}


}


} // hbnet_energy
} // scoring
} // core
