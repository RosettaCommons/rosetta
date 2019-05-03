// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/hbnet/NetworkState.hh
/// @brief specialized structs for MC HBNet to store hbond network info
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_hbnet_NetworkState_hh
#define INCLUDED_protocols_hbnet_NetworkState_hh

#include <core/conformation/Residue.hh>
#include <core/scoring/hbonds/graph/HBondInfo.hh>
#include <core/scoring/hbonds/graph/AtomInfo.hh>
#include <core/scoring/hbonds/graph/HBondGraph.hh>
#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
//#include <core/pack/rotamer_set/RotamerSets.hh>

namespace protocols {
namespace hbnet {

class NetworkState {

public:

	NetworkState(
		core::scoring::hbonds::graph::HBondEdge const * monte_carlo_seed,
		core::scoring::hbonds::graph::HBondGraphOP const & hbond_graph
	);

	~NetworkState(){}

	bool operator < ( NetworkState const & rhs) const {
		return full_twobody_energy_ < rhs.full_twobody_energy_;
	}

	///@brief NetworkState keeps track of polar atoms that need to be satisfied.
	/// This function allows NetworkState to copy the vector of polar atoms from the HBondNode
	void add_polar_atoms(
		core::scoring::hbonds::graph::HBondNode const * node
	);

	///@brief are there any buried unsats at this
	///molten residue position that are not participating in network hbonds?
	bool mres_has_unsats( unsigned int mres ) const;

	///@brief are there any heavy (not hydrogen) buried unsats at this
	///molten residue position that are not participating in network hbonds?
	bool mres_has_heavy_unsats ( unsigned int mres ) const;

public:
	boost::container::flat_map< unsigned int,
	core::scoring::hbonds::graph::AtomInfoSet >
	const & unsatisfied_sc_atoms_const() const {
		return unsatisfied_sc_atoms_;
	}

	///@brief get vector of AtomInfos for atoms that are unsatisfied
	core::scoring::hbonds::graph::AtomInfoSet *
	get_unsats_for_mres( unsigned int mres );

	///@brief get vector of AtomInfos for atoms that are unsatisfied
	core::scoring::hbonds::graph::AtomInfoSet const *
	get_unsats_for_mres( unsigned int mres ) const;

	utility::vector1< core::scoring::hbonds::graph::HBondNode const * > & nodes() {
		return nodes_;
	}

	utility::vector1< core::scoring::hbonds::graph::HBondNode const * > const & nodes() const {
		return nodes_;
	}

	void add_node( core::scoring::hbonds::graph::HBondNode const * node ) {
		nodes_.push_back( node );
	}

	utility::vector1< core::scoring::hbonds::graph::HBondEdge const * > & edges() {
		return edges_;
	}

	utility::vector1< core::scoring::hbonds::graph::HBondEdge const * > const & edges() const {
		return edges_;
	}

	void add_edge( core::scoring::hbonds::graph::HBondEdge const * edge ) {
		edges_.push_back( edge );
	}

	core::scoring::hbonds::graph::HBondEdge const * monte_carlo_seed() const {
		return monte_carlo_seed_;
	}

	core::Real full_twobody_energy() const {
		return full_twobody_energy_;
	}

	void full_twobody_energy( core::Real setting ){
		full_twobody_energy_ = setting;
	}

	void append_twobody_energy( core::Real setting ){
		full_twobody_energy_ += setting;
	}

	///@brief This holds whatever metric is used for sorting (getter)
	core::Real score() const {
		return score_;
	}

	///@brief This holds whatever metric is used for sorting (setter)
	void score( core::Real setting ) {
		score_ = setting;
	}

private:
	utility::vector1< core::scoring::hbonds::graph::HBondNode const * > nodes_;
	utility::vector1< core::scoring::hbonds::graph::HBondEdge const * > edges_;

	core::scoring::hbonds::graph::HBondEdge const * monte_carlo_seed_;//"Seed" hbond to branch off of
	core::Real full_twobody_energy_;//Sum of hbond score + clash score for all residue pairs in "residues" data object
	core::Real score_;//This holds whatever metric is used for sorting

	boost::container::flat_map <
		unsigned int /*mres*/,
		core::scoring::hbonds::graph::AtomInfoSet >
		unsatisfied_sc_atoms_;
};

inline
core::scoring::hbonds::graph::AtomInfoSet *
NetworkState::get_unsats_for_mres( unsigned int mres ) {
	auto iter = unsatisfied_sc_atoms_.find( mres );
	if ( iter == unsatisfied_sc_atoms_.end() ) {
		return 0;
	} else {
		return &( iter->second );
	}
}

inline
core::scoring::hbonds::graph::AtomInfoSet const *
NetworkState::get_unsats_for_mres( unsigned int mres ) const {
	auto iter = unsatisfied_sc_atoms_.find( mres );
	if ( iter == unsatisfied_sc_atoms_.end() ) {
		return 0;
	} else {
		return &( iter->second );
	}
}


inline
bool
NetworkState::mres_has_unsats( unsigned int mres ) const {
	return get_unsats_for_mres( mres ) != nullptr;
}

inline
bool
NetworkState::mres_has_heavy_unsats ( unsigned int mres ) const {
	core::scoring::hbonds::graph::AtomInfoSet const * unsats =
		get_unsats_for_mres( mres );

	if ( unsats == nullptr ) return false;

	if ( unsats->empty() ) return false;

	//These vectors are arranged in a way where all heavy atoms are listed before the hydrogens.
	//If the first element is a hydrogen, there must be no heavy atoms

#ifndef NDEBUG
	//Assert sorted
	bool we_have_hit_a_H = false;
	for ( core::scoring::hbonds::graph::AtomInfo const & ai : * unsats ) {
		if ( ai.is_hydrogen() ) {
			we_have_hit_a_H = true;
		} else {
			debug_assert( ! we_have_hit_a_H );
		}
	}
#endif

	return ! unsats->begin()->is_hydrogen();
}


struct NetworkStateScoreComparator{
	static bool compare( NetworkState const & a, NetworkState const & b ){
		return a.score() < b.score();
	}
};

}//hbnet
}//protocols
#endif
