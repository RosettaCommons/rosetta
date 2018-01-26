// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/hbnet/NetworkState.hh
/// @brief specialized structs for MC-HBNet to store hbond network info
/// @author Jack Maguire, jack@med.unc.edu

#ifndef INCLUDED_protocols_hbnet_NetworkState_hh
#define INCLUDED_protocols_hbnet_NetworkState_hh

#include <core/conformation/Residue.hh>
#include <core/scoring/hbonds/graph/HBondInfo.hh>
#include <core/scoring/hbonds/graph/AtomInfo.hh>
#include <core/scoring/hbonds/graph/AtomLevelHBondGraph.hh>
//#include <core/pack/rotamer_set/RotamerSets.hh>

namespace protocols {
namespace hbnet {

using mres_unsat_pair=std::pair< unsigned int /*mres*/, utility::vector1< core::scoring::hbonds::graph::AtomInfo > >;

struct compare_mres_unsat_pair : public std::binary_function< mres_unsat_pair, mres_unsat_pair, bool >{
	bool operator()( mres_unsat_pair const & a, mres_unsat_pair const & b ) const {
		return a.first < b.first;
	}

	bool operator()( mres_unsat_pair const & a, unsigned int b ) const {
		return a.first < b;
	}

	bool operator()( unsigned int a, mres_unsat_pair const & b ) const {
		return a < b.first;
	}

};


class NetworkState {

public:

	NetworkState(
		core::scoring::hbonds::graph::AtomLevelHBondEdge const * monte_carlo_seed,
		core::scoring::hbonds::graph::AtomLevelHBondGraphOP const & hbond_graph
	);

	~NetworkState(){}

	bool operator < ( NetworkState const & rhs) const {
		return full_twobody_energy_ < rhs.full_twobody_energy_;
	}

	///@brief NetworkState keeps track of polar atoms that need to be satisfied.
	/// This function allows NetworkState to copy the vector of polar atoms from the HBondNode
	void add_polar_atoms(
		core::scoring::hbonds::graph::AtomLevelHBondNode const * node
	);

	///@brief are there any buried unsats at this
	///molten residue position that are not participating in network hbonds?
	bool mres_has_unsats( unsigned int mres ) const;

	///@brief are there any heavy (not hydrogen) buried unsats at this
	///molten residue position that are not participating in network hbonds?
	bool mres_has_heavy_unsats ( unsigned int mres ) const;

public:
	utility::vector1< mres_unsat_pair > const & unsatisfied_sc_atoms_const() const {
		return unsatisfied_sc_atoms_;
	}

	///@brief get vector of AtomInfos for atoms that are unsatisfied
	utility::vector1< mres_unsat_pair >::iterator
	get_unsats_for_mres( unsigned int mres );

	///@brief get vector of AtomInfos for atoms that are unsatisfied
	utility::vector1< mres_unsat_pair >::const_iterator
	get_unsats_for_mres( unsigned int mres ) const;

	utility::vector1< core::scoring::hbonds::graph::AtomLevelHBondNode const * > & nodes() {
		return nodes_;
	}

	utility::vector1< core::scoring::hbonds::graph::AtomLevelHBondNode const * > const & nodes() const {
		return nodes_;
	}

	void add_node( core::scoring::hbonds::graph::AtomLevelHBondNode const * node ) {
		nodes_.push_back( node );
	}

	utility::vector1< core::scoring::hbonds::graph::AtomLevelHBondEdge const * > & edges() {
		return edges_;
	}

	utility::vector1< core::scoring::hbonds::graph::AtomLevelHBondEdge const * > const & edges() const {
		return edges_;
	}

	void add_edge( core::scoring::hbonds::graph::AtomLevelHBondEdge const * edge ) {
		edges_.push_back( edge );
	}

	core::scoring::hbonds::graph::AtomLevelHBondEdge const * monte_carlo_seed() const {
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
	utility::vector1< core::scoring::hbonds::graph::AtomLevelHBondNode const * > nodes_;
	utility::vector1< core::scoring::hbonds::graph::AtomLevelHBondEdge const * > edges_;

	core::scoring::hbonds::graph::AtomLevelHBondEdge const * monte_carlo_seed_;//"Seed" hbond to branch off of
	core::Real full_twobody_energy_;//Sum of hbond score + clash score for all residue pairs in "residues" data object
	core::Real score_;//This holds whatever metric is used for sorting

	utility::vector1< mres_unsat_pair > unsatisfied_sc_atoms_;
	compare_mres_unsat_pair sorter_;
};

inline utility::vector1< mres_unsat_pair >::iterator
NetworkState::get_unsats_for_mres( unsigned int mres ) {
	auto iter = std::lower_bound( unsatisfied_sc_atoms_.begin(), unsatisfied_sc_atoms_.end(), mres, sorter_ );
	debug_assert( iter != unsatisfied_sc_atoms_.end() );//if this fails, there is no element in the vector for this mres
	debug_assert( iter->first == mres );//if this fails, there is no element in the vector for this mres
	return iter;
}

inline utility::vector1< mres_unsat_pair >::const_iterator
NetworkState::get_unsats_for_mres( unsigned int mres ) const {
	auto const iter = std::lower_bound( unsatisfied_sc_atoms_.begin(), unsatisfied_sc_atoms_.end(), mres, sorter_ );
	debug_assert( iter != unsatisfied_sc_atoms_.end() );//if this fails, there is no element in the vector for this mres
	debug_assert( iter->first == mres );//if this fails, there is no element in the vector for this mres
	return iter;
}

inline bool
NetworkState::mres_has_unsats( unsigned int mres ) const {
	auto iter = std::lower_bound( unsatisfied_sc_atoms_.begin(), unsatisfied_sc_atoms_.end(), mres, sorter_ );
	if ( iter == unsatisfied_sc_atoms_.end() ) return false;//if this is false, there is no element in the vector for this mres
	if ( iter->first != mres ) return false;//if this is false, there is no element in the vector for this mres
	return ! iter->second.empty();//if the vector of unsats is empty, there are no unsats for this mres
}

inline bool
NetworkState::mres_has_heavy_unsats ( unsigned int mres ) const {
	auto iter = std::lower_bound( unsatisfied_sc_atoms_.begin(), unsatisfied_sc_atoms_.end(), mres, sorter_ );
	if ( iter == unsatisfied_sc_atoms_.end() ) return false;
	if ( iter->first != mres ) return false;
	if ( iter->second.empty() ) return false;

	//These vectors are arranged in a way where all heavy atoms are listed before the hydrogens.
	//If the first element is a hydrogen, there must be no heavy atoms
	return ! iter->second.front().is_hydrogen();
}


struct NetworkStateScoreComparator{
	static bool compare( NetworkState const & a, NetworkState const & b ){
		return a.score() < b.score();
	}
};

}//hbnet
}//protocols
#endif
