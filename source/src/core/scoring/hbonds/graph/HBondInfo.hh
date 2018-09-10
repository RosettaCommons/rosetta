// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/hbonds/graph/AtomLevelHBondGraph.hh
/// @brief class headers for AtomLevelHBondGraph, AtomLevelHBondNode, and AtomLevelHBondEdge
/// @details This class is used to store and traverse data used for HBNet's monte carlo branching protocol. Most (if not all) of the information held in this graph is from HBNet's RotamerSets and InteractionGraph. Nodes in this graph represent rotamers from the rotamer set (node id == global rotamer id) and edges respresent hydrogen bonds.
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_core_scoring_hbonds_graph_HBondInfo_HH
#define INCLUDED_core_scoring_hbonds_graph_HBondInfo_HH

#include <numeric/xyzVector.hh>
#include <core/types.hh>
#include <algorithm>

namespace core {
namespace scoring {
namespace hbonds {
namespace graph {


class LKHBondInfo {
public:
	LKHBondInfo( LKHBondInfo const & ){}

	virtual ~LKHBondInfo(){}

public:

	bool operator==( LKHBondInfo const & ) const { return true; }

private:
};


class HBondInfo {
public:
	HBondInfo() :
		first_node_is_donor_( 0 ),
		local_atom_id_A_( 0 ),
		local_atom_id_D_( 0 ),
		local_atom_id_H_( 0 ),
		lk_info_( 0 )
	{}

	HBondInfo(
		LKHBondInfo * lk_info,
		bool first_node_is_donor,
		unsigned short int local_atom_id_A,
		unsigned short int local_atom_id_D,
		unsigned short int local_atom_id_H
	) :
		first_node_is_donor_( first_node_is_donor ),
		local_atom_id_A_( local_atom_id_A ),
		local_atom_id_D_( local_atom_id_D ),
		local_atom_id_H_( local_atom_id_H ),
		lk_info_( lk_info )
	{}

	HBondInfo(
		bool first_node_is_donor,
		unsigned short int local_atom_id_A,
		unsigned short int local_atom_id_D,
		unsigned short int local_atom_id_H
	) :
		first_node_is_donor_( first_node_is_donor ),
		local_atom_id_A_( local_atom_id_A ),
		local_atom_id_D_( local_atom_id_D ),
		local_atom_id_H_( local_atom_id_H ),
		lk_info_( 0 )
	{}


	virtual ~HBondInfo(){}

public:
	inline bool first_node_is_donor() const {
		return first_node_is_donor_;
	}

	inline void first_node_is_donor( bool setting ){
		first_node_is_donor_ = setting;
	}

	inline unsigned short int local_atom_id_A() const {
		return local_atom_id_A_;
	}

	inline void local_atom_id_A( unsigned short int setting ) {
		local_atom_id_A_ = setting;
	}

	inline unsigned short int local_atom_id_D() const {
		return local_atom_id_D_;
	}

	inline void local_atom_id_D( unsigned short int setting ) {
		local_atom_id_D_ = setting;
	}

	inline unsigned short int local_atom_id_H() const {
		return local_atom_id_H_;
	}

	inline void local_atom_id_H( unsigned short int setting ) {
		local_atom_id_H_ = setting;
	}

	inline LKHBondInfo * lk_info() {
		return lk_info_;
	}

	inline LKHBondInfo const * lk_info() const {
		return lk_info_;
	}

	inline void lk_info( LKHBondInfo * lk_info_setting ){
		lk_info_ = lk_info_setting;
	}

	bool operator==( HBondInfo const & ot ) const {
		if ( lk_info_ && ot.lk_info_ ) {
			if ( !( *lk_info_ == *(ot.lk_info_) ) ) return false;
		}
		return  ( first_node_is_donor_ == ot.first_node_is_donor_ ) &&
			( local_atom_id_A_ == ot.local_atom_id_A_ ) &&
			( local_atom_id_D_ == ot.local_atom_id_D_ ) &&
			( local_atom_id_H_ == ot.local_atom_id_H_ );
	}

private:
	bool first_node_is_donor_;
	unsigned short int local_atom_id_A_;//Acceptor
	unsigned short int local_atom_id_D_;//Donor
	unsigned short int local_atom_id_H_;//Hydrogen
	LKHBondInfo * lk_info_;//Does not own!
};



} //graph
} //hbonds
} //scoring
} //core

#endif
