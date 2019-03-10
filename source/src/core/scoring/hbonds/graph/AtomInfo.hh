// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/hbonds/HBondGraph.hh
/// @brief This class holds the info needed to track potentially unsatisfied atoms. Used initially by the HBondGraph.
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_core_scoring_hbonds_graph_AtomInfo_HH
#define INCLUDED_core_scoring_hbonds_graph_AtomInfo_HH

#include <numeric/xyzVector.hh>
#include <core/types.hh>
#include <algorithm>
#include <utility/DenseBoolMap.hh>
#include <utility/pointer/owning_ptr.hh>

#include <boost/container/flat_set.hpp>

namespace core {
namespace scoring {
namespace hbonds {
namespace graph {

class AtomInfo {

public:
	AtomInfo(
		unsigned short int atomid,
		numeric::xyzVector< float > const & atom_position,
		bool is_hydrogen_setting,
		bool is_donor_setting,
		bool is_acceptor_setting,
		bool is_hydroxyl_setting,
		bool is_backbone_setting
	) :
		local_atom_id_( atomid ),
		xyz_( atom_position )
	{
		is_hydrogen( is_hydrogen_setting );
		is_donor(    is_donor_setting );
		is_acceptor( is_acceptor_setting );
		is_hydroxyl( is_hydroxyl_setting );
		is_backbone( is_backbone_setting );
	}

	virtual ~AtomInfo(){}

private:
	enum Settings
	{
		IS_HYDROGEN = 0,
		IS_DONOR,   //1
		IS_ACCEPTOR,//2
		IS_HYDROXYL,//3
		IS_BACKBONE,//4
		count, //do not use for anything other than properties_ definition!
	};

	static_assert( Settings::count - Settings::IS_HYDROGEN == 5,
		"AtomInfo's enum is not a continuous range!" );

public://setters
	void local_atom_id( unsigned short int local_atom_id ){
		local_atom_id_ = local_atom_id;
	}

	void is_hydrogen( bool setting ){
		properties_.set< IS_HYDROGEN >( setting );
	}

	void is_donor( bool setting ){
		properties_.set< IS_DONOR >( setting );
	}

	void is_acceptor( bool setting ){
		properties_.set< IS_ACCEPTOR >( setting );
	}

	void is_hydroxyl( bool setting ){
		properties_.set< IS_HYDROXYL > ( setting );
	}

	void is_backbone( bool setting ){
		properties_.set< IS_BACKBONE > ( setting );
	}

	void xyz( numeric::xyzVector< float > const & setting ){
		xyz_ = setting;
	}

public://getters
	unsigned short int local_atom_id() const {
		return local_atom_id_;
	}

	bool is_hydrogen() const {
		return properties_.get< IS_HYDROGEN >();
	}

	bool is_donor() const {
		return properties_.get< IS_DONOR >();
	}

	bool is_acceptor() const {
		return properties_.get< IS_ACCEPTOR >();
	}

	bool is_hydroxyl() const {
		return properties_.get< IS_HYDROXYL >();
	}

	bool is_backbone() const {
		return properties_.get< IS_BACKBONE >();
	}

	numeric::xyzVector< float > const & xyz() const {
		return xyz_;
	}

	bool operator<( AtomInfo const & ot ) const {
		return local_atom_id_ < ot.local_atom_id_;
	}

private://DATA
	unsigned short int local_atom_id_;
	numeric::xyzVector< float > xyz_;

	utility::DenseBoolMap< Settings::count, Settings::IS_HYDROGEN > properties_;
};

class AtomInfoSet : public boost::container::flat_set< AtomInfo > {
public:
	AtomInfoSet::const_iterator remove( unsigned short int local_atom_id ){
		for ( auto iter = begin(); iter != end(); ++iter ) {
			if ( iter->local_atom_id() == local_atom_id ) {
				return erase( iter );
			}
		}
		return end();
	}

	bool contains( unsigned short int local_atom_id ) const {
		auto predicate = [=]( AtomInfo const & ai ){
			return ai.local_atom_id() == local_atom_id;
		};
		return std::find_if( begin(), end(), predicate ) != end();
	}
};


} //graph
} //hbonds
} //scoring
} //core

#endif
