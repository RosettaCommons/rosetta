// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/hbonds/AtomLevelHBondGraph.hh
/// @brief This class holds the info needed to track potentially unsatisfied atoms. Used initially by the AtomLevelHBondGraph.
/// @author Jack Maguire, jack@med.unc.edu

#ifndef INCLUDED_core_scoring_hbonds_graph_AtomInfo_HH
#define INCLUDED_core_scoring_hbonds_graph_AtomInfo_HH

#include <numeric/xyzVector.hh>
#include <core/types.hh>
#include <algorithm>
#include <utility/DenseBoolMap.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace hbonds {
namespace graph {

class LKAtomInfo;
typedef utility::pointer::shared_ptr< LKAtomInfo > LKAtomInfoOP;
typedef utility::pointer::shared_ptr< LKAtomInfo const > LKAtomInfoCOP;

class AtomInfo {

public:
	AtomInfo(
		unsigned short int atomid,
		numeric::xyzVector< float > const & atom_position,
		bool is_hydrogen_setting,
		bool is_donor_setting,
		bool is_acceptor_setting,
		bool is_hydroxyl_setting
	) :
		local_atom_id_( atomid ),
		xyz_( atom_position ),
		lk_info_( 0 )
	{
		is_hydrogen( is_hydrogen_setting );
		is_donor(    is_donor_setting );
		is_acceptor( is_acceptor_setting );
		is_hydroxyl( is_hydroxyl_setting );
	}

	AtomInfo(
		LKAtomInfoCOP lk_info,
		unsigned short int atomid,
		numeric::xyzVector< float > const & atom_position,
		bool is_hydrogen_setting,
		bool is_donor_setting,
		bool is_acceptor_setting,
		bool is_hydroxyl_setting
	) :
		local_atom_id_( atomid ),
		xyz_( atom_position ),
		lk_info_( lk_info )
	{
		is_hydrogen( is_hydrogen_setting );
		is_donor(    is_donor_setting );
		is_acceptor( is_acceptor_setting );
		is_hydroxyl( is_hydroxyl_setting );
	}


	virtual ~AtomInfo(){}

private:
	enum Settings
	{
		IS_HYDROGEN = 0,
		IS_DONOR,   //1
		IS_ACCEPTOR,//2
		IS_HYDROXYL,//3
		count, //do not use for anything other than properties_ definition!
	};

	static_assert( Settings::count - Settings::IS_HYDROGEN == 4,
		"AtomInfo's enum is not a continuous range!" );

public://setters
	inline void local_atom_id( unsigned short int local_atom_id ){
		local_atom_id_ = local_atom_id;
	}

	inline void is_hydrogen( bool setting ){
		properties_.set< IS_HYDROGEN >( setting );
	}

	inline void is_donor( bool setting ){
		properties_.set< IS_DONOR >( setting );
	}

	inline void is_acceptor( bool setting ){
		properties_.set< IS_ACCEPTOR >( setting );
	}

	inline void is_hydroxyl( bool setting ){
		properties_.set< IS_HYDROXYL > ( setting );
	}

	inline void lk_info( LKAtomInfoCOP setting ){
		lk_info_ = setting;
	}

	inline void xyz( numeric::xyzVector< float > const & setting ){
		xyz_ = setting;
	}

public://getters
	inline unsigned short int local_atom_id() const {
		return local_atom_id_;
	}

	inline bool is_hydrogen() const {
		return properties_.get< IS_HYDROGEN >();
	}

	inline bool is_donor() const {
		return properties_.get< IS_DONOR >();
	}

	inline bool is_acceptor() const {
		return properties_.get< IS_ACCEPTOR >();
	}

	inline bool is_hydroxyl() const {
		return properties_.get< IS_HYDROXYL >();
	}

	inline LKAtomInfoCOP lk_info() const {
		return lk_info_;
	}

	inline numeric::xyzVector< float > const & xyz() const {
		return xyz_;
	}

private://DATA
	unsigned short int local_atom_id_;
	numeric::xyzVector< float > xyz_;

	utility::DenseBoolMap< Settings::count, Settings::IS_HYDROGEN > properties_;
	LKAtomInfoCOP lk_info_;//Useless at the moment. Will be implemented in the future
};

class LKAtomInfo {
public:
	LKAtomInfo(){}
	LKAtomInfo( LKAtomInfo const & ){}

	~LKAtomInfo(){}

private://DATA
	//Water info I guess?
};

} //graph
} //hbonds
} //scoring
} //core

#endif
