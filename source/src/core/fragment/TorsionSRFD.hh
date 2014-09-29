// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragData.hh
/// @brief  A fragment as list of SingleResidue Data
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///
#ifndef INCLUDED_core_fragment_TorsionSRFD_HH
#define INCLUDED_core_fragment_TorsionSRFD_HH

// Unit Headers
#include <core/fragment/SingleResidueFragData.hh>

// Package Headers
#include <core/fragment/Frame.fwd.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/types.hh>

#include <core/conformation/Residue.hh> // for ResidueSRFD

#include <core/kinematics/types.hh>
#include <core/id/TorsionID.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace fragment {

// WARNING: this class has not been maintained ... it is a sketch of how one could do non-backbone Torsions in fragments.. but
// it has to be brought up-to date ...

// that is like TorsionID without the seqpos
typedef std::pair < id::TorsionType, Size > TorsionInfo;

typedef utility::vector1 < TorsionInfo > TorsionInfoList;

class TorsionInfoSet;
	//typedef utility::pointer::owning_ptr< TorsionInfoSet const > TorsionInfoSetCOP;

class TorsionInfoSet : public utility::pointer::ReferenceCount {
public:
	id::TorsionID get_full_id( Size seq_pos, Size dof_id ) const {
		TorsionInfo const& t_info( list_[ dof_id ] );
		return id::TorsionID( seq_pos, t_info.first, t_info.second);
	};

	Size size() const
	{ return list_.size(); };

private:
	TorsionInfoList list_;
};

/*
class TorsionSRFD : public BaseTorsionSRFD {
public:
	TorsionSRFD() {
		//setup data and verify that
		assert( data_.size() == torsion_ids_->size() );
	};

	bool apply ( pose::Pose& pose, Size seq_pos ) const {
		for ( Size dof_id = 1; dof_id < data_.size(); dof_id++ ) {
			if ( ! pose.residue( seq_pos ).is_protein() )	return false;
			pose.set_torsion( torsion_ids_->get_full_id( seq_pos, dof_id ), data_[ dof_id ] );
		};
		return true;
	};

private:
	TorsionInfoSetCOP torsion_ids_; //COP--> allows to share info between TorsionSRFD but one might be able to change the INFO from outside !
	utility::vector1 < Real > data_;
};
 */
} //fragment
} //core

#endif
