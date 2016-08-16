// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/RotamerSet/RotamerSets.hh
/// @brief  RotamerSets base class declaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_rotamer_set_RotamerSetsBase_hh
#define INCLUDED_core_pack_rotamer_set_RotamerSetsBase_hh

// Unit Headers
#include <core/pack/rotamer_set/RotamerSetsBase.fwd.hh>

// Package Headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>


// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace rotamer_set {

class RotamerSetsBase : public utility::pointer::ReferenceCount
{

public:
	RotamerSetsBase();
	virtual ~RotamerSetsBase();

	virtual uint nrotamers() const = 0;
	virtual uint nrotamers_for_moltenres( uint ) const = 0;

	virtual uint total_residue() const = 0;

	virtual uint nmoltenres() const = 0;

	virtual
	uint
	moltenres_2_resid( uint ) const = 0;

	virtual
	uint
	resid_2_moltenres( uint ) const = 0;

	virtual
	uint
	moltenres_for_rotamer( uint ) const = 0;

	virtual
	uint
	res_for_rotamer( uint ) const = 0;

	virtual
	core::conformation::ResidueCOP
	rotamer( uint ) const = 0;

	virtual
	core::conformation::ResidueCOP
	rotamer_for_moltenres( uint moltenres_id, uint rotamerid ) const = 0;

	virtual
	uint
	nrotamer_offset_for_moltenres( uint ) const = 0;

	/// @brief convert rotid in full rotamer enumeration into rotamer id on its source residue
	virtual
	uint
	rotid_on_moltenresidue( uint rotid ) const = 0;

	/// @brief convert moltenres rotid to id in full rotamer enumeration
	virtual
	uint
	moltenres_rotid_2_rotid( uint moltenres, uint moltenresrotid ) const = 0;

	/// @brief Give the pose a chance to stash any data needed by the _rotset_
	///        need nonconst access to pose
	virtual
	void
	initialize_pose_for_rotsets_creation(
		pose::Pose & pose
	) const = 0;

	virtual
	void
	show( std::ostream & out ) const = 0;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

inline
std::ostream &
operator<<( std::ostream & out, RotamerSetsBase const & rs) {
	rs.show( out );
	return out;
}

} // namespace rotamer_set
} // namespace pack
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_rotamer_set_RotamerSetsBase )
#endif // SERIALIZATION


#endif // INCLUDED_core_pack_RotamerSet_RotamerSetsBase_HH
