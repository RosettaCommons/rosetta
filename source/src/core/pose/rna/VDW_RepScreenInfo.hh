// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/rna/VDW_RepScreenInfo.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_rna_VDW_RepScreenInfo_HH
#define INCLUDED_core_pose_rna_VDW_RepScreenInfo_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/pose/rna/VDW_RepScreenInfo.fwd.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <string>
#include <core/id/AtomID_Map.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace rna {

class VDW_RepScreenInfo: public utility::pointer::ReferenceCount{

public:

	VDW_RepScreenInfo();

	VDW_RepScreenInfo( VDW_RepScreenInfo const & src );

	virtual ~VDW_RepScreenInfo();

public:

	utility::vector1< core::Size > VDW_align_res;
	utility::vector1< core::Size > working_align_res;
	utility::vector1< core::Size > full_align_res;
	utility::vector1< core::Size > VDW_ignore_res;
	core::pose::PoseOP VDW_pose;
	std::string input_string;
	std::string pose_name;
	bool in_root_partition;
	core::Size import_ID;
	core::id::AtomID_Map< core::id::AtomID > align_working_to_vdw_atom_id_map;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};



} //rna
} //pose
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pose_rna_VDW_RepScreenInfo )
#endif // SERIALIZATION


#endif
