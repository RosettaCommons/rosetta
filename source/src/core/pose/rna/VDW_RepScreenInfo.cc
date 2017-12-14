// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/rna/VDW_RepScreenInfo.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pose/rna/VDW_RepScreenInfo.hh>
#include <core/pose/Pose.hh>
#include <string>


#include <basic/Tracer.hh>

static basic::Tracer TR( "core.pose.rna.VDW_RepScreenInfo" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <core/id/AtomID_Map.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace rna {

//Constructor
VDW_RepScreenInfo::VDW_RepScreenInfo():
	input_string( "" ),
	pose_name( "" ),
	in_root_partition( false ),
	import_ID( 0 ),
	align_working_to_vdw_atom_id_map( )
{
	VDW_align_res.clear();
	working_align_res.clear();
	full_align_res.clear();
	VDW_ignore_res.clear();
}

// Copy Constructor
VDW_RepScreenInfo::VDW_RepScreenInfo( VDW_RepScreenInfo const & src ):
	VDW_align_res( src.VDW_align_res ),
	working_align_res( src.working_align_res ),
	full_align_res( src.full_align_res ),
	VDW_ignore_res( src.VDW_ignore_res ),
	VDW_pose( src.VDW_pose ),
	input_string( src.input_string ),
	pose_name( src.pose_name ),
	in_root_partition( src.in_root_partition ),
	import_ID( src.import_ID ),
	align_working_to_vdw_atom_id_map( src.align_working_to_vdw_atom_id_map )
{
}


//Destructor
VDW_RepScreenInfo::~VDW_RepScreenInfo() = default;

} //rna
} //pose
} //core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::rna::VDW_RepScreenInfo::save( Archive & arc ) const {
	arc( CEREAL_NVP( VDW_align_res ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( working_align_res ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( full_align_res ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( VDW_ignore_res ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( VDW_pose ) ); // core::pose::PoseOP
	arc( CEREAL_NVP( input_string ) ); // std::string
	arc( CEREAL_NVP( pose_name ) ); // std::string
	arc( CEREAL_NVP( in_root_partition ) ); // _Bool
	arc( CEREAL_NVP( import_ID ) ); // core::Size
	arc( CEREAL_NVP( align_working_to_vdw_atom_id_map ) ); // core::id::AtomID_Map<core::id::AtomID>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::rna::VDW_RepScreenInfo::load( Archive & arc ) {
	arc( VDW_align_res ); // utility::vector1<core::Size>
	arc( working_align_res ); // utility::vector1<core::Size>
	arc( full_align_res ); // utility::vector1<core::Size>
	arc( VDW_ignore_res ); // utility::vector1<core::Size>
	arc( VDW_pose ); // core::pose::PoseOP
	arc( input_string ); // std::string
	arc( pose_name ); // std::string
	arc( in_root_partition ); // _Bool
	arc( import_ID ); // core::Size
	arc( align_working_to_vdw_atom_id_map ); // core::id::AtomID_Map<core::id::AtomID>
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::rna::VDW_RepScreenInfo );
CEREAL_REGISTER_TYPE( core::pose::rna::VDW_RepScreenInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_rna_VDW_RepScreenInfo )
#endif // SERIALIZATION
