// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/checker/VDW_CachedRepScreenInfo.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_rna_checker_VDW_CachedRepScreenInfo_HH
#define INCLUDED_protocols_stepwise_modeler_rna_checker_VDW_CachedRepScreenInfo_HH

#include <utility/pointer/ReferenceCount.hh>
#include <basic/datacache/CacheableData.hh>
#include <protocols/stepwise/modeler/rna/checker/VDW_CachedRepScreenInfo.fwd.hh>
#include <core/pose/rna/VDW_RepScreenInfo.hh>
#include <core/pose/rna/VDW_Grid.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <utility/vector1.hh>
#include <basic/datacache/BasicDataCache.hh>



#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace checker {


class VDW_CachedRepScreenInfo : public basic::datacache::CacheableData {

public:

	// @brief default constructor
	VDW_CachedRepScreenInfo();

	// @brief constructor
	VDW_CachedRepScreenInfo( core::pose::Pose const & pose );

	// @brief copy constructor
	VDW_CachedRepScreenInfo( VDW_CachedRepScreenInfo const & src );

	// @brief default destructor
	virtual ~VDW_CachedRepScreenInfo();

	basic::datacache::CacheableDataOP
	clone() const;

        void
        read_in_VDW_rep_screen_pose( core::pose::rna::VDW_RepScreenInfo & VDW_rep_screen_info ) const; 
    
        void
        read_in_VDW_rep_screen_pose_from_command_line() const;

	utility::vector1< core::pose::rna::VDW_RepScreenInfo > &
	VDW_rep_screen_info_list() const;

	core::pose::rna::VDW_GridCOP
	VDW_screen_bin() const;

private:

	mutable utility::vector1< core::pose::rna::VDW_RepScreenInfo > VDW_rep_screen_info_list_;
	core::pose::rna::VDW_GridCOP VDW_screen_bin_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


VDW_CachedRepScreenInfo const &
const_vdw_cached_rep_screen_info_from_pose( core::pose::Pose const & pose );

VDW_CachedRepScreenInfo &
nonconst_vdw_cached_rep_screen_info_from_pose( core::pose::Pose & pose );

VDW_CachedRepScreenInfo const &
make_sure_vdw_cached_rep_screen_info_is_setup( core::pose::Pose & pose );

bool
vdw_cached_rep_screen_info_is_setup( core::pose::Pose const & pose );

bool
option_vdw_rep_screen_info_user();

void
set_vdw_cached_rep_screen_info( core::pose::Pose & pose, VDW_CachedRepScreenInfoOP & vdw_cached_rep_screen_info );

void
set_vdw_cached_rep_screen_info_from_pose( core::pose::Pose & new_pose, core::pose::Pose const & pose );

void
fill_vdw_cached_rep_screen_info_from_command_line( core::pose::Pose & pose );

void
fill_vdw_cached_rep_screen_info_from_command_line( utility::vector1< core::pose::Pose * > & pose_pointers );


} //checker
} //rna
} //modeler
} //stepwise
} //protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_stepwise_modeler_rna_checker_VDW_CachedRepScreenInfo )
#endif // SERIALIZATION


#endif
