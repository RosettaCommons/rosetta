// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file .hh file for enzdes sequence recovery cache
/// @brief
/// @author sinibjelic@gmail.com

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_EnzdesSeqRecoveryCache_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_EnzdesSeqRecoveryCache_hh

//unit headers
#include <protocols/toolbox/match_enzdes_util/EnzdesSeqRecoveryCache.fwd.hh>

//package headers

//project headers
#include <core/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>
#include <set>

#include <core/id/SequenceMapping.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

class EnzdesSeqRecoveryCache : public utility::pointer::ReferenceCount {

public:

	EnzdesSeqRecoveryCache();

	virtual ~EnzdesSeqRecoveryCache();

	//copy constructor
	EnzdesSeqRecoveryCache( EnzdesSeqRecoveryCache const & other );

	void set_sequence( core::pose::Pose & native_pose);

	std::map< core::Size, char > get_sequence();

	void set_designable_residues( std::set< core::Size > des_res );

	std::set< core::Size > get_designable_residues();

	core::Real
	sequence_recovery(
		core::pose::Pose const & designed_pose
	) const;

	void remap_residues( core::id::SequenceMapping const & smap );

private:
	std::map< core::Size, char > sequence_;
	std::set< core::Size > designable_residues_;

};

} //match_enzdes_util
} //toolbox
} //protocols

#endif
