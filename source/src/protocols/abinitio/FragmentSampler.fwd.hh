// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/ShortestPathInFoldTree.fwd.hh
/// @brief  kinematics::ShortestPathInFoldTree forward declarations header
/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_FragmentSampler_fwd_hh
#define INCLUDED_protocols_abinitio_FragmentSampler_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace abinitio {

// Forward
class FragmentSampler;

// Types
typedef  utility::pointer::shared_ptr< FragmentSampler >  FragmentSamplerOP;
typedef  utility::pointer::shared_ptr< FragmentSampler const >  FragmentSamplerCOP;

enum StageID {
	ALL_STAGES = 0, //don't change!
	STAGE_1,
	STAGE_2,
	STAGE_3,
	STAGE_3a,
	STAGE_3b,
	STAGE_4,
	END_ABINITIO,
	LOOP_CLOSURE,
	SWITCH_TO_FULLATOM,
	RELAX,
	LAST_STAGE //keep last
};


} // namespace kinematics
} // namespace core

#endif
