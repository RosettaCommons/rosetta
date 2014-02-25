// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/TagDefinition.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/TagDefinition.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Util.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.TagDefinition" );

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	TagDefinition::TagDefinition( pose::Pose & pose,
												StepWiseScreenerOP first_sampler,
												bool const sampler_include_torsion_value_in_tag,
												Size const moving_res,
												bool const is_prepend,
												std::string const extra_tag /* = ""*/  ):
		pose_( pose ),
		first_sampler_( first_sampler ),
		sampler_include_torsion_value_in_tag_( sampler_include_torsion_value_in_tag ),
		moving_res_( moving_res ),
		is_prepend_( is_prepend ),
		extra_tag_( extra_tag ),
		tag_( tag_from_pose( pose ) )
	{
	}

	//Destructor
	TagDefinition::~TagDefinition()
	{}

	///////////////////////////////////////////////////////////
	// this action could also be in the update_mover, which would move forward to all poses in later screeners.
	bool
	TagDefinition::check_screen(){
		using namespace sampling::rna;
		tag_ = create_tag( "U" + extra_tag_, first_sampler_->count() );
		if ( sampler_include_torsion_value_in_tag_ ) tag_ += create_rotamer_string( pose_, moving_res_, is_prepend_ );
		return true;
	}

} //screener
} //stepwise
} //protocols
