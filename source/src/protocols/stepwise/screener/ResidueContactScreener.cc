// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/ResidueContactScreener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/ResidueContactScreener.hh>
#include <protocols/stepwise/sampling/rna/util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.ResidueContactScreener" );

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	ResidueContactScreener::ResidueContactScreener(  pose::Pose & screening_pose,
																									 Size const last_append_res,
																									 Size const last_prepend_res,
																									 Distance const atom_atom_overlap_dist_cutoff ):
		screening_pose_( screening_pose ),
		last_append_res_( last_append_res ),
		last_prepend_res_( last_prepend_res ),
		atom_atom_overlap_dist_cutoff_( atom_atom_overlap_dist_cutoff )
	{}

	//Destructor
	ResidueContactScreener::~ResidueContactScreener()
	{}

	///////////////////////////////////////////////////////////
	bool
	ResidueContactScreener::check_screen(){
		return sampling::rna::is_residues_in_contact( last_append_res_, screening_pose_,
																									 last_prepend_res_, screening_pose_,
																									 atom_atom_overlap_dist_cutoff_, 1 /*num_atom_contacts_cutoff*/ );
	}

} //screener
} //stepwise
} //protocols
