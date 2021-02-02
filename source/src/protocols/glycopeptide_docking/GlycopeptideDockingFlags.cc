// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/glycopeptide_docking/GlycopeptideDockingFlags.cc
/// @brief Holds flags, options, and pose variables for the Glycosylation Protocol.
/// @author Yashes Srinivasan (yashess@gmail.com) , Sai Pooja Mahajan (saipooja@gmail.com)

#include <protocols/glycopeptide_docking/GlycopeptideDockingFlags.hh>
#include <basic/Tracer.hh>

// Core headers

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/chains_util.hh>

// Basic/Utility headers

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/rings.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/carbohydrates.OptionKeys.gen.hh>
#include <utility/vector1.hh>

// C++ headers

#include <string>

static basic::Tracer TR( "protocols.glycopeptide_docking.GlycopeptideDockingFlags" );

using namespace core;
using namespace basic::options;
using namespace std;
using namespace utility;

namespace protocols {
namespace glycopeptide_docking {

/// @brief Constructor sets user options
GlycopeptideDockingFlags::GlycopeptideDockingFlags()
{
	setup_from_options();

}



GlycopeptideDockingFlagsOP
GlycopeptideDockingFlags::clone() const {
	return GlycopeptideDockingFlagsOP( utility::pointer::make_shared< GlycopeptideDockingFlags>( *this ) );
}


void
GlycopeptideDockingFlags::setup_from_options(){
	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::high_res ].user() ) {
		glycosylation_high_res_refinement_ = option[ OptionKeys::carbohydrates::glycopeptide_docking::high_res ];
	}

	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::low_res ].user() ) {
		glycosylation_low_res_refinement_ = option[ OptionKeys::carbohydrates::glycopeptide_docking::low_res ];
	}

	/*if ( option[ OptionKeys::rings::idealize_rings ].active() ) {
	idealize_rings_ = option[ OptionKeys::rings::idealize_rings ];
	}*/

	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::output_debug_pdbs ].active() ) {
		debug_pdbs_ = option[ OptionKeys::carbohydrates::glycopeptide_docking::output_debug_pdbs ];
	}

	/*
	if ( option[ OptionKeys::rings::lock_rings ].active() ) {
	lock_rings_ = option[ OptionKeys::rings::lock_rings ];
	}*/

	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::substrate_type ].user() ) {
		substrate_type_ = option[ OptionKeys::carbohydrates::glycopeptide_docking::substrate_type ];
	}
	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::tree_type ].user() ) {
		tree_type_ = option[ OptionKeys::carbohydrates::glycopeptide_docking::tree_type ];
	}

	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::randomize_substrate_torsions ].user() ) {
		randomize_substrate_torsions_ = option[ OptionKeys::carbohydrates::glycopeptide_docking::randomize_substrate_torsions ];
	}

	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::enable_backbone_moves_pp ].user() ) {
		enable_backbone_moves_pp_ = option[ OptionKeys::carbohydrates::glycopeptide_docking::enable_backbone_moves_pp ];
	}

	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::prevent_anchor_repacking ].user() ) {
		prevent_anchor_repacking_ = option[ OptionKeys::carbohydrates::glycopeptide_docking::prevent_anchor_repacking ];
	}

	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::sugardonor_residue ].user() ) {
		donor_ = stoi(option[ OptionKeys::carbohydrates::glycopeptide_docking::sugardonor_residue ]);
	}

	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::nevery_interface ].user() ) {
		nevery_interface_moves_ = stoi(option[ OptionKeys::carbohydrates::glycopeptide_docking::nevery_interface ]);
	}
	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::ntotal_backbone ].user() ) {
		ntotal_backbone_moves_ = stoi(option[ OptionKeys::carbohydrates::glycopeptide_docking::ntotal_backbone ]);
	}

	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::interface_distance ].user() ) {
		interface_distance_ = option[ OptionKeys::carbohydrates::glycopeptide_docking::interface_distance ];
	}

	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::allow_glycan_torsion_moves ].user() ) {
		allow_glycan_torsion_moves_ = option[ OptionKeys::carbohydrates::glycopeptide_docking::allow_glycan_torsion_moves ];
	}

	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::score_only ].user() ) {
		score_only_ = option[ OptionKeys::carbohydrates::glycopeptide_docking::score_only ];
	}

	if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
		constraints_ = option[ OptionKeys::constraints::cst_fa_file ][ 1 ];
	}

	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::output_distance_metrics ].user() ) {
		additional_metrics_ = option[ OptionKeys::carbohydrates::glycopeptide_docking::output_distance_metrics ];
	}

	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::output_distance_metrics ].user() ) {
		additional_metrics_ = option[ OptionKeys::carbohydrates::glycopeptide_docking::output_distance_metrics ];
	}
}

/// @brief Specify the enzyme chain from pose information
void
GlycopeptideDockingFlags::set_enzyme_chain( core::pose::Pose const &pose )
{
	/* TODO: Add an option to provide receptor chain */

	vector1< core::uint > const & chain_endings( pose.conformation().chain_endings() );
	first_residue_enzyme_ =  1; // default
	last_residue_enzyme_ = chain_endings[ 1 ]; // default

}


///@Get residue from pose
core::Size
GlycopeptideDockingFlags::get_resnum(core::pose::Pose const &pose,std::string const &special_residue)
{
	pose::PDBInfoCOP pdb_info( pose.pdb_info() );

	char chain = special_residue[ special_residue.size() - 1 ];
	string res_num_str = "";
	for ( core::Size pos = 0; pos < special_residue.size() - 1; pos++ ) {
		res_num_str += special_residue[ pos ];
	}
	core::Size res_num = stoi( res_num_str );
	return pdb_info->pdb2pose( chain, res_num );
}

///@Set anchor or residue_to_glycosylate residue/s from command line or script options.
void
GlycopeptideDockingFlags::set_anchor_residue(core::pose::Pose const &pose,std::string const &special_residue)
{
	anchor_residue_substrate_ = get_resnum(pose,special_residue);
}

///@Set anchor or residue_to_glycosylate residue/s from command line or script options.
void
GlycopeptideDockingFlags::set_glycosylation_residue(core::pose::Pose const &pose,std::string const &special_residue)
{
	residue_to_glycosylate_ = get_resnum(pose,special_residue);
}


///@ brief Specify the substrate chain from pose information
void
GlycopeptideDockingFlags::set_substrate_chain( core::pose::Pose const &pose)
{
	/* TODO: Add an option to provide substrate chain */

	vector1< core::uint > const & chain_endings( pose.conformation().chain_endings() );
	Size const n_chain_endings( chain_endings.size() );
	core::uint chain_num_of_last_cut_point( n_chain_endings );
	core::uint last_cut_point( chain_endings[ chain_num_of_last_cut_point ] );
	// Clear data.
	upstream_chains_ = "";
	downstream_chains_ = "";

	while ( pose.residue( last_cut_point + 1 ).has_lower_connect() && pose.residue( last_cut_point + 1 ).connected_residue_at_lower() != last_cut_point ) {
		--chain_num_of_last_cut_point;
		last_cut_point = chain_endings[ chain_num_of_last_cut_point ];
	}
	first_residue_substrate_ = last_cut_point + 1;  // default
	last_residue_substrate_ = core::pose::chain_end_res( pose,  pose.residue(first_residue_substrate_ ).chain());

	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::residue_to_glycosylate ].active() ) {
		std::string residue_to_glycosylate = option[ OptionKeys::carbohydrates::glycopeptide_docking::residue_to_glycosylate ];
		set_glycosylation_residue(pose,residue_to_glycosylate);
		//initialize anchor residue as residue to glycosylate
		set_anchor_residue(pose,residue_to_glycosylate);
	}
	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::anchor_residue ].active() ) {
		std::string anchor_residue = option[ OptionKeys::carbohydrates::glycopeptide_docking::anchor_residue ];
		//when achor residue and residue to glycosylate are different
		//this is useful for anchored docking of a glycopeptide for glycosylation at a new site
		//E.g. GalNAcT4 and GalNAcT12 - glycosylate diglycopeptides
		set_anchor_residue(pose,anchor_residue);
	}
	/*if ( option[ OptionKeys::carbohydrates::glycosylation::residue_to_randomize ].user() ) {
	std::string residue_to_randomize = option[ OptionKeys::carbohydrates::glycosylation::residue_to_randomize];
	set_residue_to_randomize(pose,residue_to_randomize);
	enable_randomize_residue(true);
	}*/
	pose::PDBInfoCOP pdb_info( pose.pdb_info() );
	jump_num_substrate_ = 1;    // default
	core::uint cut_point;
	for ( core::uint chain_number( 1 ); chain_number <= n_chain_endings; ++chain_number ) {
		cut_point = chain_endings[ chain_number ];
		if ( cut_point < last_cut_point ) {
			upstream_chains_ += pdb_info->chain( cut_point );
		} else if ( cut_point == last_cut_point ) {
			upstream_chains_ += pdb_info->chain( cut_point );
			downstream_chains_ += pdb_info->chain( cut_point + 1 );
		} else {
			downstream_chains_ += pdb_info->chain( cut_point + 1 );
		}
	}
	TR << "Upstream and downstream "<< upstream_chains_ << " " << downstream_chains_ << std::endl;
	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::preglycosylate_residues ].user() ) {
		std::string preglycosylate_residues = option[ OptionKeys::carbohydrates::glycopeptide_docking::preglycosylate_residues ];
		preglycosylate_residues_ = core::pose::get_resnum_list_ordered(preglycosylate_residues,pose);
		TR << "Positions " << preglycosylate_residues_<< " will be preglycosylated."<<endl;
	}
	if ( option[ OptionKeys::carbohydrates::glycopeptide_docking::preglycosylate_sugar_names ].user() ) {
		std::string preglycosylate_sugar_names = option[ OptionKeys::carbohydrates::glycopeptide_docking::preglycosylate_sugar_names ];
		preglycosylate_sugars_ =  utility::string_split( preglycosylate_sugar_names , ',' );
		TR << "Sugars " << preglycosylate_sugars_<< " will be preglycosylated."<<endl;
	}

}

core::Size
GlycopeptideDockingFlags::first_residue_substrate() const {
	return first_residue_substrate_;
}

core::Size
GlycopeptideDockingFlags::last_residue_substrate() const {
	return last_residue_substrate_;
}

core::Size
GlycopeptideDockingFlags::anchor_residue_substrate() const {
	return anchor_residue_substrate_;
}

core::Size
GlycopeptideDockingFlags::glycosylation_residue_substrate() const {
	return residue_to_glycosylate_;
}

core::Size
GlycopeptideDockingFlags::jump_num_substrate() const {
	return jump_num_substrate_;
}

core::Size
GlycopeptideDockingFlags::first_residue_enzyme() const {
	return first_residue_enzyme_;
}

core::Size
GlycopeptideDockingFlags::last_residue_enzyme() const {
	return last_residue_enzyme_;
}

void
GlycopeptideDockingFlags::show() const {

	TR << endl << "------------ FLAGS ------------" << endl;
	TR << "glycosylation_high_res_refinement:  " << glycosylation_high_res_refinement_ << endl;
	TR << "glycosylation_low_res_refinement:   " << glycosylation_low_res_refinement_ << endl;
	TR << "constraints:                        " << constraints_ << endl;
	TR << "first_residue_substrate:            " << first_residue_substrate() << endl;
	TR << "last_residue_substrate:             " << last_residue_substrate() << endl;
	TR << "anchor_residue_substrate:           " << anchor_residue_substrate() << endl;
	TR << "residue_to_glycosylate:             " << glycosylation_residue_substrate() << endl;
	TR << "jump_num_substrate:                 " << jump_num_substrate() << endl;
	TR << "first_residue_enzyme:               " << first_residue_enzyme() << endl;
	TR << "last_residue_enzyme:                " << last_residue_enzyme() << endl;
}

} //glycosylation
} //protocols






