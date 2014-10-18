// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/import_pose/import_pose_options.cc
///
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

//Core headers
#include <core/types.hh>

// Unit headers
#include <core/import_pose/import_pose_options.hh>
#include <core/import_pose/import_pose_options_creator.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// C++ headers

namespace core {
namespace import_pose {

static thread_local basic::Tracer tr( "core.import_pose.import_pose_options" );


///// ImportPoseOptionsCreator /////
ImportPoseOptionsCreator::ImportPoseOptionsCreator() {}

ImportPoseOptionsCreator::~ImportPoseOptionsCreator() {}

basic::resource_manager::ResourceOptionsOP
ImportPoseOptionsCreator::create_options() const {
	return basic::resource_manager::ResourceOptionsOP( new ImportPoseOptions );
}

///@detail NOTE: This creator creates an options type called
///'PoseFromPDBOptions' to make it consistent with
///'PoseFromPDBLoader'.
std::string
ImportPoseOptionsCreator::options_type() const {
	return "PoseFromPDBOptions";
}


ImportPoseOptions::ImportPoseOptions() { init_from_options(); }

ImportPoseOptions::~ImportPoseOptions() {}

std::string ImportPoseOptions::type() const { return "ImportPoseOptions"; }

void ImportPoseOptions::parse_my_tag( utility::tag::TagCOP tag )
{
	PDB_DReaderOptions::parse_my_tag( tag );

	set_centroid( tag->getOption< bool >( "centroid", 0 ));
	set_fold_tree_io( tag->getOption< bool >( "fold_tree_io", 0 ));
	set_membrane( tag->getOption< bool >( "membrane", 0 ));
	set_no_optH( tag->getOption< bool >( "no_optH", 0 ));
	set_pack_missing_sidechains( tag->getOption< bool >( "pack_missing_sidechains", 1 ));
	set_read_fold_tree( tag->getOption< bool >( "read_fold_tree", 0 ));
	set_rna( tag->getOption< bool >( "rna", 0 ));
	set_skip_set_reasonable_fold_tree( tag->getOption< bool >( "skip_set_reasonable_fold_tree", 0 ));
	set_residue_type_set( tag->getOption< std::string >( "residue_type_set", "fa_standard" ));
	set_set_up_metal_bonds( tag->getOption<bool>("auto_setup_metals", false) );
	set_metal_bond_LJ_multiplier( tag->getOption<core::Real>("metals_detection_LJ_multiplier", 1.0) );
	set_metal_bond_dist_constraint_multiplier( tag->getOption<core::Real>("metals_distance_constraint_multiplier", 1.0) );
	set_metal_bond_angle_constraint_multiplier( tag->getOption<core::Real>("metals_angle_constraint_multiplier", 1.0) );
}

// accessors
bool ImportPoseOptions::centroid() const { return centroid_; }
bool ImportPoseOptions::fold_tree_io() const { return fold_tree_io_; }
bool ImportPoseOptions::membrane() const { return membrane_; }
bool ImportPoseOptions::no_optH() const { return no_optH_; }
bool ImportPoseOptions::pack_missing_sidechains() const { return pack_missing_sidechains_; }
bool ImportPoseOptions::read_fold_tree() const { return read_fold_tree_; }
bool ImportPoseOptions::rna() const { return rna_; }
bool ImportPoseOptions::skip_set_reasonable_fold_tree() const { return skip_set_reasonable_fold_tree_; }
bool ImportPoseOptions::set_up_metal_bonds() const { return set_up_metal_bonds_;}
bool ImportPoseOptions::set_up_metal_constraints() const { return (set_up_metal_bonds_ && (metal_bond_dist_constraint_multiplier_ > 1.0e-10  || metal_bond_angle_constraint_multiplier_ > 1.0e-10 ) );}
core::Real ImportPoseOptions::metal_bond_LJ_multiplier() const { return metal_bond_LJ_multiplier_; }
core::Real ImportPoseOptions::metal_bond_dist_constraint_multiplier() const { return metal_bond_dist_constraint_multiplier_; }
core::Real ImportPoseOptions::metal_bond_angle_constraint_multiplier() const { return metal_bond_angle_constraint_multiplier_;  }

std::string const & ImportPoseOptions::residue_type_set() const { return residue_type_set_; }

// mutators
void ImportPoseOptions::set_centroid( bool centroid ) { centroid_ = centroid; }
void ImportPoseOptions::set_fold_tree_io( bool fold_tree_io ) { fold_tree_io_ = fold_tree_io; }
void ImportPoseOptions::set_membrane( bool membrane ) { membrane_ = membrane; }
void ImportPoseOptions::set_no_optH( bool no_optH ) { no_optH_ = no_optH; }
void ImportPoseOptions::set_pack_missing_sidechains( bool pack_missing_sidechains ) { pack_missing_sidechains_ = pack_missing_sidechains; }
void ImportPoseOptions::set_read_fold_tree( bool read_fold_tree ) { read_fold_tree_ = read_fold_tree; }
void ImportPoseOptions::set_rna( bool rna ) { rna_ = rna; }
void ImportPoseOptions::set_skip_set_reasonable_fold_tree( bool skip_set_reasonable_fold_tree ) { skip_set_reasonable_fold_tree_ = skip_set_reasonable_fold_tree; }
void ImportPoseOptions::set_set_up_metal_bonds( bool invalue ) {set_up_metal_bonds_ = invalue; return;}
void ImportPoseOptions::set_metal_bond_LJ_multiplier(core::Real invalue ) { metal_bond_LJ_multiplier_ = invalue; return;}
void ImportPoseOptions::set_metal_bond_dist_constraint_multiplier(core::Real invalue ) { metal_bond_dist_constraint_multiplier_ = invalue; return;}
void ImportPoseOptions::set_metal_bond_angle_constraint_multiplier(core::Real invalue ) { metal_bond_angle_constraint_multiplier_ = invalue; return;}

void ImportPoseOptions::set_residue_type_set( std::string const & residue_type_set ) { residue_type_set_ = residue_type_set; }


void ImportPoseOptions::init_from_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	set_centroid( option[ in::file::centroid_input ]()
		|| option[ in::file::centroid ]()
		|| ( option[ in::file::fullatom ].user() && !option[ in::file::fullatom ]())
		|| ( option[ in::file::residue_type_set ].user() && option[ in::file::residue_type_set ]() == "centroid" ));

	// sanity check
	if ( centroid() &&
		( option[ in::file::fullatom ]()
			|| ( option[ in::file::residue_type_set ].user() && option[ in::file::residue_type_set ]() == "fa_standard" ))) {
		tr.Warning << "conflicting command line flags for centroid/full-atom input. Choosing fullatom!" << std::endl;
		set_centroid( false );
	}
	set_fold_tree_io( option[ inout::fold_tree_io ].user());
	set_membrane( option[ in::membrane ].user() );
	set_no_optH( option[ packing::no_optH ]());
	set_pack_missing_sidechains( option[ packing::pack_missing_sidechains ].value());
	set_read_fold_tree( false ); // no option for this parameter - it can only be set to true if you call pose_from_pdd.
	set_rna(option[ in::file::residue_type_set ].user() && option[ in::file::residue_type_set]()  == "rna");
	set_skip_set_reasonable_fold_tree( option[ run::skip_set_reasonable_fold_tree ].value());

	set_set_up_metal_bonds( option[in::auto_setup_metals].user() );
	set_metal_bond_LJ_multiplier( option[in::metals_detection_LJ_multiplier]() );
	set_metal_bond_dist_constraint_multiplier( option[in::metals_distance_constraint_multiplier]() );
	set_metal_bond_angle_constraint_multiplier( option[in::metals_angle_constraint_multiplier]() );

	set_residue_type_set( option[ in::file::residue_type_set ]());
}

} // namespace import_pose
} // namespace core
