// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <utility/options/keys/OptionKeyList.hh>

// C++ headers

namespace core {
namespace import_pose {

static THREAD_LOCAL basic::Tracer tr( "core.import_pose.import_pose_options" );


///// ImportPoseOptionsCreator /////
ImportPoseOptionsCreator::ImportPoseOptionsCreator() {}

ImportPoseOptionsCreator::~ImportPoseOptionsCreator() = default;

basic::resource_manager::ResourceOptionsOP
ImportPoseOptionsCreator::create_options() const {
	return basic::resource_manager::ResourceOptionsOP( new ImportPoseOptions );
}

/// @detail NOTE: This creator creates an options type called
///'PoseFromPDBOptions' to make it consistent with
///'PoseFromPDBLoader'.
std::string
ImportPoseOptionsCreator::options_type() const {
	return "PoseFromPDBOptions";
}


ImportPoseOptions::ImportPoseOptions() { init_from_options( basic::options::option ); }

ImportPoseOptions::ImportPoseOptions( utility::options::OptionCollection const & options ) :
	StructFileReaderOptions( options )
{ init_from_options( options ); }

ImportPoseOptions::~ImportPoseOptions() = default;

std::string ImportPoseOptions::type() const { return "ImportPoseOptions"; }

void ImportPoseOptions::parse_my_tag( utility::tag::TagCOP tag )
{
	StructFileReaderOptions::parse_my_tag( tag );

	set_centroid( tag->getOption< bool >( "centroid", 0 ));
	set_fold_tree_io( tag->getOption< bool >( "fold_tree_io", 0 ));
	set_membrane( tag->getOption< bool >( "membrane", 0 ));
	set_no_optH( tag->getOption< bool >( "no_optH", 0 ));
	set_pack_missing_sidechains( tag->getOption< bool >( "pack_missing_sidechains", 1 ));
	set_read_fold_tree( tag->getOption< bool >( "read_fold_tree", 0 ));
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
void ImportPoseOptions::set_skip_set_reasonable_fold_tree( bool skip_set_reasonable_fold_tree ) { skip_set_reasonable_fold_tree_ = skip_set_reasonable_fold_tree; }
void ImportPoseOptions::set_set_up_metal_bonds( bool invalue ) {set_up_metal_bonds_ = invalue; return;}
void ImportPoseOptions::set_metal_bond_LJ_multiplier(core::Real invalue ) { metal_bond_LJ_multiplier_ = invalue; return;}
void ImportPoseOptions::set_metal_bond_dist_constraint_multiplier(core::Real invalue ) { metal_bond_dist_constraint_multiplier_ = invalue; return;}
void ImportPoseOptions::set_metal_bond_angle_constraint_multiplier(core::Real invalue ) { metal_bond_angle_constraint_multiplier_ = invalue; return;}

void ImportPoseOptions::set_residue_type_set( std::string const & residue_type_set ) { residue_type_set_ = residue_type_set; }


void ImportPoseOptions::init_from_options( utility::options::OptionCollection const & options )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	set_centroid( options[ in::file::centroid_input ]()
		|| options[ in::file::centroid ]()
		|| ( options[ in::file::fullatom ].user() && !options[ in::file::fullatom ]())
		|| ( options[ in::file::residue_type_set ].user() && options[ in::file::residue_type_set ]() == "centroid" ));

	// sanity check
	if ( centroid() &&
			( options[ in::file::fullatom ]()
			|| ( options[ in::file::residue_type_set ].user() && options[ in::file::residue_type_set ]() == "fa_standard" )) ) {
		tr.Warning << "conflicting command line flags for centroid/full-atom input. Choosing fullatom!" << std::endl;
		set_centroid( false );
	}
	set_fold_tree_io( options[ inout::fold_tree_io ].user());
	set_membrane( options[ in::membrane ].user() );
	set_no_optH( options[ packing::no_optH ]());
	set_pack_missing_sidechains( options[ packing::pack_missing_sidechains ].value());
	set_read_fold_tree( false ); // no option for this parameter - it can only be set to true if you call pose_from_pdd.
	set_skip_set_reasonable_fold_tree( options[ run::skip_set_reasonable_fold_tree ].value());

	set_set_up_metal_bonds( options[in::auto_setup_metals].user() );
	set_metal_bond_LJ_multiplier( options[in::metals_detection_LJ_multiplier]() );
	set_metal_bond_dist_constraint_multiplier( options[in::metals_distance_constraint_multiplier]() );
	set_metal_bond_angle_constraint_multiplier( options[in::metals_angle_constraint_multiplier]() );

	set_residue_type_set( options[ in::file::residue_type_set ]());
}

void ImportPoseOptions::list_options_read( utility::options::OptionKeyList & read_options )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	StructFileReaderOptions::list_options_read( read_options );
	read_options
		+ in::file::centroid_input
		+ in::file::centroid
		+ in::file::fullatom
		+ in::file::residue_type_set
		+ inout::fold_tree_io
		+ in::membrane
		+ packing::no_optH
		+ packing::pack_missing_sidechains
		+ run::skip_set_reasonable_fold_tree
		+ in::auto_setup_metals
		+ in::metals_detection_LJ_multiplier
		+ in::metals_distance_constraint_multiplier
		+ in::metals_angle_constraint_multiplier;

}

bool
ImportPoseOptions::operator == ( ImportPoseOptions const & other ) const
{
	if ( ! StructFileReaderOptions::operator == ( other ) ) return false;

	if ( centroid_                               != other.centroid_                               ) return false;
	if ( fold_tree_io_                           != other.fold_tree_io_                           ) return false;
	if ( membrane_                               != other.membrane_                               ) return false;
	if ( no_optH_                                != other.no_optH_                                ) return false;
	if ( pack_missing_sidechains_                != other.pack_missing_sidechains_                ) return false;
	if ( read_fold_tree_                         != other.read_fold_tree_                         ) return false;
	if ( skip_set_reasonable_fold_tree_          != other.skip_set_reasonable_fold_tree_          ) return false;
	if ( set_up_metal_bonds_                     != other.set_up_metal_bonds_                     ) return false;
	if ( metal_bond_LJ_multiplier_               != other.metal_bond_LJ_multiplier_               ) return false;
	if ( metal_bond_dist_constraint_multiplier_  != other.metal_bond_dist_constraint_multiplier_  ) return false;
	if ( metal_bond_angle_constraint_multiplier_ != other.metal_bond_angle_constraint_multiplier_ ) return false;
	if ( residue_type_set_                        != other.residue_type_set_                       ) return false;
	return true;
}

bool
ImportPoseOptions::operator < ( ImportPoseOptions const & other ) const
{
	if (   StructFileReaderOptions::operator <  ( other ) ) return true;
	if ( ! StructFileReaderOptions::operator == ( other ) ) return false;

	if ( centroid_                               <  other.centroid_                               ) return true;
	if ( centroid_                               != other.centroid_                               ) return false;
	if ( fold_tree_io_                           <  other.fold_tree_io_                           ) return true;
	if ( fold_tree_io_                           != other.fold_tree_io_                           ) return false;
	if ( membrane_                               <  other.membrane_                               ) return true;
	if ( membrane_                               != other.membrane_                               ) return false;
	if ( no_optH_                                <  other.no_optH_                                ) return true;
	if ( no_optH_                                != other.no_optH_                                ) return false;
	if ( pack_missing_sidechains_                <  other.pack_missing_sidechains_                ) return true;
	if ( pack_missing_sidechains_                != other.pack_missing_sidechains_                ) return false;
	if ( read_fold_tree_                         <  other.read_fold_tree_                         ) return true;
	if ( read_fold_tree_                         != other.read_fold_tree_                         ) return false;
	if ( skip_set_reasonable_fold_tree_          <  other.skip_set_reasonable_fold_tree_          ) return true;
	if ( skip_set_reasonable_fold_tree_          != other.skip_set_reasonable_fold_tree_          ) return false;
	if ( set_up_metal_bonds_                     <  other.set_up_metal_bonds_                     ) return true;
	if ( set_up_metal_bonds_                     != other.set_up_metal_bonds_                     ) return false;
	if ( metal_bond_LJ_multiplier_               <  other.metal_bond_LJ_multiplier_               ) return true;
	if ( metal_bond_LJ_multiplier_               != other.metal_bond_LJ_multiplier_               ) return false;
	if ( metal_bond_dist_constraint_multiplier_  <  other.metal_bond_dist_constraint_multiplier_  ) return true;
	if ( metal_bond_dist_constraint_multiplier_  != other.metal_bond_dist_constraint_multiplier_  ) return false;
	if ( metal_bond_angle_constraint_multiplier_ <  other.metal_bond_angle_constraint_multiplier_ ) return true;
	if ( metal_bond_angle_constraint_multiplier_ != other.metal_bond_angle_constraint_multiplier_ ) return false;
	if ( residue_type_set_                        <  other.residue_type_set_                       ) return true;
	//if ( residue_type_set_                        == other.residue_type_set_                       ) return false;
	return false;
}

} // namespace import_pose
} // namespace core
