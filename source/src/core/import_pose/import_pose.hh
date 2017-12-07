// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/import_pose.hh
///
/// @brief  various functions to construct Pose object(s) from PDB(s)
/// @author Sergey Lyskov

#ifndef INCLUDED_core_import_pose_import_pose_hh
#define INCLUDED_core_import_pose_import_pose_hh

// Package headers
#include <core/import_pose/import_pose_options.fwd.hh>

// C++ headers
#include <iosfwd>

// Utility headers
#include <basic/Tracer.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/io/StructFileRep.fwd.hh>
#include <core/io/StructFileReaderOptions.fwd.hh>
#include <core/pose/full_model_info/FullModelParameters.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <core/import_pose/FullModelPoseBuilder.fwd.hh>

#include <utility/vector1.fwd.hh>

#include <map>

class CifFile;
class CifParser;
typedef utility::pointer::shared_ptr< CifFile > CifFileOP;
typedef utility::pointer::shared_ptr< CifParser > CifParserOP;


namespace core {
namespace import_pose {

enum FileType{
	PDB_file,
	CIF_file,
	Unknown_file
};

std::ostream & operator<<( std::ostream & stream, FileType type );

typedef std::string String;

void
read_all_poses(
	utility::vector1< std::string > const & filenames,
	utility::vector1< core::pose::PoseOP > &  poses
);

void
read_additional_pdb_data(
	std::string const & file,
	pose::Pose & pose,
	io::StructFileRepCOP fd,
	bool read_fold_tree = false
);

void
read_additional_pdb_data(
	std::string const & file,
	pose::Pose & pose,
	ImportPoseOptions const & options,
	bool read_fold_tree = false
);


/// @brief Returns a PoseOP object from the Pose created from input
/// PDB  <filename>
/// @note: in PyRosetta, this will return a Pose object
///
/// example(s):
///     pose = pose_from_file("YFP.pdb")
/// See also:
///     Pose
///     PDBInfo
///     make_pose_from_sequence
///     pose_from_rcsb
///     pose_from_sequence
pose::PoseOP pose_from_file(
	std::string const & filename,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

/// @brief Returns a PoseOP object from the Pose created from input
/// PDB  <filename>, taking a set of custom ImportPoseOptions parameters.
/// @note: in PyRosetta, this will return a Pose object
///
/// example(s):
///     pose = pose_from_file("YFP.pdb")
/// See also:
///     Pose
///     PDBInfo
///     make_pose_from_sequence
///     pose_from_rcsb
///     pose_from_sequence
pose::PoseOP
pose_from_file(
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree,
	FileType type
);
/// @brief Returns a PoseOP object from the Pose created by reading the input
/// PDB  <filename>, this constructor allows for a non-default ResidueTypeSet
/// <residue_set>
pose::PoseOP pose_from_file(
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

/// @brief Determine what file type is passed to function
/// there should only be one function that calls this, pose_from_file
/// and only calls it when the filetype is unknown
FileType
determine_file_type( std::string const &contents_of_file);

void
pose_from_file(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree = false,
	FileType file_type = Unknown_file
);

/// @brief Reads in data from input PDB  <filename>  and stores it in the Pose
/// <pose>, this constructor allows for a non-default ResidueTypeSet
/// <residue_set>
void
pose_from_file(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

/// @brief Reads in data from input PDB  <filename>  and stores it in the Pose
/// <pose>, uses the FA_STANDARD ResidueTypeSet (fullatom) by default
/// @note: will use centroid if in::file::centroid_input is true
///
/// example(s):
///     pose_from_file(pose,"YFP.pdb")
/// See also:
///     Pose
///     PDBInfo
void
pose_from_file(
	pose::Pose & pose,
	std::string const & filename,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

void
pose_from_file(
	pose::Pose & pose,
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

/// @brief Reads data from an input PDB containing multiple models named
/// <filename>  and stores it in a vector of Pose objects named  <poses>  using
/// ResidueTypeSet  <residue_set>
void
pose_from_file(
	utility::vector1< pose::Pose > & poses,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

void
pose_from_file(
	utility::vector1< pose::Pose > & poses,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

utility::vector1< core::pose::PoseOP >
poseOPs_from_files(
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

utility::vector1< core::pose::PoseOP >
poseOPs_from_files(
	utility::vector1< std::string > const & filenames,
	ImportPoseOptions const & options,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

utility::vector1< core::pose::PoseOP >
poseOPs_from_files(
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< std::string > const & filenames,
	ImportPoseOptions const & options,
	bool read_fold_tree,
	FileType type = Unknown_file
);

utility::vector1< core::pose::Pose >
poses_from_files(
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

utility::vector1< core::pose::Pose >
poses_from_files(
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree,
	FileType type = Unknown_file
);


// FA_STANDARD residue set

/// @brief Reads data from an input PDB containing multiple models named
/// <filename>  and stores it in a vector of Pose objects named  <poses>
/// using the FA_STANDARD ResidueTypeSet (fullatom)
void
pose_from_file(
	utility::vector1< pose::Pose > & poses,
	std::string const & filename,
	bool read_fold_tree = false,
	FileType type = Unknown_file
);

void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	std::string const & filename = ""
);

void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	ImportPoseOptions const & options,
	std::string const & filename = ""
);

void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename
);

void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options,
	std::string const & filename
);

void pose_from_pdb_stream(
	pose::Pose & pose,
	std::istream & pdb_stream,
	std::string const & filename,
	ImportPoseOptions const & options
);

// uses the CENTROID residue_set

/// @brief Reads in data from input PDB  <filename>  and stores it in the Pose
/// <pose>  using the CENTROID ResidueTypeSet (centroid)
void
centroid_pose_from_pdb(
	pose::Pose & pose,
	std::string const & filename,
	bool read_fold_tree = false
);

/// @brief Create pose object, using given StructFileRep object.
/// If PDB cleanin specified - it will be applied first.
/// Constructs a ImportPoseOptions object from the command line
void build_pose(
	io::StructFileRepOP fd, // const?  naaaaah
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set
);

/// @brief Create pose object, using given StructFileRep object.
/// If PDB cleanin specified - it will be applied first
void build_pose(
	io::StructFileRepOP fd, // const?  naaaaah
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options
);

/// @brief Create pose object, using given StructFileRep object.
/// No PDB cleanin will be appliend.
void build_pose_as_is(
	io::StructFileRepOP fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options
);

void build_pose_as_is2(
	io::StructFileRepCOP fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	id::AtomID_Mask & missing,
	ImportPoseOptions const & options
);



// Input from command line stuff -- as much as possible.

core::pose::PoseOP
get_pdb_with_full_model_info( std::string const & input_file,
	core::chemical::ResidueTypeSetCAP rsd_set );

core::pose::PoseOP
get_pdb_and_cleanup( std::string const & input_file );

core::pose::PoseOP
get_pdb_and_cleanup( std::string const & input_file,
	core::chemical::ResidueTypeSetCAP rsd_set );

void
get_other_poses(  utility::vector1< core::pose::PoseOP > & other_poses,
	utility::vector1< std::string > const & other_files,
	core::chemical::ResidueTypeSetCAP rsd_set );

core::pose::PoseOP
initialize_pose_and_other_poses_from_command_line( core::chemical::ResidueTypeSetCAP rsd_set );

core::pose::PoseOP
initialize_pose_and_other_poses_from_options( core::chemical::ResidueTypeSetCAP rsd_set, utility::options::OptionCollection const & options );

void
cleanup( core::pose::Pose & pose,
	bool const force_cut_at_rna_chainbreak = false );


core::pose::full_model_info::FullModelParametersOP
get_sequence_information( std::string const & fasta_file,
	utility::vector1< core::Size > & cutpoint_open_in_full_model,
	bool const add_virt = false );

void
setup_for_density_scoring( core::pose::Pose & pose );




void
get_extra_cutpoints_from_names( core::Size const nres,
	utility::vector1< core::Size > & cutpoints,
	std::map< core::Size, std::string > const & non_standard_residue_map );

void
setup_fold_trees(
	utility::vector1< core::pose::Pose * > & pose_pointers,
	utility::vector1< core::Size > & cutpoint_open_in_full_model /* can be updated here*/,
	utility::vector1< core::Size > & fixed_domain_map /* can be updated here*/,
	utility::vector1< core::Size > const & cutpoint_closed,
	utility::vector1< core::Size > const & extra_minimize_res,
	utility::vector1< core::Size > const & extra_minimize_jump_res,
	utility::vector1< core::Size > const & sample_res,
	utility::vector1< core::Size > const & working_res,
	utility::vector1< core::Size > const & jump_res,
	utility::vector1< core::Size > const & preferred_root_res,
	utility::vector1< core::Size > const & virtual_sugar_res,
	core::pose::full_model_info::FullModelParameters const & full_model_parameters,
	utility::vector1< utility::vector1< core::Size > > const & pose_res_lists );

void
update_pose_fold_tree( core::pose::Pose & pose,
	utility::vector1< core::Size > const & res_list,
	utility::vector1< core::Size > const & extra_min_res,
	utility::vector1< core::Size > const & sample_res,
	utility::vector1< core::Size > const & jump_res,
	core::pose::full_model_info::FullModelParameters const & full_model_parameters );

void
define_chains( core::pose::Pose const & pose,
	utility::vector1< utility::vector1< core::Size > > & all_res_in_chain,
	utility::vector1< utility::vector1< core::Size > > & all_fixed_res_in_chain,
	utility::vector1< core::Size > const & res_list,
	utility::vector1< core::Size > const & extra_min_res );

void
setup_user_defined_jumps(
	utility::vector1< core::Size > const & jump_res,
	utility::vector1< core::Size > & jump_partners1,
	utility::vector1< core::Size > & jump_partners2,
	utility::vector1< std::pair< core::Size, core::Size > > & chain_connections,
	utility::vector1< core::Size > const & res_list,
	utility::vector1< utility::vector1< core::Size > > const & all_res_in_chain );

core::Size
get_chain( core::Size const i, utility::vector1< utility::vector1< core::Size > > const & all_res_in_chain );

void
setup_jumps( core::pose::Pose const & pose,
	utility::vector1< core::Size > & jump_partners1,
	utility::vector1< core::Size > & jump_partners2,
	utility::vector1< std::pair< core::Size, core::Size > > & chain_connections,
	utility::vector1< utility::vector1< core::Size > > const & all_res_in_chain,
	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & resnum_and_chain_and_segid_in_pose );

core::kinematics::FoldTree
get_tree( core::pose::Pose const & pose,
	utility::vector1< core::Size > const & cuts,
	utility::vector1< core::Size > const & jump_partners1,
	utility::vector1< core::Size > const & jump_partners2 );

core::kinematics::FoldTree
get_tree( core::Size const nres,
	utility::vector1< core::Size > const & cuts,
	utility::vector1< core::Size > const & jump_partners1,
	utility::vector1< core::Size > const & jump_partners2,
	utility::vector1< std::string > const & jump_atoms1,
	utility::vector1< std::string > const & jump_atoms2 );

void
update_fixed_domain_from_extra_minimize_jump_res( utility::vector1< core::Size > & fixed_domain,
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & res_list,
	utility::vector1< core::Size > const & extra_minimize_jump_res );
void
add_cutpoint_closed( core::pose::Pose & pose,
	utility::vector1< core::Size > const & res_list,
	utility::vector1< core::Size > const & cutpoint_closed );

void
put_in_cutpoint( core::pose::Pose & pose, core::Size const i );

void
add_virtual_sugar_res( core::pose::Pose & pose,
	utility::vector1< core::Size > const & res_list,
	utility::vector1< core::Size > const & virtual_sugar_res );

utility::vector1< core::Size >
figure_out_working_res( utility::vector1< core::Size > const & input_domain_map,
	utility::vector1< core::Size > const & sample_res );

utility::vector1< core::Size >
figure_out_sample_res( utility::vector1< core::Size > const & input_domain_map,
	utility::vector1< core::Size > const & working_res );

void
check_working_res( utility::vector1< core::Size > const & working_res,
	utility::vector1< core::Size > const & input_domain_map,
	utility::vector1< core::Size > const & sample_res );

void
check_extra_minimize_res_are_input( utility::vector1< core::Size > const & extra_minimize_res,
	utility::vector1< core::Size > const & input_domain_map );

void
figure_out_motif_mode( utility::vector1< core::Size > & extra_min_res,
	utility::vector1< core::Size > & terminal_res,
	utility::vector1< core::Size > const & working_res,
	utility::vector1< core::Size > const & input_domain_map,
	utility::vector1< core::Size > const & cutpoint_open_in_full_model );

void
add_block_stack_variants( utility::vector1< core::pose::Pose * > const & pose_pointers,
	utility::vector1< utility::vector1< core::Size > > const & pose_res_lists,
	utility::vector1< core::Size > const & block_stack_above_res,
	utility::vector1< core::Size > const & block_stack_below_res );

void
update_jump_res( utility::vector1< core::Size > & jump_res,
	utility::vector1< core::Size > const & extra_minimize_jump_res );

void
check_extra_minimize_res_are_input( utility::vector1< core::Size > const & extra_minimize_res,
	utility::vector1< core::Size > const & input_domain_map );

utility::vector1< core::Size >
figure_out_fixed_domain_map( utility::vector1< core::Size > const & input_domain_map,
	utility::vector1< core::Size > const & extra_minimize_res );

utility::vector1< core::Size >
figure_out_dock_domain_map( utility::vector1< core::Size > & cutpoint_open_in_full_model,
	utility::vector1< utility::vector1< core::Size > > const & pose_res_lists,
	utility::vector1< core::Size > const & working_res,
	utility::vector1< core::Size > const & sample_res,
	core::Size const nres );

void
reorder_pose( core::pose::Pose & pose, utility::vector1< core::Size > & res_list );

bool
just_modeling_RNA( utility::vector1< std::string > const & fasta_files );

void
look_for_dna( utility::vector1< core::sequence::SequenceOP > & fasta_sequences );


void
setup_water_bank_for_magnesiums( std::map< Size, std::string > & non_standard_residue_map,
	utility::vector1< core::sequence::SequenceOP > & fasta_sequences );

utility::vector1< Size >
get_cutpoints( utility::vector1< core::sequence::SequenceCOP > const & fasta_sequences,
	std::map< Size, std::string > const & non_standard_residue_map,
	utility::vector1< char > const & conventional_chains,
	utility::vector1< int  > const & conventional_numbering,
	utility::vector1< std::string > const & conventional_segids );

} // namespace import_pose
} // namespace core

#endif // INCLUDED_core_import_pose_import_pose_HH
