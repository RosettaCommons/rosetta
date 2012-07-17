// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/import_pose/import_pose_options.hh
///
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_core_import_pose_import_pose_options_HH
#define INCLUDED_core_import_pose_import_pose_options_HH


// Unit headers
#include <core/import_pose/import_pose_options.fwd.hh>
#include <core/io/pdb/pdb_dynamic_reader_options.hh>

// C++ headers
#include <string>


namespace core {
namespace import_pose {

class ImportPoseOptions : public io::pdb::PDB_DReaderOptions
{
public:
	ImportPoseOptions();
	
	virtual ~ImportPoseOptions();
	
	virtual
	void parse_my_tag( utility::tag::TagPtr tag );
	
	virtual
	std::string type() const { return "import_pose_options"; }
	
	// accessors
	bool centroid() const { return centroid_; }
	bool fold_tree_io() const { return fold_tree_io_; }
	bool no_optH() const { return no_optH_; }
	bool pack_missing_sidechains() const { return pack_missing_sidechains_; }
	bool read_fold_tree() const { return read_fold_tree_; }
	bool rna() const { return rna_; }
	bool skip_set_reasonable_fold_tree() const { return skip_set_reasonable_fold_tree_; }
	
	std::string const & residue_type_set() const { return residue_type_set_; }

	// mutators
	void set_centroid( bool centroid ) { centroid_ = centroid; }
	void set_fold_tree_io( bool fold_tree_io ) { fold_tree_io_ = fold_tree_io; }
	void set_no_optH( bool no_optH ) { no_optH_ = no_optH; }
	void set_pack_missing_sidechains( bool pack_missing_sidechains ) { pack_missing_sidechains_ = pack_missing_sidechains; }
	void set_read_fold_tree( bool read_fold_tree ) { read_fold_tree_ = read_fold_tree; }
	void set_rna( bool rna ) { rna_ = rna; }
	void set_skip_set_reasonable_fold_tree_( bool skip_set_reasonable_fold_tree ) { skip_set_reasonable_fold_tree_ = skip_set_reasonable_fold_tree; }
	
	void set_residue_type_set( std::string const & residue_type_set ) { residue_type_set_ = residue_type_set; }

private:
	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options();

private:
	bool centroid_;
	bool fold_tree_io_;
	bool no_optH_;
	bool pack_missing_sidechains_;
	bool read_fold_tree_;
	bool rna_;
	bool skip_set_reasonable_fold_tree_;
	
	std::string residue_type_set_;
	
};

} // namespace import_pose
} // namespace core

#endif // INCLUDED_core_import_pose_import_pose_options_HH
