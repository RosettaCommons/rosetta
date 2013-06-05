// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file IO-functionality for enzyme Constraints
/// @brief
/// @author Florian Richter, floric@u.washington.edu


#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_EnzCstTemplateRes_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_EnzCstTemplateRes_hh



// Unit headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>
// Package headers

// Project headers
//#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>


// Utility Headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

//Utility Headers

// C++ Headers
#include <map>

#include <core/id/AtomID.fwd.hh>
#include <utility/vector1.hh>
#include <string>

#ifdef WIN32
	#include <core/id/AtomID.hh>
	#include <core/chemical/ResidueType.hh>
#endif


namespace protocols {
namespace toolbox {
namespace match_enzdes_util {


/// @brief helper class for EnzCstTemplateRes, holds atom ids corresponding
/// @brief to the residue types in a certain pose
class EnzCstTemplateResAtoms : public utility::pointer::ReferenceCount {

public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~EnzCstTemplateResAtoms();

	friend class EnzCstTemplateRes; //::check_data_consistency(pose::Pose const & pose);
	friend class EnzConstraintParameters; //::add_constraints_to_cst_set(pose::Pose & pose);

	void
	remap_resid( core::id::SequenceMapping const & smap );

private:

	void
	remap_atomid_vector(
		core::id::SequenceMapping const & smap,
		std::vector< core::id::AtomID > & atomid_vec
	);

	std::vector< core::id::AtomID > atom1_;
	std::vector< core::id::AtomID > atom2_;
	std::vector< core::id::AtomID > atom3_;
};


/// @brief helper class for class EnzConstraintParameters, gathers information
/// @brief from cst input and pdb input
class EnzCstTemplateRes : public utility::pointer::ReferenceCount {

	//friend class EnzConstraintParameters;
public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~EnzCstTemplateRes();

	typedef std::map< core::chemical::ResidueTypeCOP, utility::vector1< utility::vector1< core::Size > > > RestypeToTemplateAtomsMap;

public:
	typedef std::map< core::chemical::ResidueTypeCOP, utility::vector1< utility::vector1< core::Size > > > AtomIndsForRestypeMap;

public:

	EnzCstTemplateRes(
		core::chemical::ResidueTypeSetCAP src_restype_set
	);

	EnzCstTemplateRes(
		core::chemical::ResidueTypeSetCAP src_restype_set,
		EnzConstraintParametersCAP src_enzio_param );

	EnzCstTemplateRes(
		EnzCstTemplateResCOP other,
		EnzConstraintParametersCAP new_ref_param
	);

	void
	read_params(std::istringstream & line_stream);

	void
	show_params() const;

	void
	get_pose_data(core::pose::Pose & pose) const;

	void
	set_enzio_param( EnzConstraintParametersCAP src_enz_io_param ) {
		enz_io_param_ = src_enz_io_param;
	}

	void
	set_param_index( core::Size index ){
		param_index_ = index; }

	core::Size
	param_index() const{
		return param_index_; }

	EnzCstTemplateResAtomsCOP
	get_template_atoms_at_pos(core::pose::Pose const & pose, core::Size seqpos ) const;

	bool
	rb_minimizable() const {
		return rb_minimizable_; }

	bool
	is_backbone() const {
		return is_backbone_; }

	core::Size
	corresponding_res_block() const {
		return corresponding_res_block_; }

	bool
	find_in_pose_if_missing_from_header( core::pose::Pose & pose);

	void
	set_external_position( core::Size resnum );


	utility::vector1< std::string > const &
	allowed_res_types() const {
		return allowed_res_types_;
	}

	RestypeToTemplateAtomsMap::const_iterator
	atom_inds_for_restype_begin() const {
		return atom_inds_for_restype_.begin(); }

	RestypeToTemplateAtomsMap::const_iterator
	atom_inds_for_restype_end() const {
		return atom_inds_for_restype_.end(); }

	void
	clear_all();

	void
	remap_resid( core::id::SequenceMapping const & smap );

	bool
	compatible_restype( core::chemical::ResidueTypeCOP restype ) const;

	void
	determine_atom_inds_for_restype( core::chemical::ResidueTypeCOP restype ) const;

	utility::vector1< core::Size > const &
	atom_inds_for_restype(
		core::Size template_atom,
		core::chemical::ResidueTypeCOP restype
  ) const;

	/// @brief checks whether the distance
	/// between any of the template res atoms
	/// of res1 and res2 is below a certain
	/// cutoff.
	bool
	residue_conformations_redundant(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2
	) const;

	void
	identical_info_consistency_check() const;

private:

	//information from cst file
	utility::vector1< std::string > atom1_, atom2_, atom3_;   //names
	std::string at1_type_, at2_type_, at3_type_;    //types

	//contains the names of allowed res types as specified by cst file
	utility::vector1< std::string > allowed_res_types_;

	//contains all allowed atoms for a certain residue type
	//represents the translation of the above cst file info (allowed_res_types_, atom1_, at1_type_, etc )
	//into workable data
	mutable RestypeToTemplateAtomsMap atom_inds_for_restype_;

	//if this is a ligand, do we want to allow rb minimization?
	//used i.e. in cases when there are two ligands, one of them should stay fixed
	bool rb_minimizable_;

	//if this residue interaction is mediated through backbone only
	//right now can only be set through file, but should be changed
	//to automatic detection in the future.
	bool is_backbone_;

	utility::vector1< core::Size > respos_from_external_;

	//if there is another residue that's identical to this
	bool identical_tag_found_;
	core::Size corresponding_res_block_, corresponding_res_num_in_block_;

	//residue set that we are working with, needed to look up allowed res types
	core::chemical::ResidueTypeSetCAP restype_set_;

	EnzConstraintParametersCAP enz_io_param_; // the params object that this residue belongs to
	core::Size param_index_; //index in CstParams, resA=1, resB =2


}; //class EnzCstTemplateRes


}
} //protocols
} //enzdes


#endif
