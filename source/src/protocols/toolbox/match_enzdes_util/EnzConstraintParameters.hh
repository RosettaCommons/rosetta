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


#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_EnzConstraintParameters_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_EnzConstraintParameters_hh


// Unit headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>

// Package headers
#include <core/scoring/func/Func.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>


// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <set>

#include <core/types.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.fwd.hh>
#include <utility/vector1.hh>
#include <iostream>



//Utility Headers

// C++ Headers

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {


////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/// @brief helper class to allow for removal of covalent constraints
class CovalentConnectionReplaceInfo : public utility::pointer::ReferenceCount {

public:

	CovalentConnectionReplaceInfo(
		std::string resA_base_in,
			std::string resB_base_in,
			std::string resA_var_in,
			std::string resB_var_in,
			core::Size Apos_in,
			core::Size Bpos_in,
			core::chemical::ResidueTypeSetCAP restype_set_in
		);

	CovalentConnectionReplaceInfo( CovalentConnectionReplaceInfo const & other );

	virtual ~CovalentConnectionReplaceInfo();

	void
	remove_covalent_connection_from_pose( core::pose::Pose & pose ) const;

	void
	remap_resid( core::id::SequenceMapping const & smap );

private:

	std::string resA_basename_, resB_basename_;
	std::string resA_varname_, resB_varname_;
	std::string resA_modname_, resB_modname_;

	core::Size resA_seqpos_, resB_seqpos_;
	core::chemical::ResidueTypeSetCAP restype_set_;

}; //class CovalentConnection


/// @brief class that holds all the parameters for one specific constraint
class EnzConstraintParameters : public utility::pointer::ReferenceCount {


public:

	void
	set_mcfi(
		toolbox::match_enzdes_util::MatchConstraintFileInfoCOP mcfi );

	EnzConstraintParameters();

	EnzConstraintParameters(
		core::Size cst_block,
		core::chemical::ResidueTypeSetCAP src_restype_set,
		EnzConstraintIOCAP src_enz_io);

	EnzConstraintParameters( EnzConstraintParameters const & other );

	virtual ~EnzConstraintParameters();

	void
	show_definitions() const;

	void
	generate_pose_specific_data(
		core::pose::Pose & pose,
		core::scoring::ScoreFunctionCOP scofx
	) const;

	bool
	is_empty() const {
		return empty_;
	}

	bool
	is_covalent() const {
		return is_covalent_;
	}

	bool
	missing_in_pose( core::pose::Pose const & pose) const;

	/// @brief function that takes all rotamers for the ResidueType(s)
	/// of the residue that's missing in the pose and places them
	/// according to the geometry specified in the mcfi
	//utility::vector1< core::conformation::ResidueCOP >
	//inverse_rotamers_for_residue_missing_in_pose(
	//	core::pose::Pose const & pose ) const;

	core::Size
	cst_block() const {
		return cst_block_;}

	void
	set_cst_block( core::Size cst_block ){ cst_block_ = cst_block;}

	/// @brief updates the pdb remarks according to what is in the
	/// EnzCstTemplateRes member data. returns false in case any
	/// error occured
	bool
	update_pdb_remarks(
		core::pose::Pose & pose
	) const;


	EnzConstraintIOCAP
	enz_io() const{
		return enz_io_; }

	EnzCstTemplateResOP
	nonconst_resA();

	EnzCstTemplateResOP
	nonconst_resB();

	EnzCstTemplateResCOP
	resA() const;

	EnzCstTemplateResCOP
	resB() const;

	EnzCstTemplateResCOP
	get_missing_template_res( core::pose::Pose const & pose ) const;

	EnzCstTemplateResCOP
	get_missing_template_other_res( core::pose::Pose const & pose ) const;

	/// @brief all residue names specified in the cstfile
	/// returns an empty set if the constraints don't apply
	/// to the specifed position
	std::set< std::string >
	allowed_res_name3_at_position( core::pose::Pose const & pose, core::Size seqpos ) const;

	void
	set_external_position_for_resA( core::Size pos );

	void
	set_external_position_for_resB( core::Size pos );

	void
	remove_covalent_connections_from_pose( core::pose::Pose & pose ) const;

	void
	remap_resid( core::id::SequenceMapping const & smap );

	static
	core::scoring::func::FuncOP
	convert_GeomSampleInfo_to_FuncOP(
		toolbox::match_enzdes_util::GeomSampleInfoCOP gsi,
		core::Real & ideal_val);

private:

	void
	generate_active_pose_constraints(
		core::pose::Pose & pose,
		core::scoring::ScoreFunctionCOP scofx
	) const;

	void
	make_constraint_covalent(
		core::pose::Pose & pose,
		Size resA_pos,
		Size resB_pos,
		Size resA_At,
		Size resB_At
	) const;

	void
	make_constraint_covalent_helper(
		core::pose::Pose & pose,
		EnzCstTemplateResOP template_res,
		core::Size res_pos,
		core::Size Atpos,
		core::Real itorsion,
		core::Real iangle,
		core::Real idis,
	 	std::string & res_varname
	) const;

	core::Size
	determine_best_constraint(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunctionCOP scofx,
		utility::vector1< core::scoring::constraints::ConstraintCOP > candidate_csts
	) const;


private: //data

	EnzCstTemplateResOP resA_, resB_;

	toolbox::match_enzdes_util::MatchConstraintFileInfoCOP mcfi_;

	core::scoring::func::FuncOP disAB_, angleA_, angleB_, torsionA_, torsionB_, torsionAB_;

	core::Real ndisAB_, nangleA_, nangleB_, ntorsionA_, ntorsionB_, ntorsionAB_;

	bool is_covalent_;

	mutable bool empty_; //safeguard if someone reads in a cstfile with all force constants set to 0
	core::chemical::ResidueTypeSetCAP restype_set_;

	EnzConstraintIOCAP enz_io_; //the EnzConstraintIO object that this template res belongs to
	core::Size cst_block_; //the block in the cst file that corresponds to this residue

}; // class EnzConstraintParameters


}
} //enzdes
} //protocols


#endif
