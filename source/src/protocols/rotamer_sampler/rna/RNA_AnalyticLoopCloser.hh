// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_AnalyticLoopCloser.hh
/// @brief Simply close the RNA loop with KIC.
/// @detailed
/// @author Rhiju Das, Fang-Chieh Chou

#ifndef INCLUDED_protocols_rotamer_sampler_rna_RNA_AnalyticLoopCloser_HH
#define INCLUDED_protocols_rotamer_sampler_rna_RNA_AnalyticLoopCloser_HH

#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <core/types.hh>

//// C++ headers
#include <string>

// AUTO-REMOVED #include <vector>

namespace protocols {
namespace rotamer_sampler {
namespace rna {

/// @brief The RNA de novo structure modeling protocol
class RNA_AnalyticLoopCloser: public RotamerSized {
public:
	/// @brief Construct the protocol object
	RNA_AnalyticLoopCloser(
			core::pose::Pose & pose,
			core::Size const moving_suite,
			core::Size const chainbreak_suite );

	/// @brief Class name
	std::string get_name() const { return "RNA_AnalyticLoopCloser"; }

	/// @brief Initialization
	void init();

	/// @brief Reset to the first (or random if is_random()) rotamer
	void reset();

	/// @brief Move to next rotamer
	void operator++();

	/// @brief Check if there are more rotamers available
	bool not_end() const;

	/// @brief Apply the current rotamer to the pose.
	void apply( core::pose::Pose & pose ) { apply( pose, id_ ); }

	/// @brief Apply the i-th rotamer to pose
	void apply( core::pose::Pose & pose, core::Size const id );

	/// @brief Get the total number of rotamers in sampler
	core::Size size() const { return nsol_; }

	/// @brief Set the reference pose
	void set_verbose( bool const setting ) { verbose_ = setting; }

	/// @brief Set the reference pose
	void set_ref_pose( core::pose::Pose const & pose ) {
		ref_pose_ = pose;
		set_init( false );
	}

private:

	void
	figure_out_dof_ids_and_offsets ( core::pose::Pose const & pose,
	                                 utility::vector1< core::Real > const & dt_ang );

	void
	figure_out_offset (
	  core::pose::Pose const & pose,
	  core::id::DOF_ID const & dof_id,
	  core::Real const & original_torsion_value,
	  utility::vector1< core::Real > & offset_save );

	void
	fill_chainTORS (
	  core::pose::Pose const & pose,
	  utility::vector1< core::id::NamedAtomID > const & atom_ids,
	  utility::vector1< utility::vector1< core::Real > > & atoms,
	  utility::vector1< core::Real > & dt_ang,
	  utility::vector1< core::Real > & db_ang,
	  utility::vector1< core::Real > & db_len ) const;

	void
	output_chainTORS (
		utility::vector1< core::Real > const & dt_ang,
		utility::vector1< core::Real > const & db_ang,
		utility::vector1< core::Real > const & db_len ) const;

private:

	bool const verbose_;
	core::Size const moving_suite_, chainbreak_suite_;
	core::pose::Pose const & ref_pose_;
	core::Size nsol_, id_;

	utility::vector1< core::id::NamedAtomID > atom_ids_;
	utility::vector1< core::Real > offset_save_;
	utility::vector1< core::id::DOF_ID > dof_ids_;

	utility::vector1< utility::vector1< core::Real > > t_ang_, b_ang_, b_len_;
};


}
}
}

#endif
