// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/RNA_KinematicCloser.hh
/// @brief Close a RNA loop with Kinematic Closer (KIC).
/// @author Rhiju Das
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_rotamer_sampler_rna_RNA_KinematicCloser_HH
#define INCLUDED_protocols_rotamer_sampler_rna_RNA_KinematicCloser_HH

// Unit headers
#include <protocols/rotamer_sampler/rna/RNA_KinematicCloser.fwd.hh>

// Package headers
#include <protocols/rotamer_sampler/RotamerSized.hh>

// Project headers
#ifdef WIN32
#include <core/id/DOF_ID.hh>
#include <core/id/NamedAtomID.hh>
#else
#include <core/id/NamedAtomID.fwd.hh>
#endif

#include <core/id/DOF_ID.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace rotamer_sampler {
namespace rna {

/// @brief The RNA de novo structure modeling protocol
class RNA_KinematicCloser: public RotamerSized {
public:
	RNA_KinematicCloser(
		core::pose::PoseOP const & ref_pose,
		core::Size const moving_suite,
		core::Size const chainbreak_suite
	);

	~RNA_KinematicCloser();

	/// @brief Class name
	std::string get_name() const { return "RNA_KinematicCloser"; }

	/// @brief Type of class (see enum in RotamerTypes.hh)
	virtual RotamerType type() const { return RNA_KINEMATIC_CLOSER; }

	/// @brief Initialization
	void init();

	using RotamerSized::apply;

	/// @brief Apply the i-th rotamer to pose
	virtual void apply( core::pose::Pose & pose, core::Size const id );

	/// @brief Get the total number of rotamers in sampler
	core::Size size() const { return nsol_; }

	/// @brief Set verbose
	void set_verbose( bool const setting ) { verbose_ = setting; }

private:

	void figure_out_dof_ids_and_offsets();

	void figure_out_offset( core::id::DOF_ID const & dof_id,
	  core::Real const original_torsion_value );

	void fill_chainTORS();

	void output_chainTORS() const;

	//Disable copy constructor and assignment
  RNA_KinematicCloser( const RNA_KinematicCloser & );
  void operator=( const RNA_KinematicCloser & );

	core::pose::PoseOP ref_pose_;
	core::Size const moving_suite_, chainbreak_suite_;
	bool verbose_;
	core::Size nsol_;

	utility::vector1< core::id::NamedAtomID > atom_ids_;
	utility::vector1< core::Real > offset_save_;
	utility::vector1< core::id::DOF_ID > dof_ids_;

	utility::vector1< utility::vector1< core::Real > > t_ang_, b_ang_, b_len_;
	utility::vector1< utility::vector1< core::Real > > atoms_;
	utility::vector1< core::Real > dt_ang_, db_len_, db_ang_;
};


}
}
}

#endif
