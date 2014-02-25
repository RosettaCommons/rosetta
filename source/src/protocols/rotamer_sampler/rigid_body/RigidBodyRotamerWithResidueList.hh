// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueList.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_rigid_body_RigidBodyRotamerWithResidueList_HH
#define INCLUDED_protocols_rotamer_sampler_rigid_body_RigidBodyRotamerWithResidueList_HH

#include <protocols/rotamer_sampler/RotamerComb.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueList.fwd.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamer.fwd.hh>
#include <protocols/rotamer_sampler/copy_dofs/ResidueListRotamer.fwd.hh>
#include <protocols/rotamer_sampler/RotamerOneValueComb.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>

namespace protocols {
namespace rotamer_sampler {
namespace rigid_body {

	class RigidBodyRotamerWithResidueList: public RotamerComb {

	public:

		//constructor
RigidBodyRotamerWithResidueList( copy_dofs::ResidueListRotamerOP copy_dofs_rotamer,
																		 RigidBodyRotamerOP rigid_body_rotamer );

		//destructor
		~RigidBodyRotamerWithResidueList();

	public:

		/// @brief Apply the current rotamer to pose
		virtual void apply( core::pose::Pose & pose );

		void
		fast_forward();

		void
		fast_forward_to_next_rigid_body();

		void
		fast_forward_to_next_translation();

		void
		fast_forward_to_next_euler_gamma();

		ValueList const & get_rigid_body_values();

		// from rigid body rotamer
		core::kinematics::Stub get_stub();

		// from residue list rotamer
		core::conformation::ResidueOP get_residue_at_origin();

copy_dofs::ResidueListRotamerOP copy_dofs_rotamer();
		RigidBodyRotamerOP rigid_body_rotamer();


		/// @brief Name of the class
		virtual std::string get_name() const { return "RigidBodyRotamerWithResidueList"; }

		/// @brief Type of class (see enum in RotamerTypes.hh)
		virtual RotamerType type() const { return RIGID_BODY_WITH_RESIDUE_LIST; }

	private:

copy_dofs::ResidueListRotamerOP copy_dofs_rotamer_;
		RigidBodyRotamerOP rigid_body_rotamer_;

	};

} //rigid_body
} //rotamer_sampler
} //protocols

#endif
