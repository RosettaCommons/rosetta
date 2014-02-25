// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueAlternatives.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_rigid_body_RigidBodyRotamerWithResidueAlternatives_HH
#define INCLUDED_protocols_rotamer_sampler_rigid_body_RigidBodyRotamerWithResidueAlternatives_HH

#include <protocols/rotamer_sampler/RotamerComb.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueAlternatives.fwd.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamer.hh>
#include <protocols/rotamer_sampler/copy_dofs/ResidueAlternativeRotamerComb.fwd.hh>

using namespace protocols::rotamer_sampler::copy_dofs;

namespace protocols {
namespace rotamer_sampler {
namespace rigid_body {

	class RigidBodyRotamerWithResidueAlternatives: public RotamerComb {

	public:

		//constructor
		RigidBodyRotamerWithResidueAlternatives( ResidueAlternativeRotamerCombOP residue_alternative_rotamer,
																						 RigidBodyRotamerOP rigid_body_rotamer );

		//destructor
		~RigidBodyRotamerWithResidueAlternatives();

	public:

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
		core::conformation::Residue const & get_residue_at_origin();

		// from residue list rotamer
		core::conformation::Residue const & get_residue_at_origin( Size const seqpos );

		ResidueAlternativeRotamerCombOP residue_alternatives_rotamer();
		RigidBodyRotamerOP rigid_body_rotamer();

		/// @brief Name of the class
		virtual std::string get_name() const { return "RigidBodyRotamerWithResidueAlternatives"; }

		/// @brief Type of class (see enum in RotamerTypes.hh)
		virtual RotamerType type() const { return RIGID_BODY_WITH_RESIDUE_ALTERNATIVES; }

		void
		apply_rigid_body_only( pose::Pose & pose );

		void
		fast_forward_to_next_residue_pair( Size const i, Size const j);

		void
		fast_forward_to_next_residue( Size const i );

Vector
get_xyz( Size const seqpos, std::string const atom_name  );

conformation::ResidueCOP
		get_residue( Size const seqpos );

	private:

		ResidueAlternativeRotamerCombOP residue_alternatives_rotamer_;
		RigidBodyRotamerOP rigid_body_rotamer_;

std::map< Size, conformation::ResidueOP > transformed_residues;


	};

} //rigid_body
} //rotamer_sampler
} //protocols

#endif
