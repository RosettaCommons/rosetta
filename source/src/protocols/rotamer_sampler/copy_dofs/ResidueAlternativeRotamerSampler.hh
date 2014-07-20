// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/copy_dofs/ResidueAlternativeRotamerSampler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_copy_dofs_ResidueAlternativeRotamerSampler_HH
#define INCLUDED_protocols_rotamer_sampler_copy_dofs_ResidueAlternativeRotamerSampler_HH

#include <protocols/rotamer_sampler/copy_dofs/CopyDofRotamerSampler.hh>
#include <protocols/rotamer_sampler/copy_dofs/ResidueAlternativeSet.fwd.hh>
#include <protocols/rotamer_sampler/copy_dofs/ResidueAlternativeRotamerSampler.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.fwd.hh>

using namespace core::conformation;
using namespace core::pose;

namespace protocols {
namespace rotamer_sampler {
namespace copy_dofs {

	class ResidueAlternativeRotamerSampler: public CopyDofRotamerSampler {

	public:

		//constructor
		ResidueAlternativeRotamerSampler( ResidueAlternativeSet const & residue_alternative_set,
															 core::pose::Pose const & starting_pose );

		//constructor
		ResidueAlternativeRotamerSampler( ResidueAlternativeSet const & residue_alternative_set );

		//constructor
		ResidueAlternativeRotamerSampler( utility::vector1< core::pose::PoseOP > const & pose_list,
															 std::map< Size, Size > const & res_map,
															 Size const representative_seqpos,
															 core::pose::Pose const & starting_pose );

		//constructor
		ResidueAlternativeRotamerSampler( utility::vector1< core::pose::PoseOP > const & pose_list,
															 std::map< Size, Size > const & res_map,
															 Size const representative_seqpos );

		//constructor
		ResidueAlternativeRotamerSampler( utility::vector1< core::pose::PoseOP > const & pose_list,
															 Size const seqpos );

		~ResidueAlternativeRotamerSampler();

	public:

		Residue const &
		get_residue_at_origin();

		Residue const &
		get_residue_at_origin_with_matching_type( Residue const & rsd_in );

		/// @brief Name of the class
		virtual std::string get_name() const { return "ResidueAlternative"; }

		/// @brief Type of class (see enum in RotamerSamplerTypes.hh)
		virtual RotamerSamplerType type() const { return RESIDUE_ALTERNATIVE; }

		Size representative_seqpos() const { return representative_seqpos_; }

	private:

		std::map< Size, Size >
		simple_res_map( Size const i );

		void
		initialize_residues();

		void
		initialize_residues_for_type( Residue const & rsd_in );

	private:

		Size const representative_seqpos_;
		std::map< std::string, utility::vector1< ResidueOP > > residues_for_each_type_;
		std::string original_type_;

	};

} //copy_dofs
} //rotamer_sampler
} //protocols

#endif
