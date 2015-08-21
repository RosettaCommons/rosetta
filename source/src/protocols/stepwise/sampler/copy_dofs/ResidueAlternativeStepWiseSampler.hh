// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSampler.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_copy_dofs_ResidueAlternativeStepWiseSampler_HH
#define INCLUDED_protocols_sampler_copy_dofs_ResidueAlternativeStepWiseSampler_HH

#include <protocols/stepwise/sampler/copy_dofs/CopyDofStepWiseSampler.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeSet.fwd.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSampler.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.fwd.hh>

namespace protocols {
namespace stepwise {
namespace sampler {
namespace copy_dofs {

class ResidueAlternativeStepWiseSampler: public CopyDofStepWiseSampler {

public:

	//constructor
	ResidueAlternativeStepWiseSampler( ResidueAlternativeSet const & residue_alternative_set,
		core::pose::Pose const & starting_pose );

	//constructor
	ResidueAlternativeStepWiseSampler( ResidueAlternativeSet const & residue_alternative_set );

	//constructor
	ResidueAlternativeStepWiseSampler( utility::vector1< core::pose::PoseOP > const & pose_list,
		std::map< Size, Size > const & res_map,
		Size const representative_seqpos,
		core::pose::Pose const & starting_pose );

	//constructor
	ResidueAlternativeStepWiseSampler( utility::vector1< core::pose::PoseOP > const & pose_list,
		std::map< Size, Size > const & res_map,
		Size const representative_seqpos );

	//constructor
	ResidueAlternativeStepWiseSampler( utility::vector1< core::pose::PoseOP > const & pose_list,
		Size const seqpos );

	~ResidueAlternativeStepWiseSampler();

public:

	core::conformation::Residue const &
	get_residue_at_origin();

	core::conformation::Residue const &
	get_residue_at_origin_with_matching_type( core::conformation::Residue const & rsd_in );

	/// @brief Name of the class
	virtual std::string get_name() const;

	/// @brief Type of class (see enum in StepWiseSamplerTypes.hh)
	virtual StepWiseSamplerType type() const { return RESIDUE_ALTERNATIVE; }

	Size representative_seqpos() const { return representative_seqpos_; }

private:

	std::map< Size, Size >
	simple_res_map( Size const i );

	void
	initialize_residues();

	void
	initialize_residues_for_type( core::conformation::Residue const & rsd_in );

private:

	Size const representative_seqpos_;
	std::map< std::string, utility::vector1< core::conformation::ResidueOP > > residues_for_each_type_;
	std::string original_type_;

};

} //copy_dofs
} //sampler
} //stepwise
} //protocols

#endif
