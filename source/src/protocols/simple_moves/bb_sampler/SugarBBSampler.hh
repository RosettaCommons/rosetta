// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/SugarBBSampler.hh
/// @brief Sample dihdrals using sugar bb data.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) and Jason William Labonte (JWLabonte@jhu.edu)


#ifndef INCLUDED_protocols_simple_moves_bb_sampler_SugarBBSampler_hh
#define INCLUDED_protocols_simple_moves_bb_sampler_SugarBBSampler_hh

#include <protocols/simple_moves/bb_sampler/SugarBBSampler.fwd.hh>
#include <protocols/simple_moves/bb_sampler/BBDihedralSampler.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace simple_moves {
namespace bb_sampler {

///@brief Sample dihdrals using sugar bb data.
class SugarBBSampler : public BBDihedralSampler {

public:

	SugarBBSampler();

	//Constructor with options.  Sampling step size is for initialization of the sampling data - here we have phi/psi probability and angles at every .1 degrees by default.
	SugarBBSampler( core::id::MainchainTorsionType torsion_type,
		BBSampleType sampling_type = probability,
		core::Real sampling_step_size = .1);

	SugarBBSampler(SugarBBSampler const & src);

	~SugarBBSampler();

	SugarBBSamplerOP
	clone() const;

public:

	core::Real
	get_torsion(core::pose::Pose const & pose, Size resnum) const;
	///Set torsions to pose


	void
	set_torsion_to_pose(core::pose::Pose & pose, Size resnum) const;

	std::string
	name() const {
		return "SugarBBSampler";
	};

private:

	core::Real sampling_step_size_;

};

} //bb_sampler
} //simple_moves
} //protocols



#endif //INCLUDED_protocols_simple_moves_bb_sampler_SugarBBSampler_hh





