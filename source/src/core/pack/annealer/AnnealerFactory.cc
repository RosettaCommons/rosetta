// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/annealer/AnnealerFactory.cc
/// @brief  Annealer Factory class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

/// Unit headers
#include <core/pack/annealer/AnnealerFactory.hh>

/// Package headers
#include <core/pack/annealer/FixbbSimAnnealer.hh>
#include <core/pack/annealer/SequenceSymmetricAnnealer.hh>
#include <core/pack/annealer/FixbbCoupledRotamerSimAnnealer.hh>
#include <core/pack/annealer/FixbbLinkingRotamerSimAnnealer.hh>
#include <core/pack/annealer/MultiCoolAnnealer.hh>

#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.fwd.hh>

#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>

namespace core {
namespace pack {
namespace annealer {

static basic::Tracer TR( "core.pack.annealer.AnnealerFactory" );

SimAnnealerBaseOP
AnnealerFactory::create_annealer(
	core::pose::Pose const & pose,
	task::PackerTaskCOP const & task,
	utility::vector0< int > & rot_to_pack,
	ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	bool start_with_current,
	interaction_graph::AnnealableGraphBaseOP const & ig,
	rotamer_set::FixbbRotamerSetsCOP const & rotamer_sets,
	ObjexxFCL::FArray1_int & current_rot_index,
	bool calc_rot_freq,
	ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq
)
{
	if ( task->keep_sequence_symmetry() ) {
		TR.Debug << "Creating SequenceSymmetricAnnealer" << std::endl;
		return utility::pointer::make_shared< SequenceSymmetricAnnealer >(
			pose, rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, ig,
			rotamer_sets, current_rot_index, calc_rot_freq, rot_freq );
	} else if ( task->rotamer_couplings_exist() ) {
		TR.Debug << "Creating FixbbCoupledRotamerSimAnnealer" << std::endl;
		return utility::pointer::make_shared< FixbbCoupledRotamerSimAnnealer >(
			rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, ig,
			rotamer_sets, current_rot_index, calc_rot_freq, rot_freq,
			task->rotamer_couplings() );
	} else if ( task->rotamer_links_exist() ) {
		TR.Debug << "Creating FixbbLinkingRotamerSimAnnealer" << std::endl;
		return utility::pointer::make_shared< FixbbLinkingRotamerSimAnnealer >(
			rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, ig,
			rotamer_sets, current_rot_index, calc_rot_freq, rot_freq,
			task->rotamer_links() );
	} else if ( task->multi_cool_annealer() ) {
		TR.Debug << "Creating MultiCoolAnnealer" << std::endl;
		return utility::pointer::make_shared< MultiCoolAnnealer >(
			task, rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, ig,
			rotamer_sets, current_rot_index, calc_rot_freq, rot_freq );
	} else {
		TR.Debug << "Creating FixbbSimAnnealer" << std::endl;
		return utility::pointer::make_shared< FixbbSimAnnealer >(
			rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, ig,
			rotamer_sets, current_rot_index, calc_rot_freq, rot_freq );
	}

	// appease compiler
	return nullptr;
}

}// namespace annealer
}// namespace pack
}// namespace core

