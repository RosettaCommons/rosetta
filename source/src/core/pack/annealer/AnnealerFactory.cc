// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/AnnealerFactory.cc
/// @brief  Annealer Factory class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

/// Unit headers
#include <core/pack/annealer/AnnealerFactory.hh>

/// Package headers
#include <core/pack/annealer/FixbbSimAnnealer.hh>
#include <core/pack/annealer/FixbbCoupledRotamerSimAnnealer.hh>
#include <core/pack/annealer/FixbbLinkingRotamerSimAnnealer.hh>
#include <core/pack/annealer/MultiCoolAnnealer.hh>

#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>

#include <core/pack/task/PackerTask.hh>

#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace annealer {

static THREAD_LOCAL basic::Tracer TR( "core.pack.annealer.AnnealerFactory" );

SimAnnealerBaseOP
AnnealerFactory::create_annealer(
	task::PackerTaskCOP task,
	utility::vector0<int> & rot_to_pack,
	ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	bool start_with_current,
	interaction_graph::InteractionGraphBaseOP ig,
	rotamer_set::FixbbRotamerSetsCOP rotamer_sets,
	ObjexxFCL::FArray1_int & current_rot_index,
	bool calc_rot_freq,
	ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq
)
{
	if ( task->rotamer_couplings_exist() ) {
		TR.Debug << "Creating FixbbCoupledRotamerSimAnnealer" << std::endl;
		return SimAnnealerBaseOP( new FixbbCoupledRotamerSimAnnealer(
			rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, ig,
			rotamer_sets, current_rot_index, calc_rot_freq, rot_freq,
			task->rotamer_couplings() ) );
	} else if ( task->rotamer_links_exist() ) {
		TR.Debug << "Creating FixbbLinkingRotamerSimAnnealer" << std::endl;
		return SimAnnealerBaseOP( new FixbbLinkingRotamerSimAnnealer(
			rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, ig,
			rotamer_sets, current_rot_index, calc_rot_freq, rot_freq,
			task->rotamer_links() ) );
	} else if ( task->multi_cool_annealer() ) {
		TR.Debug << "Creating MultiCoolAnnealer" << std::endl;
		return SimAnnealerBaseOP( new MultiCoolAnnealer(
			task, rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, ig,
			rotamer_sets, current_rot_index, calc_rot_freq, rot_freq ) );
	} else {
		TR.Debug << "Creating FixbbSimAnnealer" << std::endl;
		return SimAnnealerBaseOP( new FixbbSimAnnealer(
			rot_to_pack, bestrotamer_at_seqpos, bestenergy, start_with_current, ig,
			rotamer_sets, current_rot_index, calc_rot_freq, rot_freq ) );
	}

	// appease compiler
	return 0;
}

}// namespace annealer
}// namespace pack
}// namespace core

