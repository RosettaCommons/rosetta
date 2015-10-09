// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/AnnealerFactory.fwd.hh
/// @brief  Annealer Factory class forward declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_annealer_AnnealerFactory_hh
#define INCLUDED_core_pack_annealer_AnnealerFactory_hh

/// Unit headers
#include <core/pack/annealer/AnnealerFactory.fwd.hh>

/// Package headers
#include <core/pack/annealer/SimAnnealerBase.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/rotamer_set/FixbbRotamerSets.fwd.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>


/// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1.fwd.hh>

/// C++ headers
#include <utility/vector0.hh> // damn shame
#include <core/types.hh>

namespace core {
namespace pack {
namespace annealer {

class AnnealerFactory
{
public:

	static
	SimAnnealerBaseOP
	create_annealer(
		task::PackerTaskCOP task,
		utility::vector0<int> & rot_to_pack,
		ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
		core::PackerEnergy & bestenergy,
		bool start_with_current, // start simulation with current rotamers
		interaction_graph::AnnealableGraphBaseOP ig,
		rotamer_set::FixbbRotamerSetsCOP rotamer_sets,
		ObjexxFCL::FArray1_int & current_rot_index,
		bool calc_rot_freq,
		ObjexxFCL::FArray1D< core::PackerEnergy > & rot_freq
	);

};

}// namespace annealer
}// namespace pack
}// namespace core

#endif
