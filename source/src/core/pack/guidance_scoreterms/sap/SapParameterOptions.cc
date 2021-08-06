// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/sap/SapParameterOptions.hh
/// @brief Parameters for sap score
/// @details Contains options that change the fundamental parameters of the sap score
/// @author Brian Coventry (bcov@uw.edu)

// Unit headers
#include <core/pack/guidance_scoreterms/sap/SapParameterOptions.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/pose/Pose.hh>

// Options system
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// File I/O
#include <basic/database/open.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>

// Other Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/memory.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <cmath>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

static basic::Tracer TR("core.pack.guidance_scoreterms.sap.SapParameterOptions");


SapParameterOptions::SapParameterOptions() :
	utility::VirtualBase(),
	hydrop_lys_arg_setting( utility::get_undefined_real() ),
	hydrop_adder( 0 )
{}




} //sap
} //guidance_scoreterms
} //pack
} //core


#ifdef SERIALIZATION

template< class Archive >
void
core::pack::guidance_scoreterms::sap::SapParameterOptions::save( Archive & arc ) const {
	arc( CEREAL_NVP( hydrop_lys_arg_setting ) );
	arc( CEREAL_NVP( hydrop_adder ) );
}

template< class Archive >
void
core::pack::guidance_scoreterms::sap::SapParameterOptions::load( Archive & arc ) {
	arc( hydrop_lys_arg_setting );
	arc( hydrop_adder );
}


SAVE_AND_LOAD_SERIALIZABLE( core::pack::guidance_scoreterms::sap::SapParameterOptions );
CEREAL_REGISTER_TYPE( core::pack::guidance_scoreterms::sap::SapParameterOptions )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_guidance_scoreterms_sap_SapParameterOptions )
#endif // SERIALIZATION
