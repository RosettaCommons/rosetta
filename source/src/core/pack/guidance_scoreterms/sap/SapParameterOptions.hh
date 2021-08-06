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



#ifndef INCLUDED_core_pack_guidance_scoreterms_sap_SapParameterOptions_hh
#define INCLUDED_core_pack_guidance_scoreterms_sap_SapParameterOptions_hh

// Unit headers
#include <core/pack/guidance_scoreterms/sap/SapParameterOptions.fwd.hh>

// Package headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueProperty.hh>

// Project headers
#include <core/types.hh>
#include <map>
#include <string>
#include <utility/vector1.hh>
#include <math.h>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {


class SapParameterOptions : public utility::VirtualBase {
public:

	SapParameterOptions();

	Real hydrop_lys_arg_setting;
	Real hydrop_adder;


#ifdef    SERIALIZATION
protected:
	friend class cereal::access;

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //sap
} //guidance_scoreterms
} //pack
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_guidance_scoreterms_sap_SapParameterOptions )
#endif // SERIALIZATION

#endif // INCLUDED_core_scoring_EtableEnergy_HH
