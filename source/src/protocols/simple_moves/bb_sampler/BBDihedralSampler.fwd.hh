// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/BBDihedralSampler.fwd.hh
/// @brief This class functions to hold, access, and set independent and dependent dihedral data.
///   It can act as a base class for particular types of data.
///   It should eventually be moved out of here.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)  and Jason W. Labonte (JWLabonte@jhu.edu)


#ifndef INCLUDED_protocols_simple_moves_bb_sampler_BBDihedralSampler_fwd_hh
#define INCLUDED_protocols_simple_moves_bb_sampler_BBDihedralSampler_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace simple_moves {
namespace bb_sampler {

class BBDihedralSamplerBase;

typedef utility::pointer::shared_ptr< BBDihedralSamplerBase > BBDihedralSamplerBaseOP;
typedef utility::pointer::shared_ptr< BBDihedralSamplerBase const > BBDihedralSamplerBaseCOP;

class BBDihedralSampler;

typedef utility::pointer::shared_ptr< BBDihedralSampler > BBDihedralSamplerOP;
typedef utility::pointer::shared_ptr< BBDihedralSampler const > BBDihedralSamplerCOP;

class BBDihedralSampler2D;

typedef utility::pointer::shared_ptr< BBDihedralSampler2D > BBDihedralSampler2DOP;
typedef utility::pointer::shared_ptr< BBDihedralSampler2D const > BBDihedralSampler2DCOP;

class BBDihedralSampler3D;

typedef utility::pointer::shared_ptr< BBDihedralSampler3D > BBDihedralSampler3DOP;
typedef utility::pointer::shared_ptr< BBDihedralSampler3D const > BBDihedralSampler3DCOP;

class BBDihedralSamplerND;

typedef utility::pointer::shared_ptr< BBDihedralSamplerND > BBDihedralSamplerNDOP;
typedef utility::pointer::shared_ptr< BBDihedralSamplerND const > BBDihedralSamplerNDCOP;

} //core
} //pose
} //carbohydrates

#endif //INCLUDED_core_pose_carbohydrates_BBDihedralSampler_fwd_hh





