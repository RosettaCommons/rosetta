// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pockets/Fingerprint.hh
/// @brief  protocols::pockets::Fingerprint header
/// @author Ragul Gowthaman

#ifndef INCLUDED_protocols_pockets_DarcElectrostatics_hh
#define INCLUDED_protocols_pockets_DarcElectrostatics_hh

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/pockets/DarcElectrostatics.fwd.hh>
// AUTO-REMOVED #include <protocols/pockets/FingerprintMultifunc.fwd.hh>
#include <protocols/pockets/PocketGrid.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <numeric/constants.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1_bool.hh>
#include <list>
#include <cmath>

#include <utility/vector1.hh>

namespace protocols {
namespace pockets {

class DarcElectrostaticsBase : public utility::pointer::ReferenceCount {

	friend class FingerprintMultifunc;

public:

  DarcElectrostaticsBase();

};

class DelphiElectrostatics : public DarcElectrostaticsBase {

	friend class FingerprintMultifunc;

public:

	DelphiElectrostatics() {};

	core::Real grid_spacing_;

	void setup_from_DelphiGrid( std::string const & input_filename, Size const & esp_grid_size, core::Real const & esp_grid_spacing, core::Real const & esp_grid_midpoint_x, core::Real const & esp_grid_midpoint_y, core::Real const & esp_grid_midpoint_z);

	core::Real get_electrostatics_energy(core::pose::Pose const & ligand_pose);

	std::list< utility::vector1<core::Real> > esp_grid_point_list_;

};


}//pockets
}//protocols


#endif
