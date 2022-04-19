// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     core/conformation/membrane/MembraneGeometry.fwd.hh
/// @brief    Base class for geometry of the membrane
///
/// @details  MembraneGeometry is a container object that describes the geometry of the membrane
///
/// @note     This object is a member of Conformation and should only be accessed using
///           pose.conformation().membrane_geometry().
///
/// @author   Hope Woods (hope.woods@vanderbilt.edu)

#ifndef INCLUDED_core_conformation_membrane_MembraneGeometry_fwd_hh
#define INCLUDED_core_conformation_membrane_MembraneGeometry_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace conformation {
namespace membrane {

//enum class MP_GEOMETRY_TRANSITION;
//enum class of possible membrane geometries
//while micelles, bicelles, and nanodisc are different options they all use the same transition function
//this is used in the MembraneInfo constructor to decide which geometry to create
enum class
	MP_GEOMETRY_TRANSITION {
	SLAB = 1,
	BICELLE = 2,
	VESICLE = 3,
	DOUBLE_VESICLE = 4
};


class MembraneGeometry;
typedef utility::pointer::shared_ptr< MembraneGeometry > MembraneGeometryOP;
typedef utility::pointer::shared_ptr< MembraneGeometry const > MembraneGeometryCOP;

} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_MembraneGeometry_fwd_hh

