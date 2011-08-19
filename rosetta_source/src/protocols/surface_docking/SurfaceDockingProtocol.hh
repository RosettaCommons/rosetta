// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;
//     rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
//     under license.
// (c) The Rosetta software is developed by the contributing members of the
//     Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about
//     this can be
// (c) addressed to University of Washington UW TechTransfer,
//     email: license@u.washington.edu.

/// @file SurfaceDocking.hh
/// @brief <add a description of the class>
/// @author Robin A Thottungal (raugust1@jhu.edu)

#ifndef INCLUDED_protocols_surface_docking_SurfaceDockingProtocol_hh
#define INCLUDED_protocols_surface_docking_SurfaceDockingProtocol_hh

// Unit Headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
// Package headers

#include <protocols/moves/MoverStatus.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/surface_docking/SurfaceParameters.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <protocols/moves/DataMap.fwd.hh>
// ObjexxFCL Headers

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>
#include <map>
#include <list>

//Auto Headers
#include <sstream>


namespace protocols {
namespace surface_docking {

class SurfaceDockingProtocol : public moves::Mover {

public:

	SurfaceDockingProtocol();

	//destructor
	~SurfaceDockingProtocol();

	void apply( core::pose::Pose & );

	virtual std::string get_name() const;

	void setupFoldTree(core::pose::Pose & pose);

	//virtual void setup_list( core::pose::Pose & ) = 0;

	//virtual void set_angles( core::Real ) = 0;

	//virtual bool make_move( core::pose::Pose & ) = 0;

private:


};


} // surfaceDockingProtocols
} // protocols

#endif
