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

/// @file SurfaceOrientMover.hh
/// @brief <add a description of the class>
/// @author Robin A Thottungal (raugust1@jhu.edu)

#ifndef INCLUDED_protocols_surface_docking_SurfaceOrientMover_hh
#define INCLUDED_protocols_surface_docking_SurfaceOrientMover_hh

// Unit Headers

// Package headers

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/pose/datacache/CacheableDataType.hh>
// AUTO-REMOVED #include <basic/datacache/BasicDataCache.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

// C++ Headers
#include <string>
#include <map>
#include <list>
// AUTO-REMOVED #include <sstream>

#include <utility/vector1.hh>



namespace protocols {
namespace surface_docking {

class SurfaceOrientMover;
typedef utility::pointer::owning_ptr< SurfaceOrientMover > SurfaceOrientMoverOP;

class SurfaceOrientMover : public moves::Mover {

public:

    SurfaceOrientMover();

    //destructor
    ~SurfaceOrientMover();

    // virtual functions that get overloaded or
    //                           called from the inheriting classes
    void apply( core::pose::Pose & );

    virtual std::string get_name() const;

    core::Vector CalcTransVec(core::Real ProjectionDistance,
				core::Real VectorDistance,core::Vector Vec);

};


} // surfaceDockingProtocols
} // protocols

#endif
