// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/downstream/GeometrySecMatchRPE.fwd.hh
/// @brief
/// @author Florian Richter, floric@u.washington.edu, june 2009

#ifndef INCLUDED_protocols_match_downstream_GeometrySecMatchRPE_fwd_hh
#define INCLUDED_protocols_match_downstream_GeometrySecMatchRPE_fwd_hh

// Unit headers
// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace match {
namespace downstream {

class AtomGeometrySecMatchRPE;
typedef utility::pointer::shared_ptr< AtomGeometrySecMatchRPE > AtomGeometrySecMatchRPEOP;
typedef utility::pointer::shared_ptr< AtomGeometrySecMatchRPE const > AtomGeometrySecMatchRPECOP;

class AtomDistanceSecMatchRPE;
typedef utility::pointer::shared_ptr< AtomDistanceSecMatchRPE > AtomDistanceSecMatchRPEOP;
typedef utility::pointer::shared_ptr< AtomDistanceSecMatchRPE const > AtomDistanceSecMatchRPECOP;

class AtomAngleSecMatchRPE;
typedef utility::pointer::shared_ptr< AtomAngleSecMatchRPE > AtomAngleSecMatchRPEOP;
typedef utility::pointer::shared_ptr< AtomAngleSecMatchRPE const > AtomAngleSecMatchRPECOP;

class AtomDihedralSecMatchRPE;
typedef utility::pointer::shared_ptr< AtomDihedralSecMatchRPE > AtomDihedralSecMatchRPEOP;
typedef utility::pointer::shared_ptr< AtomDihedralSecMatchRPE const > AtomDihedralSecMatchRPECOP;

class GeometrySecMatchRPE;
typedef utility::pointer::shared_ptr< GeometrySecMatchRPE > GeometrySecMatchRPEOP;
typedef utility::pointer::shared_ptr< GeometrySecMatchRPE const > GeometrySecMatchRPECOP;


}
}
}

#endif
