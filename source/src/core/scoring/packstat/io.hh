// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packstat/io.hh
///
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_packstat_io_hh
#define INCLUDED_core_scoring_packstat_io_hh

#include <core/scoring/packstat/SimplePDB.fwd.hh>
#include <core/scoring/packstat/SimplePDB_Atom.fwd.hh>
#include <core/scoring/packstat/types.hh>

#include <iosfwd>

namespace core {
namespace scoring {
namespace packstat {

std::istream & operator>> ( std::istream & in , SimplePDB      & pdb  );
std::istream & operator>> ( std::istream & in , SimplePDB_Atom & atom );

std::ostream & operator<< ( std::ostream & out, SimplePDB      const & pdb    );
std::ostream & operator<< ( std::ostream & out, SimplePDB_Atom const & atom   );
std::ostream & operator<< ( std::ostream & out, Sphere         const & sphere );

} // namespace packstat
} // namespace scoring
} // namespace core

#endif
