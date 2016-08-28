// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/NucleotideTools.hh
/// @author Christoffer Norn (ch.norn@gmail.com)

//#include <core/types.hh>
#include <string>

namespace core {
namespace chemical {
namespace NucleotideTools {

    std::string codon2aa( std::string const & codon );
    std::string aa2randomCodon( char const & aa );
    //void add_nt_seq_to_pose( core::pose::Pose & pose ); // get the segment names for those segments that are constant in this splice function

} // NucleotideTools
} // chemical
} // core
