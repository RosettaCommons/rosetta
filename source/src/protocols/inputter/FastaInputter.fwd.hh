// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/inputter/FastaInputter.fwd.hh
/// @brief An inputter that takes a list of fasta
/// @author Ken Jung

#ifndef INCLUDED_protocols_inputter_FastaInputter_fwd_hh
#define INCLUDED_protocols_inputter_FastaInputter_fwd_hh

#include <boost/shared_ptr.hpp>

namespace protocols {
namespace inputter {

class FastaInputter;
typedef boost::shared_ptr< FastaInputter > FastaInputterSP;

} // inputter
} // protocols

#endif //INCLUDED_protocols_inputter_Inputter_fwd_HH
