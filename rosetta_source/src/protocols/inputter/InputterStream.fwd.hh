// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/inputter/InputterStream.fwd.hh
/// @brief InputterStream holds a list of of streams and does VERY basic things with them
/// like controls whether structures are duplicated across masters
/// and controls whether to take the list sequentially or round robin
/// @author Ken Jung

#ifndef INCLUDED_protocols_inputter_InputterStream_fwd_hh
#define INCLUDED_protocols_inputter_InputterStream_fwd_hh

#include <boost/shared_ptr.hpp>

namespace protocols {
namespace inputter {

class InputterStream;
typedef boost::shared_ptr< InputterStream > InputterStreamSP;

} // inputter
} // protocols

#endif //INCLUDED_protocols_inputter_InputterStream_fwd_HH
