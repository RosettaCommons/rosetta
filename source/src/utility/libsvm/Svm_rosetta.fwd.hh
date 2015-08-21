// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/libsvm/Svm_rosetta.fwd.hh
/// @brief  owning pointer for libsvm
/// @author TJ Brunette

#ifndef INCLUDED_utility_libsvm_svm_fwd_hh
#define INCLUDED_utility_libsvm_svm_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace utility {
namespace libsvm {
class Svm_node_rosetta;
typedef utility::pointer::shared_ptr< Svm_node_rosetta > Svm_node_rosettaOP;

class Svm_rosetta;
typedef utility::pointer::shared_ptr< Svm_rosetta > Svm_rosettaOP;

} //namespace libsvm
} //namespace utility

#endif
