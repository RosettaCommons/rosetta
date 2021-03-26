// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/tensorflow_manager/util.hh
/// @brief Utility functions for the tensorflow manager.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_basic_tensorflow_manager_util_hh
#define INCLUDED_basic_tensorflow_manager_util_hh

#include <string>

namespace basic {
namespace tensorflow_manager {

/// @brief Writes the following message:
/// The <module_name> requires compilation with Tensorflow support.  To compile with Tensorflow
/// support...
/// (Followed by instructions terminating in a carriage return character.)
/// @note If use_single_quotations is true, quoted strings use a single quotation (') instead of
/// a double (").  This is necessary in the context of, for example, the RosettaScripts XSD, where
/// double quotations are prohibited in descriptions.
std::string get_tensorflow_compilation_instructions( std::string const & module_name, bool const use_single_quotations = false );

} //tensorflow_manager
} //basic


#endif //basic/tensorflow_manager_util_hh

