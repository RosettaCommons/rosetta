// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedDefOptions.fwd.hh
///
/// @brief      Options for membrane protein chain embeddings
/// @details    Protein embedding is defined by a normal and center vector positioned
///             with respect to the membrane as well as a specified depth
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedDefOptions_fwd_hh
#define INCLUDED_core_membrane_io_EmbedDefOptions_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace membrane {
namespace io {

/// @brief Embedding Definition Options Class
/// @details Embedding Definition Options for protein chains
class EmbedDefOptions;
typedef utility::pointer::owning_ptr< EmbedDefOptions > EmbedDefOptionsOP;
typedef utility::pointer::owning_ptr< EmbedDefOptions const > EmbedDefOptionsCOP;

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedDefOptions_fwd_hh

