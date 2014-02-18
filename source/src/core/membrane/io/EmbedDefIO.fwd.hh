// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedDefIO.fwd.hh
///
/// @brief      Reads in membrane protein chain embeddings
/// @details    Protein embedding is defined by a normal and center vector positioned
///             with respect to the membrane as well as a specified depth
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedDefIO_fwd_hh
#define INCLUDED_core_membrane_io_EmbedDefIO_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace membrane {
namespace io {

/// @brief Embedding Definition IO Class
/// @details Read in membrane protien chain embedding data from files
class EmbedDefIO;
typedef utility::pointer::owning_ptr< EmbedDefIO > EmbedDefIOOP;
typedef utility::pointer::owning_ptr< EmbedDefIO const > EmbedDefIOCOP;

} // io
} // membrane
} // core


#endif // INCLUDED_core_membrane_io_EmbedDefIO_fwd_hh

