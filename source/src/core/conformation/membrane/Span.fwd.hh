// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/conformation/membrane/Span.fwd.hh
///
/// @brief  Object for describing start and end of a transmembrane span
/// @details The Span object stores 2 SSizes - a stard and end position of a transmembrane span.
///    Should be kept in a vector of spans toward describing the total spanning topology of a
///    membrane protein.
///    Last Modified: 7/23/14
///
/// @author  Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_Span_fwd_hh
#define INCLUDED_core_conformation_membrane_Span_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace conformation {
namespace membrane {

class Span;
typedef utility::pointer::shared_ptr< Span > SpanOP;
typedef utility::pointer::shared_ptr< Span const > SpanCOP;

} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_Span_fwd_hh
