// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/BodyPosition.fwd.hh
/// @brief  numeric::BodyPosition forward declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_BodyPosition_fwd_hh
#define INCLUDED_numeric_BodyPosition_fwd_hh


namespace numeric {


// Forward
template< typename > class BodyPosition;


// Types
typedef  BodyPosition< float >        BodyPosition_float;
typedef  BodyPosition< double >       BodyPosition_double;
typedef  BodyPosition< long double >  BodyPosition_longdouble;


} // namespace numeric


#endif // INCLUDED_numeric_BodyPosition_FWD_HH

