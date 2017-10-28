// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/option.hh
/// @brief  Program options global and initialization function
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)

#ifndef INCLUDED_basic_options_option_hh
#define INCLUDED_basic_options_option_hh

#include <utility/options/OptionCollection.hh>                 // for Option...

#include <basic/options/keys/OptionKeys.hh> // Convenience

#include <utility/options/BooleanOption.hh>        // for BooleanOption
#include <utility/options/BooleanVectorOption.hh>  // for BooleanVectorOption
#include <utility/options/FileOption.hh>           // for FileOption
#include <utility/options/FileVectorOption.hh>     // for FileVectorOption
#include <utility/options/IntegerOption.hh>        // for IntegerOption
#include <utility/options/IntegerVectorOption.hh>  // for IntegerVectorOption
#include <utility/options/Option.hh>               // for Option
#include <utility/options/OptionCollection.hh>     // for OptionCollection
#include <utility/options/PathOption.hh>           // for PathOption
#include <utility/options/PathVectorOption.hh>     // for PathVectorOption
#include <utility/options/RealOption.hh>           // for RealOption
#include <utility/options/RealVectorOption.hh>     // for RealVectorOption
#include <utility/options/StringOption.hh>         // for StringOption
#include <utility/options/StringVectorOption.hh>   // for StringVectorOption

namespace basic {
namespace options {

// Types
typedef  utility::options::OptionCollection  OptionCollection;
typedef  utility::options::Option  Option;
typedef  utility::options::BooleanOption  BooleanOption;
typedef  utility::options::IntegerOption  IntegerOption;
typedef  utility::options::RealOption  RealOption;
typedef  utility::options::StringOption  StringOption;
typedef  utility::options::FileOption  FileOption;
typedef  utility::options::PathOption  PathOption;
typedef  utility::options::BooleanVectorOption  BooleanVectorOption;
typedef  utility::options::IntegerVectorOption  IntegerVectorOption;
typedef  utility::options::RealVectorOption  RealVectorOption;
typedef  utility::options::StringVectorOption  StringVectorOption;
typedef  utility::options::FileVectorOption  FileVectorOption;
typedef  utility::options::PathVectorOption  PathVectorOption;


/// @brief OptionCollection global
extern OptionCollection option;


/// @brief Named verbosity levels
extern int const silent  ; // No messages output
extern int const quiet   ;
extern int const standard;
extern int const inform  ;
extern int const chat    ;
extern int const yap     ;
extern int const gush    ;
extern int const verbose ; // All messages output


/// @brief Initialize the options
OptionCollection &
initialize();

/// @brief Process complex option inter-dependencies, prior to Tracer initializations
OptionCollection &
pre_tracer_process();

/// @brief Process complex option inter-dependencies (after tracer system is initialized)
OptionCollection &
process();


} // namespace options
} // namespace basic


#endif // INCLUDED_basic_options_option_HH
