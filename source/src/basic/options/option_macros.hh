// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

#ifndef INCLUDED_basic_options_option_macros_hh
#define INCLUDED_basic_options_option_macros_hh

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/after_opts.hh>

#include <platform/types.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>


#define OPT(akey)							\
  basic::options::option.add_relevant( basic::options::OptionKeys::akey )

#define NEW_OPT(akey,help,adef)						\
  basic::options::option.add( basic::options::OptionKeys::akey , help ).def( adef ); \
  OPT( akey )

// option macros for Vector options with multiple default values...
#define NEW_OPT2(akey,help,adef,adef2)																				\
  basic::options::option.add( basic::options::OptionKeys::akey , help ).def( adef ).def( adef2 ); \
  OPT( akey )

#define NEW_OPT3(akey,help,adef,adef2,adef3)																	\
  basic::options::option.add( basic::options::OptionKeys::akey , help ).def( adef ).def( adef2 ).def(adef3); \
  OPT( akey )

#define NEW_OPT4(akey,help,adef,adef2,adef3,adef4)														\
  basic::options::option.add( basic::options::OptionKeys::akey , help ).def( adef ).def( adef2 ).def(adef3).def(adef4); \
  OPT( akey )

#define NEW_OPT5(akey,help,adef,adef2,adef3,adef4,adef5)											\
  basic::options::option.add( basic::options::OptionKeys::akey , help ).def( adef ).def( adef2 ).def(adef3).def(adef4).def(adef5); \
  OPT( akey )

#define NEW_OPT6(akey,help,adef,adef2,adef3,adef4,adef5,adef6)								\
  basic::options::option.add( basic::options::OptionKeys::akey , help ).def( adef ).def( adef2 ).def(adef3).def(adef4).def(adef5).def(adef6); \
  OPT( akey )

#define NEW_OPT7(akey,help,adef,adef2,adef3,adef4,adef5,adef6,adef7)					\
  basic::options::option.add( basic::options::OptionKeys::akey , help ).def( adef ).def( adef2 ).def(adef3).def(adef4).def(adef5).def(adef6).def(adef7); \
  OPT( akey )

#define NEW_OPT8(akey,help,adef,adef2,adef3,adef4,adef5,adef6,adef7,adef8)		\
  basic::options::option.add( basic::options::OptionKeys::akey , help ).def( adef ).def( adef2 ).def(adef3).def(adef4).def(adef5).def(adef6).def(adef7).def(adef8); \
  OPT( akey )

#define NEW_OPT9(akey,help,adef,adef2,adef3,adef4,adef5,adef6,adef7,adef8,adef9)	\
  basic::options::option.add( basic::options::OptionKeys::akey , help ).def( adef ).def( adef2 ).def(adef3).def(adef4).def(adef5).def(adef6).def(adef7).def(adef8).def(adef9); \
  OPT( akey )

#define NEW_OPT10(akey,help,adef,adef2,adef3,adef4,adef5,adef6,adef7,adef8,adef9,adef10) \
  basic::options::option.add( basic::options::OptionKeys::akey , help ).def( adef ).def( adef2 ).def(adef3).def(adef4).def(adef5).def(adef6).def(adef7).def(adef8).def(adef9).def(adef10); \
  OPT( akey )

#define OPT_KEY( type, key )					       \
  namespace basic { 	namespace options {	namespace OptionKeys {	\
	basic::options::type##OptionKey const key( #key );		\
      }}}

#define OPT_1GRP_KEY( type, grp, key )					\
  namespace basic { 	namespace options {	namespace OptionKeys { namespace grp { \
	  basic::options::type##OptionKey const key( #grp":"#key );	\
	}}}}


#define OPT_2GRP_KEY( type, grp1, grp2, key )				\
  namespace basic { 	namespace options {	namespace OptionKeys { namespace grp1 { namespace grp2 { \
	    basic::options::type##OptionKey const key( #grp1":"#grp2":"#key ); \
	  }}}}}

#define OPT_3GRP_KEY( type, grp1, grp2, grp3, key )				\
  namespace basic { 	namespace options {	namespace OptionKeys { namespace grp1 { namespace grp2 { namespace grp3 { \
	    basic::options::type##OptionKey const key( #grp1":"#grp2":"#grp3":"#key ); \
	  }}}}}}


#define  EXTERN_OPT_KEY( type, key )					       \
  namespace basic { 	namespace options {	namespace OptionKeys {	\
	extern basic::options::type##OptionKey const key;		\
      }}}

#define EXTERN_OPT_1GRP_KEY( type, grp, key )					\
  namespace basic { 	namespace options {	namespace OptionKeys { namespace grp { \
			extern basic::options::type##OptionKey const key;			\
	}}}}

#define EXTERN_OPT_2GRP_KEY( type, grp1, grp2, key )				\
  namespace basic { 	namespace options {	namespace OptionKeys { namespace grp1 { namespace grp2 { \
	    extern basic::options::type##OptionKey const key; \
	  }}}}}

#define EXTERN_OPT_3GRP_KEY( type, grp1, grp2, grp3, key )				\
  namespace basic { 	namespace options {	namespace OptionKeys { namespace grp1 { namespace grp2 { namespace grp3 { \
	    extern basic::options::type##OptionKey const key; \
	  }}}}}}


#endif
