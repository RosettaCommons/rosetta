// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/dunbrack/RotamerLibrary.py.hh
/// @brief   Explicit instantiation of *SingleResidueDunbrackLibrary classes for PyRosetta build
/// @author  Sergey Lyskov

#ifndef INCLUDED_core_pack_dunbrack_RotamerLibrary_py_hh
#define INCLUDED_core_pack_dunbrack_RotamerLibrary_py_hh

#ifdef PYROSETTA

#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.tmpl.hh>

namespace core {
namespace pack {
namespace dunbrack {

// We using both explicit instantiation and function call with by-reference argument so LLVM does not drop unused instantiation of template class.
// We also have to use constant directly due to template argument resolution issue

template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::ONE >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::TWO >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::THREE >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::FOUR >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::FIVE >;

template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::ONE >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::TWO >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::THREE >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::FOUR >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::FIVE >;

template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::THREE, core::pack::dunbrack::ONE >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::THREE, core::pack::dunbrack::TWO >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::THREE, core::pack::dunbrack::THREE >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::THREE, core::pack::dunbrack::FOUR >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::THREE, core::pack::dunbrack::FIVE >;

template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::FOUR, core::pack::dunbrack::ONE >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::FOUR, core::pack::dunbrack::TWO >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::FOUR, core::pack::dunbrack::THREE >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::FOUR, core::pack::dunbrack::FOUR >;
template class RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::FOUR, core::pack::dunbrack::FIVE >;

template class SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::ONE >;
template class SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::TWO >;
template class SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::THREE >;
template class SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::FOUR >;
template class SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::FIVE >;

template class SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::ONE >;
template class SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::TWO >;
template class SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::THREE >;
template class SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::FOUR >;
template class SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::FIVE >;


inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::ONE > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::TWO > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::THREE > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::FOUR > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::FIVE > & ) {};

inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::ONE > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::TWO > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::THREE > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::FOUR > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::FIVE > & ) {};

inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::THREE, core::pack::dunbrack::ONE > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::THREE, core::pack::dunbrack::TWO > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::THREE, core::pack::dunbrack::THREE > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::THREE, core::pack::dunbrack::FOUR > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::THREE, core::pack::dunbrack::FIVE > & ) {};

inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::FOUR, core::pack::dunbrack::ONE > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::FOUR, core::pack::dunbrack::TWO > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::FOUR, core::pack::dunbrack::THREE > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::FOUR, core::pack::dunbrack::FOUR > & ) {};
inline void __instantiate__(RotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::FOUR, core::pack::dunbrack::FIVE > & ) {};

inline void __instantiate__(SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::ONE > & ) {};
inline void __instantiate__(SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::TWO > & ) {};
inline void __instantiate__(SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::THREE > & ) {};
inline void __instantiate__(SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::FOUR > & ) {};
inline void __instantiate__(SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::ONE, core::pack::dunbrack::FIVE > & ) {};

inline void __instantiate__(SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::ONE > & ) {};
inline void __instantiate__(SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::TWO > & ) {};
inline void __instantiate__(SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::THREE > & ) {};
inline void __instantiate__(SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::FOUR > & ) {};
inline void __instantiate__(SemiRotamericSingleResidueDunbrackLibrary< core::pack::dunbrack::TWO, core::pack::dunbrack::FIVE > & ) {};

} // dunbrack
} // pack
} // core
#endif


#endif // INCLUDED_core_pack_dunbrack_RotamerLibrary_py_hh
