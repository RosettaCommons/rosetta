// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/types.hh
/// @brief  core::scoring package type declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_core_scoring_types_hh
#define INCLUDED_core_scoring_types_hh


// Project headers
#include <core/types.hh>
//#include <core/conformation/types.hh>

// ObjexxFCL headers
#include <ObjexxFCL/CArray.fwd.hh>
#include <ObjexxFCL/CArrayP.fwd.hh>

#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray4D.fwd.hh>
#include <ObjexxFCL/FArray5D.fwd.hh>
#include <ObjexxFCL/KeyFArray1D.fwd.hh>
#include <ObjexxFCL/KeyFArray2D.fwd.hh>
#include <ObjexxFCL/KeyFArray3D.fwd.hh>


#ifdef WIN32
#include <ObjexxFCL/FArray2D.hh> // WIN32 INCLUDE
#endif


namespace core {
namespace scoring {


// Floating point scalars
typedef  core::Real  Probability;
typedef  Real  Weight;
typedef  Real  Score;
typedef  float  TableEnergy;
typedef  float  TableProbability;

// Floating point arrays
typedef  ObjexxFCL::CArray< Energy >  CArray_Energy;
typedef  ObjexxFCL::CArrayP< Energy >  CArrayP_Energy;
typedef  ObjexxFCL::CArray< TableEnergy >  CArray_TableEnergy;
typedef  ObjexxFCL::CArrayP< TableEnergy >  CArrayP_TableEnergy;
typedef  ObjexxFCL::FArray1D< Length >  FArray1D_Length;
typedef  ObjexxFCL::FArray2D< Length >  FArray2D_Length;
typedef  ObjexxFCL::FArray3D< Length >  FArray3D_Length;
typedef  ObjexxFCL::FArray4D< Length >  FArray4D_Length;
typedef  ObjexxFCL::FArray5D< Length >  FArray5D_Length;
typedef  ObjexxFCL::FArray1D< Weight >  FArray1D_Weight;
typedef  ObjexxFCL::FArray2D< Weight >  FArray2D_Weight;
typedef  ObjexxFCL::FArray3D< Weight >  FArray3D_Weight;
typedef  ObjexxFCL::FArray4D< Weight >  FArray4D_Weight;
typedef  ObjexxFCL::FArray5D< Weight >  FArray5D_Weight;
typedef  ObjexxFCL::FArray1D< Energy >  FArray1D_Energy;
typedef  ObjexxFCL::FArray2D< Energy >  FArray2D_Energy;
typedef  ObjexxFCL::FArray3D< Energy >  FArray3D_Energy;
typedef  ObjexxFCL::FArray4D< Energy >  FArray4D_Energy;
typedef  ObjexxFCL::FArray5D< Energy >  FArray5D_Energy;
typedef  ObjexxFCL::FArray1D< TableEnergy >  FArray1D_TableEnergy;
typedef  ObjexxFCL::FArray2D< TableEnergy >  FArray2D_TableEnergy;
typedef  ObjexxFCL::FArray3D< TableEnergy >  FArray3D_TableEnergy;
typedef  ObjexxFCL::FArray4D< TableEnergy >  FArray4D_TableEnergy;
typedef  ObjexxFCL::FArray5D< TableEnergy >  FArray5D_TableEnergy;
typedef  ObjexxFCL::FArray2D< CArrayP_TableEnergy >  AtomPairEnergyTable;
typedef  ObjexxFCL::FArray1D< Probability >  FArray1D_Probability;
typedef  ObjexxFCL::FArray2D< Probability >  FArray2D_Probability;
typedef  ObjexxFCL::FArray3D< Probability >  FArray3D_Probability;
typedef  ObjexxFCL::FArray4D< Probability >  FArray4D_Probability;
typedef  ObjexxFCL::FArray5D< Probability >  FArray5D_Probability;
typedef  ObjexxFCL::FArray1D< TableProbability >  FArray1D_TableProbability;
typedef  ObjexxFCL::FArray2D< TableProbability >  FArray2D_TableProbability;
typedef  ObjexxFCL::FArray3D< TableProbability >  FArray3D_TableProbability;
typedef  ObjexxFCL::FArray4D< TableProbability >  FArray4D_TableProbability;
typedef  ObjexxFCL::FArray5D< TableProbability >  FArray5D_TableProbability;
typedef  ObjexxFCL::KeyFArray1D< Real >  KeyFArray1D_Real;
typedef  ObjexxFCL::KeyFArray2D< Real >  KeyFArray2D_Real;
typedef  ObjexxFCL::KeyFArray3D< Real >  KeyFArray3D_Real;
typedef  ObjexxFCL::KeyFArray1D< Weight >  KeyFArray1D_Weight;
typedef  ObjexxFCL::KeyFArray2D< Weight >  KeyFArray2D_Weight;
typedef  ObjexxFCL::KeyFArray3D< Weight >  KeyFArray3D_Weight;
typedef  ObjexxFCL::KeyFArray1D< Energy >  KeyFArray1D_Energy;
typedef  ObjexxFCL::KeyFArray2D< Energy >  KeyFArray2D_Energy;
typedef  ObjexxFCL::KeyFArray3D< Energy >  KeyFArray3D_Energy;
typedef  ObjexxFCL::KeyFArray1D< Probability >  KeyFArray1D_Probability;
typedef  ObjexxFCL::KeyFArray2D< Probability >  KeyFArray2D_Probability;
typedef  ObjexxFCL::KeyFArray3D< Probability >  KeyFArray3D_Probability;


} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_types_HH
