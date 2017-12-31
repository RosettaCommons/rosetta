// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file util.cc
/// @brief util functions for general rna use
/// @details See also core/pose/rna/util.hh, core/chemical/rna/util.hh
/// @author Rhiju Das

#include <protocols/rna/util.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/rna/BasePair.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

using namespace core;
using namespace core::chemical::rna;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

static basic::Tracer TR( "protocols.rna.util" );

namespace protocols {
namespace rna {



} //rna
} //protocols
