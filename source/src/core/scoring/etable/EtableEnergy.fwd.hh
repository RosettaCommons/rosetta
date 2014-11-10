// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunction.hh
/// @brief  Score function class
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_etable_EtableEnergy_fwd_hh
#define INCLUDED_core_scoring_etable_EtableEnergy_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace etable {

class EtableEnergy;

class EtableEvaluator;
class TableLookupEvaluator;
class AnalyticEtableEvaluator;

typedef utility::pointer::shared_ptr< EtableEvaluator > EtableEvaluatorOP;
typedef utility::pointer::shared_ptr< EtableEvaluator const > EtableEvaluatorCOP;
typedef utility::pointer::weak_ptr< EtableEvaluator const > EtableEvaluatorCAP;

typedef utility::pointer::shared_ptr< AnalyticEtableEvaluator > AnalyticEtableEvaluatorOP;

typedef utility::pointer::shared_ptr< EtableEnergy > EtableEnergyOP;

} // etable
} // scoring
} // core


#endif
