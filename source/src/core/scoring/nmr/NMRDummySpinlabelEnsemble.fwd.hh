// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRDummySpinlabelEnsemble.fwd.hh
/// @brief   forward declaration for NMRDummySpinlabelEnsemble.hh
/// @details last Modified: 11/28/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_NMRDummySpinlabelEnsemble_FWD_HH
#define INCLUDED_core_scoring_nmr_NMRDummySpinlabelEnsemble_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace nmr {

class NMRDummySpinlabelConformer;
class NMRDummySpinlabelEnsemble;

typedef utility::pointer::shared_ptr< NMRDummySpinlabelConformer > NMRDummySpinlabelConformerOP;
typedef utility::pointer::shared_ptr< NMRDummySpinlabelConformer const > NMRDummySpinlabelConformerCOP;
typedef utility::pointer::weak_ptr< NMRDummySpinlabelConformer > NMRDummySpinlabelConformerAP;
typedef utility::pointer::weak_ptr< NMRDummySpinlabelConformer const > NMRDummySpinlabelConformerCAP;

typedef utility::pointer::shared_ptr< NMRDummySpinlabelEnsemble > NMRDummySpinlabelEnsembleOP;
typedef utility::pointer::shared_ptr< NMRDummySpinlabelEnsemble const > NMRDummySpinlabelEnsembleCOP;
typedef utility::pointer::weak_ptr< NMRDummySpinlabelEnsemble > NMRDummySpinlabelEnsembleAP;
typedef utility::pointer::weak_ptr< NMRDummySpinlabelEnsemble const > NMRDummySpinlabelEnsembleCAP;


} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_NMRDummySpinlabelEnsemble_FWD_HH
