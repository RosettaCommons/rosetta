// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  A filter that, for each dimer in a pose, outputs the peptide which contributes most to the interface.

/// @author Nir London
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @date   Jul. 30, 2014

#ifndef INCLUDED_protocols_analysis_PeptideDeriverFilter_fwd_hh
#define INCLUDED_protocols_analysis_PeptideDeriverFilter_fwd_hh

// Project headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace peptide_deriver {

class PeptideDeriverOutputter;
class PeptideDeriverStreamOutputter;
class PeptideDeriverOutputterContainer;
class PeptideDeriverPoseOutputter;

class PeptideDeriverFilter;

class CyclizedPeptideInfo;
class DerivedPeptideEntry;

typedef utility::pointer::shared_ptr<PeptideDeriverFilter> PeptideDeriverFilterOP;
typedef utility::pointer::shared_ptr<PeptideDeriverFilter const> PeptideDeriverFilterCOP;

typedef utility::pointer::shared_ptr<PeptideDeriverOutputter> PeptideDeriverOutputterOP;
typedef utility::pointer::shared_ptr<PeptideDeriverOutputter const> PeptideDeriverOutputterCOP;

typedef utility::pointer::shared_ptr<CyclizedPeptideInfo> CyclizedPeptideInfoOP;
typedef utility::pointer::shared_ptr<CyclizedPeptideInfo const> CyclizedPeptideInfoCOP;

typedef utility::pointer::shared_ptr<DerivedPeptideEntry> DerivedPeptideEntryOP;
typedef utility::pointer::shared_ptr<DerivedPeptideEntry const> DerivedPeptideEntryCOP;
} // namespace peptide_deriver
} // namespace protocols

#endif
