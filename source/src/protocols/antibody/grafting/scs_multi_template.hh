// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/antibody/grafting/scs_multi_template.hh
/// @brief Structural Component Selector (SCS) implementation of Multi-template selection
/// @author Sergey Lyskov


#ifndef INCLUDED_protocols_antibody_grafting_scs_multi_template_hh
#define INCLUDED_protocols_antibody_grafting_scs_multi_template_hh

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/scs_blast.hh>

namespace protocols {
namespace antibody {
namespace grafting {


/// Proxy-like class for SCS_Base that alter SCS results:
///   - for results not in 'regions_' return top result obtained from source_
///   - for results in 'regions_' return N best results
class SCS_MultiTemplate : public SCS_Base
{
	typedef utility::vector0<string> Regions;

public:
	explicit SCS_MultiTemplate(SCS_BaseOP const &source, Regions const & regions=Regions(), basic::ReportOP report=basic::ReportOP()) : SCS_Base(report), source_(source), regions_(regions) {}

	/// @brief Select CDR's template without filtering or sorting. In general you probably need to call select(...) instead
	/// @throw _AE_scs_failed_ on failure
	SCS_ResultsOP raw_select(AntibodySequence const &) override;

	/// @brief Select CDR's template, filter it and sort. Try to provide at least 'n' templates if possible
	/// @throw _AE_scs_failed_ on failure
	SCS_ResultsOP select(uint n, AntibodySequence const &) override;

	/// @brief Pad results vectors for each region (if possible) by adding arbitraty but compatible templates so at least n templates for each region is avalible
	void pad_results(uint, AntibodySequence const &, SCS_Results &) override {}

	/// @brief specify mutable regions for multi-template selection mode.
	void set_multi_template_regions(Regions const & regions) { regions_=regions; }

private:
	SCS_BaseOP source_;
	Regions regions_;
};


} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__

#endif // INCLUDED_protocols_antibody_grafting_scs_multi_template_hh
