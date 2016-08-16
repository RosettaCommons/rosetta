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


#include <protocols/antibody/grafting/scs_multi_template.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <basic/Tracer.hh>

namespace protocols {
namespace antibody {
namespace grafting {

//static THREAD_LOCAL basic::Tracer TR("protocols.antibody.grafting");

/// @brief Select CDR's template without filtering or sorting. In general you probably need to call select(...) instead
/// @throw _AE_scs_failed_ on failure
SCS_ResultsOP SCS_MultiTemplate::raw_select(AntibodySequence const &A)
{
	return source_->raw_select(A);
}

/// @brief Pad results vectors for each region (if possible) by adding arbitraty but compatible templates so at least n templates for each region is avalible
// void SCS_MultiTemplate::pad_results(uint n, AntibodySequence const &, SCS_Results &)
// {
// }


/// @brief Select CDR's template, filter it and sort. Try to provide at least 'n' templates if possible
/// @throw _AE_scs_failed_ on failure
SCS_ResultsOP SCS_MultiTemplate::select(uint n, AntibodySequence const &A)
{
	//TR << "SCS_MultiTemplate::select regions: " << TR.bgGreen << TR.Black; for(auto & s : regions_) TR << s; TR << std::endl;
	SCS_ResultsOP S = source_->select(n, A);

	SCS_ResultsOP R = std::make_shared<SCS_Results>();

	struct {
		string name;
		SCS_ResultVector SCS_Results::*region;
	} J[] {
		{"h1", &SCS_Results::h1}, {"h2", &SCS_Results::h2}, {"h3", &SCS_Results::h3},
		{"l1", &SCS_Results::l1}, {"l2", &SCS_Results::l2}, {"l3", &SCS_Results::l3},
		{"frh", &SCS_Results::frh}, {"frl", &SCS_Results::frl}, {"orientation", &SCS_Results::orientation},
	};

	for(auto &j : J) {
		for(uint i=0; i<n; ++i) {
			SCS_ResultVector &r = (*S).*j.region;

			SCS_ResultOP item;
			if( std::find(regions_.begin(), regions_.end(), j.name) != regions_.end() ) {
				if( r.size() > i ) item = r[i];
			}
			else if( r.size() > 1 ) item = r[0];

			(*R.*j.region).push_back(item);
		}
	}

	*this << "SCS_MultiTemplate( multi_template_regions = [";  for(auto & s : regions_) *this << s << ", ";  *this << "] ) resutls:\n";

	report(R, n);

	return R;
}



} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__
