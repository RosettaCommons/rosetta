// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/RotamerBoltzmannWeight2Creator.hh
/// @brief Next-generation RotamerBoltzmannWeight filter
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_simple_filters_RotamerBoltzmannWeight2Creator_hh
#define INCLUDED_protocols_simple_filters_RotamerBoltzmannWeight2Creator_hh

#include <protocols/filters/FilterCreator.hh>

namespace protocols {
namespace simple_filters {

class RotamerBoltzmannWeight2Creator : public protocols::filters::FilterCreator {
public:
	virtual protocols::filters::FilterOP create_filter() const;
	virtual std::string keyname() const;
	static std::string filter_name();
};

} //protocols
} //simple_filters

#endif //INCLUDED_protocols_simple_filters_RotamerBoltzmannWeight2_fwd_hh



