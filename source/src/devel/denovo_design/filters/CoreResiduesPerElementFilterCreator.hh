// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file    src/devel/denovo_design/filters/CoreResiduesPerElementFilterCreator.fwd.hh
/// @brief   Creator for Tom's denovo protocol
/// @author  Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_devel_denovo_design_filters_CoreResiduesPerElementFilterCreator_hh
#define INCLUDED_devel_denovo_design_filters_CoreResiduesPerElementFilterCreator_hh

// Project headers
#include <protocols/filters/FilterCreator.hh>

namespace devel {
namespace denovo_design {
namespace filters {

class CoreResiduesPerElementFilterCreator : public protocols::filters::FilterCreator {

public:
	virtual protocols::filters::FilterOP create_filter() const;
	virtual std::string keyname() const;
	static  std::string filter_name();

};

}
}
}


#endif
