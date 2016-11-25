// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file    src/devel/denovo_design/filters/SSShapeComplementarityFilterCreator.fwd.hh
/// @brief   Creator for Tom's denovo protocol
/// @author  Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_devel_denovo_design_filters_SSShapeComplementarityFilterCreator_hh
#define INCLUDED_devel_denovo_design_filters_SSShapeComplementarityFilterCreator_hh

// Project headers
#include <protocols/filters/FilterCreator.hh>

namespace devel {
namespace denovo_design {
namespace filters {

class SSShapeComplementarityFilterCreator : public protocols::filters::FilterCreator {

public:
	// XRW TEMP  virtual protocols::filters::FilterOP create_filter() const;
	// XRW TEMP  virtual std::string keyname() const;
	// XRW TEMP  static  std::string filter_name();
	protocols::filters::FilterOP create_filter() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;

};

}
}
}


#endif
