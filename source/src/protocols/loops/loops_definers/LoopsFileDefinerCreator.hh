// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loops_definers/LoopsFileDefinerCreator.hh
/// @brief  Header for LoopsFileDefinerCreator for the LoosFileDefiner load-time factory registration scheme
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_loops_loops_definers_LoopsFileDefinerCreator_HH
#define INCLUDED_protocols_loops_loops_definers_LoopsFileDefinerCreator_HH

// Unit Headers
#include <protocols/loops/loops_definers/LoopsDefinerCreator.hh>
#include <protocols/loops/loops_definers/LoopsFileDefinerCreator.fwd.hh>

#include <core/types.hh>
#include <utility/vector1.hh>

#include <string>

namespace protocols {
namespace loops {
namespace loops_definers {

/// @brief creator for the LoopsFileDefiner class
class LoopsFileDefinerCreator : public LoopsDefinerCreator {
public:
	LoopsFileDefinerCreator();
	virtual ~LoopsFileDefinerCreator();

	virtual LoopsDefinerOP create_loops_definer() const;
	virtual std::string type_name() const;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;
};

} //namespace
} //namespace
} //namespace

#endif // include guard
