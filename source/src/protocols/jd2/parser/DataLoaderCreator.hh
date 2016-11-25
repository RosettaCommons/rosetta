// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/parser/DataLoaderCreator.hh
/// @brief  Creator class for the parser's DataLoader classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd2_parser_DataLoaderCreator_HH
#define INCLUDED_protocols_jd2_parser_DataLoaderCreator_HH

// Package headers
#include <protocols/jd2/parser/DataLoader.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

// Boost headers -- temp -- move to base class
#include <boost/function.hpp>

namespace protocols {
namespace jd2 {
namespace parser {

/// @brief A class for creating new instances of DataLoaders for use in the XML Parser
class DataLoaderCreator : public utility::pointer::ReferenceCount
{
public:
	typedef boost::function< std::string ( std::string const & ) > DerivedNameFunction;

public:
	DataLoaderCreator();
	virtual ~DataLoaderCreator();

	virtual DataLoaderOP create_loader() const = 0;
	virtual std::string keyname() const = 0;
	virtual DerivedNameFunction schema_ct_naming_function() const = 0;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const = 0;

};

typedef utility::pointer::shared_ptr< DataLoaderCreator > DataLoaderCreatorOP;
typedef utility::pointer::shared_ptr< DataLoaderCreator const > DataLoaderCreatorCOP;

} //namespace parser
} //namespace jd2
} //namespace protocols

#endif //INCLUDED_protocols_jd2_parser_DataLoaderCreator_HH
