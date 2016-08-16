// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/parser/StandardDataLoaderCreator.hh
/// @brief  Creator classes for the default DataLoader classes, TaskOperationLoader and ScoreFunctionLoader
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd2_parser_StandardLoaderCreators_hh
#define INCLUDED_protocols_jd2_parser_StandardLoaderCreators_hh

// Package headers
#include <protocols/jd2/parser/DataLoaderCreator.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace jd2 {
namespace parser {

class ScoreFunctionLoaderCreator : public DataLoaderCreator
{
public:
	virtual DataLoaderOP create_loader() const;
	virtual std::string keyname() const;
	virtual DerivedNameFunction schema_ct_naming_function() const;
};

class TaskOperationLoaderCreator : public DataLoaderCreator
{
public:
	virtual DataLoaderOP create_loader() const;
	virtual std::string keyname() const;
	virtual DerivedNameFunction schema_ct_naming_function() const;
};

class FragSetLoaderCreator : public DataLoaderCreator
{
public:
	virtual DataLoaderOP create_loader() const;
	virtual std::string keyname() const;
};

class MonteCarloLoaderCreator : public DataLoaderCreator
{
public:
	virtual DataLoaderOP create_loader() const;
	virtual std::string keyname() const;
};


} //namespace parser
} //namespace jd2
} //namespace protocols

#endif
