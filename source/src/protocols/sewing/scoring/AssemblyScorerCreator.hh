// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/sewing/AssemblyScorerCreator.hh
/// @brief  Base class for AssemblyScorerCreators for the AssemblyScorer load-time factory registration scheme
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_sewing_scoring_AssemblyScorerCreator_hh
#define INCLUDED_protocols_sewing_scoring_AssemblyScorerCreator_hh

// Unit Headers
#include <protocols/sewing/scoring/AssemblyScorerCreator.fwd.hh>

// Package Headers
#include <protocols/sewing/scoring/AssemblyScorer.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
// c++ headers
#include <string>

namespace protocols {
namespace sewing  {
namespace scoring {

/// @brief The Creator class is responsible for creating a particular
/// GlobalRequirement class.
class AssemblyScorerCreator : public utility::pointer::ReferenceCount
{
public:
	AssemblyScorerCreator() {}
	virtual ~AssemblyScorerCreator() {}
	virtual std::string keyname() const = 0;
	virtual AssemblyScorerOP create_assembly_scorer() const = 0;
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const =0;
	//static std::string assembly_scorer_ct_namer( std::string tag_name );
};


} //namespace
} //namespace
} //namespace

#endif
