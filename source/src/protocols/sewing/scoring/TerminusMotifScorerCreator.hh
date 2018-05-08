// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/sewing/TerminusMotifScorerCreator.hh
/// @brief  Base class for AssemblyScorerCreators for the AssemblyScorer load-time factory registration scheme
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_sewing_scoring_TerminusMotifScorerCreator_hh
#define INCLUDED_protocols_sewing_scoring_TerminusMotifScorerCreator_hh

// Unit Headers
#include <protocols/sewing/scoring/AssemblyScorerCreator.hh>
#include <protocols/sewing/scoring/TerminusMotifScorerCreator.fwd.hh>


namespace protocols {
namespace sewing  {
namespace scoring {

/// @brief The Creator class is responsible for creating a particular
/// GlobalRequirement class.
class TerminusMotifScorerCreator : public AssemblyScorerCreator
{
public:
	virtual AssemblyScorerOP create_assembly_scorer() const;
	virtual std::string keyname() const;

	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};


} //namespace
} //namespace
} //namespace

#endif
