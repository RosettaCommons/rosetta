// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/sewing/scoring/AssemblyScorer.hh
/// @brief an interface for making Requirements that deal with Assemblies
/// @author frankdt (frankdt@email.unc.edu)
/// @note   This is interface: it has no fields, and only
///         pure virtual methods.  No further constructors should
///         be defined.


#ifndef INCLUDED_protocols_sewing_scoring_AssemblyScorer_hh
#define INCLUDED_protocols_sewing_scoring_AssemblyScorer_hh


// Project forward headers
#include <protocols/sewing/scoring/AssemblyScorer.fwd.hh>
#include <core/types.hh>
#include <protocols/sewing/data_storage/SmartAssembly.fwd.hh>


// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
namespace protocols {
namespace sewing {
namespace scoring {


/// @brief an interface for making Requirements that deal with Assemblies
class AssemblyScorer  {


public: // Types


protected: // Types




public: // Constants


protected: // Constants




public: // Creation


	/// @brief Destructor
	virtual
	~AssemblyScorer() = default;


	// Creation


	/// @brief Prevent direct instantiation: No other constructors allowed.
	AssemblyScorer() = default;

	AssemblyScorer(AssemblyScorer const &) = default;


public: // Methods
	// Further subsections of methods allowed

	virtual core::Real
	score(data_storage::SmartAssemblyCOP assembly)=0;

	virtual std::string
	get_name() const=0;

	virtual void
	set_options_from_tag(
		utility::tag::TagCOP scorer_tag,
		basic::datacache::DataMap& datamap,
		protocols::filters::Filters_map const & filtermap,
		protocols::moves::Movers_map const & movermap,
		core::pose::Pose const & pose)=0;
	/*
	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & );
	*/
	virtual core::Real
	get_weight() const=0;

	virtual void
	set_weight( core::Real )=0;

	virtual core::Real
	get_last_score() const=0;

	virtual void
	set_last_score( core::Real )=0;

	virtual core::Real
	get_old_last_score() const=0;

	virtual void
	set_old_last_score( core::Real) = 0;
protected: // Methods
	// Further subsections of methods allowed




public: // Properties


protected: // Properties

	// NO FIELDS ALLOWED


}; // AssemblyScorer


} //protocols
} //sewing
} //scoring


#endif //INCLUDED_protocols_sewing_scoring_AssemblyScorer_hh
