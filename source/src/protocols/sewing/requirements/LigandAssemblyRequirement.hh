// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/sewing/requirements/LigandAssemblyRequirement.hh
/// @brief an interface for making Ligand specific Requirements for Assemblies
/// @author Minnie Langlois (minnie@email.unc.edu)
/// @note   This is interface: it has no fields, and only
///         pure virtual methods.  No further constructors should
///         be defined.


#ifndef INCLUDED_protocols_sewing_requirements_LigandAssemblyRequirement_hh
#define INCLUDED_protocols_sewing_requirements_LigandAssemblyRequirement_hh


// Project forward headers
#include <protocols/sewing/requirements/LigandAssemblyRequirement.fwd.hh>
#include <protocols/sewing/requirements/AssemblyRequirement.hh>
#include <protocols/sewing/data_storage/SmartAssembly.fwd.hh>

// Core headers

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
namespace protocols {
namespace sewing {
namespace requirements {


/// @brief an interface for making Requirements that deal with Assemblies
class LigandAssemblyRequirement : public AssemblyRequirement {


public: // Types


protected: // Types




public: // Constants


protected: // Constants




public: // Creation


	/// @brief Destructor
	~LigandAssemblyRequirement() override =default;


protected: // Creation


	/// @brief Prevent direct instantiation: No other constructors allowed.
	LigandAssemblyRequirement()=default;
	LigandAssemblyRequirement( LigandAssemblyRequirement const & )=default;

public: // Methods
	// Further subsections of methods allowed
	std::pair<bool,bool>
	test(data_storage::SmartAssemblyOP assembly) override =0;

	void
	set_options_from_tag(
		utility::tag::TagCOP requirement_tag,
		basic::datacache::DataMap& datamap
	) override =0;

	std::string
	get_name() override =0;
protected: // Methods
	// Further subsections of methods allowed




public: // Properties


protected: // Properties


	// NO FIELDS ALLOWED


}; // LigandAssemblyRequirement


} //protocols
} //sewing
} //requirements


#endif //INCLUDED_protocols_sewing_requirements_LigandAssemblyRequirement_hh



