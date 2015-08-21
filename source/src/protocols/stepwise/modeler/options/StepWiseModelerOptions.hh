// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/options/StepWiseModelerOptions.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_StepWiseModelerOptions_HH
#define INCLUDED_protocols_stepwise_modeler_StepWiseModelerOptions_HH

#include <protocols/stepwise/modeler/options/StepWiseBasicModelerOptions.hh>
#include <protocols/stepwise/modeler/options/StepWiseProteinModelerOptions.hh>
#include <protocols/stepwise/modeler/options/StepWiseRNA_ModelerOptions.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.fwd.hh>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace options {

/////////////////////////////////////////////////////////////////////////////////////////////
// Note the use of multiple inheritance. This makes the code more concise and easy to expand,
//  in the sense that any new options just need to be defined once. But its a little
//  weird, since its hard to conceptualize that StepWiseModelerOptions 'is a'  StepWiseProteinModelerOptions
//   and 'is a' StepWiseRNA_ModelerOptions object. Anyway.
//
// The other way to do this would be by composition, but that requires new set() and get() functions
// to be defined in multiple places (in the base classes and in this class as well.)
//
//   -- rhiju, 2014.
/////////////////////////////////////////////////////////////////////////////////////////////

class StepWiseModelerOptions: public StepWiseBasicModelerOptions,
	public options::StepWiseProteinModelerOptions,
	public options::StepWiseRNA_ModelerOptions {

public:

	//constructor
	StepWiseModelerOptions();

	StepWiseModelerOptions( StepWiseModelerOptions const & src );

	//destructor
	~StepWiseModelerOptions();

public:

	StepWiseModelerOptionsOP clone() const;

	/// @brief Describe this instance to a given output stream
	virtual
	void
	show( std::ostream & ) const{}

	/// @brief Initialize from the recursive "tag" structure.
	virtual
	void
	parse_my_tag( utility::tag::TagCOP ){}

	/// @brief The class name (its type) for a particular ResourceOptions instance.
	/// This function allows for better error message delivery.
	virtual
	std::string
	type() const{ return "StepWiseModelerOptions";}

public:

	void
	initialize_from_command_line();

	void
	setup_options_for_VDW_bin_checker( rna::checker::RNA_VDW_BinCheckerOP user_input_VDW_bin_checker ) const;

	StepWiseModelerOptionsOP
	get_sampler_options() const;

private:

	void
	initialize_variables();

private:


};

} //options
} //modeler
} //stepwise
} //protocols

#endif
