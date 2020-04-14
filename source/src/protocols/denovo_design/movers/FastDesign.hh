// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/movers/FastDesign.hh
/// @brief The FastDesign Protocol
/// @details
/// @author Tom Linsky
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- Added support for D-amino acids.


#ifndef INCLUDED_protocols_denovo_design_movers_FastDesign_hh
#define INCLUDED_protocols_denovo_design_movers_FastDesign_hh

// Unit headers
#include <protocols/denovo_design/movers/FastDesign.fwd.hh>
#include <protocols/relax/FastRelax.hh>

// Protocol headers
#include <protocols/constraint_generator/ConstraintGenerator.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/task_operations/DesignAroundOperation.fwd.hh>
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/io/silent/silent.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Basic/Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace protocols {
namespace denovo_design {
namespace movers {

struct RelaxScriptCommand {
	RelaxScriptCommand():
		command( "" ),
		param1( 0.0 ),
		param2( 0.0 ),
		param3( 0.0 ),
		param4( 0.0 ),
		nparams( 0 )
	{}

	std::string command;
	core::Real  param1;
	core::Real  param2;
	core::Real  param3;
	core::Real  param4;
	core::Size  nparams;
};

class FastDesign : public protocols::relax::FastRelax {
public:
	/// @brief default constructor
	FastDesign();

	/// @brief Constructor with some options
	///
	FastDesign( core::scoring::ScoreFunctionOP scorefxn_in, core::Size standard_repeats = 0 );

	FastDesign( core::scoring::ScoreFunctionOP scorefxn_in, std::string const & script_file );

	/// @brief virtual constructor to allow derivation
	~FastDesign() override;


	/// @brief Create the default task factory.  Must be called before design can occur.
	void set_up_default_task_factory();

	/// @brief Parses the FastDesignTags
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	/// @brief Return the name of this mover.

	/// @brief return a fresh instance of this class in an owning pointer
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief return a copy of this class in an owning pointer
	protocols::moves::MoverOP
	clone() const override;

	/// @brief Apply the FastDesign. Overloaded apply function from mover base class.
	void
	apply( core::pose::Pose & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public: // CitationManager fxns:

	/// @brief Does this mover provide information about how to cite it?
	/// @details Returns true.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	bool mover_provides_citation_info() const override;

	/// @brief Provide the citation.
	/// @returns A vector of citation collections.  This allows the mover to provide citations for
	/// itself and for any modules that it invokes.
	/// @details This function provides citations for the task operations in the task factory and for FastDesign.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
	utility::vector1< basic::citation_manager::CitationCollectionCOP > provide_citation_info() const override;


protected:
	/// @brief sets constraint weights -- used with constraint ramping
	void
	set_constraint_weight(
		core::scoring::ScoreFunctionOP local_scorefxn,
		core::scoring::EnergyMap const & full_weights,
		core::Real const weight,
		core::pose::Pose & pose ) const override;

private:
	core::pack::task::TaskFactoryOP
	create_default_task_factory() const;

private:   // options
	/// @brief Should we clear designable residues (sets them to ALA so that they can be rebuilt)
	bool clear_designable_residues_;

private:   // other data
	core::Size run_count_;
	protocols::constraint_generator::ConstraintGeneratorCOPs cgs_;
};

} // movers
} // denovo_design
} // protocols

#endif
