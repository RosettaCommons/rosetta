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
#include <protocols/moves/Mover.fwd.hh>

// Core headers
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Basic/Utility headers

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

	/// @brief Adds a constraint generator. Currently only used if ramp_down_constraints is set. If constraint generators are present, each time the constraint weights are changed, the constraint generators will be used to selectively remove constraint-generated constraints and add new ones
	void
	add_constraint_generator( protocols::constraint_generator::ConstraintGeneratorCOP generator );

	/// @brief Clears the list of constraint generators
	void
	clear_constraint_generators();

	/// @brief If set, residues that are designable will be mutated to alanine before design
	void
	set_clear_designable_residues( bool const clear_res ) { clear_designable_residues_ = clear_res; }

public: // CitationManager fxns:

	/// @brief Provide the citation.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & ) const override;

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
