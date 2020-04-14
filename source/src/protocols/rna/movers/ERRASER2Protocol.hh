// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/movers/ERRASER2Protocol.hh
/// @brief Run a single-threaded, checkpoint free, RosettaScripts accessible ERRASER2 job
/// @author Andy Watkins (andy.watkins2@gmail.com)

#ifndef INCLUDED_protocols_rna_movers_ERRASER2Protocol_HH
#define INCLUDED_protocols_rna_movers_ERRASER2Protocol_HH

// Unit headers
#include <protocols/rna/movers/ERRASER2Protocol.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/stepwise/monte_carlo/mover/StepWiseMasterMover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace rna {
namespace movers {

///@brief Run a single-threaded, checkpoint free, RosettaScripts accessible ERRASER2 job
class ERRASER2Protocol : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	ERRASER2Protocol();

	/// @brief Copy constructor (not needed unless you need deep copies)
	ERRASER2Protocol( ERRASER2Protocol const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~ERRASER2Protocol() override;


public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////


	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;


public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	//ERRASER2Protocol & operator=( ERRASER2Protocol const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	void score_function( core::scoring::ScoreFunctionOP const & scorefxn_in ) { scorefxn_ = scorefxn_in; };
	core::scoring::ScoreFunctionOP score_function() const { return scorefxn_; };

	void n_rounds( core::Size const n_rounds ) { n_rounds_ = n_rounds; };
	core::Size n_rounds() const { return n_rounds_; };

	void rebuild_residue_selector( core::select::residue_selector::ResidueSelectorCOP const & residue_selector ) { rebuild_residue_selector_ = residue_selector; };
	core::select::residue_selector::ResidueSelectorCOP rebuild_residue_selector() const { return rebuild_residue_selector_; };

private: // methods

	protocols::stepwise::monte_carlo::mover::StepWiseMasterMover
	configure_master_mover( Pose const & start_pose );

	void
	resample_full_model( core::Size const resample_round, core::pose::Pose & start_pose, utility::vector1< core::Size > const & definite_residues, core::Size const nstruct );

private: // data

	core::scoring::ScoreFunctionOP scorefxn_ = nullptr;
	bool minimize_protein_ = false;
	core::Size n_rounds_ = 3;
	core::select::residue_selector::ResidueSelectorCOP rebuild_residue_selector_ = nullptr;
};

std::ostream &
operator<<( std::ostream & os, ERRASER2Protocol const & mover );

} //movers
} //rna
} //protocols

#endif //protocols_rna_movers_ERRASER2Protocol_HH
