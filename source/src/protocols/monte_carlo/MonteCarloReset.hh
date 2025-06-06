// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/monte_carlo/MonteCarloReset.hh
/// @brief perform a given mover and sample structures by MonteCarlo
/// @details The "score" evaluation of pose during MC after applying mover is done by
/// ither FilterOP that can do report_sm() or ScoreFunctionOP you gave.
/// By setting sample_type_ to high, you can also sample the pose that have higher score.
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

#ifndef INCLUDED_protocols_simple_moves_MonteCarloReset_hh
#define INCLUDED_protocols_simple_moves_MonteCarloReset_hh

// Unit Headers
#include <protocols/monte_carlo/MonteCarloReset.fwd.hh>
#include <protocols/monte_carlo/GenericMonteCarloMover.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility headers

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>



namespace protocols {
namespace monte_carlo {

class MonteCarloReset : public protocols::moves::Mover {
public:

	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;

	typedef utility::tag::TagCOP TagCOP;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::MoverOP MoverOP;


public: // constructor/destructor

	/// @brief default constructor
	MonteCarloReset();

	/// @brief destructor
	~MonteCarloReset() override;

	/// @brief create copy constructor
	MoverOP clone() const override;

	/// @brief create this type of objectt
	MoverOP fresh_instance() const override;

	/// @brief apply MonteCarloReset (Mover)
	void apply( Pose & pose ) override;

	/// @brief set mover
	void set_MC( GenericMonteCarloMoverOP mover );
	GenericMonteCarloMoverOP get_MC() const;

	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	// bool recover_low() const;
	// void recover_low( bool const recover );

private: // data
	// bool recover_low_; //dflt true; if false, recovers last

	/// @brief mover
	GenericMonteCarloMoverOP MC_mover_;
};

} // namespace monte_carlo
} // namespace protocols

#endif
