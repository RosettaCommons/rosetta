// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief      Mover partners apart and relax them separately
/// @details Run quick relax on separated partners; this emulates unbound docking
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_docking_membrane_QuickRelaxPartnersSeparately_hh
#define INCLUDED_protocols_docking_membrane_QuickRelaxPartnersSeparately_hh

// Unit Headers
#include <protocols/docking/membrane/QuickRelaxPartnersSeparately.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/conformation/membrane/SpanningTopology.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/DockingPartners.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace docking {
namespace membrane {

class QuickRelaxPartnersSeparately : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details Gets the jump from docking partners
	QuickRelaxPartnersSeparately();

	/// @brief Copy Constructor
	QuickRelaxPartnersSeparately( QuickRelaxPartnersSeparately const & src );

	/// @brief Destructor
	~QuickRelaxPartnersSeparately() override;

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (QuickRelaxPartnersSeparately)
	std::string get_name() const override;

	/// @brief Moving partners apart and relax them separately
	void apply( Pose & pose ) override;

private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options with JD2
	void register_options();

	/// @brief Initialize Mover options from the commandline
	void init_from_cmd();

	/// @brief Finalize setup
	void finalize_setup( Pose & pose );

private: // data

	/// @brief Native pose
	Pose native_;

	// docking partners
	core::pose::DockingPartners partners_;

	// jump
	int jump_;
	utility::vector1< int > jumps_;

	// SpanningTopology objects
	core::conformation::membrane::SpanningTopologyOP topo_;   // full pose
	core::conformation::membrane::SpanningTopologyOP topo_up_; // upstream partner
	core::conformation::membrane::SpanningTopologyOP topo_down_;  // downstream partner

	// scorefunction
	core::scoring::ScoreFunctionOP sfxn_;

};

} // membrane
} // docking
} // protocols

#endif // INCLUDED_protocols_docking_membrane_QuickRelaxPartnersSeparately_hh
