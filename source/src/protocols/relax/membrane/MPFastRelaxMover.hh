// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/relax/membrane/MPFastRelaxMover.hh
///
/// @brief      Membrane Fast Relax Protocol - Relax with minimization of mem position
/// @details Apply the standard fast relax protocol. Enable minimization of the memrbane
///             jump and relax from the center of mass. Also use the smoothed
///             full atom membrane energy function.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 12/2/14

#ifndef INCLUDED_protocols_membrane_relax_MPFastRelaxMover_hh
#define INCLUDED_protocols_membrane_relax_MPFastRelaxMover_hh

// Unit Headers
#include <protocols/relax/membrane/MPFastRelaxMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/scoring/ScoreFunction.hh>

// Project Headers
#include <protocols/relax/RelaxProtocolBase.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace relax {
namespace membrane {

class MPFastRelaxMover : public protocols::moves::Mover {

public:

	////////////////////////////
	/// Constructors & Setup ///
	////////////////////////////

	/// @brief Default membrane fast relax constructor
	/// @details Do normal fast relax protocol with the
	/// membrane energy function and custom foldtree
	MPFastRelaxMover();

	/// @brief Allow the user to set a custom sfxn from
	/// the PyRosetta Interface
	MPFastRelaxMover( core::scoring::ScoreFunctionOP sfxn );

	/// @brief Destructor
	virtual ~MPFastRelaxMover();

	/// @brief Create a custom foldtree anchored at the COM
	/// @details Generate a foldtree where the membrane residue
	/// is anchored at the center of mass of the chain. Also recreates
	/// other jumps in the protein
	void
	setup_relax_foldtree( core::pose::Pose & pose );

	/// @brief Show protocol settings
	/// @details Show the current setup of this fast relax mover
	void show_protocol( core::pose::Pose & pose );

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	virtual protocols::moves::MoverOP clone() const;

	/// @brief Create a Fresh Instance of this Mover
	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	/////////////////////
	/// Mover methods ///
	/////////////////////

	/// @brief Apply fast relax - do the actual protocol
	virtual void apply( core::pose::Pose & pose );

	/// @brief Get name (MPFastRelaxMover)
	/// @details Get the name of this mover
	virtual std::string get_name() const;

private: // data

	// Keep a local copy of the relax protocol
	RelaxProtocolBaseOP relax_protocol_;

	// Set custon sfxn
	core::scoring::ScoreFunctionOP sfxn_;

}; // class MPFastRelaxMover

} // membrane
} //relax
} // protocols

#endif // INCLUDED_protocols_membrane_relax_MPFastRelaxMover_hh
