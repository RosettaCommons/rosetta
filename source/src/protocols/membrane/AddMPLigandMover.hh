// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/membrane/AddMPLigandMover.hh
///
/// @brief  Add "single" ligand to to membrane pose
/// @details  Accommodate membrane protein ligand in the membrane framework by
///    reorganizing the current foldtree. Resulting foldtree will
///    keep the membrane attached to the transmembrane COM and ligand to the
///    closest binding pocket residue, provided in the constructor.
///
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// @author JKLeman (julia.koehler1982@gmail.com)
/// #RosettaMPMover

#ifndef INCLUDED_protocols_membrane_AddMPLigandMover_hh
#define INCLUDED_protocols_membrane_AddMPLigandMover_hh

// Unit Headers
#include <protocols/membrane/AddMPLigandMover.fwd.hh>

#include <protocols/moves/Mover.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>

// C++ Headers
#include <cstdlib>

namespace protocols {
namespace membrane {

/// @brief Add Membrane protein ligand mover
/// @details Accommodate membrane protein ligand in the membrane foldtree
class AddMPLigandMover : public protocols::moves::Mover {

public: // constructors

	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief Add membrane protein ligand mover
	/// @details Attach ligand downstream in the foldtree
	/// for refinement at the last residue as a default.
	/// DO NOT USE
	AddMPLigandMover();

	/// @brief Add Membrane protein ligand mover (custom)
	/// @details Attach ligand downstream in the foldtree of the
	/// closest residue to the binding pocket. Also specify
	/// sequence position of the ligand
	AddMPLigandMover( core::Size closest_rsd, core::Size ligand_seqpos );

	/// @brief Copy Constructor
	/// @details Mkae a deep copy of this mover
	AddMPLigandMover( AddMPLigandMover const & src );

	/// @brief Destructor
	~AddMPLigandMover() override;

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
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Mover Apply Method
	void apply( core::pose::Pose & pose ) override;

	/// @brief Show the name of this mvoer
	std::string get_name() const override;

private:

	core::Size closest_rsd_;
	core::Size ligand_seqpos_;

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_AddMPLigandMover_hh



