// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/symmetry/SymmetricAddMembraneMover.cc
///
/// @brief      Add Membrane Representation to a Symmetric starting pose
/// @details    Given a symmetrized pose, add the membrane components,
///             ensuring all data descriptions are for the asymmetric subunit
///             (in terms of scoring) and the membrane is the root of the
///             whole system. After applying SymmetricAddMembraneMover
///             pose.conformaiton().is_membrane() AND is_symmetric( pose )
///             should both return true
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (2/9/15)

#ifndef INCLUDED_protocols_membrane_symmetry_SymmetricAddMembraneMover_hh
#define INCLUDED_protocols_membrane_symmetry_SymmetricAddMembraneMover_hh

// Unit Headers
#include <protocols/membrane/symmetry/SymmetricAddMembraneMover.fwd.hh>

// Project Headers
#include <protocols/membrane/AddMembraneMover.hh>
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
namespace symmetry {

using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::membrane;
using namespace protocols::moves;

class SymmetricAddMembraneMover : public protocols::membrane::AddMembraneMover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor of SymmetricAddMembraneMover
	/// @details Create a membrane pose, setting the membrane center
	/// at center=(0, 0, 0), normal=(0, 0, 1) and loads in spans
	/// and lips from the command line interface. Calls the
	/// default constructor of AddMembraneMover
	SymmetricAddMembraneMover();

	/// @brief Custom Constructor for SymmetricAddMembraneMover
	/// @details Creates a membrane pose, setting the membrane center
	/// at center=(0, 0, 0), normal=(0, 0, 1). Uses the user provided
	/// spanfile in this constructor to load the spanning topology.
	/// Calls the analagous parent constructor in AddMembraneMover
	SymmetricAddMembraneMover( std::string spanfile );

	/// @brief Custom Constructor for SymmetricAddMembraneMover
	/// @details Creates a membrane pose, setting the membrane center
	/// at center=(0, 0, 0), normal=(0, 0, 1). Uses the user provided
	/// SpanningTopology to be set in MembraneInfo. Calls the analagous
	/// parent constructor in AddMembraneMover
	SymmetricAddMembraneMover( SpanningTopologyOP topology );

	/// @brief Custom Constructor for SymmetricAddMembraneMover
	/// @details Creates a membrane pose, setting the membrane center
	/// at center=(0, 0, 0), normal=(0, 0, 1). Uses the user provided
	/// spanfile and lipsfile paths to load in SpanningTopology
	/// objects and LipidAccInfo objects to MembraneInfo. Calls the analagous
	/// parent constructor in AddMembraneMover
	SymmetricAddMembraneMover( std::string spanfile, std::string lips_acc );

	/// @brief Copy Constructor for SymmetricAddMembraneMover
	/// @details Create a deep copy of SymmetricAddMembraneMover
	SymmetricAddMembraneMover( SymmetricAddMembraneMover const & src );

	/// @brief Destructor
	~SymmetricAddMembraneMover();

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP
	clone() const;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void
	parse_my_tag(
		utility::tag::TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (SymmetricAddMembraneMover)
	virtual
	std::string
	get_name() const;

	/// @brief Add a membrane virtual residue to the pose by inserting by jump
	/// @details Adds a virtual residue to the pose by selecting the VRT_0_Base
	/// as the anchoring residue and nres_complex+1 as tne new sequence position
	/// Not equivalent to an append_by_jump!
	virtual
	core::Size
	add_membrane_virtual( Pose & pose );

	/// @brief Helper Method - Check for Membrane residue already in the PDB
	/// @details If there is an MEM residue in the PDB at the end of the pose
	/// with property MEMBRANE, return a vector of all of those residues.
	virtual
	utility::vector1< core::SSize >
	check_pdb_for_mem( Pose & pose );

private:

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options from Command Line
	/// @details Register mover-relevant options with JD2 - includes
	/// mp seutp options: spanfiles, lipsfile, and
	/// spans_from_structure.
	void
	register_options();

	/// @brief Initialize Mover options from the comandline
	/// @details Initialize mover settings from the commandline
	/// mainly in the mp, setup group: spans_from_structure,
	/// spanfile and lipsfiles paths
	void
	init_from_cmd();

};

} // symmetry
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_SymmetricAddMembraneMover_hh

