// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/AddMembraneMover.hh
///
/// @brief      Initialize the RosettaMP Framework by adding membrane representations to the pose
/// @details Given a pose, initialize and configure with the RosettaMP framework by taking the
///    following steps:
///     (1) Add a membrane residue to the pose (type MEM)
///      (a) Append residue to the pose or accept a new one
///      (b) Update PDB info to acknowledge the membrane residue
///      (c) Set the MEM residue at the root of the foldtree
///     (2) Initialize transmembrane spanning topology (Object: SpanningTopology)
///     (3) Initialize the MembraneInfo object either with or without the LipidAccInfo object.
///     (4) Set the membrane starting position (either default or based on user input)
///
///    This object does a massive amount of reading from CMD, RosettaScripts or Constructor. If you add
///    a new piece of data - you must modify MembraneInfo, all of these channels for input AND the apply function!
///    If and only if AddMembraneMover is applied to the pose, pose.conformation().is_membrane() MUST return true.
///
///    Last Updated: 7/23/15
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// @author  JKLeman (julia.koehler.leman@gmail.com)

#ifndef INCLUDED_protocols_membrane_AddMembraneMover_hh
#define INCLUDED_protocols_membrane_AddMembraneMover_hh

// Unit Headers
#include <protocols/membrane/AddMembraneMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/conformation/membrane/SpanningTopology.fwd.hh>
#include <core/conformation/membrane/Span.fwd.hh>
#include <core/conformation/membrane/LipidAccInfo.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

// C++ Headers
#include <cstdlib>

namespace protocols {
namespace membrane {

/// @brief Initialize the RosettaMP framework by adding membrane components to the pose
class AddMembraneMover : public protocols::moves::Mover {

public: // Constructors & General Setup

	/// @brief Create a default RosettaMP membrane setup
	/// @brief Create a membrane positioned at the origin (0, 0, 0) and aligned with
	/// the z axis. Use a defualt lipid type DOPC.
	AddMembraneMover();

	/// @brief Create a RosettaMP setup from an existing membrane residue
	/// @brief Create a membrane using the position from the existing membrane residue.
	/// Use a defualt lipid type DOPC.
	AddMembraneMover( core::Size membrane_rsd );

	/// @brief Create a RosettaMP setup from a user specified spanfile
	/// @brief Create a membrane positioned at the origin (0, 0, 0) and aligned with
	/// the z axis. Use a defualt lipid type DOPC. Load spanning topology from the user
	/// specified spanfile
	AddMembraneMover(
		std::string const & spanfile,
		core::Size membrane_rsd=0
	);

	/// @brief Create a RosettaMP setup from a user specified SpanningTopology
	/// @brief Create a membrane positioned at the origin (0, 0, 0) and aligned with
	/// the z axis. Use a defualt lipid type DOPC. Load spanning topology from the user
	/// specified spanning topology
	AddMembraneMover(
		core::conformation::membrane::SpanningTopologyOP topology,
		core::Size anchor_rsd=1,
		core::Size membrane_rsd=0
	);

	/// @brief Create a RosettaMP setup from an existing residue at a specific anchor point
	/// @brief Create a membrane using the position from the existing membrane residue. Anchor
	/// this residue at the user-specified anchor residue. Use a defualt lipid type DOPC.
	AddMembraneMover(
		core::Size anchor_rsd,
		core::Size membrane_rsd
	);

	/// @brief Create a RosettaMP setup from a user specified spanfile and lipsfile
	/// @brief Create a membrane positioned at the origin (0, 0, 0) and aligned with
	/// the z axis. Use a defualt lipid type DOPC. Load spanning topology from the user
	/// specified spanfile and lipsfile
	AddMembraneMover(
		std::string const & spanfile,
		std::string const & lipsfile,
		core::Size membrane_rsd=0
	);

	/// @brief Create a RosettaMP setup from a user specified spanfile and lipsfile
	/// @brief Create a membrane positioned at "init_center" and aligned with
	/// "init_normal". Use a defualt lipid type DOPC.
	AddMembraneMover(
		core::Vector const & init_center,
		core::Vector const & init_normal,
		std::string const & spanfile = "",
		core::Size membrane_rsd = 0
	);

	/// @brief Create a deep copy of the data in this mover
	AddMembraneMover( AddMembraneMover const & src );

	/// @brief Destructor
	~AddMembraneMover() override;

public: // Mover methods, getters & setters

	/// @brief Get the name of this Mover (AddMembraneMover)

	/// @brief Initialize the RosettaMP elements with this pose
	void apply( core::pose::Pose & pose ) override;

	/// @brief Return the current path to the spanfile held
	/// by this mover
	std::string get_spanfile() const;

	/// @brief Set Spanfile path
	/// @details Set the path to the spanfile
	void spanfile( std::string spanfile );

	/// @brief Set lipsfile path
	/// @details Set the path to the lipsfile
	void lipsfile( std::string lipsfile );

	/// @brief Set option for including lipophilicity data
	/// @details Incidate whether lipophilicity information should be read
	/// and used in MembraneInfo
	void include_lips( bool include_lips );

public: // Rosetta Scripts Methods

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
		core::pose::Pose const & pose
	) override;

	/// @brief Helper Method - Add a membrane virtual residue
	virtual Size add_membrane_virtual( core::pose::Pose & pose );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	static
	void
	attributes_for_parse_center_normal_from_tag( utility::tag::AttributeList & attributes );

private: // Setup Methods

	/// @brief Initialize Membrane Residue given pose
	virtual Size initialize_membrane_residue( core::pose::Pose & pose, core::Size membrane_rsd );

	/// @brief Helper Method - Check for Membrane residue already in the PDB
	virtual utility::vector1< core::SSize > check_pdb_for_mem( core::pose::Pose & pose );

	/// @brief Register options from JD2
	void register_options();

	/// @brief Initialize Mover options from the comandline
	void init_from_cmd();

private:

	// Pose residye typeset & include lips
	bool include_lips_;

	// SpanningTopology
	std::string spanfile_;
	core::conformation::membrane::SpanningTopologyOP topology_;

	// Lipid Accessibility Info - Lips Files
	std::string lipsfile_;

	// the membrane residue is anchored to this
	// residue position by jump
	core::Size anchor_rsd_;

	// Membrane residue number
	core::Size membrane_rsd_;

	// Initial Center/Normal Pair
	core::Vector center_;
	core::Vector normal_;

	// Set membrane thickness & steepness (FOR HIGHRES ONLY!!!)
	core::Real thickness_;
	core::Real steepness_;
	core::Real membrane_core_;

	// User defined
	bool user_defined_;

	bool got_spans_from_xml_;
};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_AddMembraneMover_hh

