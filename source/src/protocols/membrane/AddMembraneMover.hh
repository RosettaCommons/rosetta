// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/membrane/AddMembraneMover.hh
/// @brief   Initialize the RosettaMP framework by adding representations of the membrane
///          environemnt to the conformation data held by the pose
///
/// @details Given a pose, configure RosettaMP by adding the following information:
///             1. Create and append a membrane residue (MEM) to the pose
///             2. Create and store a SpanningTopology object
///             3. Setup the initial membrane coordinates (typically centered at the origin)
///             4. (Optional) Initialize per-atom lipid accessibility data
///             5. (Optional) Initialize dimensions of the aqueous pore
///    6. Initialize the ImplicitLipidMembraneInfo
///          Upon completion, the call pose.conformation().is_membrane() will return true
///
/// @note    If you add a new step, please document and ensure all data is properly
///          initialized by constructors, parse_my_tag, init_from_cmd, serialization
///          routines, and xsd routines. This class is a data loading mammoth
///
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

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// C++ Headers
#include <cstdlib>

namespace protocols {
namespace membrane {

/// @brief   Initialize the RosettaMP framework by adding representations of the membrane
class AddMembraneMover : public protocols::moves::Mover {

public: // Constructors

	/// @brief Create a default RosettaMP membrane setup
	AddMembraneMover();

	/// @brief Create a RosettaMP setup from a user specified spanfile
	AddMembraneMover(
		std::string spanfile,
		core::Size membrane_rsd = 0,
		std::string lipid_composition = "DLPC",
		core::Real temperature = 37.0
	);

	/// @brief Create a RosettaMP setup from an existing membrane residue
	AddMembraneMover(
		core::Size membrane_rsd,
		std::string lipid_composition="DLPC",
		core::Real temperature=37.0
	);

	/// @brief Create a RosettaMP setup from a user specified SpanningTopology
	AddMembraneMover(
		core::conformation::membrane::SpanningTopologyOP topology,
		core::Size anchor_rsd=1,
		core::Size membrane_rsd=0
	);

	/// @brief Create a RosettaMP setup by setting the root of the fold tree to anchor_rsd and
	/// designating MEM as membrane_rsd. The two booleans are dummy variables intended to differentiate
	// this constructor from other constructors that appear equivalent to the compiler.
	AddMembraneMover(
		core::Size anchor_rsd,
		bool anchor_format,
		bool second_anchor_format,
		core::Size membrane_rsd=0
	);

	/// @brief Create a deep copy of the data in this mover
	AddMembraneMover( AddMembraneMover const & src );

	/// @brief Destructor
	~AddMembraneMover() override;

public: // mover-specific methods & rosetta scripts support

	/// @brief Initialize the RosettaMP elements with this pose
	void apply( core::pose::Pose & pose ) override;

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Pase Rosetta Scripts Options for this Mover
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	) override;

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

public: // read & write access to data

	/// @brief Path to file with transmembrane span data
	std::string get_spanfile() const;

	/// @brief Set path to spanfile
	void spanfile( std::string spanfile );

	// Set restore variable
	void restore_lazaridis_IMM1_behavior( bool restore );

public: // functions that can be overloaded for sub-classes of this mover

	/// @brief Create and append a membrane virtual residue to the pose
	virtual
	Size
	add_membrane_virtual( core::pose::Pose & pose );

private: // helper methods for setup

	/// @brief Initialize Membrane Residue given pose
	virtual
	Size
	initialize_membrane_residue( core::pose::Pose & pose, core::Size membrane_rsd );

	/// @brief Helper Method - Check for Membrane residue already in the PDB
	virtual
	utility::vector1< core::SSize >
	check_pdb_for_mem( core::pose::Pose & pose );

	/// @brief Register options from JD2
	void register_options();

	/// @brief Initialize Mover options from the comandline
	void init_from_cmd();

	/// @brief Read path to the spanfile from an input sql database
	std::string read_spanfile_from_db();

private:

	// Restore Lazaridis IMM1 energy function behavior?
	bool restore_lazaridis_IMM1_behavior_;

	// Information about the lipid composition
	std::string lipid_composition_;
	core::Real temperature_;
	bool user_override_pore_;

	// SpanningTopology
	std::string spanfile_;
	core::conformation::membrane::SpanningTopologyOP topology_;
	bool read_spans_from_xml_;

	// the membrane residue is anchored to this residue posiiton by jump
	core::Size anchor_rsd_;

	// Membrane residue number
	core::Size membrane_rsd_;

	// Initial Center/Normal Pair
	core::Vector center_;
	core::Vector normal_;
	bool user_defined_membrane_pos_;

	// Set membrane thickness & steepness
	core::Real thickness_;
	core::Real steepness_;
	core::Real membrane_core_;

	// Optional: members for initializing memrbane data from a relational database
	std::string database_table_;
	utility::sql_database::sessionOP db_session_;

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_AddMembraneMover_hh
