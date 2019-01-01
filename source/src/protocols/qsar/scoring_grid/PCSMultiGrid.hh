// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/qsar/scoring_grid/PCSMultiGrid.hh
/// @brief   class that setups a vector of scoring grids with PseudoContactShift (PCS) values
/// @details last Modified: 05/24/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_qsar_scoring_grid_PCSMultiGrid_HH
#define INCLUDED_protocols_qsar_scoring_grid_PCSMultiGrid_HH

// Unit headers
#include <protocols/qsar/scoring_grid/PCSMultiGrid.fwd.hh>
#include <protocols/qsar/scoring_grid/GridBase.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/UltraLightResidue.fwd.hh>
#include <protocols/qsar/qsarMap.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class PCSMultiGrid : public GridBase {

public: // Typedefs

	typedef core::pose::Pose Pose;
	typedef core::conformation::Residue Residue;
	typedef core::conformation::UltraLightResidue UltraLightResidue;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::Vector Vector;

public: // Methods

	/// @brief default constructor
	PCSMultiGrid();

	/// @brief construct from PCS input file
	PCSMultiGrid(
		std::string const & filename,
		Real const weight = 1.0
	);

	/// @brief copy constructor
	PCSMultiGrid(PCSMultiGrid const & other);

	/// @brief copy assignment
	PCSMultiGrid& operator=(PCSMultiGrid const & rhs);

	/// @brief destructor
	~PCSMultiGrid() override;

	/// @brief Make a copy of the grid, respecting the subclassing.
	GridBaseOP clone() const override;

	/// @brief setup a vector of PCSSingleGrid objects
	///        initialize each PCSSingleGrid with a given center point,
	///        width and resolution (in angstroms) and set grid values to zero.
	void
	initialize(
		Vector const & center,
		Real width,
		Real resolution
	) override;

	/// @brief populate grids in the vector with PCS values based on a passed pose
	void
	refresh(
		Pose const & pose,
		Vector const & center
	) override;

	/// @brief populate grids in the vector with PCS values based on a passed pose
	void
	refresh(
		Pose const & pose,
		Vector const & center,
		Size const & ligand_chain_id_to_exclude
	) override;

	/// @brief populate grids in the vector with PCS values based on a passed pose
	void
	refresh(
		Pose const & pose,
		Vector const & center,
		utility::vector1<Size> ligand_chain_ids_to_exclude
	) override;

	/// @brief return the current score of an UltraLightResidue using the current PCSMultiGrid
	Real
	score(
		UltraLightResidue const & residue,
		Real const max_score,
		qsarMapCOP qsar_map
	) const override;

	/// @brief return the current score of an atom using the current PCSMultiGrid
	Real
	atom_score(
		UltraLightResidue const & residue,
		Size atomno,
		qsarMapCOP qsar_map
	) const override;

	/// @brief return the current score of a residue using the current PCSMultiGrid
	Real
	score(
		Residue const & residue,
		Real const max_score,
		qsarMapCOP qsar_map
	) const override;

	/// @brief return the current score of an atom using the current PCSMultiGrid
	Real
	atom_score(
		Residue const & residue,
		Size atomno,
		qsarMapCOP qsar_map
	) const override;

	/// @brief get the type of the grid
	std::string get_type() const override;

	/// @brief set the chain the grid applies to
	void set_chain(char chain) override;

	/// @brief output a BRIX formatted grid.  This really does not work well but is being left for legacy purposes
	void dump_BRIX(std::string const & prefix) const override;

	/// @brief serialize the PCSMultiGrid to a json_spirit object
	utility::json_spirit::Value serialize() const override;

	/// @brief deserialize a json_spirit object to a PCSMultiGrid
	void deserialize(utility::json_spirit::mObject data) override;

	/// @brief setup a PCSMultiGrid based on RosettaScripts input
	void parse_my_tag(utility::tag::TagCOP tag) override;

	/// @brief determine if all residue atoms are in a grid
	bool is_in_grid(UltraLightResidue const & residue) const override;

	/// @brief determine if all residue atoms are in a grid
	bool is_in_grid(Residue const & residue) const override;

	static std::string grid_name();
	static void provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd);

	std::string hash_fingerprint() const override;

	/// @brief Print a brief summary about this grid to the provided output stream
	void show(std::ostream & out) const override;

	Size get_number_pcs_grids() const { return pcs_grid_vector_.size(); }
	utility::vector1< SingleGridOP > const & get_pcs_grids() const { return pcs_grid_vector_; }
	utility::vector1< SingleGridOP > & get_pcs_grids() { return pcs_grid_vector_; }
	Real get_weight() const { return weight_; }
	void set_weight(Real w) { weight_ = w; }

private: // Methods

	void initialize_pcs_data_from_input_file(std::string const & filename);

private: // Data

	std::string type_;
	std::string pcs_input_file_;
	// Vector of PCS grid and weight
	utility::vector1< SingleGridOP > pcs_grid_vector_;
	Real weight_;
	bool pcs_data_initialized_;

};

} // namespace scoring_grid
} // namespace qsar
} // namespace protocols

#endif //INCLUDED_protocols_qsar_scoring_grid_PCSMultiGrid_hh
