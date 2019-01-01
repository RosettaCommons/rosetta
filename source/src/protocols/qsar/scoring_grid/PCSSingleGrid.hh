// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/qsar/scoring_grid/PCSSingleGrid.hh
/// @brief   class that setups a single scoring grid with PseudoContactShift (PCS) values
/// @details last Modified: 05/17/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_qsar_scoring_grid_PCSSingleGrid_HH
#define INCLUDED_protocols_qsar_scoring_grid_PCSSingleGrid_HH

// Unit headers
#include <protocols/qsar/scoring_grid/PCSSingleGrid.fwd.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>

// Package headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/UltraLightResidue.fwd.hh>
#include <protocols/qsar/qsarMap.fwd.hh>
#include <core/scoring/nmr/pcs/PCSSingle.fwd.hh>
#include <core/scoring/nmr/pcs/PCSTensor.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

// C++ headers
#include <string>

namespace protocols {
namespace qsar {
namespace scoring_grid {

class PCSSingleGrid : public SingleGrid {

public: // typedefs

	typedef core::scoring::nmr::pcs::PCSSingle PCSSingle;
	typedef core::scoring::nmr::pcs::PCSSingleOP PCSSingleOP;
	typedef core::scoring::nmr::pcs::PCSTensor PCSTensor;
	typedef core::scoring::nmr::pcs::PCSTensorOP PCSTensorOP;
	typedef core::scoring::nmr::pcs::PCSTensorCOP PCSTensorCOP;
	typedef core::pose::Pose Pose;
	typedef core::conformation::Residue Residue;
	typedef core::conformation::UltraLightResidue UltraLightResidue;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::Vector Vector;

public: // Methods

	/// @brief Default constructor
	PCSSingleGrid();

	/// @brief construct from PCS datafile and fixed PCS tensor
	PCSSingleGrid(
		std::string const & filename,
		PCSTensorOP tensor,
		Real const weight = 1.0
	);

	/// @brief construct from PCS datafile and vector of PCS tensor values
	PCSSingleGrid(
		std::string const & filename,
		utility::vector1<Real> const & tensor_vals,
		Real const weight = 1.0
	);

	/// @brief construct from PCS datafile and fixed PCS tensor
	PCSSingleGrid(
		std::string const & filename,
		PCSTensorOP tensor,
		Pose const & pose,
		Real const weight = 1.0
	);

	/// @brief construct from PCS datafile and vector of PCS tensor values
	PCSSingleGrid(
		std::string const & filename,
		utility::vector1<Real> const & tensor_vals,
		Pose const & pose,
		Real const weight = 1.0
	);

	/// @brief copy constructor
	PCSSingleGrid(PCSSingleGrid const & other);

	/// @brief copy assignment
	PCSSingleGrid& operator=(PCSSingleGrid const & rhs);

	/// @brief Make a copy of the grid, respecting the subclassing.
	GridBaseOP clone() const override;

	/// @brief populate the grid with PCS values based on a passed pose
	void
	refresh(
		Pose const & pose,
		Vector const & center
	) override;

	/// @brief populate the grid with PCS values based on a passed pose
	void
	refresh(
		Pose const & pose,
		Vector const & center,
		Size const & ligand_chain_id_to_exclude
	) override;

	/// @brief populate the grid with PCS values based on a passed pose
	void
	refresh(
		Pose const & pose,
		Vector const & center,
		utility::vector1<Size> ligand_chain_ids_to_exclude
	) override;

	/// @brief return the current score of an UltraLightResidue using the current PCSSingleGrid
	Real
	score(
		UltraLightResidue const & residue,
		Real const max_score,
		qsarMapCOP qsar_map
	) const override;

	/// @brief return the current score of an atom using the current PCSSingleGrid
	Real
	atom_score(
		UltraLightResidue const & residue,
		Size atomno,
		qsarMapCOP qsar_map
	) const override;

	/// @brief return the current score of a residue using the current PCSSingleGrid
	Real
	score(
		Residue const & residue,
		Real const max_score,
		qsarMapCOP qsar_map
	) const override;

	/// @brief return the current score of an atom using the current PCSSingleGrid
	Real
	atom_score(
		Residue const & residue,
		Size atomno,
		qsarMapCOP qsar_map
	) const override;

	/// @brief serialize the PCSSingleGrid to a json_spirit object
	utility::json_spirit::Value serialize() const override;

	/// @brief deserialize a json_spirit object to a PCSSingleGrid
	void deserialize(utility::json_spirit::mObject data) override;

	/// @brief setup a PCSSingleGrid based on RosettaScripts input
	void parse_my_tag(utility::tag::TagCOP tag) override;

	std::string hash_fingerprint() const override;

	static std::string grid_name();
	static void provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd);

	utility::vector1<PCSSingleOP> const & get_pcs_values() const { return pcs_values_; }
	PCSTensorCOP get_tensor() const { return tensor_; }
	Real get_weight() const { return weight_; }
	void set_weight(Real w) { weight_ = w; }

private: // Methods

	void init_pcs_values_from_file(std::string const & filename, Pose const & pose);

private: // Data

	std::string pcs_file_;
	utility::vector1<PCSSingleOP> pcs_values_;
	PCSTensorOP tensor_;
	Real weight_;

};

} // namespace scoring_grid
} // namespace qsar
} // namespace protocols

#endif //INCLUDED_protocols_qsar_scoring_grid_PCSSingleGrid_hh
