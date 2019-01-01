// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/qsar/scoring_grid/PCSSingleGrid.cc
/// @brief   implementation of PCSSingleGrid
/// @details last Modified: 05/17/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <protocols/qsar/scoring_grid/PCSSingleGrid.hh>
#include <protocols/qsar/scoring_grid/PCSSingleGridCreator.hh>

// Package headers
#include <protocols/qsar/scoring_grid/schema_util.hh>
#include <core/io/nmr/AtomSelection.hh>
#include <core/io/nmr/util.hh>
#include <core/scoring/nmr/util.hh>
#include <core/scoring/nmr/pcs/PCSSingle.hh>
#include <core/scoring/nmr/pcs/PCSTensor.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/tools/make_vector.hh>
#include <utility/json_spirit/json_spirit_value.h>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <utility/fixedsizearray1.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/constants.hh>

// C++ headers
#include <string>
#include <vector>

namespace protocols {
namespace qsar {
namespace scoring_grid {

std::string
PCSSingleGridCreator::keyname() const
{
	return PCSSingleGrid::grid_name();
}

GridBaseOP
PCSSingleGridCreator::create_grid(utility::tag::TagCOP tag) const
{
	GridBaseOP pcs_grid( new PCSSingleGrid() );
	pcs_grid->parse_my_tag(tag);
	return pcs_grid;
}

GridBaseOP
PCSSingleGridCreator::create_grid() const
{
	return GridBaseOP( new PCSSingleGrid() );
}

void
PCSSingleGridCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PCSSingleGrid::provide_xml_schema( xsd );
}

std::string
PCSSingleGrid::grid_name()
{
	return "PCSSingleGrid";
}

void
PCSSingleGrid::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd)
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute("grid_name", xs_string, "The name used to insert the scoring grid into the GridManager.")
		+ XMLSchemaAttribute::required_attribute("tensor", xsct_real_cslist, "Comma-separated list of PCS tensor values from which the PCSSingleGrid is calculated. Tensor values should have the following order: Xax, Xrh, xM, yM, zM, alpha, beta, gamma.")
		+ XMLSchemaAttribute::required_attribute( "pcs_file", xs_string, "Textfile with ligand PCS values.")
		+ XMLSchemaAttribute::attribute_w_default( "pcs_weight", xsct_real, "Multiply the PCS grid score by this weight value.", "1.0");

	xsd_type_definition_w_attributes( xsd, grid_name(), "A scoring grid that contains PCS values which are calculated from the input tensor and compared with the ligand PCS values.", attributes );
}

/// @brief default constructor
PCSSingleGrid::PCSSingleGrid() :
	SingleGrid("PCSSingleGrid"),
	pcs_file_(""),
	pcs_values_(),
	tensor_( new PCSTensor ),
	weight_(1.0)
{}

/// @brief construct from PCS datafile and fixed PCS tensor
PCSSingleGrid::PCSSingleGrid(
	std::string const & filename,
	PCSTensorOP tensor,
	Real const weight
) :
	SingleGrid("PCSSingleGrid"),
	pcs_file_(filename),
	pcs_values_(),
	tensor_( new PCSTensor(*tensor) ),
	weight_(weight)
{}

/// @brief construct from PCS datafile and vector of PCS tensor values
PCSSingleGrid::PCSSingleGrid(
	std::string const & filename,
	utility::vector1<Real> const & tensor_vals,
	Real const weight
) :
	SingleGrid("PCSSingleGrid"),
	pcs_file_(filename),
	pcs_values_(),
	tensor_( new PCSTensor ),
	weight_(weight)
{
	tensor_->set_tensor_in_pas(tensor_vals);
	tensor_->reorder_tensor();
}

/// @brief construct from PCS datafile and fixed PCS tensor
PCSSingleGrid::PCSSingleGrid(
	std::string const & filename,
	PCSTensorOP tensor,
	Pose const & pose,
	Real const weight
) :
	PCSSingleGrid(filename, tensor, weight)
{
	init_pcs_values_from_file(filename, pose);
}

/// @brief construct from PCS datafile and vector of PCS tensor values
PCSSingleGrid::PCSSingleGrid(
	std::string const & filename,
	utility::vector1<Real> const & tensor_vals,
	Pose const & pose,
	Real const weight
) :
	PCSSingleGrid(filename, tensor_vals, weight)
{
	init_pcs_values_from_file(filename, pose);
}

/// @brief copy constructor
PCSSingleGrid::PCSSingleGrid(PCSSingleGrid const & other) :
	SingleGrid(other),
	pcs_file_(other.pcs_file_),
	tensor_( new PCSTensor( *(other.tensor_) ) ),
	weight_(other.weight_)
{
	pcs_values_.clear();
	pcs_values_.resize(other.pcs_values_.size());
	for ( Size i = 1, i_end = other.pcs_values_.size(); i <= i_end; ++i ) {
		pcs_values_[i] = PCSSingleOP( new PCSSingle( *(other.pcs_values_[i]) ) );
	}
}

/// @brief copy assignment
PCSSingleGrid&
PCSSingleGrid::operator=(PCSSingleGrid const & rhs) {
	if ( this != &rhs ) {
		SingleGrid::operator=(rhs);
		pcs_file_ = rhs.pcs_file_;
		pcs_values_.clear();
		pcs_values_.resize(rhs.pcs_values_.size());
		for ( Size i = 1, i_end = rhs.pcs_values_.size(); i <= i_end; ++i ) {
			pcs_values_[i] = PCSSingleOP( new PCSSingle( *(rhs.pcs_values_[i]) ) );
		}
		tensor_ = PCSTensorOP( new PCSTensor( *(rhs.tensor_) ) );
		weight_ = rhs.weight_;
	}
	return *this;
}

/// @brief Make a copy of the grid, respecting the subclassing.
GridBaseOP
PCSSingleGrid::clone() const {
	return GridBaseOP( new PCSSingleGrid( *this ) );
}

/// @brief populate the grid with PCS values based on a passed pose
void
PCSSingleGrid::refresh(
	Pose const & pose,
	Vector const & /*center*/
)
{
	if ( pcs_values_.empty() ) {
		init_pcs_values_from_file(pcs_file_, pose);
	}

	// rotation matrix from euler angles
	Vector euler_angles(tensor_->get_alpha(), tensor_->get_beta(), tensor_->get_gamma());
	numeric::xyzMatrix<core::Real> rotM = core::scoring::nmr::rotation_matrix_from_euler_angles(euler_angles, tensor_->get_euler_convention());
	utility::fixedsizearray1<Real,5> params = { tensor_->get_metal_center().x(),
		tensor_->get_metal_center().y(),
		tensor_->get_metal_center().z(),
		tensor_->get_ax(),
		tensor_->get_rh() };

	numeric::xyzVector<Size> dimensions = get_dimensions();
	for ( int x_index(0), x_stop(static_cast<int>(dimensions.x()-1)); x_index <= x_stop; ++x_index ) {

		for ( int y_index(0), y_stop(static_cast<int>(dimensions.y()-1)); y_index <= y_stop; ++y_index ) {

			for ( int z_index(0), z_stop(static_cast<int>(dimensions.z()-1)); z_index <= z_stop; ++z_index ) {

				numeric::xyzVector<int> grdPt(x_index, y_index, z_index);
				Vector box_center = get_grid().coords(grdPt);

				// The PCS at this grid point
				Real pcs = core::scoring::nmr::pcs_func(params, rotM, box_center);
				set_point(box_center, pcs);
			}
		}
	}
}

/// @brief populate the grid with PCS values based on a passed pose
void
PCSSingleGrid::refresh(
	Pose const & pose,
	Vector const & center,
	Size const & /*ligand_chain_id_to_exclude*/
)
{
	refresh(pose, center);
}

/// @brief populate the grid with PCS values based on a passed pose
void
PCSSingleGrid::refresh(
	Pose const & pose,
	Vector const & center,
	utility::vector1<Size> /*ligand_chain_ids_to_exclude*/
)
{
	refresh(pose, center);
}

/// @brief return the current score of an UltraLightResidue using the current PCSSingleGrid
core::Real
PCSSingleGrid::score(
	UltraLightResidue const & residue,
	Real const max_score,
	qsarMapCOP qsar_map
) const
{
	core::conformation::ResidueCOP base_residue = residue.residue();
	return score(*base_residue, max_score, qsar_map);
}

/// @brief return the current score of an atom using the current PCSSingleGrid
core::Real
PCSSingleGrid::atom_score(
	UltraLightResidue const & residue,
	Size atomno,
	qsarMapCOP qsar_map
) const
{
	core::conformation::ResidueCOP base_residue = residue.residue();
	return atom_score(*base_residue, atomno, qsar_map);
}

/// @brief return the current score of a residue using the current PCSSingleGrid
core::Real
PCSSingleGrid::score(
	Residue const & residue,
	Real const /*max_score*/,
	qsarMapCOP /*qsar_map*/
) const
{
	Real total_score(0.0);

	for ( Size i(1); i <= pcs_values_.size(); ++i ) {
		utility::vector1<core::id::AtomID> const & eq_spins = pcs_values_[i]->get_protein_spins();
		Size count(0);
		Real av_pcs(0.0);
		for ( Size j(1); j <= eq_spins.size(); ++j ) {

			// Is residue ID correct and atomno in residue?
			if ( residue.seqpos() == eq_spins[j].rsd() &&
					residue.natoms() >= eq_spins[j].atomno() ) {

				// Get residue atom xyz from AtomID
				Vector const & spin_xyz = residue.xyz(eq_spins[j].atomno());

				// Get PCS value found at this grid point
				if ( get_grid().is_in_grid(spin_xyz.x(), spin_xyz.y(), spin_xyz.z()) ) {
					// Average PCS in case of multiple equivalent spins
					av_pcs += get_grid().getValue(spin_xyz.x(), spin_xyz.y(), spin_xyz.z());
					++count;
				}
			}
		}
		if ( count ) {
			av_pcs /= count;
			Real dev = av_pcs - pcs_values_[i]->get_pcs_exp();
			total_score += dev*dev;
		}
	}
	return total_score * weight_;
}

/// @brief return the current score of an atom using the current PCSSingleGrid
core::Real
PCSSingleGrid::atom_score(
	Residue const & residue,
	Size atomno,
	qsarMapCOP /*qsar_map*/
) const
{
	Real atom_score(0.0);
	bool found_atom_in_pcs_data(false);

	for ( Size i(1); i <= pcs_values_.size(); ++i ) {
		if ( found_atom_in_pcs_data ) break;
		utility::vector1<core::id::AtomID> const & eq_spins = pcs_values_[i]->get_protein_spins();
		for ( Size j(1); j <= eq_spins.size(); ++j ) {

			// Is residue ID correct and found atom?
			if ( residue.seqpos() == eq_spins[j].rsd() && atomno == eq_spins[j].atomno() ) {
				found_atom_in_pcs_data=true;
				Vector const & spin_xyz = residue.xyz(eq_spins[j].atomno());

				// Get PCS value found at this grid point
				if ( get_grid().is_in_grid(spin_xyz.x(), spin_xyz.y(), spin_xyz.z()) ) {
					Real pcs = get_grid().getValue(spin_xyz.x(), spin_xyz.y(), spin_xyz.z());
					Real dev = pcs - pcs_values_[i]->get_pcs_exp();
					atom_score += dev*dev;
				}
				break;
			}
		}
	}
	return atom_score * weight_;
}

/// @brief serialize the PCSSingleGrid to a json_spirit object
utility::json_spirit::Value
PCSSingleGrid::serialize() const
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair base_data("base_data",SingleGrid::serialize());
	std::vector<Value> pcs_single_vector;
	for ( Size i = 1; i <= pcs_values_.size(); ++i ) {
		pcs_single_vector.push_back(pcs_values_[i]->serialize());
	}
	Pair pcs_value_record("pcs_values", pcs_single_vector);
	Pair file_record("pcs_file",Value(pcs_file_));
	Pair tensor_record("pcs_tensor",tensor_->serialize());
	Pair weight_record("pcs_weight",Value(weight_));

	return Value(utility::tools::make_vector(base_data,pcs_value_record,file_record,tensor_record,weight_record));
}

/// @brief deserialize a json_spirit object to a PCSSingleGrid
void
PCSSingleGrid::deserialize(utility::json_spirit::mObject data)
{
	SingleGrid::deserialize(data["base_data"].get_obj());

	utility::json_spirit::mArray pcs_data(data["pcs_values"].get_array());
	pcs_values_.resize(pcs_data.size());
	Size i(1);
	for ( utility::json_spirit::mArray::iterator it = pcs_data.begin(); it != pcs_data.end(); ++it ) {
		PCSSingleOP pcs( new PCSSingle( *(pcs_values_[i]) ) );
		pcs->deserialize(it->get_obj());
		pcs_values_[i] = pcs;
		++i;
	}
	pcs_file_ = data["pcs_file"].get_str();
	PCSTensorOP tensor( new PCSTensor(*tensor_) );
	tensor->deserialize(data["pcs_tensor"].get_obj());
	tensor_ = tensor;
	weight_ = data["pcs_weight"].get_real();
}

std::string
PCSSingleGrid::hash_fingerprint() const {
	std::stringstream ss;
	const char sep('\t');
	ss << grid_name();
	ss << sep << get_type(); // Only thing of interest from parent class
	ss << sep << pcs_file_ << sep << weight_;
	ss << sep << pcs_values_.size();
	ss << sep << tensor_->get_metal_center().x() << sep << tensor_->get_metal_center().y() << sep << tensor_->get_metal_center().z();
	ss << sep << tensor_->get_ax() << sep << tensor_->get_rh();
	ss << sep << tensor_->get_alpha() << sep << tensor_->get_beta() << sep << tensor_->get_gamma();
	return ss.str();
}

/// @brief setup a PCSSingleGrid based on RosettaScripts input
void
PCSSingleGrid::parse_my_tag(utility::tag::TagCOP tag)
{
	using utility::string2Real;

	if ( tag->hasOption("pcs_file") ) {
		pcs_file_= tag->getOption< std::string >("pcs_file");
	}

	if ( tag->hasOption("tensor") ) {
		std::string tensor_string = tag->getOption< std::string >("tensor");
		utility::vector1< std::string > tensor_str_vals = utility::string_split(tensor_string, ',');
		if ( tensor_str_vals.size() != 8 ) {
			utility_exit_with_message("ERROR during setup of PCSSingleGrid. Number of PCS tensor values must be 8.");
		}
		utility::vector1< Real > tensor_vals(8);
		tensor_vals[1] = string2Real(tensor_str_vals[6]); // alpha
		tensor_vals[2] = string2Real(tensor_str_vals[7]); // beta
		tensor_vals[3] = string2Real(tensor_str_vals[8]); // gamma
		tensor_vals[4] = string2Real(tensor_str_vals[3]); // xM
		tensor_vals[5] = string2Real(tensor_str_vals[4]); // yM
		tensor_vals[6] = string2Real(tensor_str_vals[5]); // zM
		tensor_vals[7] = string2Real(tensor_str_vals[1]); // Xax
		tensor_vals[8] = string2Real(tensor_str_vals[2]); // Xrh
		PCSTensorOP tensor( new PCSTensor );
		tensor->set_tensor_in_pas(tensor_vals);
		tensor->reorder_tensor();
		tensor_ = tensor;
	}

	if ( tag->hasOption("pcs_weight") ) {
		weight_ = tag->getOption< Real >("pcs_weight");
	}
}

void
PCSSingleGrid::init_pcs_values_from_file(
	std::string const & filename,
	Pose const & pose
)
{
	using namespace core::io::nmr;
	using namespace core::scoring::nmr::pcs;
	utility::vector1< utility::vector1< AtomSelection > > spins;
	utility::vector1< Real > values;
	utility::vector1< Real > errors;
	core::io::nmr::read_pcs_datafile(filename, spins, values, errors);
	runtime_assert_msg(spins.size() == values.size() && spins.size() == errors.size(), "Vector of spin selections and PCS values and errors must have the same size.");
	pcs_values_.reserve(spins.size());
	for ( Size i=1, i_end=spins.size(); i<=i_end; ++i ) {
		PCSSingleOP single_pcs( new PCSSingle(spins[i], pose, values[i], errors[i]) );
		pcs_values_.push_back(single_pcs);
	}
}

} // namespace scoring_grid
} // namespace qsar
} // namespace protocols
