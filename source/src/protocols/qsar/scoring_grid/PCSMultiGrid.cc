// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/qsar/scoring_grid/PCSMultiGrid.cc
/// @brief   implementation of PCSMultiGrid
/// @details last Modified: 05/24/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <protocols/qsar/scoring_grid/PCSMultiGrid.hh>
#include <protocols/qsar/scoring_grid/PCSMultiGridCreator.hh>
#include <protocols/qsar/scoring_grid/PCSSingleGrid.hh>
//#include <protocols/qsar/scoring_grid/SingleGrid.hh>

// Package headers
#include <protocols/qsar/scoring_grid/schema_util.hh>
#include <core/io/nmr/util.hh>
#include <core/scoring/nmr/util.hh>
#include <core/scoring/nmr/pcs/PCSTensor.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/tools/make_vector.hh>
#include <utility/json_spirit/json_spirit_value.h>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/string.functions.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/nmr.OptionKeys.gen.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// C++ headers
#include <string>
#include <vector>
#include <fstream>

namespace protocols {
namespace qsar {
namespace scoring_grid {

std::string
PCSMultiGridCreator::keyname() const
{
	return PCSMultiGrid::grid_name();
}

GridBaseOP
PCSMultiGridCreator::create_grid(utility::tag::TagCOP tag) const
{
	GridBaseOP pcs_grid( new PCSMultiGrid() );
	pcs_grid->parse_my_tag(tag);
	return pcs_grid;
}

GridBaseOP
PCSMultiGridCreator::create_grid() const
{
	return GridBaseOP( new PCSMultiGrid() );
}

void
PCSMultiGridCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PCSMultiGrid::provide_xml_schema( xsd );
}

std::string
PCSMultiGrid::grid_name()
{
	return "PCSMultiGrid";
}

void
PCSMultiGrid::provide_xml_schema(utility::tag::XMLSchemaDefinition & xsd)
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute("grid_name", xs_string, "The name used to insert the scoring grid into the GridManager")
		+ XMLSchemaAttribute::required_attribute( "pcs_input_file", xs_string, "Main PCS input file.")
		+ XMLSchemaAttribute::attribute_w_default( "pcs_weight", xsct_real, "Multiply the PCS grid score by this weight value.", "1.0");

	xsd_type_definition_w_attributes( xsd, grid_name(), "A collection of PCS scoring grids. For each grid, PCS values are calculated from fixed input tensor parameters and compared with the ligand PCS values.", attributes );
}

/// @brief default constructor
PCSMultiGrid::PCSMultiGrid() :
	type_("PCSMultiGrid"),
	pcs_input_file_(""),
	pcs_grid_vector_(),
	weight_(1.0),
	pcs_data_initialized_(false)
{}

/// @brief construct from PCS input file
PCSMultiGrid::PCSMultiGrid(
	std::string const & filename,
	Real const weight
) :
	type_("PCSMultiGrid"),
	pcs_input_file_(filename),
	weight_(weight)
{
	initialize_pcs_data_from_input_file(filename);
	pcs_data_initialized_=true;
}

/// @brief copy constructor
PCSMultiGrid::PCSMultiGrid(PCSMultiGrid const & other) :
	GridBase(other),
	type_(other.type_),
	pcs_input_file_(other.pcs_input_file_),
	weight_(other.weight_),
	pcs_data_initialized_(other.pcs_data_initialized_)
{
	pcs_grid_vector_.clear();
	for ( Size i = 1, i_end = other.pcs_grid_vector_.size(); i <= i_end; ++i ) {
		pcs_grid_vector_.push_back( utility::pointer::dynamic_pointer_cast< SingleGrid >( (other.pcs_grid_vector_[i])->clone() ) );
	}
}

/// @brief copy assignment
PCSMultiGrid&
PCSMultiGrid::operator=(PCSMultiGrid const & rhs) {
	if ( this != &rhs ) {
		GridBase::operator=(rhs);
		type_ = rhs.type_;
		pcs_input_file_ = rhs.pcs_input_file_;
		pcs_grid_vector_.clear();
		for ( Size i = 1, i_end = rhs.pcs_grid_vector_.size(); i <= i_end; ++i ) {
			pcs_grid_vector_.push_back( utility::pointer::dynamic_pointer_cast< SingleGrid >( (rhs.pcs_grid_vector_[i])->clone() ) );
		}
		weight_ = rhs.weight_;
		pcs_data_initialized_ = rhs.pcs_data_initialized_;
	}
	return *this;
}

/// @brief destructor
PCSMultiGrid::~PCSMultiGrid() {}

/// @brief Make a copy of the grid, respecting the subclassing.
GridBaseOP
PCSMultiGrid::clone() const {
	return GridBaseOP( new PCSMultiGrid( *this ) );
}

/// @brief get the type of the grid
std::string
PCSMultiGrid::get_type() const {
	return grid_name();
}

/// @brief set the chain the grid applies to
void
PCSMultiGrid::set_chain(char chain) {
	for ( Size i = 1, i_end = pcs_grid_vector_.size(); i <= i_end; ++i ) {
		pcs_grid_vector_[i]->set_chain(chain);
	}
}

/// @brief determine if all residue atoms are in a grid
bool
PCSMultiGrid::is_in_grid(UltraLightResidue const & residue) const {
	for ( Size i = 1, i_end = pcs_grid_vector_.size(); i <= i_end; ++i ) {
		if ( !(pcs_grid_vector_[i]->is_in_grid(residue)) ) {
			return false;
		}
	}
	return true;
}

/// @brief determine if all residue atoms are in a grid
bool
PCSMultiGrid::is_in_grid(Residue const & residue) const {
	for ( Size i = 1, i_end = pcs_grid_vector_.size(); i <= i_end; ++i ) {
		if ( !(pcs_grid_vector_[i]->is_in_grid(residue)) ) {
			return false;
		}
	}
	return true;
}

/// @brief Print a brief summary about this grid to the provided output stream
void
PCSMultiGrid::show(std::ostream & out) const {
	out << "Combined PCS scoring grid vector contains " << pcs_grid_vector_.size() << " single grids." << std::endl;
}

/// @brief output a BRIX formatted grid. This really does not work well but is being left for legacy purposes
void
PCSMultiGrid::dump_BRIX(std::string const & /*prefix*/) const {
	utility_exit_with_message("ERROR: PCSMultiGrid is currently unable to output a single BRIX grid.");
}

/// @brief serialize the PCSMultiGrid to a json_spirit object
utility::json_spirit::Value
PCSMultiGrid::serialize() const {
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair type_record("type", Value(type_));
	Pair file_record("pcs_input_file", Value(pcs_input_file_));
	std::vector<Value> pcs_grid_vector;
	for ( Size i = 1, i_end = pcs_grid_vector_.size(); i <= i_end; ++i ) {
		pcs_grid_vector.push_back(pcs_grid_vector_[i]->serialize());
	}
	Pair grid_vector_record("pcs_grids", pcs_grid_vector);
	Pair weight_record("pcs_weight", Value(weight_));

	return Value(utility::tools::make_vector(type_record,file_record,grid_vector_record,weight_record));
}

/// @brief deserialize a json_spirit object to a PCSMultiGrid
void
PCSMultiGrid::deserialize(utility::json_spirit::mObject data) {
	type_ = data["type"].get_str();
	pcs_input_file_ = data["pcs_input_file"].get_str();
	utility::json_spirit::mArray grid_data(data["pcs_grids"].get_array());
	pcs_grid_vector_.resize(grid_data.size());
	Size i(1);
	for ( utility::json_spirit::mArray::iterator it = grid_data.begin(); it != grid_data.end(); ++it ) {
		SingleGridOP grid( new PCSSingleGrid() );
		grid->deserialize(it->get_obj());
		pcs_grid_vector_[i] = grid;
		++i;
	}
	weight_ = data["pcs_weight"].get_real();
}

/// @brief setup a PCSMultiGrid based on RosettaScripts input
void
PCSMultiGrid::parse_my_tag(utility::tag::TagCOP tag) {
	if ( tag->hasOption("pcs_input_file") ) {
		pcs_input_file_= tag->getOption< std::string >("pcs_input_file");
		initialize_pcs_data_from_input_file(pcs_input_file_);
		pcs_data_initialized_=true;
	}
	if ( tag->hasOption("pcs_weight") ) {
		weight_ = tag->getOption< Real >("pcs_weight");
	}
}

/// @brief setup a vector of PCSSingleGrid objects
///        initialize each PCSSingleGrid with a given center point,
///        width and resolution (in angstroms) and set grid values to zero.
void
PCSMultiGrid::initialize(
	Vector const & center,
	Real width,
	Real resolution
)
{
	if ( !pcs_data_initialized_ ) {
		initialize_pcs_data_from_input_file(pcs_input_file_);
		pcs_data_initialized_=true;
	}
	for ( Size i = 1, i_end = pcs_grid_vector_.size(); i <= i_end; ++i ) {
		pcs_grid_vector_[i]->initialize(center, width, resolution);
	}
}

/// @brief populate grids in the vector with PCS values based on a passed pose
void
PCSMultiGrid::refresh(
	Pose const & pose,
	Vector const & center
)
{
	for ( Size i = 1, i_end = pcs_grid_vector_.size(); i <= i_end; ++i ) {
		pcs_grid_vector_[i]->refresh(pose, center);
	}
}

/// @brief populate grids in the vector with PCS values based on a passed pose
void
PCSMultiGrid::refresh(
	Pose const & pose,
	Vector const & center,
	Size const & /*ligand_chain_id_to_exclude*/
)
{
	refresh(pose, center);
}

/// @brief populate grids in the vector with PCS values based on a passed pose
void
PCSMultiGrid::refresh(
	Pose const & pose,
	Vector const & center,
	utility::vector1<Size> /*ligand_chain_ids_to_exclude*/
)
{
	refresh(pose, center);
}

/// @brief return the current score of an UltraLightResidue using the current PCSMultiGrid
core::Real
PCSMultiGrid::score(
	UltraLightResidue const & residue,
	Real const max_score,
	qsarMapCOP qsar_map
) const
{
	Real total_score(0.0);
	for ( Size i = 1, i_end = pcs_grid_vector_.size(); i <= i_end; ++i ) {
		total_score += pcs_grid_vector_[i]->score(residue, max_score, qsar_map);
	}
	return total_score * weight_;
}

/// @brief return the current score of an atom using the current PCSMultiGrid
core::Real
PCSMultiGrid::atom_score(
	UltraLightResidue const & residue,
	Size atomno,
	qsarMapCOP qsar_map
) const
{
	Real total_score(0.0);
	for ( Size i = 1, i_end = pcs_grid_vector_.size(); i <= i_end; ++i ) {
		total_score += pcs_grid_vector_[i]->atom_score(residue, atomno, qsar_map);
	}
	return total_score * weight_;
}

/// @brief return the current score of a residue using the current PCSMultiGrid
core::Real
PCSMultiGrid::score(
	Residue const & residue,
	Real const max_score,
	qsarMapCOP qsar_map
) const
{
	Real total_score(0.0);
	for ( Size i = 1, i_end = pcs_grid_vector_.size(); i <= i_end; ++i ) {
		total_score += pcs_grid_vector_[i]->score(residue, max_score, qsar_map);
	}
	return total_score * weight_;
}

/// @brief return the current score of an atom using the current PCSMultiGrid
core::Real
PCSMultiGrid::atom_score(
	Residue const & residue,
	Size atomno,
	qsarMapCOP qsar_map
) const
{
	Real total_score(0.0);
	for ( Size i = 1, i_end = pcs_grid_vector_.size(); i <= i_end; ++i ) {
		total_score += pcs_grid_vector_[i]->atom_score(residue, atomno, qsar_map);
	}
	return total_score * weight_;
}

void
PCSMultiGrid::initialize_pcs_data_from_input_file(std::string const & filename) {
	using namespace core::io::nmr;
	using namespace core::scoring::nmr::pcs;
	using utility::to_string;
	using utility::string2Real;
	using ObjexxFCL::uppercased;

	std::ifstream infile;
	std::string line;
	Size line_number(0);

	infile.open(filename.c_str(), std::ios::in);
	if ( !infile.is_open() ) {
		utility_exit_with_message( "Unable to open PCS input file " + filename );
	}

	pcs_grid_vector_.clear();

	while ( std::getline(infile, line) ) {
		++line_number;

		// Remove leading and trailing whitespace and skip comments and blank lines.
		utility::trim( line, " \t\n" );
		if ( ( line.size() < 1 ) || ( line[ 0 ] == '#' ) ) {
			continue;
		} else if ( line.find("dataset") != std::string::npos ) {
			// Read parameter of one PCS dataset
			Size idx_key = line.find("dataset");
			Size idx_equal_sign = line.find("=", idx_key + 7);
			if ( idx_equal_sign == std::string::npos ) {
				utility_exit_with_message("ERROR: No equal sign found after parameter \"dataset\" in PCS input file. Please provide input parameter \" key = value \".");
			} else {
				utility::vector1<std::string>  dataset_params = read_pcs_dataset_params_list(utility::trim(line.substr(idx_equal_sign+1)));
				// The PCS dataset string contains 14 input parameter but we need only a few, i.e. the filename, weight and tensor values
				std::string pcs_datafile     = dataset_params[1];
				Real wt                      = string2Real(dataset_params[3]);
				utility::vector1<Real> tensor_vals(8);
				tensor_vals[1] = string2Real(dataset_params[12]);  // alpha
				tensor_vals[2] = string2Real(dataset_params[13]);  // beta
				tensor_vals[3] = string2Real(dataset_params[14]);  // gamma
				tensor_vals[4] = string2Real(dataset_params[7]);  // xM
				tensor_vals[5] = string2Real(dataset_params[8]);  // yM
				tensor_vals[6] = string2Real(dataset_params[9]);  // zM
				tensor_vals[7] = string2Real(dataset_params[10]);  // Xax
				tensor_vals[8] = string2Real(dataset_params[11]);  // Xrh

				PCSTensorOP tensor( new PCSTensor );
				tensor->set_tensor_in_pas(tensor_vals);
				SingleGridOP grid( new PCSSingleGrid(pcs_datafile, tensor, wt) );
				pcs_grid_vector_.push_back(grid);
			}
		}
	}
	pcs_data_initialized_=true;
}

std::string
PCSMultiGrid::hash_fingerprint() const {
	std::stringstream ss;
	const char sep('\t');
	ss << grid_name();
	ss << sep << type_;
	ss << sep << pcs_input_file_ << sep << pcs_data_initialized_ << sep << weight_;
	ss << sep << pcs_grid_vector_.size();
	for ( auto const & entry: pcs_grid_vector_ ) {
		ss << sep << entry->hash_fingerprint();
	}
	return ss.str();
}

} // namespace scoring_grid
} // namespace qsar
} // namespace protocols
