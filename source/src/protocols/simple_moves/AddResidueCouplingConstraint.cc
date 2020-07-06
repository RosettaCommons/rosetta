// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/AddResidueCouplingConstraint.cc
/// @brief
/// @author Moritz Ertelt

/*
* The mover is designed to improve the design process by using information from co-evolutionary
* correlation matrices. These matrices can be created by different methods, e.g. GREMLIN or plmDCA.
* The matrices get input as tensor of the shape [N, A, N, A] where N is the length of the protein
* and A is the length of the alphabet (typically 21 amino acids + 1 for gaps).
*
* GRELIN alphabet order: ARNDCQEGHILKMFPSTWYV-
* plmDCA alphabet order: ACDEFGHIKLMNPQRSTVWY-
*
* The order of the alphabet used can be defined by the user. To access the measure of co-evolution
* between, e.g., residue 1 (ALA) and residue 2 (ARG) the tensor lookup would be [1, 0, 2, 1]
* (using the GREMLIN alphabet order). The co-evolution strength is multiplied by a user defined
* strength factor (default=1) and then added as new term to the scoring function.
*/

#include <protocols/simple_moves/AddResidueCouplingConstraint.hh>
#include <protocols/simple_moves/AddResidueCouplingConstraintCreator.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/constraints/ResidueCouplingConstraint.hh>
#include <core/scoring/constraints/KofNConstraint.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/mover_schemas.hh>
#include <utility/io/izstream.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <numeric/MathNTensor_io.hh>
#include <numeric/MathNTensor.hh>
#include <utility/json_spirit/json_spirit.h>
#include <utility/json_spirit/json_spirit_tools.hh>

namespace protocols {
namespace simple_moves {

static basic::Tracer TR("protocols.simple_moves.AddResidueCouplingConstraint");

std::string AddResidueCouplingConstraintCreator::keyname() const {
	return AddResidueCouplingConstraint::mover_name();
}
protocols::moves::MoverOP AddResidueCouplingConstraintCreator::create_mover() const {
	return protocols::moves::MoverOP(new AddResidueCouplingConstraint);
}

void AddResidueCouplingConstraintCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	AddResidueCouplingConstraint::provide_xml_schema( xsd );
}

// setting up function and associated variables+typ of variable
AddResidueCouplingConstraint::AddResidueCouplingConstraint() :
	tensor_file_(""),
	index_file_(""),
	alphabet_("ARNDCQEGHILKMFPSTWYV-"),
	strength_(1.0){
}

AddResidueCouplingConstraint::AddResidueCouplingConstraint(
	std::string const& tensor_file,
	std::string const& index_file,
	core::Real strength,
	std::string const& alphabet) :
	tensor_file_(tensor_file),
	index_file_(index_file),
	alphabet_(alphabet),
	strength_(strength){
}


protocols::moves::MoverOP AddResidueCouplingConstraint::clone() const {
	return std::make_shared<AddResidueCouplingConstraint>(*this);
}


// mover function, read in the tensor and pass it to constraints; create PosList
void AddResidueCouplingConstraint::apply(core::pose::Pose & pose) {
	using namespace core::scoring::constraints;


	//reading in the tensor file together with the JSON information
	CouplingTensorOP tensor(std::make_shared< CouplingTensor >());
	utility::json_spirit::mObject     mObject_json;
	numeric::read_tensor_from_file(tensor_file_, *tensor, mObject_json);

	CouplingTensorOP tensor_const(utility::pointer::dynamic_pointer_cast< CouplingTensor >(tensor));

	// reading the index list, if an index file was specified
	std::vector<int> index_list;

	// Trying to read optional index file, used by GREMLIN, not others (like plmDCA)
	if ( !index_file_.empty() ) {
		std::ifstream input_file(index_file_);
		if ( !input_file ) {
			utility_exit_with_message("Unable to open index file: '" + index_file_ + "'");
		}

		int value;
		while ( input_file >> value ) {
			index_list.push_back(value);
		}
	}

	for ( auto it = index_list.begin(); it != index_list.end(); ++it ) {
		for ( auto jt = it + 1; jt != index_list.end(); ++jt ) {
			core::Size seqpos1 = *it + 1;
			core::Size seqpos2 = *jt + 1;

			core::Size tensor_index1 =  it - index_list.begin();
			core::Size tensor_index2 = jt - index_list.begin();

			TR << seqpos1 << " " << seqpos2 << std::endl;
			TR << tensor_index1 << " " << tensor_index2 << std::endl;
			TR << alphabet_ << std::endl;

			auto linking_constraint = std::make_shared<core::scoring::constraints::ResidueCouplingConstraint>(
				pose, seqpos1, seqpos2, tensor_index1, tensor_index2, tensor_const, strength_, alphabet_);

			pose.add_constraint(linking_constraint);
		}
	}
}


std::string AddResidueCouplingConstraint::get_name() const {
	return AddResidueCouplingConstraint::mover_name();
}

std::string AddResidueCouplingConstraint::mover_name() {
	return "AddResidueCouplingConstraint";
}

void AddResidueCouplingConstraint::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/) {
	tensor_file_ = tag->getOption< std::string >("tensor_file");
	index_file_ = tag->getOption< std::string >("index_file", "");
	strength_ = tag->getOption< core::Real >("strength", 1.0);
	alphabet_ = tag->getOption< std::string >("alphabet", "ARNDCQEGHILKMFPSTWYV-");
}

void AddResidueCouplingConstraint::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;

	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute(
		"tensor_file", xs_string,
		"Path to the co-evolution tensor file");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"index_file", xs_string,
		"Index file indicating residue positions in the tensor (used by Gremlin)", "");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"strength", xsct_real,
		"Strength factor, gets multiplied with the eo-evolution value (tensor cell)", "1.0");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"alphabet", xs_string,
		"Alphabet order of the tensor. The default is the Gremlin alphabet.","ARNDCQEGHILKMFPSTWYV-");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Mover which adds ResidueCouplingConstraints which favors residue types such, that co-evolving residue pairs are retained during resign."
		"To enable the constraint type the res_type_linking_constraint term must be reweighted."
		"This mover was developed for compatibility with the Gremlin tensor format consisting of the tensor (bin) and an index file",
		attlist);
}

} // namespace simple_moves
} // namespace protocols
