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
/// @author Moritz Ertelt, moritz.ertelt@vanderbilt.edu

#ifndef INCLUDED_protocols_simple_moves_AddResidueCouplingConstraint_hh
#define INCLUDED_protocols_simple_moves_AddResidueCouplingConstraint_hh

#include <protocols/simple_moves/AddResidueCouplingConstraint.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <string>

namespace protocols {
namespace simple_moves {

/*!
The mover is designed to improve the design process by using information from co-evolutionary correlation matrices. These matrices can be created by different methods, e.g. GREMLIN or plmDCA. The matrices get input as tensor of the shape [N, A, N, A] where N is the length of the protein and A is the length of the alphabet (typically 21 amino acids + 1 for gaps).

GRELIN alphabet order: ARNDCQEGHILKMFPSTWYV-
plmDCA alphabet order: ACDEFGHIKLMNPQRSTVWY-

The order of the alphabet used can be defined by the user. To access the measure of co-evolution between, e.g., residue 1 (ALA) and residue 2 (ARG) the tensor lookup would be [1, 0, 2, 1] (using the GREMLIN alphabet order). The co-evolution strength is multiplied by a user defined strength factor (default=1) and then added as new term to the scoring function.
*/
class AddResidueCouplingConstraint :
	public moves::Mover {
public:
	AddResidueCouplingConstraint();
	AddResidueCouplingConstraint(std::string const& tensor_file, std::string const& index_file, core::Real strength, std::string const& alphabet);
	AddResidueCouplingConstraint(AddResidueCouplingConstraint const &) = default;
	~AddResidueCouplingConstraint() override {}

	protocols::moves::MoverOP clone() const override;

	void apply(core::pose::Pose & pose) override;
	std::string get_name() const override;

	void parse_my_tag (
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data) override;

	static std::string mover_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	std::string tensor_file_;
	std::string index_file_;
	std::string alphabet_;
	core::Real strength_;
};

}  // namespace simple_moves
}  // namespace protocols

#endif
