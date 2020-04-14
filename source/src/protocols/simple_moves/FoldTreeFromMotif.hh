// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/FoldTreeFromMotif.hh
/// @author longxing

#ifndef INCLUDED_protocols_simple_moves_FoldTreeFromMotif_HH
#define INCLUDED_protocols_simple_moves_FoldTreeFromMotif_HH

// Unit headers
#include <protocols/simple_moves/FoldTreeFromMotif.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

namespace protocols {
namespace simple_moves {

class FoldTreeFromMotif : public protocols::moves::Mover
{
public:
	FoldTreeFromMotif();
	~FoldTreeFromMotif() override;

	void apply( Pose & pose ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap &
	) override;


	void
	residue_selector(
		core::select::residue_selector::ResidueSelectorCOP const & selector_in
	);

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:
	core::select::residue_selector::ResidueSelectorCOP motif_selector_;
};

} // namespace simple_moves
} // namespace protocols

#endif //INCLUDED_protocols_simple_moves_FoldTreeFromMotif_HH
