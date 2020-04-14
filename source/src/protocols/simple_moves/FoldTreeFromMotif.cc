// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/FoldTreeFromMotif.cc
/// @brief
/// @author longxing

// Unit headers
#include <protocols/simple_moves/FoldTreeFromMotif.hh>
#include <protocols/simple_moves/FoldTreeFromMotifCreator.hh>

// Project headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>

#include <core/select/util.hh>
#include <protocols/rosetta_scripts/util.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace simple_moves {

static basic::Tracer TR( "protocols.loops.FoldTreeFromMotif" );

FoldTreeFromMotif::FoldTreeFromMotif() :
	Mover( FoldTreeFromMotif::mover_name() )
{
}


FoldTreeFromMotif::~FoldTreeFromMotif() = default;

protocols::moves::MoverOP FoldTreeFromMotif::clone() const
{
	return utility::pointer::make_shared< FoldTreeFromMotif >( *this );
}

protocols::moves::MoverOP FoldTreeFromMotif::fresh_instance() const
{
	return utility::pointer::make_shared< FoldTreeFromMotif >();
}

void
FoldTreeFromMotif::apply( core::pose::Pose & pose )
{

	utility::vector1<bool> subset = motif_selector_->apply( pose );

	utility::vector1<core::Size> motif_residues = core::select::get_residues_from_subset( subset );
	if ( pose.num_chains() != 2 || motif_residues.size() == 0 ) {
		TR << "number of chains: " << pose.num_chains() << std::endl;
		TR << "total number of motif residues: " << motif_residues.size() << std::endl;
		utility_exit_with_message("APPLY ERROR!!");
	}
	core::Size middle_residue = motif_residues[core::Size( motif_residues.size()/2)+1];


	core::Size binder_len = pose.conformation().chain_end(1);
	core::kinematics::FoldTree ft;
	ft.add_edge(middle_residue, 1, -1);
	ft.add_edge(middle_residue, binder_len, -1);
	ft.add_edge(binder_len+1, pose.size(), -1);
	ft.add_edge(binder_len+1, binder_len, -1);
	ft.new_jump(middle_residue, binder_len+1, binder_len);
	ft.reorder(middle_residue);
	TR << ft << std::endl;
	TR << "Check fold tree: " << ft.check_fold_tree() << std::endl;
	pose.fold_tree( ft );

}

void
FoldTreeFromMotif::residue_selector(
	core::select::residue_selector::ResidueSelectorCOP const & selector_in
) {
	runtime_assert_string_msg( selector_in, "Error in protocols::simple_moves::FoldTreeFromMotif::residue_selector(): You must pass a valid residue selector." );
	motif_selector_ = selector_in->clone();
}


void
FoldTreeFromMotif::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data
) {

	residue_selector( protocols::rosetta_scripts::parse_residue_selector( tag, data ) );

}


std::string FoldTreeFromMotif::get_name() const {
	return mover_name();
}

std::string FoldTreeFromMotif::mover_name() {
	return "FoldTreeFromMotif";
}

void FoldTreeFromMotif::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector", "ResidueSelector that defines the motif region" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Helper mover that looks for motif definitions and sets up the fold tree.", attlist );
}

std::string FoldTreeFromMotifCreator::keyname() const {
	return FoldTreeFromMotif::mover_name();
}

protocols::moves::MoverOP
FoldTreeFromMotifCreator::create_mover() const {
	return utility::pointer::make_shared< FoldTreeFromMotif >();
}

void FoldTreeFromMotifCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FoldTreeFromMotif::provide_xml_schema( xsd );
}


} // namespace loops
} // namespace protocols
