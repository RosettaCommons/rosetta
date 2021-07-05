// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file   protocols/ncbb/oop/OopCreatorMover.hh
///
/// @brief
/// @author Andrew Watkins


#ifndef INCLUDED_protocols_ncbb_oop_OopCreatorMover_hh
#define INCLUDED_protocols_ncbb_oop_OopCreatorMover_hh

#include <protocols/ncbb/oop/OopCreatorMover.fwd.hh>

#include <core/pose/Pose.fwd.hh>




#include <core/select/residue_selector/ResidueSelector.fwd.hh>


#include <utility/pointer/owning_ptr.hh>

#include <protocols/moves/Mover.hh>

namespace protocols {
namespace ncbb {
namespace oop {

class OopCreatorMover : public moves::Mover
{

public:

	//default ctor
	//OopCreatorMover(): Mover("OopCreatorMover"){}
	OopCreatorMover();
	OopCreatorMover(
		core::select::residue_selector::ResidueSelectorCOP  oop_plus_positions,
		core::select::residue_selector::ResidueSelectorCOP  oop_minus_positions,
		core::select::residue_selector::ResidueSelectorCOP  oop_d_plus_positions,
		core::select::residue_selector::ResidueSelectorCOP  oop_d_minus_positions,
		core::select::residue_selector::ResidueSelectorCOP  oop_low_e_puck_positions,
		core::Size prepend_n_residues,
		core::Size append_n_residues,
		bool final_repack,
		bool final_minimize,
		bool final_mc,
		bool final_correct_oop_post
	);

	//default dtor
	~OopCreatorMover() override{}

	//methods
	void apply( core::pose::Pose & pose ) override;
	protocols::moves::MoverOP fresh_instance() const override { return utility::pointer::make_shared< OopCreatorMover >(); }
	protocols::moves::MoverOP clone() const override;
	void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::select::residue_selector::ResidueSelectorCOP oop_plus_positions_;
	core::select::residue_selector::ResidueSelectorCOP oop_minus_positions_;
	core::select::residue_selector::ResidueSelectorCOP oop_d_plus_positions_;
	core::select::residue_selector::ResidueSelectorCOP oop_d_minus_positions_;
	core::select::residue_selector::ResidueSelectorCOP oop_low_e_puck_positions_;
	core::Size prepend_n_residues_;
	core::Size append_n_residues_;
	bool final_repack_;
	bool final_minimize_;
	bool final_mc_;
	bool final_correct_oop_post_;

};


} // namespace oop
} // namespace ncbb
} // namespace protocols

#endif // INCLUDED_protocols_ncbb_oop_OopCreator_hh
