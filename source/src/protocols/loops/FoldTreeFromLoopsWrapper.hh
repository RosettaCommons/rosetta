// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/FoldTreeFromLoops.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_loops_FoldTreeFromLoopsWrapper_HH
#define INCLUDED_protocols_loops_FoldTreeFromLoopsWrapper_HH

// Unit headers
#include <protocols/loops/FoldTreeFromLoopsWrapper.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <protocols/loops/Loops.fwd.hh>

namespace protocols {
namespace loops {

class FoldTreeFromLoops : public protocols::moves::Mover
{
public:
	FoldTreeFromLoops();
	~FoldTreeFromLoops() override;

	void apply( Pose & pose ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const &
	) override;

	void loop_str( std::string const & str );

	std::string loop_str() const;

	void loops( LoopsOP const l );

	LoopsOP loops() const;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	///@brief Add cutpoint variants that are used for scoring with chainbreak term?
	/// default false
	void
	add_cutpoint_variants( bool add_cp_variants );


private:
	std::string loop_str_; // loaded at parsetime but only realized at apply
	LoopsOP loops_; // a different interface into FoldTreeFromLoops, which takes precedence over loop_str_;
	bool add_cp_variants_ = false;
};

} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_loops_FoldTreeFromLoops_HH
