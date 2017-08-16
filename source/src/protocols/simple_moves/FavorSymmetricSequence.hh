// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/FavorSymmetricSequence.hh
/// @brief apply constraints to enforce a symmetric sequence
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_simple_moves_FavorSymmetricSequence_HH
#define INCLUDED_protocols_simple_moves_FavorSymmetricSequence_HH

#include <protocols/simple_moves/FavorSymmetricSequence.fwd.hh>
#include <protocols/moves/Mover.hh>
namespace protocols {
namespace simple_moves {

class FavorSymmetricSequence : public moves::Mover {
public:
	FavorSymmetricSequence();
	FavorSymmetricSequence(core::Size symmetric_units, core::Real penalty);
	FavorSymmetricSequence(FavorSymmetricSequence const & src);
	~FavorSymmetricSequence() override = default;

	protocols::moves::MoverOP clone() const override;

	void apply(core::pose::Pose & pose) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::Size symmetric_units_;
	core::Real penalty_;
};

}
}


#endif /* FAVORSYMMETRICSEQUENCE_HH_ */
