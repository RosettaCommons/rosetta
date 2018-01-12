// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Ingemar Andre

#ifndef INCLUDED_protocols_minimization_packing_symmetry_SymPackRotamersMover_hh
#define INCLUDED_protocols_minimization_packing_symmetry_SymPackRotamersMover_hh

// Unit headers
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <core/pack/task/PackerTask.hh>

// Project headers
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.fwd.hh>

#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector0.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace minimization_packing {
namespace symmetry {

class SymPackRotamersMover : public protocols::minimization_packing::PackRotamersMover {

public:
	// default constructor
	SymPackRotamersMover();

	SymPackRotamersMover(
		core::scoring::ScoreFunctionCOP scorefxn,
		core::pack::task::PackerTaskCOP task = 0,
		core::Size nloop = 1
	);

	// destructor (important for properly forward-declaring smart-pointer members)
	~SymPackRotamersMover();

	// copy constructor
	SymPackRotamersMover( PackRotamersMover const & other );

	// virtual void apply( core::pose::Pose & pose );

	core::pack::task::PackerTaskOP
	make_symmetric_task(
		core::pose::Pose & pose,
		core::pack::task::PackerTaskOP task
	);

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		filters::Filters_map const &filters,
		moves::Movers_map const &movers,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	// to be used/redefined by derived classes
	void setup( core::pose::Pose & pose ) override;
	// need a more elegant rot_to_pack implementation than this
	core::PackerEnergy run(
		core::pose::Pose & pose,
		utility::vector0< int > rot_to_pack = utility::vector0<int>()
	) const override;

private:

	// pointers to data that are passed in
	core::pack::rotamer_set::symmetry::SymmetricRotamerSetsOP sym_rotamer_sets_;
	core::pack::task::PackerTaskOP symmetric_task_;
	AnnealableGraphBaseOP symmetric_ig_;
};

} // symmetry
} // moves
} // protocols

#endif
