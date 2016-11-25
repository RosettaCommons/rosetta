// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/matdes/StoreQuasiSymmetricTaskMover.hh
/// @brief Headers for StoreQuasiSymmetricTaskMover class
/// @author Neil King (neilking@uw.edu)

#ifndef INCLUDED_devel_matdes_StoreQuasiSymmetricTaskMover_hh
#define INCLUDED_devel_matdes_StoreQuasiSymmetricTaskMover_hh

//unit headers
#include <devel/matdes/StoreQuasiSymmetricTaskMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/TaskFactory.hh>

namespace devel {
namespace matdes {

/// @brief mover that prepares a quasisymmetric task based upon a combination
/// of input tasks.
class StoreQuasiSymmetricTaskMover : public protocols::moves::Mover {

public:

	StoreQuasiSymmetricTaskMover();
	~StoreQuasiSymmetricTaskMover() override;

	void apply( core::pose::Pose & pose  ) override;
	// XRW TEMP  std::string get_name() const override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data_map,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	// setters
	void quasi_symm_comp( std::string const quasi_symm_comp );
	void num_quasi_repeats( core::Size const num_quasi_repeats );
	void offset_resis( core::Size const offset_resis );

	// getters
	char quasi_symm_comp() const;
	core::Size num_quasi_repeats() const;
	core::Size offset_resis() const;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::pack::task::TaskFactoryOP task_factory_;
	std::string task_name_;
	bool overwrite_;
	std::string quasi_symm_comp_;
	core::Size num_quasi_repeats_;
	core::Size offset_resis_;
};


} // matdes
} // devel

#endif

