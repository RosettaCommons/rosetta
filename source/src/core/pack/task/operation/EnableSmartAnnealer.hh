// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/operation/EnableSmartAnnealer.hh
/// @brief The smart annealer uses tensorflow to decrease the sample space of a typical packing run
/// @author Jack Maguire, jackmaguire1444@gmail.com


#ifndef INCLUDED_core_pack_task_operation_EnableSmartAnnealer_hh
#define INCLUDED_core_pack_task_operation_EnableSmartAnnealer_hh

#include <core/pack/task/operation/EnableSmartAnnealer.fwd.hh>

#include <core/pack/task/operation/TaskOperation.hh>

#include <core/types.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {

///@brief Task Operation for turning on the multi-cool annealer
class EnableSmartAnnealer: public TaskOperation {
public:

	EnableSmartAnnealer();

	EnableSmartAnnealer(EnableSmartAnnealer const & src);

	~EnableSmartAnnealer() override;

	TaskOperationOP clone() const override;

	/// @brief Configure from a RosettaScripts XML tag.
	void
	parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & ) override;

	//////////////////////

	void
	apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const override;

	/// @brief Return the name used to construct this TaskOperation from an XML file
	static std::string keyname();

	/// @brief Describe the format of XML file used to initialize this TaskOperation
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void
	set_smart_annealer_model( std::string const & smart_annealer_model ){
		smart_annealer_model_ = smart_annealer_model;
	}

	std::string const & smart_annealer_model() const {
		return smart_annealer_model_;
	}

	void
	set_smart_annealer_cutoff( core::Real const smart_annealer_cutoff ){
		smart_annealer_cutoff_ = smart_annealer_cutoff;
	}

	core::Real smart_annealer_cutoff() const {
		return smart_annealer_cutoff_;
	}

	void
	set_smart_annealer_pick_again( bool const smart_annealer_pick_again ){
		smart_annealer_pick_again_ = smart_annealer_pick_again;
	}

	bool smart_annealer_pick_again() const {
		return smart_annealer_pick_again_;
	}

	void
	set_smart_annealer_disable_during_quench( bool const smart_annealer_disable_during_quench ){
		smart_annealer_disable_during_quench_ = smart_annealer_disable_during_quench;
	}

	bool smart_annealer_disable_during_quench() const {
		return smart_annealer_disable_during_quench_;
	}

private:
	std::string smart_annealer_model_ = "generation2";
	core::Real smart_annealer_cutoff_ = 0.25;
	bool smart_annealer_pick_again_ = true;
	bool smart_annealer_disable_during_quench_ = true;

};

} //operation
} //task
} //pack
} //core

#endif
