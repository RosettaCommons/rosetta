// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/protein_interface_design/filters/SpecificResiduesNearInterfaceFilterCreator.hh
/// @brief  Pass if a specific set of residues are near the interface
/// @author Arpit Tandon
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#include <protocols/protein_interface_design/filters/SpecificResiduesNearInterfaceFilter.hh>
#include <protocols/protein_interface_design/filters/SpecificResiduesNearInterfaceFilterCreator.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/scoring/Interface.hh>

#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/vector1.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

// C++ Headers
#include <utility/excn/Exceptions.hh>
#include <sstream>

namespace protocols {
namespace protein_interface_design {
namespace filters {

using std::endl;
using std::string;
using std::stringstream;
using std::ostream;
using core::pose::Pose;
using core::Real;
using core::Size;
using protocols::rosetta_scripts::parse_task_operations;
using protocols::scoring::Interface;
using protocols::filters::FilterOP;
using core::pack::task::TaskFactoryOP;
using core::pack::task::PackerTaskCOP;
using utility::vector1;

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.filters.SpecificResiduesNearInterfaceFilter" );

////////////  Creator ////////////////////////
FilterOP
SpecificResiduesNearInterfaceFilterCreator::create_filter() const {
	return FilterOP( new SpecificResiduesNearInterfaceFilter );
}

string
SpecificResiduesNearInterfaceFilterCreator::keyname() const {
	return "SpecificResiduesNearInterface";
}
/////////// End Creator /////////////////////


/// @brief default ctor
SpecificResiduesNearInterfaceFilter::SpecificResiduesNearInterfaceFilter() :
	parent( "SpecificResiduesNearInterfaceFilter" ),
	task_factory_(/* NULL */),
	rb_jump_(1)
{}

SpecificResiduesNearInterfaceFilter::SpecificResiduesNearInterfaceFilter(
	SpecificResiduesNearInterfaceFilter const & src) :
	parent( "SpecificResiduesNearInterfaceFilter" ),
	task_factory_(src.task_factory_),
	rb_jump_(src.rb_jump_)
{}

SpecificResiduesNearInterfaceFilter::~SpecificResiduesNearInterfaceFilter(){}

FilterOP
SpecificResiduesNearInterfaceFilter::fresh_instance() const{
	return FilterOP( new SpecificResiduesNearInterfaceFilter() );
}

FilterOP
SpecificResiduesNearInterfaceFilter::clone() const{
	return FilterOP( new SpecificResiduesNearInterfaceFilter( *this ) );
}


TaskFactoryOP
SpecificResiduesNearInterfaceFilter::task_factory() const{
	return task_factory_;
}

void
SpecificResiduesNearInterfaceFilter::task_factory(
	TaskFactoryOP tf
){
	task_factory_ = tf;
}

void
SpecificResiduesNearInterfaceFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & )
{
	if ( !tag->hasOption("task_operations" ) ) {
		throw utility::excn::EXCN_RosettaScriptsOption(
			"Please specify a task operation for the residues that should be near the other chain.");
	}
	task_factory(parse_task_operations( tag, data ));

	rb_jump_ = tag->getOption<core::Size>( "jump_number", 1 );

	// TODO make this work
	// interface_distance_threshold_ = tag->getOption<core::Size>( "interface_distance_threshold", 8 );
	//
	// TR
	//  << "Specific Residues Near Interface over jump number " << rb_jump_
	//  << " with threshold " << residues_in_interface_threshold_ << endl;
}


Size
SpecificResiduesNearInterfaceFilter::compute(
	Pose const & pose
) const {

	if ( !task_factory_ ) {
		utility_exit_with_message(
			"Please specify a task operation for the residues that should be near the other chain.");
	}

	if ( pose.conformation().num_chains() <= rb_jump_ ) {
		stringstream error_msg;
		error_msg
			<< "You have specied jump number '" << rb_jump_ << "', "
			<< "however the pose only has '" << pose.conformation().num_chains() << "' "
			<< "chains." << endl;
		utility_exit_with_message(error_msg.str());
	}

	Interface interface(rb_jump_);
	interface.calculate(pose);

	PackerTaskCOP task(task_factory_->create_task_and_apply_taskoperations(pose));
	vector1< bool > relevant_residues(task->repacking_residues());

	bool all_at_interface(true);
	for ( Size residue_number=1;
			residue_number <= pose.total_residue(); ++residue_number ) {
		if ( !relevant_residues[residue_number] ) continue;

		if ( !interface.is_interface(residue_number) ) {
			all_at_interface = false;
			break;
		}
	}

	return all_at_interface;
}

bool
SpecificResiduesNearInterfaceFilter::apply(
	Pose const & pose
) const {
	return compute(pose);
}


void
SpecificResiduesNearInterfaceFilter::report(
	ostream & out,
	Pose const & pose
) const {
	out
		<< "SpecificResiduesNearInterfaceFilter returns "
		<< (compute(pose) ? "true" : "false") << endl;
}


Real
SpecificResiduesNearInterfaceFilter::report_sm(
	Pose const & pose
) const {
	return(static_cast<Real>(compute(pose)));
}


} // namespace
} // namespace
} // namespace
