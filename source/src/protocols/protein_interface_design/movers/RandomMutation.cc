// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/RandomMutation.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/RandomMutation.hh>
#include <protocols/protein_interface_design/movers/RandomMutationCreator.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
// Package headers
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>
#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/vector0.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace std;
using namespace core::scoring;

static basic::Tracer TR( "protocols.protein_interface_design.movers.RandomMutation" );

// XRW TEMP std::string
// XRW TEMP RandomMutationCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return RandomMutation::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP RandomMutationCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new RandomMutation );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP RandomMutation::mover_name()
// XRW TEMP {
// XRW TEMP  return "RandomMutation";
// XRW TEMP }

RandomMutation::RandomMutation() :
	Mover( RandomMutation::mover_name() ),
	task_factory_( /* NULL */ ),
	scorefxn_( /* NULL */ )
{
}


RandomMutation::~RandomMutation() = default;

void
RandomMutation::apply( core::pose::Pose & pose )
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::chemical;

	PackerTaskCOP task;
	if ( cache_task_ && task_ ) {
		if ( pose.size() == task_->total_residue() ) {
			task = task_;
		} else {
			task_.reset(); // Invalidate cached task.
		}
	}
	if ( ! task ) {
		task = task_factory()->create_task_and_apply_taskoperations( pose );
	}
	if ( cache_task_ && !task_ ) {
		task_ = task;
	}

	utility::vector1< core::Size > being_designed;
	being_designed.clear();

	for ( core::Size resi = 1; resi <= pose.size(); ++resi ) {
		if ( task->residue_task( resi ).being_designed() && pose.residue(resi).is_protein() ) {
			being_designed.push_back( resi );
		}
	}
	if ( being_designed.empty() ) {
		TR.Warning << "No residues are listed as designable." << std::endl;
		return;
	}
	core::Size const random_entry = being_designed[ (core::Size) floor( numeric::random::rg().uniform() * being_designed.size() )+1 ];
	using ResidueTypeCOPList = list<ResidueTypeCOP>;
	ResidueTypeCOPList const & allowed( task->residue_task( random_entry ).allowed_residue_types() );
	utility::vector1< AA > allow_temp;
	allow_temp.clear();
	for ( ResidueTypeCOP restype : allowed ) {
		if ( restype->aa() != pose.residue( random_entry ).aa() ) {
			allow_temp.push_back( restype->aa() );
		}
	}

	AA const target_aa( allow_temp[ (core::Size) floor( numeric::random::rg().uniform() * allow_temp.size() ) + 1 ] );
	utility::vector1< bool > allowed_aas;
	allowed_aas.clear();
	allowed_aas.assign( num_canonical_aas, false );
	allowed_aas[ target_aa ] = true;
	//PackerTaskOP mutate_residue = task_factory()->create_task_and_apply_taskoperations( pose );
	PackerTaskOP mutate_residue( task->clone() );
	mutate_residue->initialize_from_command_line().or_include_current( true );
	for ( core::Size resi = 1; resi <= pose.size(); ++resi ) {
		if ( resi != random_entry ) {
			mutate_residue->nonconst_residue_task( resi ).restrict_to_repacking();
		} else {
			mutate_residue->nonconst_residue_task( resi ).restrict_absent_canonical_aas( allowed_aas );
		}
	}
	TR<<"Mutating residue "<<pose.residue( random_entry ).name3()<<random_entry<<" to ";
	protocols::minimization_packing::PackRotamersMoverOP pack;
	pack = protocols::minimization_packing::PackRotamersMoverOP( new protocols::minimization_packing::PackRotamersMover( scorefxn(), mutate_residue ) );
	pack->apply( pose );
	TR<<pose.residue( random_entry ).name3()<<std::endl;
	(*scorefxn())(pose);
}

// XRW TEMP std::string
// XRW TEMP RandomMutation::get_name() const {
// XRW TEMP  return RandomMutation::mover_name();
// XRW TEMP }

void
RandomMutation::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	cache_task_ = tag->getOption< bool >( "cache_task", false );
}

protocols::moves::MoverOP
RandomMutation::clone() const {
	return( protocols::moves::MoverOP( new RandomMutation( *this ) ));
}

core::scoring::ScoreFunctionOP
RandomMutation::scorefxn() const{
	return scorefxn_;
}

void
RandomMutation::scorefxn( core::scoring::ScoreFunctionOP scorefxn )
{
	scorefxn_ = scorefxn;
}

core::pack::task::TaskFactoryOP
RandomMutation::task_factory() const{
	return( task_factory_ );
}

void
RandomMutation::task_factory( core::pack::task::TaskFactoryOP task_factory){
	task_factory_ = task_factory;
}

bool RandomMutation::cache_task() const {
	return( cache_task_ );
}

void RandomMutation::cache_task( bool cache ) {
	cache_task_ = cache;
}

std::string RandomMutation::get_name() const {
	return mover_name();
}

std::string RandomMutation::mover_name() {
	return "RandomMutation";
}

void RandomMutation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "cache_task", xsct_rosetta_bool, "Cache the packer task every time", "false" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string RandomMutationCreator::keyname() const {
	return RandomMutation::mover_name();
}

protocols::moves::MoverOP
RandomMutationCreator::create_mover() const {
	return protocols::moves::MoverOP( new RandomMutation );
}

void RandomMutationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RandomMutation::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
