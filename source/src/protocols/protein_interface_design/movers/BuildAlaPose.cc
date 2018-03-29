// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/BuildAlaPose.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/BuildAlaPose.hh>
#include <protocols/protein_interface_design/movers/BuildAlaPoseCreator.hh>

// Package headers
#include <protocols/rosetta_scripts/util.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/id/types.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/calc_taskop_movers/DesignRepackMover.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.BuildAlaPose" );

// XRW TEMP std::string
// XRW TEMP BuildAlaPoseCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return BuildAlaPose::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP BuildAlaPoseCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new BuildAlaPose );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP BuildAlaPose::mover_name()
// XRW TEMP {
// XRW TEMP  return "build_Ala_pose";
// XRW TEMP }

BuildAlaPose::BuildAlaPose() : calc_taskop_movers::DesignRepackMover( BuildAlaPose::mover_name() ),
	AA_("ALA")
{}

BuildAlaPose::BuildAlaPose(
	bool const partner1,
	bool const partner2,
	core::Real interface_distance_cutoff,
	std::string AA
) :
	calc_taskop_movers::DesignRepackMover( BuildAlaPose::mover_name() )
{
	repack_partner1_=design_partner1_=partner1;
	repack_partner2_=design_partner2_=partner2;
	interface_distance_cutoff_ = interface_distance_cutoff;
	AA_=AA;
	runtime_assert( interface_distance_cutoff_ >= 0 );
}

BuildAlaPose::~BuildAlaPose() = default;

void
BuildAlaPose::apply( pose::Pose & pose )
{
	using namespace core::scoring;
	using namespace core::pack::task::operation;

	allowed_aas_.assign( core::chemical::num_canonical_aas, false );
	allowed_aas_[ core::chemical::aa_from_name(AA_) ] = true;
	//allowed_aas_[ core::chemical::aa_ala ] = true;
	//allowed_aas_[ core::chemical::aa_gly ] = true;
	/* if( repack_partner1_ ^ repack_partner2_ ){
	bool const prevent_chain1( !repack_partner1_ );
	bool const prevent_chain2( !repack_partner2_ );
	core::Size const prevent_chain_begin( pose.conformation().chain_begin( prevent_chain1 ? 1 : 2 ) );
	core::Size const prevent_chain_end( pose.conformation().chain_end( prevent_chain2 ? 2 : 1 ) );
	for( core::Size res=prevent_chain_begin; res<=prevent_chain_end; ++res )
	prevent_repacking_.push_back( res );
	} */
	setup_packer_and_movemap( pose );

	using namespace core::scoring;
	core::scoring::ScoreFunctionOP scorefxn( get_score_function() );
	pack::pack_rotamers( pose, *scorefxn, task_ );
	(*scorefxn)( pose );
	/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( pose );
}

// XRW TEMP std::string
// XRW TEMP BuildAlaPose::get_name() const {
// XRW TEMP  return BuildAlaPose::mover_name();
// XRW TEMP }

void
BuildAlaPose::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	design_partner1_ = tag->getOption<bool>( "partner1", false );
	design_partner2_ = tag->getOption<bool>( "partner2", true );
	repack_partner1_ = design_partner1_;
	repack_partner2_ = design_partner2_;
	interface_distance_cutoff_ = tag->getOption<core::Real>( "interface_cutoff_distance", 20.0 );
	AA_ = tag->getOption<std::string>( "AA", "ALA" );

	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	TR<<"defined BuildAlaPose mover "<<" for partners "<<( repack_partner1_ ? "1" : "" )<<( repack_partner2_ ? "2": "" )<<" with distance cutoff "<< interface_distance_cutoff_ << " and convert to type " << core::chemical::aa_from_name(AA_) << " NOT WORK for GLY"<< std::endl;
}

protocols::moves::MoverOP
BuildAlaPose::clone() const {
	return( protocols::moves::MoverOP( new BuildAlaPose( *this ) ));
}

std::string BuildAlaPose::get_name() const {
	return mover_name();
}

std::string BuildAlaPose::mover_name() {
	return "build_Ala_pose";
}

void BuildAlaPose::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default( "partner1", xsct_rosetta_bool, "Design/repack the first chain", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "partner2", xsct_rosetta_bool, "Design/repack the second chain", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "interface_cutoff_distance", xsct_real, "Distance from the interface that counts for backrubbing", "8.0" )
		+ XMLSchemaAttribute::attribute_w_default( "AA", xs_string, "Amino acid, by name, from which to build the pose", "ALA" );

	rosetta_scripts::attributes_for_parse_task_operations(attlist);

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string BuildAlaPoseCreator::keyname() const {
	return BuildAlaPose::mover_name();
}

protocols::moves::MoverOP
BuildAlaPoseCreator::create_mover() const {
	return protocols::moves::MoverOP( new BuildAlaPose );
}

void BuildAlaPoseCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BuildAlaPose::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
