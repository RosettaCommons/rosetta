// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/enzdes/RemoveLigandMover.cc
/// @brief
/// @author Javier Castellanos (javiercv@uw.edu)

// unit headers
#include <protocols/enzdes/RemoveLigandFilter.hh>
#include <protocols/enzdes/RemoveLigandFilterCreator.hh>

// package headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/types.hh>

#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/protein_interface_design/filters/RmsdFilter.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>

//utility headers
#include <utility/tag/Tag.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

using namespace core;
using namespace core::scoring;

static basic::Tracer TR( "protocols.enzdes.RemoveLigandFilter" );
namespace protocols {
namespace enzdes {

using namespace protocols::filters;
using namespace protocols::moves;
using namespace utility::tag;
using protocols::simple_moves::MinMover;
using protocols::protein_interface_design::filters::RmsdFilter;
using core::kinematics::MoveMapOP;
using core::pose::PoseOP;
using core::pose::Pose;
using core::Real;

// XRW TEMP std::string
// XRW TEMP RemoveLigandFilterCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return "RemoveLigandFilter";
// XRW TEMP }

// XRW TEMP FilterOP
// XRW TEMP RemoveLigandFilterCreator::create_filter() const
// XRW TEMP {
// XRW TEMP  return FilterOP( new RemoveLigandFilter );
// XRW TEMP }



RemoveLigandFilter::RemoveLigandFilter( ):
	Filter( "RemoveLigandFilter" ),
	threshold_( 99.99 ),
	mover_( MoverOP( new MinMover ) ),
	filter_( FilterOP( new RmsdFilter ) )
{
	MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_bb(true);
	movemap->set_chi(true);

	auto* min_mover = dynamic_cast< MinMover* >( mover_.get() );
	min_mover->movemap( movemap );
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	min_mover->score_function( scorefxn );
}


RemoveLigandFilter::RemoveLigandFilter( Real threshold ):
	Filter( "RemoveLigandFilter" ),
	threshold_( threshold ),
	mover_( MoverOP( new MinMover ) ),
	filter_( FilterOP( new RmsdFilter ) )
{ }

RemoveLigandFilter::RemoveLigandFilter( RemoveLigandFilter const & rval ):
	Filter( "RemoveLigandFilter" ),
	threshold_( rval.threshold_ ),
	mover_( rval.mover_ ),
	filter_( rval.filter_ )
{

}

Real
RemoveLigandFilter::report_sm( Pose const & pose ) const
{

	Pose no_lig_pose( pose );
	protocols::toolbox::pose_manipulation::remove_non_protein_residues( no_lig_pose );
	if ( filter_->get_type() == "Rmsd" ) {
		PoseOP init_pose( new Pose( no_lig_pose ) );
		auto* rmsd_filter = dynamic_cast< RmsdFilter* >( filter_.get() );
		rmsd_filter->reference_pose( init_pose );
		rmsd_filter->superimpose( true );
		using namespace core::select::residue_selector;
		rmsd_filter->set_selection( ResidueSelectorOP( new TrueResidueSelector ) ); // All residues

		mover_->apply( no_lig_pose );
		return rmsd_filter->report_sm( no_lig_pose );

	} else {
		Real start_score = filter_->report_sm( no_lig_pose );
		mover_->apply( no_lig_pose );
		Real end_score = filter_->report_sm( no_lig_pose );
		return end_score - start_score;
	}
}

bool
RemoveLigandFilter::apply(Pose const & pose) const
{
	return report_sm( pose ) < threshold_;
}

void
RemoveLigandFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, Filters_map const & filters, Movers_map const & movers, Pose const &)
{
	threshold_ = tag->getOption< Real >( "threshold", 3.0 );
	const std::string mover_name = tag->getOption< std::string >( "mover", "");
	const std::string filter_name = tag->getOption< std::string >( "filter", "");

	auto  find_mover ( movers.find( mover_name ));
	auto find_filter( filters.find( filter_name ));
	if ( find_mover == movers.end() && mover_name != "" ) {
		TR.Error << "mover not found in map: \n" << tag << std::endl;
		runtime_assert( find_mover != movers.end() );
	}
	if ( find_filter == filters.end() && filter_name != "" ) {
		TR.Error << "filter not found in map: \n" << tag << std::endl;
		runtime_assert( find_filter != filters.end() );
	}

	if ( mover_name != "" ) {
		mover_ = find_mover->second;
	}
	if ( filter_name != "" ) {
		filter_ = find_filter->second;
	}
}

std::string RemoveLigandFilter::name() const {
	return class_name();
}

std::string RemoveLigandFilter::class_name() {
	return "RemoveLigandFilter";
}

void RemoveLigandFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(
		"threshold", xsct_real,
		"XSD XRW: TO DO score threshold",
		"3.0");

	attlist + XMLSchemaAttribute(
		"mover", xs_string,
		"XSD XRW: TO DO");

	attlist + XMLSchemaAttribute(
		"filter", xs_string,
		"XSD XRW: TO DO rmsd type filter");

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Check if the ligand's pocket is stable by removing the ligand, "
		"relaxing the structure and calculating rms to the starting structure.",
		attlist );
}

std::string RemoveLigandFilterCreator::keyname() const {
	return RemoveLigandFilter::class_name();
}

protocols::filters::FilterOP
RemoveLigandFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new RemoveLigandFilter );
}

void RemoveLigandFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RemoveLigandFilter::provide_xml_schema( xsd );
}


} // enzdes
} // protocols
