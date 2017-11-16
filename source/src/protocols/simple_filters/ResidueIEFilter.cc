// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ResidueIEFilter.cc
/// @brief
/// @author Sagar Khare (khares@u.washington.edu)

#include <protocols/simple_filters/ResidueIEFilter.hh>
#include <protocols/simple_filters/ResidueIEFilterCreator.hh>

#include <protocols/filters/Filter.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <utility/vector1.hh>
#include <core/pose/selection.hh>
#include <protocols/simple_filters/EnergyPerResidueFilter.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer tr( "protocols.simple_filters.ResidueIEFilter" );

ResidueIEFilter::ResidueIEFilter():
	filters::Filter( "ResidueIE" ),
	resnum_str_( "" ),
	restype_( "TRP" ),
	score_type_( core::scoring::total_score ),
	threshold_( 0.0 ),
	whole_pose_ ( false ),
	whole_interface_ ( false ),
	rb_jump_( 1 ),
	interface_distance_cutoff_( 8.0 ),
	max_penalty_( 1000.0 ),
	penalty_factor_( 1.0 ),
	report_energy_( false ),
	selector_()
{
}

ResidueIEFilter::ResidueIEFilter(
	std::string  resnums,
	std::string  restype,
	core::scoring::ScoreFunctionCOP scorefxn,
	core::scoring::ScoreType const score_type,
	core::Real const threshold,
	bool const whole_pose,
	bool const whole_interface,
	core::Size const rb_jump,
	core::Real const interface_distance_cutoff,
	core::Real max_penalty,
	core::Real penalty_factor,
	bool const report_energy
) :
	filters::Filter( "ResidueIE" ),
	resnum_str_(std::move( resnums )),
	restype_(std::move( restype )),
	score_type_( score_type ),
	threshold_( threshold ),
	whole_pose_( whole_pose ),
	whole_interface_( whole_interface ),
	rb_jump_( rb_jump ),
	interface_distance_cutoff_( interface_distance_cutoff ),
	max_penalty_( max_penalty ),
	penalty_factor_( penalty_factor ),
	report_energy_( report_energy ),
	selector_()
{
	using namespace core::scoring;

	if ( scorefxn ) scorefxn_ = scorefxn->clone();
	if ( score_type_ != total_score ) {
		core::Real const old_weight( scorefxn_->get_weight( score_type_ ) );
		scorefxn_->reset();
		scorefxn_->set_weight( score_type_, old_weight );

	}
	use_resE_ = false;
}

ResidueIEFilter::ResidueIEFilter( ResidueIEFilter const &init ) :
	Filter( init ), resnum_str_( init.resnum_str_ ),
	restype_( init.restype_ ),
	score_type_( init.score_type_ ),
	threshold_( init.threshold_ ),
	whole_pose_ (init.whole_pose_),
	whole_interface_ (init.whole_interface_),
	rb_jump_ (init.rb_jump_),
	interface_distance_cutoff_ ( init.interface_distance_cutoff_),
	max_penalty_ ( init.max_penalty_),
	penalty_factor_ ( init.penalty_factor_),
	report_energy_ ( init.report_energy_ ),
	use_resE_ ( init.use_resE_ ),
	selector_ ( init.selector_ )
{
	using namespace core::scoring;
	if ( init.scorefxn_ ) scorefxn_ = init.scorefxn_->clone();
}

core::scoring::ScoreFunctionOP
ResidueIEFilter::scorefxn() const{
	return scorefxn_;
}

void
ResidueIEFilter::scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
	scorefxn_ = scorefxn;
}

core::scoring::ScoreType
ResidueIEFilter::score_type() const{
	return score_type_;
}

void
ResidueIEFilter::score_type( core::scoring::ScoreType score_type ){
	score_type_ = score_type;
}

core::Real
ResidueIEFilter::threshold() const{
	return threshold_;
}

void
ResidueIEFilter::threshold( core::Real const th ){
	threshold_ = th;
}

std::string const &
ResidueIEFilter::resnums() const{
	return resnum_str_;
}

void
ResidueIEFilter::resnums( std::string const & resnums_str ){
	resnum_str_ = resnums_str;
}

ResidueIEFilter::~ResidueIEFilter() = default;

void
ResidueIEFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	using namespace core::scoring;
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	score_type_ = core::scoring::score_type_from_name( tag->getOption<std::string>( "score_type", "total_score" ) );
	threshold_ = tag->getOption<core::Real>( "energy_cutoff", threshold_ );
	whole_pose_ = tag->getOption<bool>( "whole_pose" , whole_pose_ );
	whole_interface_ = tag->getOption<bool>( "interface" , whole_interface_ );
	rb_jump_ = tag->getOption<core::Size>( "jump_number", rb_jump_ );
	interface_distance_cutoff_ = tag->getOption<core::Real>( "interface_distance_cutoff", interface_distance_cutoff_ );
	restype_ =  tag->getOption<std::string>( "restype3", restype_ );
	penalty_factor_ =  tag->getOption<core::Real>( "penalty_factor", penalty_factor_ );
	max_penalty_ =  tag->getOption<core::Real>( "max_penalty", max_penalty_ );
	report_energy_ = tag->getOption<bool>( "report_energy", report_energy_ );
	resnum_str_ = tag->getOption<std::string>( "resnums", resnum_str_ );

	if ( tag->hasOption("selector") ) {
		std::string const selector_name = tag->getOption< std::string >( "selector" );
		try {
			selector_ = data.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find ResidueSelector named '" << selector_name << "' from the Datamap from DisulfidizeMover.\n";
			error_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}
		debug_assert( selector_ );
		tr << "Using residue selector " << selector_name << std::endl;
	}

	runtime_assert( tag->hasOption("residues") || whole_pose_ || whole_interface_ || selector_ );

	use_resE_ = tag->getOption<bool>( "use_resE" , 0 );
}

bool
ResidueIEFilter::apply( core::pose::Pose const & pose ) const
{
	core::Real const penalty( compute( pose ) );
	if ( penalty>max_penalty_ ) return false;
	return true;
}


void
ResidueIEFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	core::Real const penalty( compute( pose ) );
	out << "Total penalty for restype "<< restype_ << "is "<< penalty << std::endl;
}


core::Real
ResidueIEFilter::report_sm( core::pose::Pose const & pose ) const
{

	core::Real const penalty( compute( pose ) );
	return( penalty );
}


core::Real
ResidueIEFilter::compute( core::pose::Pose const & pose ) const
{
	using namespace core::scoring;
	using namespace utility::graph;

	std::set< core::Size > const resnums = compute_resnums( pose );

	tr << "The following residues will be considered for interaction energy calculation:"<< std::endl;
	for ( core::Size const res : resnums ) {
		tr << pose.residue_type( res ).name3() << res <<" + ";
		if ( !(pose.residue_type( res ).name3() == restype_) ) {
			tr << "Residue " << res << " in pose is of type "<< pose.residue_type( res ).name3() << ". Requested restype3 is "<< restype_<<". Skipping!"<<std::endl;
			return false;
		}
	}

	core::pose::Pose in_pose = pose;
	(*scorefxn_)(in_pose);
	core::Real penalty (0.0);
	EnergyMap const curr_weights = in_pose.energies().weights();

	if ( resnums.size() == 0 ) {
		tr << "No residues found. Skipping calculation."<< std::endl;

		return (0.0);
	}

	for ( core::Size const res : resnums ) {

		core::Real res_intE (0.0);

		if ( use_resE_ ) {
			protocols::simple_filters::EnergyPerResidueFilter const eprf(res, scorefxn_, score_type_, 100000.0/*dummy threshold*/);
			res_intE = eprf.compute( pose );
		} else {
			core::Real res_intE_samechain (0.0);
			core::Real res_intE_differentchain (0.0);

			core::Size res_chain = in_pose.chain(res);

			//Fill residue energies by traversing energy graph
			for ( EdgeListConstIterator egraph_it = in_pose.energies().energy_graph().get_node( res )->const_edge_list_begin(); egraph_it != in_pose.energies().energy_graph().get_node( res )->const_edge_list_end(); ++egraph_it ) {
				EnergyEdge const * Eedge = static_cast< EnergyEdge const * > (*egraph_it);
				res_intE += Eedge->dot( curr_weights );


				if ( in_pose.chain(Eedge->get_other_ind(res)) == res_chain ) {
					res_intE_samechain += Eedge->dot( curr_weights );
				} else {
					res_intE_differentchain += Eedge->dot( curr_weights );
				}
			}//for each egraph_it
			tr << "Residue "<< pose.residue_type( res ).name3()<<res<< " has a intra-chain interaction energy "<< res_intE_samechain << std::endl;
			tr << "Residue "<< pose.residue_type( res ).name3()<<res<< " has a inter-chain interaction energy "<< res_intE_differentchain << std::endl;
		}

		tr << "Residue "<< pose.residue_type( res ).name3()<<res<< " has an (interaction) energy "<< res_intE <<", threshold is "<< threshold_<<" and penalty is ";

		core::Real const res_penalty = penalty_from_score( res_intE );
		penalty += res_penalty;
		tr << res_penalty << std::endl;
	} //foreach res

	return( penalty * penalty_factor_ );
}

core::Real
ResidueIEFilter::penalty_from_score( core::Real const res_intE ) const
{
	if ( report_energy_ ) return res_intE;
	if ( res_intE > threshold_ ) return (res_intE - threshold_);
	return 0.0;
}

std::string ResidueIEFilter::name() const {
	return class_name();
}

std::string ResidueIEFilter::class_name() {
	return "ResidueIE";
}

void ResidueIEFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::rosetta_scripts::attributes_for_parse_score_function(attlist);

	attlist + XMLSchemaAttribute::attribute_w_default(
		"score_type", xs_string,
		"XSD XRW: TO DO",
		"total_score");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"energy_cutoff", xsct_real,
		"XSD XRW: TO DO",
		"0.0");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"whole_pose", xsct_rosetta_bool,
		"XSD XRW: TO DO",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"interface", xsct_rosetta_bool,
		"XSD XRW: TO DO",
		"true");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"jump_number", xsct_non_negative_integer,
		"XSD XRW: TO DO",
		"1");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"interface_distance_cutoff", xsct_real,
		"XSD XRW: TO DO",
		"8.0");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"restype3", xs_string,
		"XSD XRW: TO DO",
		"TRP");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"penalty_factor", xsct_real,
		"XSD XRW: TO DO",
		"1.0");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"max_penalty", xsct_real,
		"XSD XRW: TO DO",
		"1000.0");

	attlist + XMLSchemaAttribute(
		"report_energy", xsct_rosetta_bool,
		"XSD XRW: TO DO");

	attlist + XMLSchemaAttribute(
		"resnums", xs_string,
		"XSD XRW: TO DO");

	attlist + XMLSchemaAttribute(
		"selector", xs_string,
		"XSD XRW: TO DO");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"use_resE", xsct_rosetta_bool,
		"XSD XRW: TO DO",
		"false");

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Finds the ineraction energy (IE) for the specified residue over the interface",
		attlist );
}

std::string ResidueIEFilterCreator::keyname() const {
	return ResidueIEFilter::class_name();
}

protocols::filters::FilterOP
ResidueIEFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ResidueIEFilter );
}

void ResidueIEFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueIEFilter::provide_xml_schema( xsd );
}


std::set< core::Size >
ResidueIEFilter::compute_resnums( core::pose::Pose const & pose ) const
{
	std::set< core::Size > resnums;

	if ( !resnum_str_.empty() ) {
		resnums = core::pose::get_resnum_list( resnum_str_, pose );
	}

	if ( selector_ ) {
		tr << "Applying residue selector to determine resnums" << std::endl;

		core::select::residue_selector::ResidueSubset subset = selector_->apply(pose);
		//sanity check
		debug_assert( subset.size() == pose.size() );
		for ( core::Size i=1, endi=pose.size(); i<=endi; ++i ) {
			if ( subset[i] ) {
				resnums.insert( i );
			}
		}
	}

	bool whole_pose = whole_pose_;
	bool whole_interface = whole_interface_;
	if ( resnums.empty() ) {
		whole_pose = true;
		whole_interface = false;
		tr << "Failed to parse residues: " << resnum_str_ << ". Using whole pose." << std::endl;
	} else {
		whole_pose = false;
		whole_interface = false;
	}

	if ( whole_interface ) {
		tr << "Detecting target resnums from interface." << std::endl;
		if ( pose.conformation().num_chains() < 2 ) {
			std::stringstream msg;
			msg << "ResidueIEFilter: pose must contain at least two chains! The given pose has " << pose.size() << " and " << pose.conformation().num_chains() << " chains" << std::endl;
			throw utility::excn::EXCN_BadInput( msg.str() );
		}

		core::pose::Pose in_pose = pose;

		(*scorefxn_)(in_pose);
		protocols::scoring::Interface interface_obj(rb_jump_);

		in_pose.update_residue_neighbors();

		interface_obj.distance( interface_distance_cutoff_ );
		interface_obj.calculate( in_pose );
		for ( core::Size resnum = 1; resnum <= pose.size(); ++resnum ) {
			if ( in_pose.residue(resnum).is_protein()  &&  interface_obj.is_interface( resnum ) && (in_pose.residue_type(resnum).name3() == restype_) ) {
				resnums.insert( resnum );
			}
		}
	} else if ( whole_pose ) { // whole_interface_
		tr << "Detecting target resnums from whole pose." << std::endl;

		core::pose::Pose in_pose = pose;
		for ( core::Size resnum = 1; resnum <= in_pose.size(); ++resnum ) {
			if ( in_pose.residue(resnum).is_protein()  && (in_pose.residue_type(resnum).name3() == restype_) ) resnums.insert( resnum );
		}
	}

	return resnums;
}

}
}
