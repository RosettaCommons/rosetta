// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/DdgFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#include <protocols/simple_filters/DdgFilter.hh>
#include <protocols/simple_filters/DdgFilterCreator.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/simple_moves/ddG.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/scoring/Interface.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/Energies.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <utility/excn/Exceptions.hh>
#include <ObjexxFCL/format.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.DdgFilter" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP DdgFilterCreator::create_filter() const { return protocols::filters::FilterOP( new DdgFilter ); }

// XRW TEMP std::string
// XRW TEMP DdgFilterCreator::keyname() const { return "Ddg"; }


DdgFilter::DdgFilter() :
	filters::Filter( "Ddg" ),
	ddg_threshold_( -15.0 ),
	scorefxn_( /* NULL */ ),
	rb_jump_( 1 ),
	use_custom_task_(false),
	repack_bound_(true),
	relax_bound_(false),
	repeats_( 1 ),
	repack_( true ),
	relax_mover_( /* NULL */ ),
	pb_enabled_(false),
	translate_by_(1000),
	extreme_value_removal_( false )
{
	scorename_ = "ddg";
}


DdgFilter::DdgFilter( core::Real const ddg_threshold,
	core::scoring::ScoreFunctionCOP scorefxn,
	core::Size const rb_jump/*=1*/,
	core::Size const repeats/*=1*/) :
	Filter("Ddg" ),
	ddg_threshold_(ddg_threshold),
	scorefxn_(scorefxn->clone()),
	rb_jump_(rb_jump),
	use_custom_task_( false ),
	repack_bound_( true ),
	relax_bound_( false ),
	repeats_(repeats),
	repack_( true ),
	relax_mover_( /* NULL */ ),
	pb_enabled_(false),
	translate_by_(1000),
	extreme_value_removal_( false )
{
	// Determine if this PB enabled.
	if ( scorefxn_->get_weight(core::scoring::PB_elec) != 0. ) {
		// Set this to PB enabled
		pb_enabled_ = true;
		TR << "PB enabled" << std::endl;
	} else {
		pb_enabled_ = false;
	}
}

DdgFilter::~DdgFilter() = default;

filters::FilterOP DdgFilter::clone() const {
	return filters::FilterOP( new DdgFilter( *this ) );
}
filters::FilterOP DdgFilter::fresh_instance() const{
	return filters::FilterOP( new DdgFilter() );
}

void
DdgFilter::repack( bool const repack )
{
	repack_ = repack;
}

bool
DdgFilter::repack() const
{
	return repack_;
}

void
DdgFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	filters::Filters_map const & filters,
	moves::Movers_map const & movers,
	core::pose::Pose const & )
{
	using namespace core::scoring;

	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	ddg_threshold_ = tag->getOption<core::Real>( "threshold", -15 );
	rb_jump_ = tag->getOption< core::Size >( "jump", 1 );
	repeats( tag->getOption< core::Size >( "repeats", 1 ) );
	repack( tag->getOption< bool >( "repack", 1 ) );
	if ( tag->hasOption( "symmetry" ) ) {
		TR << "DdgFilter autodetermines symmetry from input pose - symmetry option has no effect." << std::endl;
	}
	use_custom_task( tag->hasOption("task_operations") );
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	repack_bound( tag->getOption<bool>( "repack_bound", 1 ) );
	relax_bound( tag->getOption<bool>( "relax_bound", 0 ) );
	translate_by_ = tag->getOption<core::Real>( "translate_by", 1000 );

	if ( tag->hasOption( "relax_mover" ) ) {
		relax_mover( protocols::rosetta_scripts::parse_mover( tag->getOption< std::string >( "relax_mover" ), movers ) );
	}
	if ( tag->hasOption( "filter" ) ) {
		filter( protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "filter" ), filters ) );
	}

	if ( tag->hasOption("chain_num") ) {
		chain_ids_ = utility::string_split(tag->getOption<std::string>("chain_num"),',',core::Size());
	}
	extreme_value_removal( tag->getOption< bool >( "extreme_value_removal", false ) );
	if ( extreme_value_removal() ) {
		runtime_assert( repeats() >= 3 ); // otherwise makes no sense...
	}

	if ( repeats() > 1 && !repack() ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "ERROR: it doesn't make sense to have repeats if repack is false, since the values converge very well." );
	}

	TR
		<< "ddg filter with threshold " << ddg_threshold_
		<< " repeats=" << repeats()
		<< " and scorefxn " << rosetta_scripts::get_score_function_name(tag)
		<< " over jump " << rb_jump_
		<< "extreme_value_removal: " << extreme_value_removal()
		<< " and repack " << repack() << std::endl;

	// Determine if this PB enabled.
	if ( scorefxn_->get_weight(core::scoring::PB_elec) != 0. ) {
		// Set this to PB enabled
		pb_enabled_ = true;
		TR << "PB enabled.  Translation distance = " << translate_by_ << " A" << std::endl;
		if ( tag->hasOption("translate_by") && translate_by_ > 100 ) {
			TR.Warning << "Translation distance may be too large for PB-enabled scoring.  Consider 100 (default for PB enabled runs) if you run out of memory."<< std::endl;
			TR.Warning.flush();
		} else if ( !tag->hasOption("translate_by") ) {
			translate_by_ = 100;
			TR.Warning << "Translation distance set to 100 in order to save memory for the PB calculations."<<std::endl;
			TR.Warning.flush();
		}
	} else {
		pb_enabled_ = false;
	}
	TR.flush();
}

bool
DdgFilter::apply( core::pose::Pose const & pose ) const
{
	core::Real const pose_ddg( compute( pose ) );
	TR<<"ddg is "<<pose_ddg<<" ";
	if ( pose_ddg <= ddg_threshold_ ) {
		TR<<"passing"<<std::endl;
		return true;
	}
	TR<<"failing"<<std::endl;
	TR.flush();
	return false;
}

void
DdgFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const pose_ddg( compute( pose ) );
	out<<"ddg "<<pose_ddg<<'\n';
}

core::Real
DdgFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const pose_ddg( compute( pose ) );
	return( pose_ddg );
}

core::Size
DdgFilter::repeats() const
{
	return( repeats_ );
}

void
DdgFilter::repeats( core::Size const repeats )
{
	repeats_ = repeats;
}

core::Real
DdgFilter::translate_by() const
{
	return( translate_by_ );
}

void
DdgFilter::translate_by( core::Real const translate_by )
{
	translate_by_ = translate_by;
}

core::Real
DdgFilter::compute( core::pose::Pose const & pose_in ) const {
	core::pose::Pose pose(pose_in);
	if ( repack() ) {
		protocols::simple_moves::ddG ddg( scorefxn_, rb_jump_, chain_ids_ );
		if ( use_custom_task() ) {
			ddg.use_custom_task( use_custom_task() );
			ddg.task_factory( task_factory() );
		}
		if ( repack_bound() ) {
			ddg.repack_bound( repack_bound() );
		}
		if ( relax_bound() ) {
			ddg.relax_bound( relax_bound() );
		}
		ddg.translate_by( translate_by() );
		ddg.relax_mover( relax_mover() );
		ddg.filter( filter() );
		core::Real average( 0.0 );
		utility::vector1< core::Real > repeat_values;
		repeat_values.clear();
		for ( core::Size i = 1; i<=repeats_; ++i ) {
			ddg.calculate( pose );
			repeat_values.push_back( ddg.sum_ddG() );
			ddg.report_ddG( TR );
		}
		if ( extreme_value_removal() ) {
			runtime_assert( repeat_values.size() >= 3 );
			utility::vector1< core::Real > non_extreme_vals( repeat_values.size() - 2 );
			TR<<"removing extreme values. Considering values: ";
			std::sort( repeat_values.begin(), repeat_values.end() );
			std::copy( repeat_values.begin() + 1, repeat_values.end() - 1, non_extreme_vals.begin() );
			for ( core::Real const val : non_extreme_vals ) {
				average += val;
				TR<<val<<", ";
			}
			average /= (core::Real) non_extreme_vals.size();
			TR<<'\n'<<" average value: "<<average<<std::endl;
		} else { // fi extreme_value_removal
			for ( core::Real const val : repeat_values ) {
				average += val;
				TR<<val<<", ";
			}
			average /= (core::Real) repeat_values.size();
			TR<<'\n'<< " average value: "<< average<<std::endl;
		}

		return average;
	} else {
		if ( repeats() > 1 && !repack() ) {
			utility_exit_with_message( "ERROR: it doesn't make sense to have repeats if repack is false, since the values converge very well." );
		}
		using namespace protocols::moves;

		filters::FilterCOP scoring_filter;
		if ( filter_ ) {
			scoring_filter = filter_;
		} else {
			scoring_filter = filters::FilterCOP( filters::FilterOP( new simple_filters::ScoreTypeFilter( scorefxn_, core::scoring::total_score, 10000/*threshold*/ ) ) );
		}
		//JBB This is not handling symmetric poses properly. This needs to be fixed.
		core::pose::Pose split_pose( pose );
		if ( chain_ids_.size() > 0 ) {
			//We want to translate each chain the same direction, though it doesnt matter much which one
			core::Vector translation_axis(1,0,0);
			for ( core::Size current_chain_id : chain_ids_ ) {
				core::Size current_jump_id = core::pose::get_jump_id_from_chain_id(current_chain_id,split_pose);
				rigid::RigidBodyTransMoverOP translate( new rigid::RigidBodyTransMover( split_pose, current_jump_id) );
				translate->step_size( translate_by_ );
				translate->trans_axis(translation_axis);
				translate->apply( split_pose );
			}
		} else {
			rigid::RigidBodyTransMoverOP translate( new rigid::RigidBodyTransMover( split_pose, rb_jump_ ) );
			translate->step_size( translate_by_ );
			translate->apply( split_pose );
		}

		core::Real const bound_energy( scoring_filter->report_sm( pose ));
		core::Real const unbound_energy( scoring_filter->report_sm( split_pose ));
		core::Real const dG( bound_energy - unbound_energy );
		return( dG );
	}
}

void
DdgFilter::relax_mover( protocols::moves::MoverOP m ){
	relax_mover_ = m;
}

protocols::moves::MoverOP
DdgFilter::relax_mover() const{ return relax_mover_; }

void
DdgFilter::filter( protocols::filters::FilterOP f ){
	filter_ = f;
}

protocols::filters::FilterOP
DdgFilter::filter() const{ return filter_; }

std::string DdgFilter::name() const {
	return class_name();
}

std::string DdgFilter::class_name() {
	return "Ddg";
}

void DdgFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "threshold" , xsct_real , "Threshold below which the filter computes the binding energy." , "-15" )
		+ XMLSchemaAttribute::attribute_w_default( "jump" , xsct_positive_integer , "Specifies which chains to separate. Jump=1 would separate the chains interacting across the first chain termination, jump=2, second etc." , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "repeats" , xsct_positive_integer , "Averages the calculation over the number of repeats. Note that ddg calculations show noise of about 1-1.5 energy units, so averaging over 3-5 repeats is recommended for many applications." , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "repack" , xsct_rosetta_bool , "Should the complex be repacked in the bound and unbound states prior to taking the energy difference? If false, the filter turns to a dG evaluator. If repack=false repeats should be turned to 1, b/c the energy evaluations converge very well with repack=false." , "1" )
		+ XMLSchemaAttribute( "symmetry" , xs_string , "Note: DdgFilter autodetermines symmetry from input pose - symmetry option has no effect." ) //seems like this could be removed
		+ XMLSchemaAttribute::attribute_w_default( "repack_bound" , xsct_rosetta_bool , "Should the complex be repacked in the bound state? Note: If repack=true, then the complex will be repacked in the bound and unbound state by default. However, if the complex has already been repacked in the bound state prior to calling the DdgFilter then setting repack_bound=false allows one to avoid unnecessary repetition." , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "relax_bound" , xsct_rosetta_bool , "Should the relax mover (if specified) be applied to the bound as well as the unbound state? Note: the bound state is not relaxed by default." , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "translate_by" , xsct_real , "How far to translate the unbound pose. Note: Default is now 100 Angstroms rather than 1000." , "100" )
		+ XMLSchemaAttribute( "relax_mover" , xs_string , "Optionally define a mover which will be applied prior to computing the system energy in the unbound state." )
		+ XMLSchemaAttribute( "filter" , xs_string , "If specified, the given filter will be calculated in the bound and unbound state for the score, rather than the given scorefunction. Repacking, if any, will be done with the provided scorefunction." )
		+ XMLSchemaAttribute( "chain_num" , xs_string , "Allows you to specify a list of chain numbers to use to calculate the ddg, rather than a single jump. You cannot move chain 1, moving all the other chains is the same thing as moving chain 1, so do that instead. Use independently of jump." )
		+ XMLSchemaAttribute::attribute_w_default( "extreme_value_removal" , xsct_rosetta_bool , "Compute ddg value times, sort and remove the top and bottom evaluation. This should reduce the noise levels in trajectories involving 1000s of evaluations. If set to true, repeats must be set to at least 3." , "false" ) ;

	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist ) ;
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist ) ;
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Computes the binding energy for the complex and if it is below the threshold returns true. o/w false. Useful for identifying complexes that have poor binding energy and killing their trajectory.", attlist );
}

std::string DdgFilterCreator::keyname() const {
	return DdgFilter::class_name();
}

protocols::filters::FilterOP
DdgFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new DdgFilter );
}

void DdgFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DdgFilter::provide_xml_schema( xsd );
}

}

}
