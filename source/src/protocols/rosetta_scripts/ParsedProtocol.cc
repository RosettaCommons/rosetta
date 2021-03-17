// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rosetta_scripts/ParsedProtocol.cc
/// @brief  A mover that applies a protocol parsed from a ROSETTASCRIPTS script
/// @author Sarel Fleishman (sarelf@u.washington.edu)
/// @author Luki Goldschmidt (lugo@uw.edu)
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- Modified this to facilitate use of ParsedProtocols to combine movers and filters in code outside of a RosettaScripts context.

// Unit Headers
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/rosetta_scripts/ParsedProtocolCreator.hh>

// Project Headers
#include <basic/datacache/DataMapObj.hh>
#include <basic/datacache/DataMap.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/simple_metrics/SimpleMetric.hh>
#include <core/simple_metrics/util.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/moves/ResId.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/NullMover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/filters/FilterFactory.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/rosetta_scripts/MultiplePoseMover.hh>

// Package Headers
#include <basic/Tracer.hh>
#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

// JD2 headers
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>

// Utility Headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

// C++ headers
#include <map>
#include <string>
#include <algorithm>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace rosetta_scripts {

static basic::Tracer TR( "protocols.rosetta_scripts.ParsedProtocol" );
static basic::Tracer TR_call_order( "protocols.rosetta_scripts.ParsedProtocol_call_order" );
static basic::Tracer TR_report( "protocols.rosetta_scripts.ParsedProtocol.REPORT" );

using Real = core::Real;
using Pose = core::pose::Pose;

using namespace core;
using namespace std;

ParsedProtocol::ParsedProtocol() :
	protocols::moves::Mover( "ParsedProtocol" ),
	final_scorefxn_( /* 0 */ ), // By default, don't rescore with any scorefunction.
	mode_("sequence"),
	last_attempted_mover_idx_( 0 ),
	report_call_order_( false ),
	last_mover_(/* NULL */),
	resume_support_(false)
{
}

ParsedProtocol::~ParsedProtocol() = default;

/// @details Takes care of the docking, design and filtering moves. pre_cycle and pose_cycle can
/// be setup in derived classes to setup variables before and after these cycles.
void
ParsedProtocol::apply( Pose & pose )
{
	n_steps_passed_in_previous_run_ = 0;

	ParsedProtocolAP this_weak_ptr(
		utility::pointer::dynamic_pointer_cast< ParsedProtocol >( get_self_ptr() )
	);

	if ( protocols::jd2::jd2_used() ) {
		protocols::jd2::JobDistributor::get_instance()->current_job()->add_output_observer( this_weak_ptr );
	}

	try {
		last_mover_=protocols::moves::MoverOP(); //Reset this

		protocols::moves::Mover::set_last_move_status( protocols::moves::FAIL_RETRY );
		//  pose.update_residue_neighbors();

		//fpd search the mover-filter pairs backwards for movers that have remaining poses
		ParsedProtocolStepVector::const_reverse_iterator rmover_it = steps_.rbegin();
		const ParsedProtocolStepVector::const_reverse_iterator movers_crend = steps_.rend();

		if ( resume_support_ && mode_ == "sequence" ) {
			for ( rmover_it=steps_.rbegin() ; rmover_it != movers_crend; ++rmover_it ) {
				if ( (*rmover_it).mover == nullptr ) { continue; }
				core::pose::PoseOP checkpoint = (*rmover_it).mover->get_additional_output();

				// otherwise continue where we left off
				if ( checkpoint ) {
					std::string const mover_name( rmover_it->mover->get_name() );

					// if mode_ is not 'sequence' then checkpointing is unsupported
					if ( checkpoint && mode_ != "sequence" ) {
						utility_exit_with_message("Mover "+mover_name+" returned multiple poses in a ParsedProtocol with mode!=sequence");
					}

					TR<<"=======================RESUMING FROM "<<mover_name<<"======================="<<std::endl;
					pose = *checkpoint;

					if ( ! apply_step( pose, *rmover_it, true) ) {
						//      final_score(pose);
						if ( protocols::jd2::jd2_used() ) {
							protocols::jd2::JobDistributor::get_instance()->current_job()->remove_output_observer( this_weak_ptr );
						}
						return;
					} else {
						break;
					}
				}
			}
		} else {
			rmover_it = steps_.rend();
		}

		if ( mode_ == "sequence" ) {
			sequence_protocol(pose, rmover_it.base());
		} else if ( mode_ =="random_order" ) {
			random_order_protocol(pose);
		} else if ( mode_ =="single_random" ) {
			random_single_protocol(pose);
		} else {
			TR.Warning << "mode is " << mode_ << " .This is not a valid ParsedProtocol Mode, your pose is being ignored" <<std::endl;
		}

		if ( get_last_move_status() == protocols::moves::MS_SUCCESS ) { // no point scoring a failed trajectory (and sometimes you get etable vs. pose atomset mismatches
			final_score(pose);
			runtime_assert( n_steps_passed_in_previous_run_ == steps_.size() || mode_ == "single_random" );
		} else {
			runtime_assert( n_steps_passed_in_previous_run_ < steps_.size() );
		}

		if ( protocols::jd2::jd2_used() ) {
			protocols::jd2::JobDistributor::get_instance()->current_job()->remove_output_observer( this_weak_ptr );
		}

	} catch( utility::excn::Exception & excn ) {
		TR.Error << "Exception while processing procotol: " << excn.msg() << std::endl;
		if ( protocols::jd2::jd2_used() ) {
			protocols::jd2::JobDistributor::get_instance()->current_job()->remove_output_observer( this_weak_ptr );
		}

		throw(excn);

	} catch( ... ) { /*To handle other exceptions*/
		TR.Error << "Exception while processing procotol:" << std::endl;
		if ( protocols::jd2::jd2_used() ) {
			protocols::jd2::JobDistributor::get_instance()->current_job()->remove_output_observer( this_weak_ptr );
		}

		throw;
	}
}


void ParsedProtocol::final_scorefxn( core::scoring::ScoreFunctionCOP scorefxn )
{
	final_scorefxn_ = scorefxn;
}

core::scoring::ScoreFunctionCOP ParsedProtocol::final_scorefxn() const
{
	return final_scorefxn_;
}

void
ParsedProtocol::final_score(core::pose::Pose & pose) const {
	core::scoring::ScoreFunctionCOP scorefxn( final_scorefxn() );

	if ( ! scorefxn ) { return; }

	(*scorefxn)(pose);
}

void
ParsedProtocol::report_all( Pose const & pose ) const {
	for ( ParsedProtocolStep const & step : steps_ ) {
		if ( step.report_filter_at_end_ && step.filter != nullptr ) {
			TR_report<<"============Begin report for "<<step.filter->get_user_defined_name()<<"=================="<<std::endl;
			step.filter->report( TR_report, pose );
			TR_report<<"============End report for "<<step.filter->get_user_defined_name()<<"=================="<<std::endl;
		}
	}
	TR_report.flush();
}

void
ParsedProtocol::add_values_to_job( Pose const & pose, protocols::jd2::Job & job ) const {
	for ( auto const & step : steps_ ) {
		if ( step.report_filter_at_end_ && step.filter != nullptr ) {
			core::Real const filter_value( step.filter->report_sm( pose ) );
			if ( filter_value > -9999 ) {
				job.add_string_real_pair(step.filter->get_user_defined_name(), filter_value);
			}
		}
	}
}

void
ParsedProtocol::report_filters_to_pose( Pose & pose ) const {
	for ( auto const & step : steps_ ) {
		if ( step.filter == nullptr ) { continue; }
		protocols::filters::Filter const & filter( *step.filter );
		if ( step.report_filter_at_end_ ) {
			core::Real const filter_value( filter.report_sm( pose ) );
			if ( filter_value > -9999 ) {
				setPoseExtraScore(pose, filter.get_user_defined_name(), (float)filter_value);
			}
		}
	}
}

ParsedProtocol::iterator
ParsedProtocol::begin(){
	return steps_.begin();
}

ParsedProtocol::const_iterator
ParsedProtocol::begin() const{
	return steps_.begin();
}

ParsedProtocol::iterator
ParsedProtocol::end(){
	return steps_.end();
}

ParsedProtocol::const_iterator
ParsedProtocol::end() const{
	return steps_.end();
}

/// @brief Add a mover-filter pair.
/// @details Indended for use OUTSIDE of a RosettaScripts context.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
ParsedProtocol::add_step(
	protocols::moves::MoverOP mover,
	std::string const &mover_name,
	protocols::filters::FilterOP filter,
	bool const report_filter_at_end
) {
	protocols::moves::MoverOP mover_to_add = nullptr;
	if ( mover ) mover_to_add = mover->clone();
	protocols::filters::FilterOP filter_to_add = nullptr;
	if ( filter ) filter_to_add = filter; //I don't know why the mover is cloned while the filter is not, but I'm keeping this consistent with the parse_my_tag function. VKM 17 Sept 2015.
	steps_.push_back( ParsedProtocolStep( mover_to_add, mover_name, filter_to_add, report_filter_at_end ) );
	return;
}

void
ParsedProtocol::add_step(
	ParsedProtocolStep const & step
) {
	steps_.push_back( step );
}

/// @details sets resid for the constituent filters and movers
void
ParsedProtocol::set_resid( core::Size const resid ){
	for ( auto & step : steps_ ) {
		using namespace protocols::moves;
		if ( step.mover != nullptr ) {
			modify_ResId_based_object( step.mover, resid );
		}
		if ( step.filter != nullptr ) {
			modify_ResId_based_object( step.filter, resid );
		}
	}
}

/// @details sets resid for the constituent filters and movers
void
ParsedProtocol::set_resid( core::pose::ResidueIndexDescriptionCOP r ){
	for ( auto & step : steps_ ) {
		using namespace protocols::moves;
		if ( step.mover != nullptr ) {
			modify_ResId_based_object( step.mover, r );
		}
		if ( step.filter != nullptr ) {
			modify_ResId_based_object( step.filter, r );
		}
	}
}

protocols::moves::MoverOP ParsedProtocol::clone() const
{
	return utility::pointer::make_shared< protocols::rosetta_scripts::ParsedProtocol >( *this );
}

/// @details May return nullptr for the Mover if one isn't present in the tag.
std::pair< moves::MoverOP, std::string >
parse_mover_subtag( utility::tag::TagCOP const tag_ptr,
	basic::datacache::DataMap& data
) {
	using namespace protocols::moves;

	MoverOP mover_to_add = nullptr;
	std::string mover_name; // user must specify a mover name. there is no valid default.

	runtime_assert( !( tag_ptr->hasOption("mover_name") && tag_ptr->hasOption("mover") ) );
	if ( tag_ptr->hasOption( "mover_name" ) ) {
		mover_name = tag_ptr->getOption<string>( "mover_name" );
		mover_to_add = protocols::rosetta_scripts::parse_mover_or_null( mover_name, data );
		if ( ! mover_to_add ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Mover " + mover_name + " not found in map");
		}
	} else if ( tag_ptr->hasOption( "mover" ) ) {
		mover_name = tag_ptr->getOption<string>( "mover" );
		mover_to_add = protocols::rosetta_scripts::parse_mover_or_null( mover_name, data );
		if ( ! mover_to_add ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Mover " + mover_name + " not found in map");
		}
	} else if ( tag_ptr->getName() != "Add" ) {
		MoverOP new_mover( MoverFactory::get_instance()->newMover( tag_ptr, data ) );
		debug_assert( new_mover );

		if ( tag_ptr->hasOption("name") ) {
			TR.Warning << "The mover named '" << tag_ptr->getOption<std::string>( "name" )
				<< " will not be accessible, since it was defined on-the-fly in the PROTOCOLS"
				<< " section." << std::endl;
		}

		TR.Info << "Defined on-the-fly '" << new_mover->name() << "' mover of type "
			<< tag_ptr->getName() << std::endl;
		mover_to_add = new_mover;
		mover_name = "OnTheFly("+tag_ptr->getName()+")";
	} else {
		mover_to_add = nullptr;
		mover_name = "NONE";
	}

	utility::tag::TagCOP tag_parent( tag_ptr->getParent() );
	if ( data.has( "stopping_condition", mover_name ) && tag_parent->hasOption( "name" ) ) {
		TR.Info << "ParsedProtocol's mover " << mover_name
			<< " requests its own stopping condition. This ParsedProtocol's stopping_condition will point at the mover's"
			<< std::endl;
		data.add( "stopping_condition", tag_parent->getOption< std::string >( "name" ), data.get_ptr< basic::datacache::DataMapObj< bool > >( "stopping_condition", mover_name ) );
	}

	return std::make_pair( mover_to_add, mover_name );
}

void
ParsedProtocol::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &data
)
{
	using namespace protocols::moves;
	using namespace utility::tag;

	TR<<"ParsedProtocol mover with the following settings" << std::endl;

	mode_=tag->getOption<string>("mode", "sequence");
	if ( mode_ != "sequence" && mode_ != "random_order" && mode_ != "single_random" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Error: mode must be sequence, random_order, or single_random");
	}

	utility::vector0< TagCOP > const & dd_tags( tag->getTags() );
	utility::vector1< core::Real > a_probability( dd_tags.size(), 1.0/dd_tags.size() );
	core::Size count( 1 );
	for ( auto const & tag_ptr: dd_tags ) {
		/////// Mover
		std::pair< MoverOP, std::string > mover_add_pair = parse_mover_subtag( tag_ptr, data );
		std::string const& mover_name( mover_add_pair.second );
		MoverOP mover_to_add( mover_add_pair.first );

		/////// Filter
		runtime_assert( !( tag_ptr->hasOption("filter_name") && tag_ptr->hasOption( "filter" ) ) );
		std::string filter_name;
		if ( tag_ptr->hasOption( "filter_name" ) ) {
			filter_name = tag_ptr->getOption<string>( "filter_name", "true_filter" );
		} else if ( tag_ptr->hasOption( "filter" ) ) {
			filter_name = tag_ptr->getOption<string>( "filter", "true_filter" );
		}

		protocols::filters::FilterOP filter_to_add;
		if ( ! filter_name.empty() ) {
			filter_to_add = protocols::rosetta_scripts::parse_filter_or_null( filter_name, data );
			if ( ! filter_to_add ) {
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Filter " + filter_name + " not found in map");
			}
		}

		////// Metrics
		utility::vector1< core::simple_metrics::SimpleMetricCOP > metrics_to_add;
		utility::vector1< std::string > metric_labels;
		if ( tag_ptr->hasOption( "metrics" ) ) {
			metrics_to_add = core::simple_metrics::get_metrics_from_datamap_and_subtags(tag_ptr, data);
			utility::vector1< std::string > metric_names = utility::string_split( tag_ptr->getOption<string>( "metrics" ), ',' );
			runtime_assert( metric_names.size() == metrics_to_add.size() );
			if ( tag_ptr->hasOption( "labels" ) ) {
				metric_labels = utility::string_split( tag_ptr->getOption<string>( "labels" ), ',' );
				if ( metric_labels.size() > metric_names.size() ) {
					TR.Error << "For metrics=\""<< tag_ptr->getOption<string>( "metrics" ) << "\" there are "
						<< metric_labels.size() << " labels and only " << metric_names.size() << " metrics." << std::endl;
					throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Too many labels for the number of metrics.");
				}
				TR.Debug << "Resizing metric label length from " << metric_labels.size() << " to " << metric_names.size() << std::endl;
				metric_labels.resize( metric_names.size() ); // Fill extra with empty
				for ( core::Size ii(1); ii <= metric_labels.size(); ++ii ) {
					if ( metric_labels[ii].empty() ) {
						TR.Debug << "Metric label " << ii << " is empty, replacing with " << metric_names[ii] << std::endl;
						metric_labels[ii] = metric_names[ii]; // Then use the names.
					}
				}
			} else {
				TR.Debug << "No metric labels specified, using metric names" << std::endl;
				metric_labels = metric_names;
			}
		}

		////// Report

		TR << "Added";
		if ( mover_to_add != nullptr ) {
			TR << " mover \"" << mover_name << "\"";
		}
		if ( filter_to_add != nullptr ) {
			TR << " filter \"" << filter_name << "\"";
		}
		if ( ! metric_labels.empty() ) {
			TR << " metrics:";
			for ( auto const & ml: metric_labels ) {
				TR << " \"" << ml << "\"";
			}
		}
		TR << std::endl;

		if ( mode_ == "single_random" ) {
			a_probability[ count ] = tag_ptr->getOption< core::Real >( "apply_probability", 1.0/dd_tags.size() );
			TR<<"and execution probability of "<<a_probability[ count ]<<std::endl;
		}
		count++;

		bool const report_at_end( tag_ptr->getOption< bool >( "report_at_end", true ) );
		if ( mover_to_add != nullptr ) {
			mover_to_add = mover_to_add->clone();
		}
		ParsedProtocolStep step( mover_to_add, mover_name, filter_to_add, report_at_end );
		step.metrics = metrics_to_add;
		step.metric_labels = metric_labels;
		steps_.push_back( step );
	}
	if ( mode_ == "single_random" ) {
		apply_probability( a_probability );
	}
	report_call_order( tag->getOption< bool >( "report_call_order", false ) );

	resume_support_ = tag->getOption<bool>("resume_support", false);
	if ( resume_support_ ) {
		TR << "Legacy protocol resume support enabled." << std::endl;
	}

	TR.flush();
}

/// @details Looks for any submovers that have additional output poses to process.
/// If any are found, run remainder of parsed protocol.
core::pose::PoseOP
ParsedProtocol::get_additional_output( )
{
	if ( !resume_support_ ) {
		// Get output from last specified mover
		const ParsedProtocolStepVector::const_reverse_iterator movers_crend = steps_.rend();
		for ( auto rmover_it=steps_.rbegin() ; rmover_it != movers_crend; ++rmover_it ) {
			protocols::moves::MoverOP mover = (*rmover_it).mover;
			if ( mover && mover->get_name() != "NullMover" ) {
				core::pose::PoseOP additional_pose = mover->get_additional_output();

				if ( additional_pose ) {
					final_score(*additional_pose);
				}

				return additional_pose;
			}
		}
		return nullptr;
	}

	// Legacy Protocol resume suppor
	core::pose::PoseOP pose = nullptr;
	auto rmover_it = steps_.rbegin();

	//fpd search the mover-filter pairs backwards; look for movers that have remaining poses
	for ( auto movers_crend = steps_.rend(); rmover_it != movers_crend; ++rmover_it ) {
		if ( (*rmover_it).mover == nullptr ) { continue; }
		core::pose::PoseOP checkpoint = (*rmover_it).mover->get_additional_output();
		if ( checkpoint ) {
			//std::string const mover_name( rmover_it->first.first->get_name() );
			pose = checkpoint;

			if ( ! apply_step( *pose, *rmover_it, true ) ) {
				return pose;
			} else {
				break;
			}
		}
	}

	// no saved poses?  return now
	if ( !pose ) return nullptr;

	// if mode_ is not 'sequence' then checkpointing is unsupported
	if ( mode_ != "sequence" ) {
		utility_exit_with_message("ParsedProtocol returned multiple poses in a ParsedProtocol with mode!=sequence");
	}

	// otherwise pick up from the checkpoint
	for ( auto step_it = rmover_it.base();
			step_it!=steps_.end(); ++step_it ) {
		if ( ! apply_step( *pose, *step_it ) ) {
			return pose;
		}
		// This used to not increment n_steps_passed_in_previous_run_
	}
	protocols::moves::Mover::set_last_move_status( protocols::moves::MS_SUCCESS ); // tell jobdistributor to save pose
	TR<<"setting status to success"<<std::endl;

	// report filter values to pose DataCache
	report_filters_to_pose( *pose );

	// report filter values to tracer output
	if ( TR_report.visible() ) {
		report_all( *pose );
	}

	final_score( *pose );


	return pose;
}

bool
ParsedProtocol::apply_step( Pose & pose, ParsedProtocolStep const & step, bool skip_mover /*=false*/ ) {
	if ( ! skip_mover ) {
		if ( ! apply_mover( pose, step ) ) {
			return false;
		}
	}
	if ( !apply_filter( pose, step ) ) {
		return false;
	}
	apply_metrics( pose, step );
	if ( ! skip_mover ) {
		++n_steps_passed_in_previous_run_;
	}
	return true;
}

bool
ParsedProtocol::apply_mover( Pose & pose, ParsedProtocolStep const & step ) {
	if ( step.mover == nullptr ) { return true; } // No mover no failure

	std::string const & mover_name( step.mover->get_name() );
	std::string const & mover_user_name( step.mover_user_name );

	// If the mover about to be applied is a MultiplePoserMover,
	// tell it about the previous mover where it should pull poses from
	TR.Debug << "apply_mover: Last mover: " << ( last_mover_ ? last_mover_->get_name() : "(None)" ) << std::endl;
	TR.Debug << "apply_mover: This mover: " << mover_name << std::endl;

	if ( last_mover_ ) {
		protocols::rosetta_scripts::MultiplePoseMoverOP mp_mover(
			utility::pointer::dynamic_pointer_cast< protocols::rosetta_scripts::MultiplePoseMover >( step.mover )
		);
		if ( mp_mover ) {
			mp_mover->set_previous_mover(last_mover_);
		}
	}

	TR<<"=======================BEGIN MOVER "<<mover_name<<" - "<<mover_user_name<<"======================="<<std::endl;
	step.mover->set_native_pose( get_native_pose() );
	step.mover->apply( pose );

	last_mover_ = step.mover;

	info().insert( info().end(), step.mover->info().begin(), step.mover->info().end() );

	moves::MoverStatus status( step.mover->get_last_move_status() );
	if ( status != protocols::moves::MS_SUCCESS ) {
		TR << "Mover " << step.mover->get_name() << " reports failure!" << std::endl;
		protocols::moves::Mover::set_last_move_status( status ); // Set status for this ParsedProtocol mover
		return false;
	}

	return true;
}

bool
ParsedProtocol::apply_filter( Pose & pose, ParsedProtocolStep const & step ) {
	if ( step.filter == nullptr ) { return true; } // No filter means we pass

	std::string const filter_name( step.filter->get_user_defined_name() );

	TR << "=======================BEGIN FILTER " << filter_name << "=======================" << std::endl;
	// Since filters get const poses, they don't necessarily have an opportunity to update neighbors themselves
	pose.update_residue_neighbors();
	bool pass = step.filter->apply( pose );
	if ( !step.report_filter_at_end_ ) { //report filter now
		core::Real const filter_value( step.filter->report_sm( pose ) );
		setPoseExtraScore(pose, step.filter->get_user_defined_name(), (float)filter_value);
		TR_report << "============Begin report for " << step.filter->get_user_defined_name() << "==================" << std::endl;
		step.filter->report( TR_report, pose );
		TR_report << "============End report for " << step.filter->get_user_defined_name() << "==================" << std::endl;
	}
	TR << "=======================END FILTER " << filter_name << "=======================" << std::endl;

	//filter failed -- set status in mover and return false
	if ( !pass ) {
		TR << "Filter " << filter_name << " reports failure!" << std::endl;
		protocols::moves::Mover::set_last_move_status( protocols::moves::FAIL_RETRY ); // Set status for this ParsedProtocol mover
		return false;
	}
	return true;
}

void
ParsedProtocol::apply_metrics( Pose & pose, ParsedProtocolStep const & step ) {
	for ( core::Size ii(1); ii <= step.metrics.size(); ++ii ) {
		if ( step.metrics[ii] == nullptr ) { continue; }

		std::string metric_label;
		if ( ii <= step.metric_labels.size() ) {
			metric_label = step.metric_labels[ii];
		}
		if ( metric_label.empty() || metric_label == "-" ) {
			metric_label = step.metrics[ii]->get_final_sm_type();
		}
		TR << "=======================BEGIN METRIC " << metric_label << "=======================" << std::endl;
		step.metrics[ii]->apply( metric_label, pose );
		TR << "=======================END METRIC " << metric_label << "=======================" << std::endl;
	}
}

/// @brief Provide the citation.
void
ParsedProtocol::provide_citation_info( basic::citation_manager::CitationCollectionList & citations ) const {
	using namespace basic::citation_manager;

	// Add citations for all movers and filters:
	for ( ParsedProtocolStep const & step : steps_ ) {
		citations.add( step.mover );
		citations.add( step.filter );
		for ( core::simple_metrics::SimpleMetricCOP const & metric: step.metrics ) {
			citations.add( metric );
		}
	}

	citations.add( last_mover_ );
}


void ParsedProtocol::finish_protocol(Pose & pose) {
	protocols::moves::Mover::set_last_move_status( protocols::moves::MS_SUCCESS ); // tell jobdistributor to save pose
	TR.Info << "setting status to success" << std::endl; // JRP/OFL that sounds like a debug statement to me...

	// report filter values to pose DataCache
	report_filters_to_pose( pose );
	// report filter values to tracer output
	if ( TR_report.visible() ) {
		report_all( pose );
	}

	if ( report_call_order() ) {
		std::string job_name ( protocols::jd2::jd2_used() ? protocols::jd2::current_output_name() : "CURRENT_JOB" );
		TR_call_order << job_name << " ";
		for ( ParsedProtocolStep const & p : steps_ ) {
			TR_call_order<<p.mover_user_name<<" ";
		}
		TR_call_order<<std::endl;
	}
}

void
ParsedProtocol::sequence_protocol( Pose & pose,
	ParsedProtocolStepVector::const_iterator step_it_in ) {

	for ( auto step_it = step_it_in;
			step_it!=steps_.end(); ++step_it ) {
		if ( !apply_step( pose, *step_it ) ) {
			return;
		}
	}
	// we're done! mark as success
	finish_protocol( pose );
}

void ParsedProtocol::random_order_protocol(Pose & pose){
	numeric::random::random_permutation(steps_.begin(),steps_.end(),numeric::random::rg());

	for ( ParsedProtocolStepVector::const_iterator it = steps_.begin(); it != steps_.end(); ++it ) {
		if ( ! apply_step( pose, *it ) ) {
			return;
		}
	}
	// we're done! mark as success
	finish_protocol( pose );
}

utility::vector1< core::Real >
ParsedProtocol::apply_probability() {
	core::Real sum( 0 );
	for ( core::Real const prob : apply_probability_ ) {
		sum += prob;
	}
	runtime_assert( sum >= 0.999 && sum <= 1.001 );
	return apply_probability_;
}

void
ParsedProtocol::apply_probability( utility::vector1< core::Real > const & a ){
	apply_probability_ = a;
	runtime_assert( apply_probability_.size() == steps_.size() );
	core::Real sum( 0 );
	for ( core::Real const prob : apply_probability_ ) {
		sum += prob;
	}
	runtime_assert( sum >= 0.999 && sum <= 1.001 );
}

void ParsedProtocol::random_single_protocol(Pose & pose){
	core::Real const random_num( numeric::random::rg().uniform() );
	core::Real sum( 0.0 );
	core::Size mover_index( 0 );
	for ( core::Real const probability : apply_probability() ) {
		sum += probability; mover_index++;
		if ( sum >= random_num ) {
			break;
		}
	}
	last_attempted_mover_idx( mover_index );
	if ( !apply_step( pose, steps_[ mover_index ] ) ) {
		return;
	}

	// we're done! mark as success
	finish_protocol( pose );
}

std::string ParsedProtocol::get_name() const {
	return mover_name();
}

std::string ParsedProtocol::mover_name() {
	return "ParsedProtocol";
}

std::string complex_type_name_for_parsed_protocol_subelement( std::string const & foo ) {
	return "parsed_protocol_subelement_" + foo + "_complex_type";
}

void ParsedProtocol::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	XMLSchemaRestriction mode_enum;
	mode_enum.name("mode_types");
	mode_enum.base_type(xs_string);
	mode_enum.add_restriction(xsr_enumeration, "sequence");
	mode_enum.add_restriction(xsr_enumeration, "random_order");
	mode_enum.add_restriction(xsr_enumeration, "single_random");
	xsd.add_top_level_element(mode_enum);

	attlist + XMLSchemaAttribute::attribute_w_default(
		"mode", "mode_types",
		"\"sequence\" (default) - perform the Mover/Filter pair in the specified sequence; "
		"\"random_order\" - perform EACH of the defined Mover/Filter pairs one time in a random order; "
		"\"single_random\" - randomly pick a SINGLE Mover/Filter pair from the list",
		"sequence");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"filter_name", xs_string,
		"XSD XRW: TO DO",
		"true_filter");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"filter", xs_string,
		"XSD XRW: TO DO",
		"true_filter");

	attlist + XMLSchemaAttribute(
		"apply_probability", xsct_real,
		"by default equal probability for all tags");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"resume_support", xsct_rosetta_bool,
		"XSD XRW: TO DO",
		"false");

	XMLSchemaSimpleSubelementList ssl;
	AttributeList add_subattlist;

	// Either mover_name or mover.
	add_subattlist + XMLSchemaAttribute( "mover_name", xs_string, "The mover whose execution is desired" )
		+ XMLSchemaAttribute( "mover", xs_string, "The mover whose execution is desired" );
	add_subattlist + XMLSchemaAttribute( "filter_name", xs_string, "The filter whose execution is desired" )
		+ XMLSchemaAttribute( "filter", xs_string, "The filter whose execution is desired" );
	add_subattlist + XMLSchemaAttribute( "metrics", xs_string, "A comma-separated list of metrics to run at this point." )
		+ XMLSchemaAttribute( "labels", xs_string, "A comma-separated list of labels to use for the provided metrics in the output. If empty/missing, use the metric names from the metrics setting. If '-', use the metric's default." );
	add_subattlist + XMLSchemaAttribute("apply_probability", xsct_real,"by default equal probability for all tags");
	add_subattlist + XMLSchemaAttribute::attribute_w_default(
		"report_at_end", xsct_rosetta_bool,
		"Report filter value via filter re-evaluation on final pose after "
		"conclusion of protocol. Otherwise report filter value as evaluated "
		"mid-protocol.",
		"true");


	ssl.add_simple_subelement( "Add", add_subattlist, "The steps to be applied."/*, 0 minoccurs*/ )
		.complex_type_naming_func( & complex_type_name_for_parsed_protocol_subelement );

	ssl.add_group_subelement( & protocols::filters::FilterFactory::filter_xml_schema_group_name );
	ssl.add_group_subelement( & protocols::moves::MoverFactory::mover_xml_schema_group_name );


	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, mover_name(),
		"This is a special mover that allows making a single compound mover and filter vector "
		"(just like protocols). The optional option mode changes the order of operations within "
		"the protocol, as defined by the option. If undefined, mode defaults to the historical "
		"functionality, which is operation of the Mover/Filter pairs in the defined order.",
		attlist,
		ssl );
}

std::string ParsedProtocolCreator::keyname() const {
	return ParsedProtocol::mover_name();
}

protocols::moves::MoverOP
ParsedProtocolCreator::create_mover() const {
	return utility::pointer::make_shared< ParsedProtocol >();
}

void ParsedProtocolCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ParsedProtocol::provide_xml_schema( xsd );
}


} //rosetta_scripts
} //protocols
