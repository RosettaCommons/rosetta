// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file CCDLoopCloser.cc
///
/// @brief
/// @author Tim Jacobs

//unit
#include <devel/loop_creation/CCDLoopCloser.hh>
#include <devel/loop_creation/CCDLoopCloserCreator.hh>

//core
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Edge.hh>

//protocols
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/loop_closure/ccd/RamaCheck.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/pose_metric_calculators/ClashCountCalculator.hh>

//basic
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>

//numeric
#include <numeric/random/random.hh>

//utility
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace devel {
namespace loop_creation {

static basic::Tracer TR( "devel.loop_creation.CCDLoopCloser" );

//****CREATOR METHODS****//
// XRW TEMP std::string
// XRW TEMP CCDLoopCloserCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return CCDLoopCloser::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP CCDLoopCloserCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new CCDLoopCloser );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CCDLoopCloser::mover_name()
// XRW TEMP {
// XRW TEMP  return "CCDLoopCloser";
// XRW TEMP }
//****END CREATOR METHODS****//


/// @brief default constructor
CCDLoopCloser::CCDLoopCloser():
	max_closure_attempts_(10),
	max_ccd_moves_per_closure_attempt_(10000),
	tolerance_( 0.1 ),
	max_rama_score_increase_( 2.0 ),
	max_total_delta_helix_( 15 ),
	max_total_delta_strand_( 15 ),
	max_total_delta_loop_( 15 ),
	early_exit_cutoff_( 0.01 )
{
	init();
}

/// @brief explicit constructor
CCDLoopCloser::CCDLoopCloser(
	core::Size max_closure_attempts,
	bool prevent_nonloop_modifications,
	core::Size max_ccd_moves_per_closure_attempt,
	core::Real tolerance,
	core::Real max_rama_score_increase,
	core::Real max_total_delta_helix,
	core::Real max_total_delta_strand,
	core::Real max_total_delta_loop,
	core::Real early_exit_cutoff
):
	LoopCloser(prevent_nonloop_modifications),
	max_closure_attempts_(max_closure_attempts),
	max_ccd_moves_per_closure_attempt_(max_ccd_moves_per_closure_attempt),
	tolerance_(tolerance),
	max_rama_score_increase_(max_rama_score_increase),
	max_total_delta_helix_(max_total_delta_helix),
	max_total_delta_strand_(max_total_delta_strand),
	max_total_delta_loop_(max_total_delta_loop),
	early_exit_cutoff_(early_exit_cutoff)
{
	init();
}

protocols::moves::MoverOP
CCDLoopCloser::clone() const {
	return( protocols::moves::MoverOP( new CCDLoopCloser( *this ) ) );
}
protocols::moves::MoverOP
CCDLoopCloser::fresh_instance() const {
	return protocols::moves::MoverOP( new CCDLoopCloser );
}

// XRW TEMP std::string
// XRW TEMP CCDLoopCloser::get_name() const {
// XRW TEMP  return "CCDLoopCloser";
// XRW TEMP }

void
CCDLoopCloser::init(){
}

void
CCDLoopCloser::apply(
	core::pose::Pose & pose
){
	using namespace core;
	success_=false;

	core::kinematics::FoldTree saved_ft = pose.fold_tree();
	if ( prevent_nonloop_modifications() ) {
		//prepare_fold_tree(pose);
		protocols::loops::set_single_loop_fold_tree( pose, loop() );
	}


	protocols::loops::add_single_cutpoint_variant(pose, loop());

	// setup movemap
	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	for ( Size ii=loop().start(); ii<=loop().stop(); ++ii ) {
		mm->set_bb( ii, true );

		//don't change phi for prolines
		if ( pose.residue(ii).aa() == chemical::aa_pro ) {
			mm->set( id::TorsionID( id::phi_torsion, id::BB, ii ), false );
		}
	}
	TR.Debug << "CCD residues: " << loop().start() << " " << loop().stop() << std::endl;
	TR.Debug << "CCD cutpoint: " << loop().cut() << std::endl;

	//output variable for fast_ccd_closure
	core::Real deviation, torsion_delta, rama_delta;

	core::pose::Pose saved_pose = pose;
	protocols::loops::loop_closure::ccd::CCDLoopClosureMover ccd_loop_closure_mover( loop(), mm );
	ccd_loop_closure_mover.max_cycles( max_ccd_moves_per_closure_attempt_ );
	ccd_loop_closure_mover.tolerance( early_exit_cutoff_ );
	ccd_loop_closure_mover.rama()->max_rama_score_increase( max_rama_score_increase_ );
	ccd_loop_closure_mover.max_total_torsion_delta_per_residue(
		max_total_delta_helix_,
		max_total_delta_strand_,
		max_total_delta_loop_ );
	for ( core::Size i=1; i<=max_closure_attempts_; ++i ) {
		ccd_loop_closure_mover.apply( pose );

		deviation = ccd_loop_closure_mover.deviation();
		torsion_delta = ccd_loop_closure_mover.torsion_delta();
		rama_delta = ccd_loop_closure_mover.rama_delta();
		core::Size num_moves_needed = ccd_loop_closure_mover.actual_cycles();

		TR.Debug << "Loop closure attempt " << i << " returned in " << num_moves_needed << " iterations (max of "
			<< max_ccd_moves_per_closure_attempt_ << ")" << std::endl;

		TR.Debug << "torsion rmsd: " << torsion_delta << std::endl;
		TR.Debug << "rama delta: " << rama_delta << std::endl;
		TR.Debug << "ccd deviation: " << deviation << std::endl;

		if ( deviation < tolerance_ ) {
			//Calculator for backbone clash detection
			core::pose::metrics::PoseMetricCalculatorOP clash_calculator( new protocols::pose_metric_calculators::ClashCountCalculator(2.0) );
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( "clash_calculator", clash_calculator );

			basic::MetricValue<core::Size> bb_clash_metric;
			clash_calculator->get("bb", bb_clash_metric, pose);
			core::Size starting_clashes=bb_clash_metric.value();
			TR << "Number of BB clashes in input pose: " << starting_clashes << std::endl;

			basic::MetricValue<core::Size> loop_clash_metric;
			clash_calculator->get("bb", loop_clash_metric, pose);
			core::Size loop_clashes = loop_clash_metric.value() - starting_clashes;
			TR << "Number of BB clashes from loop: " << loop_clashes << std::endl;

			core::pose::metrics::CalculatorFactory::Instance().remove_calculator("clash_calculator");
			if ( loop_clashes == 0 ) {
				pose.fold_tree(saved_ft);
				TR << "Loop " << loop().start() << " " << loop().stop() << " successfully closed in " << i
					<< " attempts (max of " << max_closure_attempts_ << ")" << std::endl;
				success_=true;
				break;
			}
		}
		//restore pose to original and try closure again
		pose=saved_pose;
	}
	pose.fold_tree(saved_ft);
}

/// @brief Create a fold tree that prevents downstream propogation
/// of loop insertions
void
CCDLoopCloser::prepare_fold_tree(
	core::pose::Pose & pose
){
	using namespace core;

	//prepare special foldtree
	core::kinematics::FoldTree new_ft;
	new_ft.add_edge(1, loop().start()-1, kinematics::Edge::PEPTIDE);
	new_ft.add_edge(loop().start()-1, loop().cut(), kinematics::Edge::PEPTIDE);
	new_ft.add_edge(loop().start()-1, loop().stop()+1, 1);
	new_ft.add_edge(loop().stop()+1, loop().cut()+1, kinematics::Edge::PEPTIDE);
	new_ft.add_edge(loop().stop()+1, pose.size(), kinematics::Edge::PEPTIDE);

	pose.fold_tree(new_ft);
}

/// @brief parse tag for use in RosettaScripts
void
CCDLoopCloser::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
){
	using namespace core;

	if ( tag->hasOption("prevent_nonloop_modification") ) {
		prevent_nonloop_modifications(tag->getOption< bool >("prevent_nonloop_modifications"));
	}

	if ( tag->hasOption("max_ccd_moves_per_closure_attempt") ) {
		max_ccd_moves_per_closure_attempt_ = tag->getOption< Size >("max_ccd_moves_per_closure_attempt");
	}

	if ( tag->hasOption("tolerance") ) {
		tolerance_ = tag->getOption< Real >("tolerance");
	}

	if ( tag->hasOption("max_closure_attempts") ) {
		max_closure_attempts_ = tag->getOption< Size >("max_closure_attempts");
	}

	if ( tag->hasOption("max_rama_score_increase") ) {
		max_rama_score_increase_ = tag->getOption< Real >("max_rama_score_increase");
	}

	if ( tag->hasOption("max_total_delta_helix") ) {
		max_total_delta_helix_ = tag->getOption< Real >("max_total_delta_helix");
	}

	if ( tag->hasOption("max_total_delta_strand") ) {
		max_total_delta_strand_ = tag->getOption< Real >("max_total_delta_strand");
	}

	if ( tag->hasOption("max_total_delta_loop") ) {
		max_total_delta_loop_ = tag->getOption< Real >("max_total_delta_loop");
	}

	if ( tag->hasOption("early_exit_cutoff") ) {
		early_exit_cutoff_ = tag->getOption< Real >("early_exit_cutoff");
	}
}

std::string CCDLoopCloser::get_name() const {
	return mover_name();
}

std::string CCDLoopCloser::mover_name() {
	return "CCDLoopCloser";
}

void CCDLoopCloser::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "prevent_nonloop_modification", xsct_rosetta_bool, "XRW TO DO" )
		+ XMLSchemaAttribute( "max_ccd_moves_per_closure_attempt", xsct_non_negative_integer, "XRW TO DO" )
		+ XMLSchemaAttribute( "tolerance", xsct_real, "XRW TO DO" )
		+ XMLSchemaAttribute( "max_closure_attempts", xsct_non_negative_integer, "XRW TO DO" )
		+ XMLSchemaAttribute( "max_rama_score_increase", xsct_real, "XRW TO DO" )
		+ XMLSchemaAttribute( "max_total_delta_helix", xsct_real, "XRW TO DO" )
		+ XMLSchemaAttribute( "max_total_delta_strand", xsct_real, "XRW TO DO" )
		+ XMLSchemaAttribute( "max_total_delta_loop", xsct_real, "XRW TO DO" )
		+ XMLSchemaAttribute( "early_exit_cutoff", xsct_real, "XRW TO DO" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string CCDLoopCloserCreator::keyname() const {
	return CCDLoopCloser::mover_name();
}

protocols::moves::MoverOP
CCDLoopCloserCreator::create_mover() const {
	return protocols::moves::MoverOP( new CCDLoopCloser );
}

void CCDLoopCloserCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CCDLoopCloser::provide_xml_schema( xsd );
}


} //loop creation
} //devel
