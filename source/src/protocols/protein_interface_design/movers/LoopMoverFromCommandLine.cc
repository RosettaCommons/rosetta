// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/Loops/LoopMoverFromCommandLine.cc
/// @brief Parseable class to do full loop remodeling with input fragment files from command line
/// @author Jordan Willis (jordan.r.willis@vanderbilt.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/LoopMoverFromCommandLine.hh>
#include <protocols/protein_interface_design/movers/LoopMoverFromCommandLineCreator.hh>

// Package headers

// Project headers
#include <protocols/loops/loop_mover/perturb/LoopMover_CCD.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/loops/loop_mover/perturb/LoopMover_KIC.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/loops/loops_main.hh> // for various loop utility fxns
#include <protocols/loops/Loops.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <core/conformation/Conformation.hh>
#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/forge/methods/util.hh>

#include <protocols/protein_interface_design/movers/SaveAndRetrieveSidechains.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh> //getting SS from frag files
#include <core/scoring/ScoreFunction.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>

//create option keys for loop movers

#include <core/fragment/FragSet.hh>
#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

static THREAD_LOCAL basic::Tracer TR( "protocols.moves.LoopRemodelFromCommandLine" );
static THREAD_LOCAL basic::Tracer TR_report( "protocols.moves.LoopRemodelFromCommandLine.REPORT" );

// XRW TEMP std::string
// XRW TEMP LoopMoverFromCommandLineCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return LoopMoverFromCommandLine::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP LoopMoverFromCommandLineCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new LoopMoverFromCommandLine );
// XRW TEMP }

// XRW TEMP std::string LoopMoverFromCommandLine::mover_name()
// XRW TEMP {
// XRW TEMP  return "LoopMoverFromCommandLine";
// XRW TEMP }

LoopMoverFromCommandLine::~LoopMoverFromCommandLine() {}


protocols::moves::MoverOP
LoopMoverFromCommandLine::clone() const
{
	return( protocols::moves::MoverOP( new LoopMoverFromCommandLine( *this ) ) );
}


//call on empty constructor
LoopMoverFromCommandLine::LoopMoverFromCommandLine() :
	simple_moves::DesignRepackMover( LoopMoverFromCommandLine::mover_name() ),
	intermedrelax_( "no" ),
	remodel_( "no" ),
	relax_( "no" ),
	string_refine_( "no" )
{
	design(false);
}

//full member variables defined in constructor

LoopMoverFromCommandLine::LoopMoverFromCommandLine(
	std::string const & protocol,
	bool const perturb,
	bool const refine,
	core::scoring::ScoreFunctionOP & hires_score,
	core::scoring::ScoreFunctionOP & lores_score,
	std::string const & loop_file_name,
	protocols::loops::LoopsCOP loops
) :
	simple_moves::DesignRepackMover ( LoopMoverFromCommandLine::mover_name()),
	protocol_ ( protocol ),
	perturb_( perturb),
	refine_(refine),
	intermedrelax_( "no" ),
	remodel_( "no" ),
	relax_( "no" ),
	string_refine_( "no" )
{
	hires_score_ = hires_score;
	lores_score = lores_score->clone();
	loop_file_name_= loop_file_name;
	loops_ = protocols::loops::LoopsOP( new protocols::loops::Loops( *loops ) );
	design(false);
}


//apply to pose
void
LoopMoverFromCommandLine::apply ( core::pose::Pose & pose)
{
	using namespace protocols::loops;
	using core::pack::task::operation::TaskOperationCOP;
	pose.conformation().detect_disulfides(); // I don't think that this is important but just in case
	core::pose::Pose native_pose = pose;
	loops::set_secstruct_from_psipred_ss2( pose );
	LoopsOP loops( new protocols::loops::Loops( loop_file_name_ ) );
	loops->verify_against(pose);
	loops->auto_choose_cutpoints(pose);
	if ( loops->size() == 0 )  {
		TR << "No loops found!" << std::endl;
		return; // bounce out if we didn't define any loops
	} else {
		TR << *loops << std::endl;
	}
	if ( loops->size() > 0 ) {
		core::pack::task::TaskFactoryOP task_factory( new core::pack::task::TaskFactory );
		task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
		task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::IncludeCurrent ) );
		task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::NoRepackDisulfides ) );
		// set up temporary fold tree for loop closure
		TR.Debug << "Original FoldTree " << pose.fold_tree() << std::endl;
		core::kinematics::FoldTree old_ft( pose.fold_tree() );
		for ( Loops::iterator it = loops->v_begin(); it != loops->v_end(); ++it ) {
			it->set_extended( true ); // set all loops to extended (needed for CCD mover to really perturb)
			protocols::loops::LoopsOP single_loop( new protocols::loops::Loops() );
			single_loop->add_loop(*it);
			core::kinematics::FoldTree new_ft = protocols::forge::methods::fold_tree_from_loops( pose, *single_loop );
			pose.fold_tree( new_ft );
			add_cutpoint_variants( pose );
			core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
			//pose will always start full atom
			//protocols::moves::MonteCarlo mc( pose, *scorefxn_repack_, mc_kt );
			protocols::protein_interface_design::movers::SaveAndRetrieveSidechains retrieve_sc( pose );
			retrieve_sc.allsc( true );

			if ( protocol_ == "automatic" ) {
				utility::vector1< core::fragment::FragSetOP > frag_libs;
				protocols::loops::read_loop_fragments( frag_libs );

				protocols::comparative_modeling::LoopRelaxMover lrm;
				lrm.frag_libs( frag_libs );
				lrm.loops( single_loop );
				lrm.relax( relax() );
				lrm.refine( string_refine() );
				lrm.remodel( remodel() );
				lrm.intermedrelax( intermedrelax() );
				lrm.scorefxns( lores_score_, hires_score_ );
				lrm.apply( pose );
				return;
			}
			core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID_t );
			if ( protocol_ == "kinematic" ) {
				if ( perturb_ ) {
					protocols::loops::loop_mover::perturb::LoopMover_Perturb_KIC perturb(single_loop, lores_score_ );
					perturb.set_native_pose( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose ( native_pose ) ) ) );
					perturb.apply( pose );
				}
				core::util::switch_to_residue_type_set( pose, core::chemical::FULL_ATOM_t );
				retrieve_sc.apply( pose ); // recover sidechains from pre-centroid pose
				if ( refine_ ) {
					protocols::loops::loop_mover::refine::LoopMover_Refine_KIC refine( single_loop, hires_score_ );
					refine.set_redesign_loop(false); // design?
					refine.set_native_pose( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose ( native_pose ) ) ) );
					pose.update_residue_neighbors();
					refine.apply( pose );
				}
			} else if ( protocol_ == "ccd" ) { // protocol == kinematic
				TR << "Task Factory =" << task_factory;
				TR << "ccd protocol" << std::endl;
				pose.update_residue_neighbors();
				core::scoring::dssp::Dssp dssp( pose );
				dssp.insert_ss_into_pose( pose );
				//std::string const full_ss = pose.secstruct();
				//std::string const full_sequence = pose.sequence();
				utility::vector1< core::fragment::FragSetOP > frag_libs;
				protocols::loops::read_loop_fragments( frag_libs );
				if ( perturb_ ) {
					protocols::loops::loop_mover::perturb::LoopMover_Perturb_CCD perturb(single_loop, lores_score_ );
					for ( core::Size i = 1; i <= frag_libs.size(); ++i ) {
						perturb.add_fragments( frag_libs[i] );
					}
					perturb.set_strict_loops( true );
					perturb.set_native_pose( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose ( native_pose ) ) ) );
					perturb.apply( pose );
				}
				core::util::switch_to_residue_type_set( pose, core::chemical::FULL_ATOM_t );
				retrieve_sc.apply( pose ); // recover sidechains from pre-centroid pose
				if ( refine_ ) {
					protocols::loops::loop_mover::refine::LoopMover_Refine_CCD refine(single_loop, hires_score_ );
					for ( core::Size i = 1; i <= frag_libs.size(); ++i ) {
						refine.add_fragments( frag_libs[i] );
					}
					//core::pack::task::PackerTaskOP task = task_factory->create_task_and_apply_taskoperations( pose );
					refine.set_redesign_loop( false );
					refine.set_native_pose( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose ( native_pose ) ) ) );
					refine.apply( pose );
				}//refine
			}//ccd
		}//end single loop
	}//loops>0
}
// XRW TEMP std::string
// XRW TEMP LoopMoverFromCommandLine::get_name() const {
// XRW TEMP  return LoopMoverFromCommandLine::mover_name();
// XRW TEMP }
void
LoopMoverFromCommandLine::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	protocol_ = tag->getOption<std::string>( "protocol", "ccd" );
	perturb_ = tag->getOption<bool>( "perturb", 1 );
	if ( protocol_ == "automatic" ) { // ugly, but LoopRemodelMover accepts string whereas the other movers accept bool
		refine( tag->getOption< std::string >( "refine", "no" ) );
	} else {
		refine( tag->getOption<bool>( "refine", 1 ) );
	}
	intermedrelax( tag->getOption< std::string >( "intermedrelax", "no" ) );
	remodel( tag->getOption< std::string >( "remodel", "no" ) );
	relax( tag->getOption< std::string > ("relax", "no" ) );

	hires_score_ = protocols::rosetta_scripts::parse_score_function( tag, "refine_score", data )->clone();
	lores_score_ = protocols::rosetta_scripts::parse_score_function( tag, "perturb_score", data, "score4L" )->clone();

	loop_file_name_ = tag->getOption<std::string>("loop_file", "loops.loops");
	// task_factory(protocols::rosetta_scripts::parse_task_operations( tag, data ));

}//parsemytags

std::string LoopMoverFromCommandLine::get_name() const {
	return mover_name();
}

std::string LoopMoverFromCommandLine::mover_name() {
	return "LoopMoverFromCommandLine";
}

void LoopMoverFromCommandLine::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	XMLSchemaRestriction closuremethod;
	closuremethod.name( "closuremethod" );
	closuremethod.base_type( xs_string );
	closuremethod.add_restriction( xsr_enumeration, "ccd" );
	closuremethod.add_restriction( xsr_enumeration, "automatic" );
	closuremethod.add_restriction( xsr_enumeration, "kinematic" );
	xsd.add_top_level_element( closuremethod );

	attlist + XMLSchemaAttribute::attribute_w_default( "protocol", "closuremethod", "Which loop closure protocol to use?", "ccd" )
		+ XMLSchemaAttribute::attribute_w_default( "perturb", xsct_rosetta_bool, "Perturb before closure?", "1" )
		// In parse_my_tag, refine has two default values, contingent on the value of protocol_.
		// Since the default value must be universal, we have no default.
		+ XMLSchemaAttribute( "refine", xsct_rosetta_bool, "Use a refinement protocol?" )
		+ XMLSchemaAttribute::attribute_w_default( "intermedrelax", xsct_rosetta_bool, "Use an intermediate relax protocol?", "no" )
		+ XMLSchemaAttribute::attribute_w_default( "remodel", xsct_rosetta_bool, "Perform a remodeling protocol?", "no" )
		+ XMLSchemaAttribute::attribute_w_default( "relax", xsct_rosetta_bool, "A final round of fullatom relaxation?", "no" );

	rosetta_scripts::attributes_for_parse_score_function( attlist, "refine_score" );
	rosetta_scripts::attributes_for_parse_score_function( attlist, "perturb_score" );

	attlist + XMLSchemaAttribute::attribute_w_default( "loop_file", xs_string, "File from which loops definitions might be read", "loops.loops" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string LoopMoverFromCommandLineCreator::keyname() const {
	return LoopMoverFromCommandLine::mover_name();
}

protocols::moves::MoverOP
LoopMoverFromCommandLineCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoopMoverFromCommandLine );
}

void LoopMoverFromCommandLineCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoopMoverFromCommandLine::provide_xml_schema( xsd );
}

}//movers
}//protein_interface_design
}//protocols
