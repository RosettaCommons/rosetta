// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/enzdes/EnzRepackMinimize.cc
/// @brief Repack/minimize class for rosetta_scripts-compatible enzdes
/// @author Sagar Khare (khares@uw.edu)

// Unit headers
#include <protocols/enzdes/EnzRepackMinimize.hh>
#include <protocols/enzdes/EnzRepackMinimizeCreator.hh>
#include <utility/string_util.hh>
//package headers
#include <protocols/enzdes/enzdes_util.hh>
#include <protocols/enzdes/EnzdesBaseProtocol.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>

// Project Headers
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/ScoreFunction.hh>

// AMW DEBUG
#include <core/conformation/Residue.hh>

#include <numeric/random/random.hh>

#include <basic/datacache/DataMap.hh>
#include <protocols/protein_interface_design/movers/BackrubDDMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>

#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <boost/functional/hash.hpp>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using namespace core;
using namespace core::scoring;

static basic::Tracer TR( "protocols.enzdes.EnzRepackMinimize" );


namespace protocols {
namespace enzdes {

// XRW TEMP std::string
// XRW TEMP EnzRepackMinimizeCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return EnzRepackMinimize::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP EnzRepackMinimizeCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new EnzRepackMinimize );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP EnzRepackMinimize::mover_name()
// XRW TEMP {
// XRW TEMP  return "EnzRepackMinimize";
// XRW TEMP }

EnzRepackMinimize::EnzRepackMinimize() :
	protocols::moves::Mover( "EnzRepackMinimize" ),
	cst_opt_( false ),
	design_( false ), repack_( false ), fix_catalytic_( false ),
	minimize_in_stages_( false ), min_rb_( true ), min_sc_( false ), min_bb_( false ),
	min_lig_( false ), minimize_prot_jumps_( false ), backrub_( false ),
	task_factory_(/* NULL */),
	n_cycles_( 1 )
{}

EnzRepackMinimize::EnzRepackMinimize(std::string const & name) :
	protocols::moves::Mover ( name ),
	cst_opt_( false ),
	design_( false ), repack_( false ), fix_catalytic_( false ),
	minimize_in_stages_( false ), min_rb_( true ), min_sc_( false ), min_bb_( false ),
	min_lig_( false ), minimize_prot_jumps_( false ), backrub_( false ),
	task_factory_(/* NULL */),
	n_cycles_( 1 )
{}

EnzRepackMinimize::~EnzRepackMinimize() = default;

protocols::moves::MoverOP
EnzRepackMinimize::clone() const
{
	return protocols::moves::MoverOP( new EnzRepackMinimize( *this ) );
}

protocols::moves::MoverOP
EnzRepackMinimize::fresh_instance() const
{
	return protocols::moves::MoverOP( new EnzRepackMinimize );
}

void
EnzRepackMinimize::apply( pose::Pose & pose )
{
	ensure_scofx_cstfication( pose );
	protocols::enzdes::EnzdesBaseProtocolOP enzprot( new protocols::enzdes::EnzdesBaseProtocol() );
	enzprot->set_fix_cataa( fix_catalytic_ );
	enzprot->set_minimize_options(min_sc_, min_bb_,min_rb_,min_lig_);
	enzprot->rb_min_jumps( rb_min_jumps() ); /// override the min_rb_ option, if multiple jumps were set

	for ( core::Size i=1; i<=n_cycles_; ++i ) {
		(*scorefxn_repack_)(pose);
		if ( cst_opt_ ) design_=true; // to enable poly-ala conversion
		if ( task_factory_ ==nullptr ) task_ = enzprot->create_enzdes_pack_task( pose, design_ );
		else task_ = create_ptask( pose );
		if ( cst_opt_ ) {
			TR<<"Starting PolyAla CstOptimization..."<<std::endl;
			if ( minimize_in_stages_ ) minimize_in_stages( pose, task_, true/*cst_opt*/, min_sc_, min_rb_, min_lig_ );
			enzprot->cst_minimize(pose, task_, true/*cst_opt*/);
		} else {
			if ( design_ || repack_ ) {
				TR<<"Starting Design/Repack..."<<std::endl;
				enzprot->enzdes_pack( pose, task_, scorefxn_repack_, 1/*cycles*/, false /*minimize_after_packing*/,false /*pack_unconstrained*/, false /*favor_native*/);
			}
			if ( min_sc_ || min_bb_ || min_lig_ || min_rb_ ) {
				TR<<"Starting Minimization..."<<std::endl;
				if ( minimize_in_stages_ ) minimize_in_stages( pose, task_, false/*cst_opt*/, min_sc_, min_rb_, min_lig_ );
				enzprot->set_scorefxn(scorefxn_minimize_);
				enzprot->set_minimize_options(min_sc_, min_bb_,min_rb_,min_lig_);
				enzprot->cst_minimize(pose, task_, false );
			}
		}
		if ( backrub_ ) {
			core::pose::Pose old_Pose = pose;
			design_ = false;
			if ( task_factory_ ==nullptr ) task_ = enzprot->create_enzdes_pack_task( pose, design_ );
			else task_ = create_ptask( pose );
			enzprot->set_minimize_options(min_sc_, false/*min_bb*/,min_rb_,min_lig_,true/*backrub*/);
			core::kinematics::MoveMapOP movemap = enzprot->create_enzdes_movemap( pose, task_, minimize_prot_jumps_  );
			core::scoring::ScoreFunctionCOP br_scorefxn = scorefxn_minimize_;
			utility::vector1<core::Size> residues;
			for ( core::Size i =1; i<pose.size(); ++i ) {
				if ( movemap->get_bb(i) ) residues.push_back(i);
			}
			TR<<"Now Backrub minimizing:  min_sc "<<min_sc_<<" min_bb "<< min_bb_<<std::endl;
			protocols::protein_interface_design::movers::BackrubDDMover br( br_scorefxn, true/*br_partner1?*/, false/*brpartner2=lig*/, 7.5/*interface distance*/, 1000/*num moves*/, 0.6/*kT*/,0.25/*sidechain move probability*/, residues/*vector of residues to backrub*/ ); //default params; residues vector will over-ride other interface detection stuff specified by other params
			br.apply(pose);
			pose.constraint_set( old_Pose.constraint_set()->clone() );
			pose.fold_tree( old_Pose.fold_tree() );
			pose.update_residue_neighbors();
		}
		pose.update_residue_neighbors();
		(*scorefxn_minimize_)(pose);
		TR<<"Finished Cyle#"<< i <<" of EnzRepackMinimize"<<std::endl;
	}
	TR<<"Finished EnzRepackMinimize"<<std::endl;
}

core::pack::task::PackerTaskOP
EnzRepackMinimize::create_ptask( core::pose::Pose & pose )
{

	using namespace core::pack::task;
	TR<<"Creating packer task based on specified task operations..."<< std::endl;
	task_factory_->push_back( operation::TaskOperationCOP( new operation::InitializeFromCommandline ) );
	PackerTaskOP task = task_factory_->create_task_and_apply_taskoperations( pose );
	return task;
}

void
EnzRepackMinimize::minimize_in_stages(
	core::pose::Pose & pose,
	core::pack::task::PackerTaskCOP task,
	bool const & cstopt,
	bool const & min_sc,
	bool const & min_rb,
	bool const & min_lig
)
{
	protocols::enzdes::EnzdesBaseProtocolOP enzprot( new protocols::enzdes::EnzdesBaseProtocol() );
	enzprot->set_scorefxn( scorefxn_minimize() );

	enzprot->set_minimize_options( min_sc, false/*min_bb_*/,min_rb,min_lig );
	enzprot->cst_minimize(pose, task, cstopt);

	TR<<"Finished non-bb dof minimization"<<std::endl;

	enzprot->set_minimize_options(false/*min_sc_*/, true/*min_bb_*/,min_rb_,false/*min_lig_*/);
	enzprot->cst_minimize(pose, task, cstopt);

	TR<<"Finished bb dof minimization"<<std::endl;

}

void
EnzRepackMinimize::ensure_scofx_cstfication(core::pose::Pose const & pose )
{
	toolbox::match_enzdes_util::EnzdesCacheableObserverCOP enzobs(toolbox::match_enzdes_util::get_enzdes_observer( pose ) );
	if ( !enzobs ) return;

	if ( enzobs->cst_cache() ) { // Make sure scorefunction has cst terms if constraints (seem to ) exist
		if ( !( enzutil::is_scofx_cstfied(scorefxn_repack_) && enzutil::is_scofx_cstfied(scorefxn_minimize_) ) ) {
			TR<<"Setting up Scorefunction to include constraints..."<<std::endl;
			enzutil::enable_constraint_scoreterms(scorefxn_minimize_);
			enzutil::enable_constraint_scoreterms(scorefxn_repack_);
			enzutil::scorefxn_update_from_options(scorefxn_minimize_);
			enzutil::scorefxn_update_from_options(scorefxn_repack_);
		}
	}
}

void
EnzRepackMinimize::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ){

	if ( tag->hasOption("task_operations") ) task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	else task_factory_ = nullptr;

	n_cycles_ = tag->getOption<core::Size>( "cycles", 1 );

	minimize_in_stages_ = tag->getOption<bool>( "min_in_stages", 0 );
	design_ = tag->getOption<bool>( "design", 0 );
	repack_ = tag->getOption<bool>( "repack_only", 0 );
	fix_catalytic_ = tag->getOption<bool>( "fix_catalytic", 0 );
	cst_opt_ = tag->getOption<bool>( "cst_opt", 0 );
	backrub_ = tag->getOption<bool>( "backrub", 0 );

	if ( tag->hasOption( "minimize_rb" ) ) {
		set_min_rb( tag->getOption<bool>( "minimize_rb", 1 ) );
	}
	if ( tag->hasOption( "minimize_bb" ) ) {
		set_min_bb(  tag->getOption<bool>( "minimize_bb", 0 ) );
	}
	if ( tag->hasOption( "minimize_sc" ) ) {
		set_min_sc(  tag->getOption<bool>( "minimize_sc", 1 ) );
	}
	if ( tag->hasOption( "minimize_lig" ) ) {
		set_min_lig(  tag->getOption<bool>( "minimize_lig", 0 ) );
	}
	if ( tag->hasOption( "minimize_prot_jumps" ) ) {
		minimize_prot_jumps_ = tag->getOption<bool>( "minimize_prot_jumps", 0 );
	}

	scorefxn_repack_ = protocols::rosetta_scripts::parse_score_function( tag, "scorefxn_repack", data )->clone();
	scorefxn_minimize_ = protocols::rosetta_scripts::parse_score_function( tag, "scorefxn_minimize", data )->clone();

	rb_min_jumps( utility::string_split< core::Size >( tag->getOption< std::string >( "rb_min_jumps", "" ), ',', core::Size() ) );
	if ( rb_min_jumps().size() > 0 ) {
		TR<<"rb_min_jumps set to "<<tag->getOption< std::string >( "rb_min_jumps" )<<" superseding the minimize_rb option"<<std::endl;
	}
	if ( design_ && repack_ ) utility_exit_with_message("Can't both Design and Repack_Only. Check xml file");
	if ( minimize_in_stages_ && (!min_bb_) ) utility_exit_with_message( "EnzRepackMinimize cant minimize in stages without minimize_bb set to 1. Check xml file." );

	//task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, data );

	TR
		<< "design="<<design_<< ","
		<< " with repack scorefxn " << protocols::rosetta_scripts::get_score_function_name(tag, "scorefxn_repack" )
		<< " and minimize scorefxn " << protocols::rosetta_scripts::get_score_function_name(tag, "scorefxn_minimize" )
		<< std::endl;

}

void
EnzRepackMinimize::set_scorefxn_repack( core::scoring::ScoreFunctionCOP scorefxn ) {
	scorefxn_repack_ = scorefxn->clone();
}

void
EnzRepackMinimize::set_scorefxn_minimize( core::scoring::ScoreFunctionCOP scorefxn ) {
	scorefxn_minimize_ = scorefxn->clone();
}

core::scoring::ScoreFunctionOP
EnzRepackMinimize::scorefxn_repack() const {
	return scorefxn_repack_;
}

core::scoring::ScoreFunctionOP
EnzRepackMinimize::scorefxn_minimize() const {
	return scorefxn_minimize_;
}

void
EnzRepackMinimize::task_factory( core::pack::task::TaskFactoryOP p ) {
	task_factory_ = p;
}

utility::vector1< core::Size >
EnzRepackMinimize::rb_min_jumps() const{ return rb_min_jumps_; }

void
EnzRepackMinimize::rb_min_jumps( utility::vector1< core::Size > const v ){
	rb_min_jumps_ = v; }

std::string EnzRepackMinimize::get_name() const {
	return mover_name();
}

std::string EnzRepackMinimize::mover_name() {
	return "EnzRepackMinimize";
}

void EnzRepackMinimize::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	rosetta_scripts::attributes_for_parse_task_operations_w_factory(attlist);

	attlist + XMLSchemaAttribute::attribute_w_default(
		"cycles", xsct_non_negative_integer,
		"number of cycles of repack-minimize (default=1 cycle) "
		"(Note: In contrast to the enzyme_design application, all cycles use the provided scorefunction.)",
		"1");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"min_in_stages", xsct_rosetta_bool,
		"first minimize non-backbone dofs, followed by backbone dofs only, "
		"and then everything together (default=0)",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"design", xsct_rosetta_bool,
		"optimize sequence of residues spatially around the ligand "
		"(detection of neighbors need to be specified in the flagfile or "
		"resfile, default=0)",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"repack_only", xsct_rosetta_bool,
		"if true, only repack sidechains without changing sequence. "
		"(default =0) If both design and repack_only are false, don't repack at all, only minimize",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"fix_catalytic", xsct_rosetta_bool,
		"fix catalytic residues during repack/minimization (default =0)",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"cst_opt", xsct_rosetta_bool,
		"perform minimization of enzdes constraints with a reduced scorefunction and in a polyAla background. (default= 0)",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"backrub", xsct_rosetta_bool,
		"use backrub to minimize (default=0).",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"minimize_rb", xsct_rosetta_bool,
		"minimize rigid body orientation of ligand (default=1)",
		"true");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"minimize_bb", xsct_rosetta_bool,
		"minimize back bone conformation of backbone segments that surround the ligand "
		"(contiguous neighbor segments of more than 3 residues are automatically chosen, default=0)",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"minimize_sc", xsct_rosetta_bool,
		"minimize sidechains (default=1)",
		"true");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"minimize_lig", xsct_rosetta_bool,
		"minimize ligand internal torsion degrees of freedom (allowed deviation needs to be specified by flag, default =0)",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"minimize_prot_jumps", xsct_rosetta_bool,
		"Minimize all jumps",
		"false");

	rosetta_scripts::attributes_for_parse_score_function(attlist, "scorefxn_repack");
	rosetta_scripts::attributes_for_parse_score_function(attlist, "scorefxn_minimize");

	attlist + XMLSchemaAttribute(
		"rb_min_jumps", xsct_nnegative_int_cslist,
		"specify which jumps to minimize. If this is specified it takes precedence "
		"over minimize_rb above. Useful if you have more than one ligand in the system "
		"and you only want to optimize one of the ligands, e.g., rb_min_jumps=1,2 "
		"would minimize only across jumps 1 and 2.");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"EnzRepackMinimize, similar in spirit to RepackMinimizeMover, does the "
		"design/repack followed by minimization of a protein-ligand (or TS "
		"model) interface with enzyme design style constraints (if present, "
		"see AddOrRemoveMatchCstsMover) using specified score functions and "
		"minimization dofs. Only design/repack or minimization can be done by "
		"setting appropriate tags. A shell of residues around the ligand are "
		"repacked/designed and/or minimized. If constrained optimization or "
		"cst_opt is specified, ligand neighbors are converted to Ala, "
		"minimization performed, and original neighbor sidechains are "
		"placed back.",
		attlist );
}

std::string EnzRepackMinimizeCreator::keyname() const {
	return EnzRepackMinimize::mover_name();
}

protocols::moves::MoverOP
EnzRepackMinimizeCreator::create_mover() const {
	return protocols::moves::MoverOP( new EnzRepackMinimize );
}

void EnzRepackMinimizeCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	EnzRepackMinimize::provide_xml_schema( xsd );
}



} //enzdes
} //protocols
