// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Kale Kundert
/// @author Roland A. Pache, PhD

// Headers
#include <protocols/loop_modeling/LoopBuilder.hh>
#include <protocols/loop_modeling/LoopBuilderCreator.hh>

// Core headers
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Basic headers
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/option.hh>

// Protocol headers
#include <protocols/kinematic_closure/KicMover.hh>
#include <protocols/kinematic_closure/perturbers/IdealizeNonPhiPsi.hh>
#include <protocols/kinematic_closure/perturbers/Rama2bPerturber.hh>
#include <protocols/kinematic_closure/perturbers/FragmentPerturber.hh>
#include <protocols/kinematic_closure/pivot_pickers/LoopPivots.hh>
#include <protocols/kinematic_closure/solution_pickers/FilteredSolutions.hh>
#include <protocols/loop_modeling/refiners/MinimizationRefiner.hh>
#include <protocols/loop_modeling/utilities/rosetta_scripts.hh>
#include <protocols/loop_modeling/utilities/TrajectoryLogger.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.hh>

// Namespaces
using namespace std;
using namespace basic::options;

using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;
using utility::tag::TagCOP;
using basic::datacache::DataMap;

namespace protocols {
namespace loop_modeling {

static basic::Tracer TR( "protocols.loop_modeling.LoopBuilder" );

LoopBuilder::LoopBuilder() {
	using protocols::kinematic_closure::KicMover;
	using protocols::kinematic_closure::KicMoverOP;
	using protocols::kinematic_closure::perturbers::IdealizeNonPhiPsi;
	using protocols::kinematic_closure::perturbers::Rama2bPerturber;
	using protocols::kinematic_closure::perturbers::PerturberOP;
	using protocols::kinematic_closure::pivot_pickers::LoopPivots;
	using protocols::kinematic_closure::pivot_pickers::PivotPickerOP;
	using protocols::kinematic_closure::solution_pickers::FilteredSolutions;
	using protocols::kinematic_closure::solution_pickers::FilteredSolutionsOP;
	using protocols::loop_modeling::refiners::MinimizationRefiner;
	using protocols::loop_modeling::refiners::MinimizationRefinerOP;
	using protocols::loop_modeling::utilities::TrajectoryLogger;
	using protocols::loop_modeling::utilities::TrajectoryLoggerOP;

	max_attempts_ = option[OptionKeys::loops::max_kic_build_attempts]();

	// The rama check currently works by comparing the generated torsions to the
	// input torsions.  Since the purpose of loop rebuilding is to forget
	// everything about the input structure, the rama check shouldn't be used.
	// One motivation for changing the rama check to use a static threshold (in
	// addition to simplicity) is that it could be used here.

	FilteredSolutionsOP solution_picker( new FilteredSolutions );
	solution_picker->dont_check_rama();
	solution_picker->be_lenient();

	kic_mover_ = add_child( utility::pointer::make_shared< KicMover >() );
	kic_mover_->add_perturber( utility::pointer::make_shared< IdealizeNonPhiPsi >() );
	kic_mover_->add_perturber( utility::pointer::make_shared< Rama2bPerturber >() );
	kic_mover_->set_pivot_picker( utility::pointer::make_shared< LoopPivots >() );
	kic_mover_->set_solution_picker( solution_picker );

	minimizer_ = add_child( utility::pointer::make_shared< MinimizationRefiner >() );

	logger_ = utility::pointer::make_shared< TrajectoryLogger >();
}

LoopBuilder::~LoopBuilder() = default;

void LoopBuilder::parse_my_tag(
	TagCOP tag,
	DataMap & data
) {

	LoopMover::parse_my_tag(tag, data);
	utilities::set_scorefxn_from_tag(*this, tag, data);
	max_attempts_ = tag->getOption<core::Size>("max_attempts", max_attempts_);
}

void LoopBuilder::use_fragments(
	utility::vector1<core::fragment::FragSetCOP> const & frag_libs) {

	using protocols::kinematic_closure::perturbers::IdealizeNonPhiPsi;
	using protocols::kinematic_closure::perturbers::Rama2bPerturber;
	using protocols::kinematic_closure::perturbers::FragmentPerturber;
	using protocols::kinematic_closure::perturbers::PerturberOP;

	// Note that a Rama2bPerturber is added just before the FragmentPerturber.
	// This is very important for benchmark runs seeking to recover the input
	// structure.  The FragmentPerturber will use a fragment even if it only
	// overlaps with the region being sampled by one residue.  When this happens,
	// KIC will often generate loops that are quite similar to the input loop.
	// For the benchmarks mentioned above, where the input loop is also the
	// target loop, this is a subtle but effective form of cheating.
	//
	// In production runs, the Rama2bPerturber doesn't really need to be here.
	// But there's also no reason for it not to be here, since it takes a
	// negligible amount of time to run.  And it's probably best to use the same
	// algorithm in the production runs as in the benchmark runs.

	kic_mover_->clear_perturbers();
	kic_mover_->add_perturber(utility::pointer::make_shared< IdealizeNonPhiPsi >());
	kic_mover_->add_perturber(utility::pointer::make_shared< Rama2bPerturber >());
	kic_mover_->add_perturber(utility::pointer::make_shared< FragmentPerturber >(frag_libs));
}

bool LoopBuilder::do_apply(Pose & pose, Loop const & loop) {

	// Only attempt to rebuild loops that are marked as "extended".

	if ( ! loop.is_extended() ) return true;

	// Setup the loop movers.

	kic_mover_->set_loop(loop);
	minimizer_->set_loop(loop);
	logger_->init(this);

	// Idealize the loop. The fold tree is changed temporarily
	// to preserve the positions of the start and stop residues.

	core::kinematics::FoldTree new_ft;
	new_ft.add_edge(1, loop.stop() - 1, core::kinematics::Edge::PEPTIDE);
	new_ft.add_edge(loop.start() - 1, loop.stop() + 1, 1);
	new_ft.add_edge(loop.stop() + 1, loop.stop(), core::kinematics::Edge::PEPTIDE);
	new_ft.add_edge(loop.stop() + 1, pose.size(), core::kinematics::Edge::PEPTIDE);

	idealize_loop(pose, loop, new_ft);

	// Make a strong effort to rebuild the loop with KIC.

	kic_mover_->was_successful(false);
	for ( core::Size i = 1; i <= max_attempts_ && ! kic_mover_->was_successful(); i++ ) {
		kic_mover_->apply(pose);

		bool success = kic_mover_->was_successful();
		logger_->record_move(TR, pose, i, success, success);
	}

	if ( ! kic_mover_->was_successful() ) return false;

	// Idealize loop again to correct the bond angles of O and H

	core::kinematics::FoldTree new_ft2(pose.size());
	idealize_loop(pose, loop, new_ft2);

	// Minimize the loop if it was successfully built.

	minimizer_->apply(pose);
	logger_->record_endpoint(TR, pose);

	return minimizer_->was_successful();
}

Size LoopBuilder::get_max_attempts() const {
	return max_attempts_;
}

void LoopBuilder::set_max_attempts(core::Size attempts) {
	max_attempts_ = attempts;
}

ScoreFunctionOP LoopBuilder::get_score_function() const {
	return get_tool<ScoreFunctionOP>(ToolboxKeys::SCOREFXN);
}

void LoopBuilder::set_score_function(ScoreFunctionOP score_function) {
	set_tool(ToolboxKeys::SCOREFXN, score_function);
}

utilities::TrajectoryLoggerOP LoopBuilder::get_logger() const {
	return logger_;
}

string LoopBuilder::get_name() const {
	return mover_name();
}

string LoopBuilder::mover_name() {
	return "LoopBuilder";
}

// handle parse_score_function for derived classes
void LoopBuilder::get_score_function_attributes( utility::tag::AttributeList &attlist ) {
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );
}

void LoopBuilder::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	// Create a complex type and  get the LoopMover attributes, as parse_my_tag calls LoopMover::parse_my_tag
	XMLSchemaComplexTypeGenerator ct_gen;
	// subelement_list for the LoopMover subelements
	XMLSchemaSimpleSubelementList subelement_list;
	LoopMover::define_composition_schema( xsd, ct_gen, subelement_list );

	ct_gen.element_name( mover_name() )//complex_type_naming_func( & moves::complex_type_name_for_mover )
		.description( "Builds in backbone atoms for loop regions where they are missing. "
		"The backbones created by LoopBuilder will have ideal bond lengths, ideal bond angles, "
		"and torsions picked from a Ramachandran distribution."  )
		.add_attribute( XMLSchemaAttribute("max_attempts", xsct_non_negative_integer, "Stop after n cycles irrespective of convergence.")  )
		.write_complex_type_to_schema( xsd );
}

std::string LoopBuilderCreator::keyname() const {
	return LoopBuilder::mover_name();
}

protocols::moves::MoverOP
LoopBuilderCreator::create_mover() const {
	return utility::pointer::make_shared< LoopBuilder >();
}

void LoopBuilderCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoopBuilder::provide_xml_schema( xsd );
}

void LoopBuilder::idealize_loop(Pose & pose, Loop const & loop, core::kinematics::FoldTree &ft) const{
	core::kinematics::FoldTree original_ft(pose.fold_tree());
	pose.fold_tree(ft);
	loops::idealize_loop(pose, loop);
	pose.fold_tree(original_ft);
}



}
}
