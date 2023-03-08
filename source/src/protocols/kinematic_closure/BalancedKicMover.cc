// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/kinematic_closure/BalancedKicMover.cc
/// @brief implementation of KIC with detailed balance
/// @author Kale Kundert, Ameya Harmalkar (aharmal1@jhu.edu)

// Headers {{{1
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/internal.hh>
#include <protocols/kinematic_closure/BalancedKicMover.hh>
#include <protocols/kinematic_closure/BalancedKicMoverCreator.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.fwd.hh>
#include <protocols/kinematic_closure/perturbers/PerturberSet.hh>
#include <protocols/kinematic_closure/perturbers/RamaPerturber.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.fwd.hh>
#include <protocols/kinematic_closure/pivot_pickers/StandardPivots.hh>
#include <protocols/rosetta_scripts/util.hh>

// Core headers
#include <core/id/types.hh>
#include <core/id/TorsionID.hh>
#include <core/id/TorsionID_Range.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>

// Utility headers
#include <utility/exit.hh>
#include <numeric/constants.hh>
#include <numeric/random/random.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

// Option key headers
#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// XSD headers
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ headers
#include <iostream>
#include <limits>

// Global Names {{{1

using namespace std;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using core::id::TorsionID;
using core::id::TorsionID_Range;
using protocols::kinematic_closure::pivot_pickers::PivotPickerOP;
using utility::vector1;
static basic::Tracer TR( "protocols.kinematic_closure.BalancedKicMover" );
// }}}1

namespace protocols {
namespace kinematic_closure {

// Member functions:

BalancedKicMover::BalancedKicMover() { // {{{1
	init_from_options();
	perturbers_ = utility::pointer::make_shared< perturbers::PerturberSet >();
	perturbers_->add(utility::pointer::make_shared< perturbers::RamaPerturber >());
	perturbers_->mark_as_default();

	loop_ = Loop(0, 0);
	pivot_picker_ = utility::pointer::make_shared< pivot_pickers::StandardPivots >();
	is_fold_tree_stale_ = true;
	selector_on_ = false;
}

protocols::moves::MoverOP BalancedKicMover::clone() const {
	return utility::pointer::make_shared< BalancedKicMover >( *this );
}

protocols::moves::MoverOP BalancedKicMover::fresh_instance() const {
	return utility::pointer::make_shared< BalancedKicMover >();
}

BalancedKicMover::~BalancedKicMover() = default; // {{{1

void BalancedKicMover::init_from_options(){
	if ( loops_file_.empty() ) {
		if ( option[basic::options::OptionKeys::loops::loop_file].user() ) {
			loops_file_ = option[basic::options::OptionKeys::loops::loop_file]()[1];
		}
	}
}

void BalancedKicMover::apply(Pose & pose) { // {{{1
	ClosureProblemOP problem( new ClosureProblem() );
	ClosureSolutionCOP solution;
	SolutionList unperturbed_solutions, perturbed_solutions;
	TR.Debug << "Commence loop/s set-up" << std::endl;

	if ( !(loops_file_.empty()) ) {
		TR.Debug << "Loading loops from loops_file" << std::endl;
		if ( loops_ == nullptr ) {
			protocols::loops::LoopsOP loops( utility::pointer::make_shared< protocols::loops::Loops >(loops_file_) );
			TR.Info << "Initial loops " << *loops << std::endl;
			set_loops(loops);
		}
		TR.Debug << "Warning: Residue selectors will override the flag inputs." << std::endl;
	}

	if ( selector_on_ ) {
		utility::vector1< bool > selection = residue_selector_->apply(pose);
		protocols::loops::LoopsOP loops( utility::pointer::make_shared< protocols::loops::Loops >(selection) );
		TR.Warning << "Residue selector is overwriting over the loops_file input (if any)." << std::endl;
		TR.Info << "Initial loops " << *loops << std::endl;
		set_loops(loops);
	}
	// To do: check that cutpoint is not outside the loop and if necessary notify user that only first loop is used

	if ( loops_ == nullptr ) {
		utility_exit_with_message(
			"Before calling BalancedKicMover.apply(), you must provide a loop either "
			"via BalancedKicMover.set_loop(), a loops file or residue_selector tag in XML.");
	}

	TR.Debug << "Loop object has been set-up successfully." << std::endl;
	TR.Debug << "Original Foldtree : " << pose.fold_tree() << std::endl;
	core::kinematics::FoldTree og_ft( pose.fold_tree() );

	for ( auto it = loops_->v_begin(); it != loops_->v_end(); ++it ) {
		// setting all loops to be extended from loop related loop modeling code. Is it really necessary?
		it->set_extended( true );
		protocols::loops::set_single_loop_fold_tree(pose, *it);

		// Solve the unperturbed problem.
		problem->frame(pose, *it, pivot_picker_);
		unperturbed_solutions = problem->solve();

		// Solve the perturbed problem.
		perturbers_->perturb_with_balance(pose, problem);
		perturbed_solutions = problem->solve();

		// Pick a solution to apply.
		solution = pick_solution(unperturbed_solutions, perturbed_solutions);
		solution->apply(pose);

		// Decide if a real move was made.
		bool is_trivial = is_solution_trivial(
			problem, solution, unperturbed_solutions, perturbed_solutions);
		type(is_trivial ? "balanced-kic-no-op" : "balanced-kic");
	}
}

void BalancedKicMover::set_loop(Loop const & loop) { // {{{1
	loop_ = loop;
	is_fold_tree_stale_ = true;
}

void BalancedKicMover::set_loops(protocols::loops::LoopsOP const loops){
	loops_ = loops;
	protocols::loops::Loop extend_loop( *(loops_->v_begin()) ); // make a local copy of the loop
	set_loop( extend_loop );  // Reset the loop_ object to the first loop in the list
	is_fold_tree_stale_ = true;
}

void BalancedKicMover::set_residue_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	runtime_assert_string_msg( selector, "ERROR: protocols::kinematic_closure::BalancedKicMover::set_residue_selector(): A null pointer was passed to this function." );
	residue_selector_ = selector;
}

void BalancedKicMover::add_perturber(perturbers::PerturberOP perturber) { // {{{1
	perturbers_->add(perturber) ;
}

void BalancedKicMover::set_pivot_picker(PivotPickerOP picker) { // {{{1
	pivot_picker_ = picker;
}

vector1<TorsionID_Range> BalancedKicMover::torsion_id_ranges(Pose &) { // {{{1
	using core::id::BB;

	vector1<TorsionID_Range> results;
	Real static const pi = numeric::constants::r::pi;



	// For all the loop residues
	for ( auto it = loops_->v_begin(); it != loops_->v_end(); ++it ) {
		//protocols::loops::Loop const & single_loop(new protocols::loops::Loop() );
		protocols::loops::Loop const & single_loop(*it);
		std::cout << "Printing out the loop here : " << single_loop << std::endl;
		for ( core::Size residue = single_loop.start(); residue <= single_loop.stop(); residue++ ) {
			TorsionID phi(residue, BB, 1);
			TorsionID psi(residue, BB, 2);
			TorsionID omega(residue, BB, 3);
			results.push_back(TorsionID_Range(phi, -pi, pi));
		} // all residues of single loop
	} // add residues of all the loops
	return results;
}
// }}}1

// Static member functions:

// {{{1
/// @details Note that this is a static method, which means that it can be used
/// outside the context of this class.  The inputs are two sets of closure
/// solutions.  The first set should contain the starting conformation, and the
/// second should contain novel conformations.  The first set can be generated
/// by running the KIC algorithm without perturbing the non-pivot torsions.
/// For each solution in these two sets, a jacobian will have to be calculated.
/// The solution that is returned will be free of bias.  See apply() for an
/// example of how this method is used.

ClosureSolutionCOP BalancedKicMover::pick_solution(
	SolutionList const & unperturbed_solutions,
	SolutionList const & perturbed_solutions) {

	Real total_jacobian = 0;
	Real selection_chance = 0;
	Real random_threshold = numeric::random::uniform();
	ChainedSolutionList all_solutions(
		unperturbed_solutions, perturbed_solutions);

	runtime_assert(! all_solutions.empty());

	// Calculate the Jacobian for each solution and keep the sum.  The Jacobian
	// relates to the probability that a certain set of pivot torsions will lead
	// to a closed solution.

	for ( ClosureSolutionCOP solution : all_solutions ) {
		total_jacobian += solution->get_jacobian();
	}

	// Use the Jacobian weights to pick a balanced solution.  The pick must be
	// made from both the perturbed and unperturbed pools of solutions, otherwise
	// the forward and reverse move probabilities won't be equivalent.

	for ( ClosureSolutionCOP solution : all_solutions ) {
		selection_chance += solution->get_jacobian() / total_jacobian;
		if ( selection_chance >= random_threshold ) return solution;
	}

	// Execution will only get this far if random_threshold is very close to one
	// and floating point error causes the sum of all the selection_chance terms
	// to be slightly less than that.  In this case, the right course of action
	// is clearly to return the last solution.

	return all_solutions.back();
}

// {{{1
/// @details In order to obey detailed balance, the set of unperturbed
/// solutions must contain one solution that is identical to the input pose.
/// When this solution is picked and applied, the resulting move will pass the
/// Metropolis criterion and lead to an artificially inflated acceptance rate.
/// This method provides a way to report what really happened in the move, so
/// that an accurate acceptance rate can be conveyed.  Note that this is a
/// static method, so it can be used outside the context of this class.  The
/// input pose conformation is inferred from the given problem.

bool BalancedKicMover::is_solution_trivial(
	ClosureProblemCOP problem,
	ClosureSolutionCOP picked_solution,
	SolutionList const & unperturbed_solutions,
	SolutionList const & perturbed_solutions)
{
	Real distance, closest_distance = numeric_limits<Real>::infinity();
	ClosureSolutionCOP closest_solution = nullptr;
	ChainedSolutionList all_solutions(
		unperturbed_solutions, perturbed_solutions);

	for ( ClosureSolutionCOP solution : all_solutions ) {
		distance = solution->get_distance(problem.get());
		if ( distance < closest_distance ) {
			closest_distance = distance;
			closest_solution = solution;
		}
	}

	runtime_assert(closest_solution.get() != nullptr);
	return picked_solution.get() == closest_solution.get();
}
// }}}1

std::string BalancedKicMover::get_name() const {
	return mover_name();
}

std::string BalancedKicMover::mover_name() {
	return "BalancedKicMover";
}

void BalancedKicMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( tag->hasOption("residue_selector") ) {
		set_residue_selector( protocols::rosetta_scripts::parse_residue_selector( tag, data ) );
		TR.Debug << "INFO: The residue selectors option has been set." << std::endl;
		// Make sure the loops file is ignored if residue_selectors are provided.
		selector_on_ = true;
	}

	if ( tag->hasOption("pivot_residues") ) {
		// pivot_residue_selector_ = core::pose::get_resnum_selector(tag, "pivot_residues");
		// Clear the command-line input movemap
		utility_exit_with_message(
			"INFO: The pivot residues option has not been set yet. Please use the residue_selector tag instead.");
	}

	if ( tag->hasOption("preserve_detailed_balance") ) {
		TR.Warning << "INFO: The detailed balanced option doesn't require to be set." << std::endl;
	}
	//set_preserve_detailed_balance( tag->getOption<bool>( "preserve_detailed_balance", preserve_detailed_balance() ) );
}

void BalancedKicMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	auto ct_gen = complex_type_generator_for_balancedKIC_mover( xsd );
	ct_gen->element_name( mover_name() )
		.description( "Performs kinematic closure moves with detailed balance along two pivot atoms" )
		.write_complex_type_to_schema( xsd );
}

utility::tag::XMLSchemaComplexTypeGeneratorOP
BalancedKicMover::complex_type_generator_for_balancedKIC_mover( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute(
		"pivot_residues", xs_string,
		"residues for which contiguous stretches can contain segments "
		"(comma separated) can use PDB numbers ([resnum][chain]) or "
		"absolute Rosetta numbers (integer). Note that this feature has not been implemented yet.");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"preserve_detailed_balance", xsct_rosetta_bool,
		"if set to true, does not change branching atom angles during apply and "
		"sets ideal branch angles during initialization if used with MetropolisHastings. "
		"Note that this feature has not been implemented yet.",
		"false");

	// get the attributes for parsing residue_selector
	// core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "residue_selector", "Pre-defined residue_selector" );
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector", "Pre-defined residue_selector" );

	XMLSchemaSimpleSubelementList subelements;
	rosetta_scripts::append_subelement_for_parse_movemap_factory_legacy( xsd, subelements );

	//XMLSchemaComplexTypeGeneratorOP ct_gen( utility::pointer::make_shared< XMLSchemaComplexTypeGenerator >() );
	XMLSchemaComplexTypeGeneratorOP ct_gen( utility::pointer::make_shared< XMLSchemaComplexTypeGenerator >() );
	ct_gen->complex_type_naming_func( & moves::complex_type_name_for_mover )
		.add_attributes( attlist )
		.set_subelements_repeatable( subelements )
		.add_optional_name_attribute();
	return ct_gen;
}

std::string BalancedKicMoverCreator::keyname() const {
	return BalancedKicMover::mover_name();
}

protocols::moves::MoverOP
BalancedKicMoverCreator::create_mover() const {

	return utility::pointer::make_shared< BalancedKicMover >();
}

void BalancedKicMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BalancedKicMover::provide_xml_schema( xsd );
}

}
}
