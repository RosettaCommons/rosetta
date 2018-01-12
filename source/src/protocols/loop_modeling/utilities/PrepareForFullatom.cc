// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Headers {{{1
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/utilities/PrepareForFullatom.hh>
#include <protocols/loop_modeling/utilities/PrepareForFullatomCreator.hh>
#include <protocols/loop_modeling/utilities/rosetta_scripts.hh>

// Protocol headers
#include <protocols/moves/Mover.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>

// Core headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// Utility headers
#include <utility/exit.hh>

// RosettaScripts headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Namespaces {{{1
using namespace std;

using core::Size;
using core::chemical::DISULFIDE;
using core::chemical::FA_STANDARD;
using core::kinematics::MoveMap;
using core::optimization::AtomTreeMinimizer;
using core::optimization::AtomTreeMinimizerOP;
using core::optimization::MinimizerOptions;
using core::optimization::symmetry::SymAtomTreeMinimizer;
using core::pose::Pose;
using core::pose::symmetry::is_symmetric;
using core::pose::symmetry::make_residue_mask_symmetric;
using core::pose::symmetry::make_symmetric_movemap;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;
using core::util::switch_to_residue_type_set;

using protocols::minimization_packing::PackRotamersMover;
using protocols::minimization_packing::PackRotamersMoverOP;
using protocols::minimization_packing::symmetry::SymPackRotamersMover;
// }}}1

namespace protocols {
namespace loop_modeling {
namespace utilities {

// XRW TEMP moves::MoverOP PrepareForFullatomCreator::create_mover() const { // {{{1
// XRW TEMP  return moves::MoverOP( new PrepareForFullatom );
// XRW TEMP }

// XRW TEMP string PrepareForFullatomCreator::keyname() const { // {{{1
// XRW TEMP  return "PrepareForFullatom";
// XRW TEMP }
// }}}1

PrepareForFullatom::PrepareForFullatom() // {{{1
: force_repack_(false) {}

void PrepareForFullatom::parse_my_tag( // {{{1
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose) {

	LoopMover::parse_my_tag(tag, data, filters, movers, pose);
	utilities::set_scorefxn_from_tag(*this, tag, data);

	force_repack_ = tag->getOption<bool>("force_repack", force_repack_);
}

bool PrepareForFullatom::do_apply(Pose & pose) { // {{{1
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using protocols::simple_moves::ReturnSidechainMover;

	if ( original_pose_.size() == 0 && ! force_repack_ ) {
		utility_exit_with_message("<PrepareForFullatom> cannot copy sidechains because no original pose was specified.");
	}

	// If the pose is already fullatom, this means that the centroid stages were
	// skipped and that the pose has been fullatom the whole time.  So we don't
	// need to do anything unless a full repack has been requested.

	if ( pose.is_fullatom() && ! force_repack_ ) { return true; }

	// Convert the pose to fullatom mode.

	switch_to_residue_type_set(pose, core::chemical::FULL_ATOM_t );

	// Copy sidechains from the original pose.

	if ( original_pose_.is_fullatom() ) {
		ReturnSidechainMover return_sidechains(original_pose_);
		return_sidechains.apply(pose);
	}

	// Decide which sidechains to repack.  Every sidechain will be repacked if
	// the original pose is in centroid mode or if a full repack was requested.
	// Otherwise only sidechains within the loop are repacked.

	vector1<bool> residues_to_repack(pose.size(), false);

	for ( Size i = 1; i < pose.size(); i++ ) {
		residues_to_repack[i] =
			force_repack_ ||
			original_pose_.is_centroid() ||
			get_loops()->is_loop_residue(i);
	}

	pose.conformation().detect_disulfides();

	ExtraRotamersGenericOP extra_rotamers( new ExtraRotamersGeneric );
	extra_rotamers->ex1(true);
	extra_rotamers->ex2(true);
	extra_rotamers->extrachi_cutoff(0);

	TaskFactoryOP task_factory( new TaskFactory );
	task_factory->push_back(TaskOperationCOP( new RestrictToRepacking ));
	task_factory->push_back(TaskOperationCOP( new NoRepackDisulfides ));
	task_factory->push_back(extra_rotamers);

	PackerTaskOP task = task_factory->create_task_and_apply_taskoperations(pose);
	make_residue_mask_symmetric(pose, residues_to_repack);
	task->restrict_to_residues(residues_to_repack);

	// Setup a move map for the minimizer that will allow only sidechain DOFs to
	// move.  Again, disulfides are explicitly left in place.

	MoveMap move_map;
	move_map.set_bb(false);
	move_map.set_chi(true);

	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( pose.residue(i).has_variant_type(DISULFIDE) ) {
			move_map.set_chi(i, false);
		}
	}

	// Apply the packer and the minimizer, accounting for symmetry.

	PackRotamersMoverOP packer;
	AtomTreeMinimizerOP minimizer;
	MinimizerOptions min_options("lbfgs_armijo_nonmonotone", 1e-5, true, false);
	ScoreFunctionCOP fa_score_function = get_score_function();

	if ( is_symmetric(pose) ) {
		packer = PackRotamersMoverOP( new SymPackRotamersMover(fa_score_function, task) );
		minimizer = AtomTreeMinimizerOP( new SymAtomTreeMinimizer );
		make_symmetric_movemap(pose, move_map);
	} else {
		packer = PackRotamersMoverOP( new PackRotamersMover(fa_score_function, task) );
		minimizer = AtomTreeMinimizerOP( new AtomTreeMinimizer );
	}

	packer->apply(pose);
	minimizer->run(pose, move_map, *fa_score_function, min_options);

	return true;
}

void PrepareForFullatom::set_original_pose(Pose const & pose) { // {{{1
	original_pose_ = pose;
}

ScoreFunctionOP PrepareForFullatom::get_score_function() { // {{{1
	return get_tool<ScoreFunctionOP>(ToolboxKeys::SCOREFXN);
}

void PrepareForFullatom::set_score_function(ScoreFunctionOP score_function) { // {{{1
	set_tool(ToolboxKeys::SCOREFXN, score_function);
}
void PrepareForFullatom::set_force_repack(bool value) { // {{{1
	force_repack_ = value;
}

std::string PrepareForFullatom::get_name() const {
	return mover_name();
}

std::string PrepareForFullatom::mover_name() {
	return "PrepareForFullatom";
}

void PrepareForFullatom::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	// add score function tags
	utilities::attributes_for_set_scorefxn_from_tag( attlist );
	// Create a complex type and  get the LoopMover attributes, as parse_my_tag calls LoopMover::parse_my_tag
	XMLSchemaComplexTypeGenerator ct_gen;
	// subelement_list for the LoopMover subelements
	XMLSchemaSimpleSubelementList subelement_list;
	LoopMover::define_composition_schema( xsd, ct_gen, subelement_list );
	ct_gen.element_name( mover_name() )
		.description( "Switch from centroid to full atom representation." )
		.add_attribute( XMLSchemaAttribute("force_repack", xsct_rosetta_bool, "Enforce side chain repackig")  )
		.write_complex_type_to_schema( xsd );
}

std::string PrepareForFullatomCreator::keyname() const {
	return PrepareForFullatom::mover_name();
}

protocols::moves::MoverOP
PrepareForFullatomCreator::create_mover() const {
	return protocols::moves::MoverOP( new PrepareForFullatom );
}

void PrepareForFullatomCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PrepareForFullatom::provide_xml_schema( xsd );
}

// }}}1

}
}
}
