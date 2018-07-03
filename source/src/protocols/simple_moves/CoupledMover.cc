// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file CoupledMover.cc
/// @brief implementation of CoupledMover class and functions
/// @author Noah Ollikainen (nollikai@gmail.com)
/// @author Anum Glasgow (anumazam@gmail.com)
/// @author Amanda Loshbaugh (aloshbau@gmail.com)

// Unit headers
#include <protocols/simple_moves/CoupledMover.hh>
#include <protocols/simple_moves/CoupledMoverCreator.hh>
#include <protocols/kinematic_closure/KicMover.hh>


// Project headers
#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>

#include <protocols/backrub/BackrubMover.hh>
#include <core/conformation/Residue.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <protocols/minimization_packing/BoltzmannRotamerMover.hh>

// For fixing up the neighborhood
#include <protocols/minimization_packing/MinPackMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>

// Parser headers
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/angle.functions.hh>
#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/xyz.functions.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.simple_moves.CoupledMover");

namespace protocols {
namespace simple_moves {

// default constructor
CoupledMover::CoupledMover() : protocols::moves::Mover()
{
	protocols::moves::Mover::type( "Coupled" );
	resnum_ = 0;
	randomize_resnum_ = false;
	fix_backbone_ = false;
	rotation_std_dev_ = 4.572016;
	uniform_backrub_ = false;
	temperature_ = 1.0;
	bias_sampling_ = true;
	bump_check_ = true;
	ligand_resnum_ = 0;
	ligand_jump_id_ = 0;
	ligand_weight_ = 1.0;
	rotation_magnitude_ = 1;
	translation_magnitude_ = 0.1;
	short_backrub_mover_ = protocols::simple_moves::ShortBackrubMoverOP( new protocols::simple_moves::ShortBackrubMover() );
	short_backrub_mover_->set_rotation_std_dev( rotation_std_dev_ );
	loop_size_ = 4;
	boltzmann_rotamer_mover_ = protocols::minimization_packing::BoltzmannRotamerMoverOP( new protocols::minimization_packing::BoltzmannRotamerMover() );
	boltzmann_rotamer_mover_->set_temperature( temperature_ );
	boltzmann_rotamer_mover_->set_bias_sampling( bias_sampling_ );
	boltzmann_rotamer_mover_->set_ligand_resnum( ligand_resnum_ );
	boltzmann_rotamer_mover_->set_ligand_weight( ligand_weight_ );
	repack_neighborhood_ = false; // Default legacy behavior
}

/// @brief constructor that sets input pose, score function and packer task
// This is the constructer used by CoupledMovesProtocol during NOT ligand mode
CoupledMover::CoupledMover(
	core::pose::PoseOP pose,
	core::scoring::ScoreFunctionOP score_fxn,
	core::pack::task::PackerTaskOP packer_task ) : protocols::simple_moves::CoupledMover()
{
	short_backrub_mover_ = protocols::simple_moves::ShortBackrubMoverOP( new protocols::simple_moves::ShortBackrubMover( pose ) );
	set_score_fxn( score_fxn );
	set_packer_task( packer_task );
}

/// @brief constructor that sets input pose, score function, packer task, and ligand residue number
// This is the constructor used by CoupledMovesProtocol during ligand mode
CoupledMover::CoupledMover(
	core::pose::PoseOP pose,
	core::scoring::ScoreFunctionOP score_fxn,
	core::pack::task::PackerTaskOP packer_task,
	core::Size ligand_resnum ) : protocols::simple_moves::CoupledMover( pose, score_fxn, packer_task )
{
	ligand_resnum_ = ligand_resnum;
	ligand_jump_id_ = pose->fold_tree().get_jump_that_builds_residue( ligand_resnum );
	rigid_body_mover_ = protocols::rigid::RigidBodyPerturbMoverOP( new protocols::rigid::RigidBodyPerturbMover( ligand_jump_id_, rotation_magnitude_, translation_magnitude_ ) );
	set_rigid_body_mover( rigid_body_mover_ );
}

// copy constructor
CoupledMover::CoupledMover( CoupledMover const & )= default;

// destructor
CoupledMover::~CoupledMover()= default;

// clone this object
CoupledMover::MoverOP
CoupledMover::clone() const
{
	return CoupledMover::MoverOP( new protocols::simple_moves::CoupledMover( *this ) );
}

// create this type of object
CoupledMover::MoverOP
CoupledMover::fresh_instance() const
{
	return CoupledMover::MoverOP( new protocols::simple_moves::CoupledMover() );
}

void
CoupledMover::apply( core::pose::Pose & pose )
{

	if ( resnum_ == 0 ) {
		randomize_resnum_ = true;
	}

	if ( randomize_resnum_ ) {
		utility::vector1< core::Size > move_positions;
		for ( core::Size i = 1; i <= packer_task_->total_residue(); i++ ) {
			if ( packer_task_->pack_residue( i ) || packer_task_->design_residue( i ) ) {
				move_positions.push_back( i );
			}
		}
		core::Size random_index = numeric::random::rg().random_range( 1, move_positions.size() );
		resnum_ = move_positions[ random_index ];
	}

	if ( pose.residue( resnum_ ).is_protein() ) {
		if ( fix_backbone_ == false ) {
			// Backbone move is backrub
			if ( backbone_mover_ == "backrub" ) {
				short_backrub_mover_->set_resnum( resnum_ );
				short_backrub_mover_->apply( pose );
			} else if ( backbone_mover_ == "kic" ) {
				// Backbone move is KIC
				protocols::kinematic_closure::KicMoverOP fullatom_kic_mover( new protocols::kinematic_closure::KicMover() );
				fullatom_kic_mover->add_perturber( perturber_ );
				using protocols::loops::Loop;
				if ( loop_size_ == 0 ) {
					loop_size_ = numeric::random::rg().random_range( 3, 7 );
				}

				// Make sure you give KIC residues that are on the same chain and in the pose
				core::Size chainID = pose.residue( resnum_ ).chain();
				core::Size start_index = 0;
				if ( ( resnum_-loop_size_ >= 1 ) && ( resnum_-loop_size_ <= pose.size() ) ) {
					start_index = (resnum_-loop_size_);
					if ( pose.residue( start_index ).chain() != chainID ) {
						start_index = 0;
					}
				}
				core::Size stop_index = 0;
				if ( ( resnum_+loop_size_ >= 1 ) && ( resnum_+loop_size_ <= (pose.size()-1) ) ) {
					stop_index = (resnum_+loop_size_);
					if ( pose.residue( stop_index ).chain() != chainID ) {
						stop_index = 0;
					}
				}
				if ( start_index && stop_index ) {
					// KIC segfaults if it's fed a loop containing an incompatible residue type, e.g. ligand. This loop catches and prevents that.
					core::Size loop_is_protein = 1;
					for ( core::Size i = start_index; i <= stop_index; i++ ) {
						if ( pose.residue_type(i).is_protein() || pose.residue_type(i).is_peptoid() || pose.residue_type(i).is_carbohydrate() ) {
							continue;
						} else { loop_is_protein = 0; }
					}
					if ( loop_is_protein ) {
						Loop loop( start_index, stop_index, resnum_, 0.0, false ); // start_in, stop_in, cut_in, skip_rate, extended_in
						fullatom_kic_mover->set_loop( loop );
						TR << TR.White << "Applying KIC to loop " << start_index << "-" << stop_index << TR.Reset << std::endl;
						fullatom_kic_mover->apply( pose );
					} else if ( ! loop_is_protein ) {
						TR << TR.Red << "[ WARNING ] Didn't apply KIC move to loop around residue number " << resnum_ << " (Rosetta numbering) because the loop contains a residue type not compatible with KIC (e.g. a ligand)." << TR.Reset << std::endl;
					}
				} else if ( ! start_index || ! stop_index ) {
					TR << TR.Red << "[ WARNING ] Didn't apply KIC move to loop around residue number " << resnum_ << " (Rosetta numbering) because it's too close to end of chain or pose." << TR.Reset << std::endl;
				}
			}
		}
	} else {
		runtime_assert_string_msg( (get_ligand_resnum() != 0), "In CoupledMover, ligand mode not active but a ligand is present: please make CoupledMoves aware of the ligands.");
		rigid_body_mover_->apply( pose );
	}

	boltzmann_rotamer_mover_ = get_boltzmann_rotamer_mover();
	boltzmann_rotamer_mover_->set_task( packer_task_ );
	boltzmann_rotamer_mover_->set_scorefxn( score_fxn_ );
	// Apply boltzmann_rotamer_mover to design residue
	boltzmann_rotamer_mover_->set_resnum( resnum_ );
	TR << "Applying boltzmann_rotamer_mover to residue number " << resnum_  << TR.Reset << std::endl;
	boltzmann_rotamer_mover_->apply( pose );

	if ( repack_neighborhood_ ) {
		// Repack shell around design position
		TR << "Repacking neighborhood around design position "  << resnum_ << TR.Reset << std::endl;
		core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
		// Set all positions to repack only
		tf->push_back( core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );
		// Prevent non-neighbor positions from repacking
		core::select::residue_selector::NeighborhoodResidueSelector nbr_selector = core::select::residue_selector::NeighborhoodResidueSelector();
		utility::vector1< bool > focus(pose.size(), false);
		focus[ resnum_ ] = true;
		nbr_selector.set_focus( focus );
		nbr_selector.set_distance( 5 ); // distance of 4 results in no neighbors found. distance of 6 results in drastically lowered acceptance rate.
		nbr_selector.set_include_focus_in_subset( true );
		core::pack::task::operation::PreventRepackingRLTOP prevent_repacking = core::pack::task::operation::PreventRepackingRLTOP( new core::pack::task::operation::PreventRepackingRLT() );
		core::select::residue_selector::ResidueSelectorOP nbr_selectorOP = nbr_selector.clone();
		core::pack::task::operation::OperateOnResidueSubsetOP operate_on_neighbors = core::pack::task::operation::OperateOnResidueSubsetOP( new core::pack::task::operation::OperateOnResidueSubset() );
		operate_on_neighbors->op( prevent_repacking );
		operate_on_neighbors->selector( nbr_selectorOP );

		operate_on_neighbors->flip_subset( true );
		tf->push_back( operate_on_neighbors );
		// Apply to pose
		core::pack::task::PackerTaskOP pt( tf->create_packer_task( pose ) );
		protocols::minimization_packing::PackRotamersMover prm = protocols::minimization_packing::PackRotamersMover( score_fxn_ );

		//protocols::simple_moves::MinPackMoverOP minpack(new protocols::simple_moves::MinPackMover( score_fxn_ ));
		//minpack->task_factory( tf );
		//minpack->apply( pose );

		prm.task_factory( tf );
		prm.apply( pose );
	}

} // apply

// setters
void CoupledMover::set_resnum( core::Size resnum ) { resnum_ = resnum; }
void CoupledMover::set_randomize_resnum( bool randomize_resnum ) { randomize_resnum_ = randomize_resnum; }
void CoupledMover::set_fix_backbone( bool fix_backbone ) { fix_backbone_ = fix_backbone; }
void CoupledMover::set_rotation_std_dev( core::Real rotation_std_dev ) {
	rotation_std_dev_ = rotation_std_dev;
	short_backrub_mover_->set_rotation_std_dev(rotation_std_dev);
}
void CoupledMover::set_uniform_backrub( bool uniform_backrub ) {
	uniform_backrub_ = uniform_backrub;
	short_backrub_mover_->set_uniform_backrub( uniform_backrub );
}
void CoupledMover::set_input_pose( core::pose::PoseCOP pose ) {
	short_backrub_mover_->set_input_pose( pose );
}
void CoupledMover::set_temperature( core::Real temperature ) {
	temperature_ = temperature;
	boltzmann_rotamer_mover_->set_temperature( temperature );
}
void CoupledMover::set_bias_sampling( bool bias_sampling ) {
	bias_sampling_ = bias_sampling;
	boltzmann_rotamer_mover_->set_bias_sampling( bias_sampling );
}
void CoupledMover::set_bump_check( bool bump_check ) {
	bump_check_ = bump_check;
	boltzmann_rotamer_mover_->set_bump_check( bump_check );
}
void CoupledMover::set_ligand_resnum( core::Size ligand_resnum, core::pose::PoseCOP pose ) {
	ligand_resnum_ = ligand_resnum;
	boltzmann_rotamer_mover_->set_ligand_resnum( ligand_resnum );
	set_ligand_jump_id( pose->fold_tree().get_jump_that_builds_residue( ligand_resnum ) );
}
void CoupledMover::set_ligand_jump_id( core::Size ligand_jump_id ) {
	ligand_jump_id_ = ligand_jump_id;
	rigid_body_mover_ = protocols::rigid::RigidBodyPerturbMoverOP( new protocols::rigid::RigidBodyPerturbMover( ligand_jump_id_, rotation_magnitude_, translation_magnitude_ ) );
}
void CoupledMover::set_ligand_weight( core::Real ligand_weight ) {
	ligand_weight_ = ligand_weight;
	boltzmann_rotamer_mover_->set_ligand_weight( ligand_weight );
}
void CoupledMover::set_rotation_magnitude( core::Real rotation_magnitude ) {
	rotation_magnitude_ = rotation_magnitude;
	rigid_body_mover_->rot_magnitude( rotation_magnitude );
}
void CoupledMover::set_translation_magnitude( core::Real translation_magnitude ) {
	translation_magnitude_ = translation_magnitude;
	rigid_body_mover_->trans_magnitude( translation_magnitude );
}
void CoupledMover::set_short_backrub_mover( protocols::simple_moves::ShortBackrubMoverOP short_backrub_mover ) {
	short_backrub_mover_ = short_backrub_mover;
}
void CoupledMover::set_boltzmann_rotamer_mover( protocols::minimization_packing::BoltzmannRotamerMoverOP boltzmann_rotamer_mover ) { boltzmann_rotamer_mover_ = boltzmann_rotamer_mover;
}
void CoupledMover::set_rigid_body_mover( protocols::rigid::RigidBodyPerturbMoverOP rigid_body_mover ) {
	rigid_body_mover_ = rigid_body_mover;
}
void CoupledMover::set_score_fxn( core::scoring::ScoreFunctionOP score_fxn ) {
	score_fxn_ = score_fxn;
}
void CoupledMover::set_packer_task( core::pack::task::PackerTaskOP packer_task ) {
	packer_task_ = packer_task;
}
void CoupledMover::set_loop_size( core::Size const loop_size ) {
	loop_size_ = loop_size;
}
void CoupledMover::set_perturber( kinematic_closure::perturbers::PerturberOP perturber ) {
	perturber_ = perturber;
}
void CoupledMover::set_backbone_mover ( std::string const & backbone_mover ) {
	backbone_mover_ = backbone_mover;
}
void CoupledMover::set_repack_neighborhood ( bool repack_neighborhood ) {
	repack_neighborhood_ = repack_neighborhood;
}

// getters
core::Size
CoupledMover::get_resnum() const {
	return resnum_;
}
bool
CoupledMover::get_randomize_resnum() const {
	return randomize_resnum_;
}
bool
CoupledMover::get_fix_backbone() const {
	return fix_backbone_;
}
core::Real
CoupledMover::get_rotation_std_dev() const {
	return rotation_std_dev_;
}
bool
CoupledMover::get_uniform_backrub() const {
	return uniform_backrub_;
}
core::Real
CoupledMover::get_temperature() const {
	return temperature_;
}
bool
CoupledMover::get_bias_sampling() const {
	return bias_sampling_;
}
bool
CoupledMover::get_bump_check() const {
	return bump_check_;
}
core::Size
CoupledMover::get_ligand_resnum() const {
	return ligand_resnum_;
}
core::Size
CoupledMover::get_ligand_jump_id() const {
	return ligand_jump_id_;
}
core::Real
CoupledMover::get_ligand_weight() const {
	return ligand_weight_;
}
core::Real
CoupledMover::get_rotation_magnitude() const {
	return rotation_magnitude_;
}
core::Real
CoupledMover::get_translation_magnitude() const {
	return translation_magnitude_;
}
protocols::simple_moves::ShortBackrubMoverOP
CoupledMover::get_short_backrub_mover() const {
	return short_backrub_mover_;
}
protocols::minimization_packing::BoltzmannRotamerMoverOP
CoupledMover::get_boltzmann_rotamer_mover() const {
	return boltzmann_rotamer_mover_;
}
protocols::rigid::RigidBodyPerturbMoverOP
CoupledMover::get_rigid_body_mover() const {
	return rigid_body_mover_;
}
core::scoring::ScoreFunctionOP
CoupledMover::get_score_fxn() const {
	return score_fxn_;
}
core::pack::task::PackerTaskOP
CoupledMover::get_packer_task() const {
	return packer_task_;
}
core::Size CoupledMover::get_loop_size() const {
	return loop_size_;
}
kinematic_closure::perturbers::PerturberOP CoupledMover::get_perturber() const {
	return perturber_;
}
const std::string & CoupledMover::get_backbone_mover() const {
	return backbone_mover_;
}
const bool & CoupledMover::get_repack_neighborhood() const {
	return repack_neighborhood_;
}

// Additional functions
std::string CoupledMover::get_name() const {
	return mover_name();
}
std::string CoupledMover::mover_name() {
	return "CoupledMover";
}
std::string CoupledMoverCreator::keyname() const {
	return CoupledMover::mover_name();
}
protocols::moves::MoverOP
CoupledMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new CoupledMover );
}

// XML functions
void CoupledMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Runs a Coupled Move. There are no XML attributes.", attlist );
}
void CoupledMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CoupledMover::provide_xml_schema( xsd );
}

} // moves
} // protocols
