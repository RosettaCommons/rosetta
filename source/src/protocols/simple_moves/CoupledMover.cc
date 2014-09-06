// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CoupledMover.cc
/// @brief implementation of CoupledMover class and functions
/// @author Noah Ollikainen (nollikai@gmail.com)

// Unit headers
#include <protocols/simple_moves/CoupledMover.hh>
#include <protocols/simple_moves/CoupledMoverCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <core/conformation/Residue.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/kinematics/FoldTree.hh>

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


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.simple_moves.CoupledMover");
static numeric::random::RandomGenerator RG(415610);

namespace protocols {
namespace simple_moves {

std::string
CoupledMoverCreator::keyname() const {
	return CoupledMoverCreator::mover_name();
}

protocols::moves::MoverOP
CoupledMoverCreator::create_mover() const {
	return new CoupledMover;
}

std::string
CoupledMoverCreator::mover_name() {
	return "CoupledMover";
}

// default constructor
CoupledMover::CoupledMover() : protocols::moves::Mover()
{
	protocols::moves::Mover::type( "Coupled" );
	resnum_ = 0;
	randomize_resnum_ = false;
	fix_backbone_ = false;
	rotation_std_dev_ = 4.572016;
	temperature_ = 1.0;
	bias_sampling_ = true;
	ligand_resnum_ = 0;
	rotation_magnitude_ = 1;
	translation_magnitude_ = 0.1;
	short_backrub_mover_ = new protocols::simple_moves::ShortBackrubMover();
	short_backrub_mover_->set_rotation_std_dev(rotation_std_dev_);
	boltzmann_rotamer_mover_ = new protocols::simple_moves::BoltzmannRotamerMover();
	boltzmann_rotamer_mover_->set_temperature(temperature_);
	boltzmann_rotamer_mover_->set_bias_sampling(bias_sampling_);
	rigid_body_mover_ = new protocols::rigid::RigidBodyPerturbMover();
	rigid_body_mover_->rot_magnitude(rotation_magnitude_);
	rigid_body_mover_->trans_magnitude(translation_magnitude_);
}

/// @brief constructor that sets input pose, score function and packer task
CoupledMover::CoupledMover( 
	core::pose::PoseOP pose,
	core::scoring::ScoreFunctionOP score_fxn,
	core::pack::task::PackerTaskOP packer_task ) : protocols::moves::Mover()
{
	protocols::moves::Mover::type( "Coupled" );
	resnum_ = 0;
	randomize_resnum_ = false;
	fix_backbone_ = false;
	rotation_std_dev_ = 4.572016;
	temperature_ = 1.0;
	bias_sampling_ = true;
	ligand_resnum_ = 0;
	ligand_jump_id_ = 0;
	rotation_magnitude_ = 1;
	translation_magnitude_ = 0.1;
	short_backrub_mover_ = new protocols::simple_moves::ShortBackrubMover(pose);
	short_backrub_mover_->set_rotation_std_dev(rotation_std_dev_);
	boltzmann_rotamer_mover_ = new protocols::simple_moves::BoltzmannRotamerMover(score_fxn, packer_task);
	boltzmann_rotamer_mover_->set_temperature(temperature_);
	boltzmann_rotamer_mover_->set_bias_sampling(bias_sampling_);
	rigid_body_mover_ = new protocols::rigid::RigidBodyPerturbMover();
	rigid_body_mover_->rot_magnitude(rotation_magnitude_);
	rigid_body_mover_->trans_magnitude(translation_magnitude_);
	score_fxn_ = score_fxn;
	packer_task_ = packer_task;
}

/// @brief constructor that sets input pose, score function and packer task and ligand residue number
CoupledMover::CoupledMover( 
	core::pose::PoseOP pose,
	core::scoring::ScoreFunctionOP score_fxn,
	core::pack::task::PackerTaskOP packer_task,
	core::Size ligand_resnum ) : protocols::moves::Mover()
{
	protocols::moves::Mover::type( "Coupled" );
	resnum_ = 0;
	randomize_resnum_ = false;
	fix_backbone_ = false;
	rotation_std_dev_ = 4.572016;
	temperature_ = 1.0;
	bias_sampling_ = true;
	ligand_resnum_ = ligand_resnum;
	ligand_jump_id_ = pose->fold_tree().get_jump_that_builds_residue(ligand_resnum);
	rotation_magnitude_ = 1;
	translation_magnitude_ = 0.1;
	short_backrub_mover_ = new protocols::simple_moves::ShortBackrubMover(pose);
	short_backrub_mover_->set_rotation_std_dev(rotation_std_dev_);
	boltzmann_rotamer_mover_ = new protocols::simple_moves::BoltzmannRotamerMover(score_fxn, packer_task);
	boltzmann_rotamer_mover_->set_temperature(temperature_);
	boltzmann_rotamer_mover_->set_bias_sampling(bias_sampling_);
	rigid_body_mover_ = new protocols::rigid::RigidBodyPerturbMover(ligand_jump_id_, rotation_magnitude_, translation_magnitude_);
	score_fxn_ = score_fxn;
	packer_task_ = packer_task;
}

// copy constructor
CoupledMover::CoupledMover( CoupledMover const & rval ):
	protocols::moves::Mover( rval ),
	resnum_( rval.resnum_ ),
	randomize_resnum_( rval.randomize_resnum_ ),
	fix_backbone_( rval.fix_backbone_ ),
	rotation_std_dev_( rval.rotation_std_dev_ ),
	temperature_( rval.temperature_ ),
	bias_sampling_( rval.bias_sampling_ ),
	ligand_resnum_( rval.ligand_resnum_ ),
	ligand_jump_id_( rval.ligand_jump_id_ ),
	rotation_magnitude_( rval.rotation_magnitude_ ),
	translation_magnitude_( rval.translation_magnitude_ ),
	short_backrub_mover_( rval.short_backrub_mover_ ),
	boltzmann_rotamer_mover_( rval.boltzmann_rotamer_mover_ ),
	rigid_body_mover_( rval.rigid_body_mover_ ),
	score_fxn_( rval.score_fxn_ ),
	packer_task_( rval.packer_task_ )
{}

// destructor
CoupledMover::~CoupledMover(){}

// clone this object
CoupledMover::MoverOP
CoupledMover::clone() const
{
	return new protocols::simple_moves::CoupledMover( *this );
}

// create this type of object
CoupledMover::MoverOP
CoupledMover::fresh_instance() const
{
	return new protocols::simple_moves::CoupledMover();
}

void
CoupledMover::apply( core::pose::Pose & pose )
{
	if (resnum_ == 0) 
		randomize_resnum_ = true;
		
	if (randomize_resnum_) {
		utility::vector1< core::Size > move_positions;
		for(core::Size i = 1; i <= packer_task_->total_residue(); i++) {
			if ( packer_task_->pack_residue(i) || packer_task_->design_residue(i) )
				move_positions.push_back(i);
		}
		core::Size random_index = RG.random_range(1, move_positions.size());
		resnum_ = move_positions[random_index];
	}
	
	if (pose.residue(resnum_).is_protein()) {
		if (fix_backbone_ == false) {
			short_backrub_mover_->set_resnum(resnum_);
			short_backrub_mover_->apply(pose);
		}
	} else {
		rigid_body_mover_->apply(pose);
	}
	
	boltzmann_rotamer_mover_->set_resnum(resnum_);
	boltzmann_rotamer_mover_->apply(pose);
	
}

std::string
CoupledMover::get_name() const {
	return "CoupledMover";
}
	
// setters
void CoupledMover::set_resnum( core::Size resnum ) { resnum_ = resnum; }
void CoupledMover::set_randomize_resnum( bool randomize_resnum ) { randomize_resnum_ = randomize_resnum; }
void CoupledMover::set_fix_backbone( bool fix_backbone ) { fix_backbone_ = fix_backbone; }
void CoupledMover::set_rotation_std_dev( core::Real rotation_std_dev ) {
	rotation_std_dev_ = rotation_std_dev; 
	short_backrub_mover_->set_rotation_std_dev(rotation_std_dev);
}
void CoupledMover::set_input_pose( core::pose::PoseOP pose ) {
	short_backrub_mover_->set_input_pose(pose);
}
void CoupledMover::set_temperature( core::Real temperature ) {
	temperature_ = temperature;
	boltzmann_rotamer_mover_->set_temperature(temperature);
}
void CoupledMover::set_bias_sampling( bool bias_sampling ) {
	bias_sampling_ = bias_sampling;
	boltzmann_rotamer_mover_->set_bias_sampling(bias_sampling);
}
void CoupledMover::set_ligand_resnum( core::Size ligand_resnum ) {
	ligand_resnum_ = ligand_resnum;
}
void CoupledMover::set_ligand_jump_id( core::Size ligand_jump_id ) {
	ligand_jump_id_ = ligand_jump_id;
	rigid_body_mover_ = new protocols::rigid::RigidBodyPerturbMover(ligand_jump_id_, rotation_magnitude_, translation_magnitude_);
}
void CoupledMover::set_rotation_magnitude( core::Real rotation_magnitude ) {
	rotation_magnitude_ = rotation_magnitude;
	rigid_body_mover_->rot_magnitude(rotation_magnitude);
}
void CoupledMover::set_translation_magnitude( core::Real translation_magnitude ) {
	translation_magnitude_ = translation_magnitude;
	rigid_body_mover_->trans_magnitude(translation_magnitude);
}
void CoupledMover::set_short_backrub_mover( protocols::simple_moves::ShortBackrubMoverOP short_backrub_mover ) { short_backrub_mover_ = short_backrub_mover; }
void CoupledMover::set_boltzmann_rotamer_mover( protocols::simple_moves::BoltzmannRotamerMoverOP boltzmann_rotamer_mover ) { boltzmann_rotamer_mover_ = boltzmann_rotamer_mover; }
void CoupledMover::set_rigid_body_mover( protocols::rigid::RigidBodyPerturbMoverOP rigid_body_mover ) { rigid_body_mover_ = rigid_body_mover; }
void CoupledMover::set_score_fxn( core::scoring::ScoreFunctionOP score_fxn ) { score_fxn_ = score_fxn; }
void CoupledMover::set_packer_task( core::pack::task::PackerTaskOP packer_task ) { packer_task_ = packer_task; }

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
core::Real
CoupledMover::get_temperature() const {
	return temperature_;
}
bool
CoupledMover::get_bias_sampling() const {
	return bias_sampling_;
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
protocols::simple_moves::BoltzmannRotamerMoverOP
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

/// @brief parse xml
void
CoupledMover::parse_my_tag(
	TagCOP const,
	basic::datacache::DataMap &,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & pose)
{
	short_backrub_mover_ = new protocols::simple_moves::ShortBackrubMover(new core::pose::Pose(pose));
}

} // moves
} // protocols
