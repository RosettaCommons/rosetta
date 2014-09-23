// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/enzdes/EnzdesMovers.hh
/// @brief a collection of movers that are used at different stages in enzyme design
/// @author Sinisa Bjelic sinibjelic@gmail.com, Florian Richter, floric@u.washington.edu

#include <protocols/enzdes/EnzdesMovers.hh>
#include <protocols/enzdes/EnzdesMoversCreator.hh>

// AUTO-REMOVED #include <protocols/enzdes/EnzdesFixBBProtocol.hh>
#include <protocols/enzdes/EnzdesBaseProtocol.hh>
#include <protocols/enzdes/enzdes_util.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/enzdes/EnzdesTaskOperations.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <protocols/ligand_docking/ligand_functions.hh>
// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <utility/tag/Tag.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
// AUTO-REMOVED
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <core/id/AtomID.hh>

#include <utility/string_util.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
#include <protocols/rosetta_scripts/util.hh>

// option key includes
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace enzdes {

static thread_local basic::Tracer mv_tr( "protocols.enzdes.PredesignPerturbMover" );
//PredesignPerturbMoverCreator

EnzdesConstraintReporter::EnzdesConstraintReporter() : utility::pointer::ReferenceCount(),
	ligand_seqpos_( 0 )
{}

EnzdesConstraintReporter::~EnzdesConstraintReporter() {}

EnzdesConstraintReporter::EnzdesConstraintReporter( EnzdesConstraintReporter const & src ) :
		utility::pointer::ReferenceCount(),
		constrained_lig_atoms_( src.constrained_lig_atoms_ ),
		constrained_nonligand_atoms_( src.constrained_nonligand_atoms_ ),
		ligand_seqpos_( src.ligand_seqpos_ )
{}


void
EnzdesConstraintReporter::find_constraints_to_ligand(
	core::pose::Pose const & pose
)
{
	using namespace core::scoring::constraints;
	
	//Get a set of constraints for a ligand.
	//It is a vector of ResidueConstraints.
	ConstraintSetCOP constraint_set;
	constraint_set=pose.constraint_set();
	
	// Each residue constraint is a map container
	// of Size and ConstraintsOP, defining a residue number
	//and constraints object respectively.
	for ( ConstraintSet::ResiduePairConstraintsIterator
			rpc_start = constraint_set->residue_pair_constraints_begin(ligand_seqpos_),
			rpc_end=constraint_set->residue_pair_constraints_end(ligand_seqpos_);
			rpc_start != rpc_end; ++rpc_start ) {
		//Constraints object is a pointer vector of constraint objects.
		// Each contraint object contains type of constraints.
		ConstraintsOP constraints;
		constraints=rpc_start->second;

		for ( Constraints::const_iterator
				iter     = constraints->begin(),
				iter_end = constraints->end();
				iter != iter_end; ++iter ) {
			mv_tr.Info <<(*iter)->type() << std::endl;
			if ( (*iter)->type() == "MultiConstraint" || (*iter)->type() == "AmbiguousConstraint" ) {
				add_constrained_atoms_from_multiconstraint( utility::pointer::dynamic_pointer_cast <MultiConstraint const > (*iter) );
			} else if ( ((*iter)->type()) == "AtomPair") {
				add_constrained_atoms_from_atom_pair_constraint( utility::pointer::dynamic_pointer_cast <AtomPairConstraint const > (*iter) );
			} // else, ignore this constraint
		}
	}
}

void
EnzdesConstraintReporter::add_constrained_atoms_from_multiconstraint(
  core::scoring::constraints::MultiConstraintCOP real_multi_constraint )
{
	core::scoring::constraints::ConstraintCOPs multi_constraint_members;
	multi_constraint_members=real_multi_constraint->member_constraints();
	
	for (core::scoring::constraints::ConstraintCOPs::const_iterator
			MC_it=multi_constraint_members.begin();
			MC_it!=multi_constraint_members.end(); MC_it++) {
		if ( ((*MC_it)->type()) == "AtomPair") {
			assert( utility::pointer::dynamic_pointer_cast <core::scoring::constraints::AtomPairConstraint const > ((*MC_it)) );
			add_constrained_atoms_from_atom_pair_constraint( utility::pointer::dynamic_pointer_cast <core::scoring::constraints::AtomPairConstraint const > ((*MC_it)) );
		} else if ( ((*MC_it)->type()) == "MultiConstraint" || ((*MC_it)->type()) == "AmbiguousConstraint" ) {
			assert( utility::pointer::dynamic_pointer_cast <core::scoring::constraints::MultiConstraint const > ((*MC_it)) );
			add_constrained_atoms_from_multiconstraint( (utility::pointer::dynamic_pointer_cast <core::scoring::constraints::MultiConstraint const > ((*MC_it))) );
		} // else, ignore this constraint
	}
}

void
EnzdesConstraintReporter::add_constrained_atoms_from_atom_pair_constraint(
	core::scoring::constraints::AtomPairConstraintCOP atom_pair_constraint
)
{
	// if ligand equals residue of atom1 then assign atom1
	// otherwise asign the other atom
	if ( atom_pair_constraint->atom(1).rsd() == ligand_seqpos_ ) {
		add_constrained_lig_atom( atom_pair_constraint->atom(1).atomno() );
		add_constrained_nonligand_atom( atom_pair_constraint->atom(2)    );
	} else {
		add_constrained_lig_atom( atom_pair_constraint->atom(2).atomno() );
		add_constrained_nonligand_atom( atom_pair_constraint->atom(1)    );
	}
}


void
EnzdesConstraintReporter::add_constrained_lig_atom(
	core::Size atom_no 
)
{
	for(utility::vector1< core::Size >::const_iterator
			it  = constrained_lig_atoms_.begin();
			it != constrained_lig_atoms_.end(); it++){
		if ( (*it) == atom_no) return;
	}
	constrained_lig_atoms_.push_back(atom_no);
	if ( mv_tr.Info.visible()) {
		mv_tr.Info << "Constrained ligand atom: " << std::endl;
		mv_tr.Info << atom_no << std::endl;
	}
}

void
EnzdesConstraintReporter::add_constrained_nonligand_atom(
	core::id::AtomID const & atid
)
{
	for(utility::vector1< core::id::AtomID > ::const_iterator
			it  = constrained_nonligand_atoms_.begin();
			it != constrained_nonligand_atoms_.end(); it++){
		if ( (*it) == atid ) return;
	}
	constrained_nonligand_atoms_.push_back( atid );
	if ( mv_tr.Info.visible()) {
		mv_tr.Info << "Constrained non-ligand atom: res " << std::endl;
		mv_tr.Info << atid.rsd() << " atom " << atid.atomno() << std::endl;
	}
}


moves::MoverOP
PredesignPerturbMoverCreator::create_mover() const
{
	return new PredesignPerturbMover;
}

std::string
PredesignPerturbMoverCreator::keyname() const
{
	return PredesignPerturbMoverCreator::mover_name();
}

std::string
PredesignPerturbMoverCreator::mover_name()
{
	return "PredesignPerturbMover";
}

//-------------PredesignPerturbMover-----------------//

PredesignPerturbMover::PredesignPerturbMover():
	protocols::rigid::RigidBodyPerturbMover()
{
	trans_magnitude(basic::options::option[basic::options::OptionKeys::enzdes::trans_magnitude]);
	rot_magnitude(basic::options::option[basic::options::OptionKeys::enzdes::rot_magnitude]);
	dock_trials_ = basic::options::option[basic::options::OptionKeys::enzdes::dock_trials];
}

PredesignPerturbMover::~PredesignPerturbMover(){}

void
PredesignPerturbMover::set_docking_pose(
	core::pose::Pose &pose,
	core::pack::task::PackerTaskCOP task )
{
	toolbox::match_enzdes_util::EnzdesCstCacheOP cst_cache( toolbox::match_enzdes_util::get_enzdes_observer( pose )->cst_cache() );

	for(core::Size i = 1, i_end = pose.total_residue(); i <= i_end; ++i) {
		if( task -> design_residue(i) && !(cst_cache && cst_cache->contains_position(i)) )
			positions_to_replace_.push_back(i);
	}

	protocols::toolbox::pose_manipulation::construct_poly_ala_pose(
		pose, positions_to_replace_, true, true, true );
}

void
PredesignPerturbMover::reinstate_pose(
	core::pose::Pose &pose,
	core::pose::Pose const &old_Pose )
{
	core::Size ires;
  for(core::Size i=1, i_end=positions_to_replace_.size(); i<=i_end; ++i) {
		ires=positions_to_replace_[i];
 		pose.replace_residue( ires, old_Pose.residue(ires), true);
	}
}

void
PredesignPerturbMover::find_constraints_to_ligand(
	core::pose::Pose const & pose
)
{
	constraint_reporter_.find_constraints_to_ligand( pose ); // delegate
}

core::Vector
PredesignPerturbMover::find_rotation_center( core::pose::Pose const &pose )
{
	if ( constraint_reporter_.constrained_lig_atoms().size() != 0 ) {
		return find_geometric_center_for_constrained_lig_atoms( pose );
	} else { // Use geometric center of atoms
		core::Vector geometric_center( 0.0 );
		core::conformation::Residue const & res( pose.residue( constraint_reporter_.ligand_resno() ) );
		for ( core::Size ii(1); ii <= res.natoms(); ++ii) {
			geometric_center += res.xyz(ii);
		}
		geometric_center /= res.natoms();
		return geometric_center;
	}
}

core::Vector
PredesignPerturbMover::find_geometric_center_for_constrained_lig_atoms(
	core::pose::Pose const &pose)
{
	assert( constraint_reporter_.constrained_lig_atoms().size() != 0 );

	core::Vector geometric_center( 0.0 );
	for( utility::vector1< core::Size >::const_iterator
			it = constraint_reporter_.constrained_lig_atoms().begin();
			it != constraint_reporter_.constrained_lig_atoms().end(); ++it ){
		geometric_center+=pose.residue( constraint_reporter_.ligand_resno() ).xyz(*it);
	}

	geometric_center /= constraint_reporter_.constrained_lig_atoms().size();

	return geometric_center;
}

void
PredesignPerturbMover::apply(
	core::pose::Pose & pose
)
{
  //make a poly ala of the designable
 	protocols::enzdes::EnzdesBaseProtocolOP enzprot = new protocols::enzdes::EnzdesBaseProtocol();
  core::pose::Pose org_Pose(pose);
  core::pack::task::PackerTaskOP task;
	if ( task_factory_ !=0 ) task = task_factory_->create_task_and_apply_taskoperations( pose );
  else
	task	= enzprot -> create_enzdes_pack_task( pose, true );
  set_docking_pose( pose, task );

  protocols::moves::MonteCarloOP MCpredock = new protocols::moves::MonteCarlo(
		pose,
		*(core::scoring::ScoreFunctionFactory::create_score_function( "enzdes_polyA_min" ) ),
		2.0 /* temperature, from RosettaLigand paper */);
  MCpredock->reset( pose );
  MCpredock->reset_counters();

  find_constraints_to_ligand(pose);

  //itereate through all constraints in the pose and check for the constrained atoms
	//ligand is always connected through the last jump
	mv_tr.Info << "starting predocking ... " << std::endl;
  for( core::Size i=1; i <= dock_trials_; ++i ){
    rot_center_= find_rotation_center(pose);

    core::kinematics::Jump flexible_jump = pose.jump( pose.num_jump() );
    core::kinematics::Stub downstream_stub = pose.conformation().downstream_jump_stub( pose.num_jump() );
    flexible_jump.set_rb_center( dir_, downstream_stub, rot_center_ );

    flexible_jump.gaussian_move( dir_, trans_mag_, rot_mag_ );
    pose.set_jump( pose.num_jump(), flexible_jump );
    MCpredock->boltzmann(pose);

  }
  MCpredock->show_counters();
  MCpredock->recover_low( pose );
  mv_tr.Info << "... done predocking" << std::endl;
  //put back the old pose
  reinstate_pose( pose, org_Pose );
}

void PredesignPerturbMover::set_ligand(core::Size res_no)
{
	constraint_reporter_.ligand_resno( res_no );
}


void
PredesignPerturbMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & pose)
{
	trans_magnitude( tag -> getOption< core::Real >( "trans_magnitude", 0.1 ) );
	rot_magnitude( tag -> getOption< core::Real >( "rot_magnitude", 2.0 ) );
	dock_trials_ = tag -> getOption< core::Size >( "dock_trials", 100 );
	constraint_reporter_.ligand_resno( (core::Size) pose.fold_tree().downstream_jump_residue( pose.num_jump() ));
  if( tag->hasOption("task_operations") ) task_factory_ = ( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
  else task_factory_ = NULL;

}

protocols::moves::MoverOP
PredesignPerturbMover::clone() const
{
	return protocols::moves::MoverOP( new PredesignPerturbMover( *this ) );
}

protocols::moves::MoverOP
PredesignPerturbMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new PredesignPerturbMover );
}

std::string
PredesignPerturbMover::get_name() const
{
	return "PredesignPerturbMover";
}

//-------RepackLigandSiteWithoutLigandMover----------//

RepackLigandSiteWithoutLigandMover::RepackLigandSiteWithoutLigandMover()
	: sfxn_(NULL), lig_seqpos_(0), enzcst_io_(NULL), calculate_silent_Es_(false)
{
	silent_Es_.clear();
}

RepackLigandSiteWithoutLigandMover::RepackLigandSiteWithoutLigandMover(
	core::scoring::ScoreFunctionCOP sfxn,
	bool calculate_silent_Es
	) : sfxn_(sfxn), lig_seqpos_(0), enzcst_io_(NULL), calculate_silent_Es_(calculate_silent_Es)
{
	silent_Es_.clear();
}

RepackLigandSiteWithoutLigandMover::~RepackLigandSiteWithoutLigandMover(){}

void
RepackLigandSiteWithoutLigandMover::set_sfxn(
	core::scoring::ScoreFunctionCOP sfxn )
{
	sfxn_ = sfxn;
}

void
RepackLigandSiteWithoutLigandMover::set_cstio(
	toolbox::match_enzdes_util::EnzConstraintIOCOP enzcst_io )
{
	enzcst_io_ = enzcst_io;
}

void
RepackLigandSiteWithoutLigandMover::set_calculate_silent_Es(
	bool calculate )
{
	calculate_silent_Es_ = calculate;
}

    
void
RepackLigandSiteWithoutLigandMover::set_separate_prt_ligand( bool separate_prt_ligand )
{
    separate_prt_ligand_ = separate_prt_ligand;
}
    
void
RepackLigandSiteWithoutLigandMover::apply(
	core::pose::Pose & pose )
{
	runtime_assert( sfxn_ != 0 );
	//tmp hack
	//the constraints can be a headache in this situation, so well completely take them out for now
	// the problem is that some constrained interactions can be covalent, and the EnzConstraintIO
	// object at the moment can't separately take out the constraints and the covalent connections
	core::scoring::ScoreFunctionOP tmpsfxn = sfxn_->clone();
	tmpsfxn->set_weight(core::scoring::coordinate_constraint, 0.0 );
	tmpsfxn->set_weight(core::scoring::atom_pair_constraint, 0.0 );
	tmpsfxn->set_weight(core::scoring::angle_constraint, 0.0 );
	tmpsfxn->set_weight(core::scoring::dihedral_constraint, 0.0 );
	sfxn_ = tmpsfxn;
	//tmp hack over

	if( lig_seqpos_ == 0 ){
		utility::vector1< core::Size > all_ligands( ligand_docking::get_ligand_seqpos( pose ) );
		if( all_ligands.size() != 1 ) utility_exit_with_message( "Pose has more or less than one ligand. This mover atm can only hadndle poses with one ligand.");
		lig_seqpos_ = all_ligands[1];
	}
	core::pose::PoseOP startpose;
	utility::vector1< core::Size > special_res;
	if( calculate_silent_Es_ ) startpose = new core::pose::Pose( pose );

	//1. if there are constraints between protein and ligand, we should take them out.
	if( enzcst_io_ ){
		if( calculate_silent_Es_ ) special_res = enzcst_io_->ordered_constrained_positions( pose );
		enzcst_io_->remove_constraints_from_pose( pose, true /* keep covalent*/, false /*fail on missing*/ );
	}
	core::Real start_score( (*sfxn_)( pose ) );
	if( enzcst_io_ ){
		enzcst_io_->remove_constraints_from_pose( pose, false /* keep covalent*/, false /*fail on missing*/ );
	}

	//2. construct the proper task
	DetectProteinLigandInterfaceOP detect_enzdes_interface = new DetectProteinLigandInterface();
  detect_enzdes_interface->set_design(false);
	core::pack::task::TaskFactory taskfactory;
  taskfactory.push_back( core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline() ) );
  taskfactory.push_back( detect_enzdes_interface);
	ptask_ = taskfactory.create_task_and_apply_taskoperations( pose );

	//3. shoot the ligand into space
	separate_protein_and_ligand( pose );
	(*sfxn_)( pose );
	//pose.dump_pdb( "rlswlm_after_rigid.pdb");

	//4. repack
	protocols::simple_moves::PackRotamersMoverOP packer = new protocols::simple_moves::PackRotamersMover(sfxn_, ptask_);
	packer->apply( pose );
	//pose.dump_pdb( "rlswlm_after_repack.pdb");

	//5. if requested, do more shit
	if( calculate_silent_Es_ ){
		//5a. Ediff
		core::Real end_score( (*sfxn_)( pose ) );
		//std::cerr << "start score=" << start_score <<", int score=" << int_score << ", end_score=" << end_score << std::endl;
		silent_Es_.push_back( core::io::silent::SilentEnergy("nlr_dE", end_score - start_score, 1, 12 ) );

		//5b. rmsd of repackable region
		ObjexxFCL::FArray1D_bool pack_region( ptask_->total_residue(), false );
		for( core::Size i = 1; i <= ptask_->total_residue(); ++i ){
			if( ptask_->residue_task( i ).being_packed() && pose.residue( i ).is_protein() ) pack_region( i ) = true;
		}
		core::Real pack_region_rmsd( core::scoring::rmsd_no_super_subset( *startpose, pose, pack_region, core::scoring::is_protein_sidechain_heavyatom ) );
		silent_Es_.push_back( core::io::silent::SilentEnergy("nlr_totrms", pack_region_rmsd, 1, 12 ) );

		//5c. rmsds of any eventual special residues
		for( core::Size i =1; i <= special_res.size(); ++i){
			if( pose.residue_type( special_res[i] ).is_ligand() ) continue;
			ObjexxFCL::FArray1D_bool pack_region( ptask_->total_residue(), false );
			pack_region( special_res[i] ) = true;
			core::Real spec_res_rmsd( core::scoring::rmsd_no_super_subset( *startpose, pose, pack_region, core::scoring::is_protein_sidechain_heavyatom ) );
			std::string title( "nlr_SR"+utility::to_string( i )+"_rms");
			silent_Es_.push_back( core::io::silent::SilentEnergy(title, spec_res_rmsd, 1, 14 ) );
		}
	}
	//6. finally let's remove the ligand from the pose for completeness
	if( enzcst_io_ ){
		if( enzcst_io_->contains_position( pose, lig_seqpos_ ) ){
			enzcst_io_->remove_position_from_template_res( pose, lig_seqpos_ );
			//for now we'll wipe out the cst cache
			toolbox::match_enzdes_util::get_enzdes_observer( pose )->set_cst_cache( NULL );
		}
	}
    // PG 21-05-2013
    // EnzFilters relies on this code but continues with the
    // pose without ligand - this is probably??? fine but
    // to run in debug the rms_util.tmpl.hh asserts the two
    // poses to be of equil size
    if( separate_prt_ligand_ ){
	pose.conformation().delete_residue_slow( lig_seqpos_ );
	(*sfxn_)( pose ); //make sure energies are up to date
	//pose.dump_pdb( "rlswlm_after_repackdel.pdb");
    }
}

std::string
RepackLigandSiteWithoutLigandMover::get_name() const {
	return "RepackLigandSiteWithoutLigandMover";
}


void
RepackLigandSiteWithoutLigandMover::separate_protein_and_ligand(
	core::pose::Pose & pose ) const
{
	protocols::rigid::RigidBodyTransMover trans_mover( pose, pose.fold_tree().get_jump_that_builds_residue( lig_seqpos_ ) );
	trans_mover.step_size( 666 );
	trans_mover.apply( pose );
}

core::pack::task::PackerTaskCOP
RepackLigandSiteWithoutLigandMover::get_ptask() const {
	return ptask_;
}

//-------UpdateEnzdesHeaderMover----------//

moves::MoverOP
UpdateEnzdesHeaderMoverCreator::create_mover() const
{
	return new UpdateEnzdesHeaderMover();
}

std::string
UpdateEnzdesHeaderMoverCreator::keyname() const
{
	return UpdateEnzdesHeaderMoverCreator::mover_name();
}

std::string
UpdateEnzdesHeaderMoverCreator::mover_name()
{
	return "UpdateEnzdesHeader";
}

UpdateEnzdesHeaderMover::UpdateEnzdesHeaderMover(){}

UpdateEnzdesHeaderMover::~UpdateEnzdesHeaderMover(){}

void
UpdateEnzdesHeaderMover::apply(
	core::pose::Pose & pose )
{
	enzutil::create_remark_headers_from_cstcache( pose );
}

std::string
UpdateEnzdesHeaderMover::get_name() const {
	return UpdateEnzdesHeaderMoverCreator::mover_name();
}

protocols::moves::MoverOP
UpdateEnzdesHeaderMover::clone() const
{
	return protocols::moves::MoverOP( new UpdateEnzdesHeaderMover( *this ) );
}

} //enzdes
} //protocols

