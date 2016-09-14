// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/farna/RNAIdealizeMover.cc
/// @brief Slowly accomodate movement from non-ideal to ideal bond lengths and angles by repeated minimization
/// @author Andy Watkins (amw579@nyu.edu)

// Unit headers
#include <protocols/farna/RNAIdealizeMover.hh>
#include <protocols/farna/RNAIdealizeMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/AtomID.hh>


#include <protocols/stepwise/modeler/rna/o2prime/O2PrimePacker.hh>

#include <protocols/simple_moves/MinMover.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

#include <iostream> 

static THREAD_LOCAL basic::Tracer TR( "protocols.farna.RNAIdealizeMover" );

namespace protocols {
namespace farna {

using namespace core;
using namespace core::id;
using namespace core::scoring;
using namespace core::conformation;
using namespace scoring::constraints;
using namespace scoring::func;	

RNAIdealizeMover::RNAIdealizeMover():
	protocols::moves::Mover( RNAIdealizeMover::class_name() ),
	iterations_( 100 ),
	noise_( false ),
	final_minimize_( false )
{}

RNAIdealizeMover::~RNAIdealizeMover(){}

void
RNAIdealizeMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	pose::Pose const & )
{
	iterations_ = tag->getOption< Size >( "iterations", iterations_ );
	noise_ = tag->getOption< bool >( "noise", noise_ );
	final_minimize_ = tag->getOption< bool >( "final_minimize", final_minimize_ );
}

protocols::moves::MoverOP
RNAIdealizeMover::clone() const
{
	return protocols::moves::MoverOP( new RNAIdealizeMover( *this ) );
}

protocols::moves::MoverOP
RNAIdealizeMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new RNAIdealizeMover );
}

std::string
RNAIdealizeMover::get_name() const
{
	return RNAIdealizeMover::class_name();
}

std::string
RNAIdealizeMover::class_name()
{
	return "RNAIdealizeMover";
}

void
RNAIdealizeMover::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );
}

std::ostream &
operator<<( std::ostream & os, RNAIdealizeMover const & mover )
{
	mover.show(os);
	return os;
}

/// @brief Add Gaussian noise to the xyz coordinates of each atom.
/// @details Intentionally uses set_xyz because we don't want changes to propagate
void
RNAIdealizeMover::perturb_pose( pose::Pose & pose ) const
{
	using namespace numeric;
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		for ( Size jj = 1; jj <= pose.residue_type( ii ).natoms(); ++jj ) {
			pose.conformation().set_xyz(
				AtomID( jj, ii ), 
				pose.conformation().xyz( AtomID( jj, ii ) ) 
					+ xyzVector< Real >( 
						random::rg().gaussian() * 0.02, 
						random::rg().gaussian() * 0.02, 
						random::rg().gaussian() * 0.02 ) );
		}
	}
}

bool	
connected_in_atom_tree(
	core::pose::Pose const & pose,
	AtomID const & atom_id1,
	AtomID const & atom_id2
) {
	core::kinematics::tree::AtomCOP atom1( pose.atom_tree().atom( atom_id1 ).get_self_ptr() );
	core::kinematics::tree::AtomCOP atom2( pose.atom_tree().atom( atom_id2 ).get_self_ptr() );
		
	if ( atom1->parent() == atom2 ) return true;
	if ( atom2->parent() == atom1 ) return true;
		
	return false;
}

void
RNAIdealizeMover::add_bond_constraint(
	AtomID const & atom_id1,
	AtomID const & atom_id2,
	core::pose::Pose & pose
) const {
	std::string const & atom_name1 = pose.residue_type( atom_id1.rsd() ).atom_name( atom_id1.atomno() );
	std::string const & atom_name2 = pose.residue_type( atom_id2.rsd() ).atom_name( atom_id2.atomno() );
	
	if ( !ref_pose_.residue_type( atom_id1.rsd() ).has( atom_name1 ) ) return;
	if ( !ref_pose_.residue_type( atom_id2.rsd() ).has( atom_name2 ) ) return;
	
	Real const bond_length_sd( 0.05 );
	Real const bond_length = ( ref_pose_.residue( atom_id1.rsd() ).xyz( atom_name1 ) -
		ref_pose_.residue( atom_id2.rsd() ).xyz( atom_name2 ) ).length();
	func::FuncOP dist_harm_func( new core::scoring::func::HarmonicFunc( bond_length, bond_length_sd ) );
	pose.add_constraint( ConstraintCOP( new AtomPairConstraint( atom_id1, atom_id2, dist_harm_func, atom_pair_constraint ) ) );
	
	TR.Trace << "PUTTING CONSTRAINT ON DISTANCE: " <<
	atom_id2.rsd() << " " << atom_name1 << "; "  <<
	atom_id1.rsd() << " " << atom_name2 << " "  <<
	bond_length << std::endl;
}
	
void
RNAIdealizeMover::add_bond_angle_constraint(
	AtomID const & atom_id1,
	AtomID const & atom_id2,
	AtomID const & atom_id3,
	core::pose::Pose & pose
) const {
	using namespace numeric::conversions;
	
	if ( atom_id2 == atom_id3 ) return;
	if ( atom_id1 == atom_id3 ) return;
	if ( atom_id1 == atom_id2 ) return;
	
	std::string const & atom_name1 = pose.residue_type( atom_id1.rsd() ).atom_name( atom_id1.atomno() );
	std::string const & atom_name2 = pose.residue_type( atom_id2.rsd() ).atom_name( atom_id2.atomno() );
	std::string const & atom_name3 = pose.residue_type( atom_id3.rsd() ).atom_name( atom_id3.atomno() );
	
	if ( !ref_pose_.residue_type( atom_id1.rsd() ).has( atom_name1 ) ) return;
	if ( !ref_pose_.residue_type( atom_id2.rsd() ).has( atom_name2 ) ) return;
	if ( !ref_pose_.residue_type( atom_id3.rsd() ).has( atom_name3 ) ) return;
	
	Real const bond_angle_sd_( radians ( 3.0 ) );
	Real const bond_angle = angle_radians(
		ref_pose_.residue( atom_id2.rsd() ).xyz( atom_name2 ),
		ref_pose_.residue( atom_id1.rsd() ).xyz( atom_name1 ),
		ref_pose_.residue( atom_id3.rsd() ).xyz( atom_name3 ) );
	
	if ( bond_angle < 0.001 ) TR.Warning << "WHAT THE HELL????????? " << std::endl;
	
	pose.add_constraint( ConstraintCOP( new AngleConstraint(
		atom_id2, atom_id1, atom_id3,
		FuncOP( new CircularHarmonicFunc( bond_angle, bond_angle_sd_ ) ),
		angle_constraint ) ) );
	
	TR.Trace << "PUTTING CONSTRAINT ON ANGLE: "
	<< atom_id2.rsd() << " " << atom_name2 << "; "
	<< atom_id1.rsd() << " " << atom_name1 << "; "
	<< atom_id3.rsd() << " " << atom_name3 << " ==> "
	<< degrees( bond_angle ) << " " << degrees( bond_angle_sd_ ) << std::endl;
}
	
/// @brief For dofs not in the AtomTree, constrain to ideal
/// @details If you don't do this, all the error ends up in the dofs not in the
/// AtomTree
void
RNAIdealizeMover::constrain_to_ideal( pose::Pose & pose ) const {
	
	for ( Size ii = 1; ii <= pose.size(); ++ii )  {
		if ( pose.residue_type( ii ).aa() == core::chemical::aa_vrt ) continue;
		
		chemical::ResidueType const & residue_type( pose.residue_type( ii ) );
		
		for ( Size jj = 1; jj <= residue_type.natoms(); ++jj ) {
			
			AtomID jj_atomid( jj, ii );
			utility::vector1< AtomID > nbrs( pose.conformation().bonded_neighbor_all_res( jj_atomid ) );
			
			for ( auto const & nbr : nbrs ) {
				if ( nbr.rsd() > pose.size() || nbr.rsd() < 1 ) continue;
				if ( connected_in_atom_tree( pose, jj_atomid, nbr ) ) continue;
				
				add_bond_constraint( jj_atomid, nbr, pose );
			}
			
			// Bond angles
			for ( auto const & nbr : nbrs ) {
				if ( nbr.rsd() > pose.size() || nbr.rsd() < 1 ) continue;				
				
				for ( auto const & ang_nbr : nbrs ) {
					if ( ang_nbr.rsd() > pose.size() || ang_nbr.rsd() < 1 ) continue;										
					if ( connected_in_atom_tree( pose, jj_atomid, nbr ) 
							&& connected_in_atom_tree( pose, jj_atomid, ang_nbr ) ) continue;
					
					add_bond_angle_constraint( jj_atomid, nbr, ang_nbr, pose );
				}
				
				utility::vector1< AtomID > nbrs2( pose.conformation().bonded_neighbor_all_res( nbr ) );
				
				for ( auto const & nbr2 : nbrs2 ) {
					if ( nbr2.rsd() > pose.size() || nbr2.rsd() < 1 ) continue;
					
					if ( connected_in_atom_tree( pose, jj_atomid, nbr ) 
							&& connected_in_atom_tree( pose, nbr, nbr2 ) ) continue;
										
					add_bond_angle_constraint( jj_atomid, nbr, nbr2, pose );
				}
			}
		}
	}		
}

void
RNAIdealizeMover::apply( pose::Pose & pose )
{
	TR << "Idealizing pose with " << iterations_ << " iterations." << std::endl;
	TR << "A final round of minimization " << ( final_minimize_ ? "will" : "won't" ) << " follow." << std::endl;
	TR << "Noise " << ( noise_ ? "is" : "isn't" ) << " added to starting coords." << std::endl;
	
	ScoreFunctionOP scorefxn = get_score_function();
	
	// Perturb input pose with Gaussian xyz noise
	if ( noise_ ) perturb_pose( pose );
	
	// Store original cst set and foldtree
	auto const cst_set = pose.constraint_set()->clone();
	auto const orig_ft = pose.fold_tree();
	
	// For every ba/bl in the pose, identify the ideal value. 
	// Easiest implementation first: create a pose of the same sequence
	// as a reference.
	pose::make_pose_from_sequence( ref_pose_, pose.annotated_sequence(), core::chemical::FA_STANDARD );
	
	scorefxn->set_weight( rna_bond_geometry, 1 );
	scorefxn->set_weight( atom_pair_constraint, 1 );
	//scorefxn->set_weight( angle_constraint, 1 );
	constrain_to_ideal( pose );
	
	// Add high-value bounded func coord constraints to every P/CA
	scorefxn->set_weight( coordinate_constraint, 10 );

	// 1. Add virt res
	bool we_added_the_virt = false;
	if ( pose.residue_type( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
		we_added_the_virt = true;
		// attach virt res there
		core::chemical::ResidueTypeSet const & rts = *pose.residue_type( 1 ).residue_type_set();
		ResidueOP new_res( ResidueFactory::create_residue( *( rts.get_representative_type_name3( "VRT" ) ) ) );
		pose.append_residue_by_jump( *new_res, pose.size() );
		// make the virt atom the root
		kinematics::FoldTree newF( pose.fold_tree() );
		newF.reorder( pose.size() );
		pose.fold_tree( newF );
	}
	
	utility::vector1< Size > all_res( pose.size() );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) all_res[ ii ] = ii;
	protocols::stepwise::modeler::rna::o2prime::O2PrimePacker o2prime_packer( pose, scorefxn, all_res );
	
	// 2. Add funcs
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( pose.residue_type( ii ).aa() == core::chemical::aa_vrt ) continue;
		
		Real const coord_sdev( 0.3 );
		Real const coord_tol(  0.3 );
		Size const my_anchor( pose.size() ); //anchor on virtual residue
		Residue const & rsd( pose.residue( ii ) );
		Size const atm_indexP = rsd.has( "P" ) ? rsd.atom_index( "P" ) : rsd.atom_index( "CA" );
		pose.add_constraint( ConstraintCOP( new CoordinateConstraint(
			AtomID( atm_indexP, ii ),
			AtomID( 1, my_anchor ), rsd.xyz( atm_indexP ),
			FuncOP( new FlatHarmonicFunc( 0.0, coord_sdev, coord_tol ) ) ) ) );
	}
	
	
	std::map< DOF_ID, Real > ref_dofs;
	std::map< DOF_ID, Real > start_dofs;
	std::set< Size > residues_with_idealizing_dofs;

	for ( Size jj = 1; jj <= ref_pose_.size(); ++jj ) {
		for ( Size kk = 1; kk <= ref_pose_.residue_type( jj ).natoms(); ++kk ) {
			
			// Broken for C5' of upper terminus
			if ( jj == ref_pose_.size() && kk == 5 ) continue;
			
			auto const atom_id  = AtomID( kk, jj );
			
			TR.Trace << ref_pose_.residue_type( jj ).name() << " " << atom_id << ref_pose_.residue_type( jj ).atom_name( kk ) << std::endl;
			
			kinematics::tree::AtomCOP current_atom( ref_pose_.atom_tree().atom_dont_do_update( atom_id ).get_self_ptr() );
			if ( !current_atom ) continue;
			if ( current_atom->is_jump() ) continue;
			if ( !current_atom->parent() ) continue;

			kinematics::tree::AtomCOP input_stub_atom1( current_atom->input_stub_atom1() );
			if ( !input_stub_atom1 ) continue;
			if ( input_stub_atom1->is_jump() ) continue;

			auto const dis_id   = DOF_ID( atom_id, D     );
			// if distance is within 0.05A we don't care.
			if ( std::abs( ref_pose_.conformation().dof( dis_id ) - pose.conformation().dof( dis_id ) ) > 0.05 ) {
				ref_dofs[ dis_id ] = ref_pose_.conformation().dof( dis_id );
				start_dofs[ dis_id ] = pose.conformation().dof( dis_id );
				residues_with_idealizing_dofs.insert( jj );
			}
			
 			kinematics::tree::AtomCOP input_stub_atom2( current_atom->input_stub_atom2() );
			if ( !input_stub_atom2 ) continue;
			if ( input_stub_atom2->is_jump() ) continue;
			if ( input_stub_atom2 == current_atom ) continue;
			
			auto const ang_id   = DOF_ID( atom_id, THETA );
			// if distance is within 5 degrees we don't care.
			if ( std::abs( ref_pose_.conformation().dof( ang_id ) - 
					numeric::nearest_angle_degrees( pose.conformation().dof( ang_id ), ref_pose_.conformation().dof( ang_id ) ) ) > 5 ) {
				ref_dofs[ ang_id ] = ref_pose_.conformation().dof( ang_id );
				start_dofs[ ang_id ] = pose.conformation().dof( ang_id );
				residues_with_idealizing_dofs.insert( jj );
			}
		}
	}
	

	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	//mm->set_bb( true );
	//mm->set_chi( true );
	mm->set_bb( false );
	mm->set_chi( false );
	for ( Size const res : residues_with_idealizing_dofs ) {
		mm->set_bb( true, res );
		mm->set_chi( true, res );
		
		if ( res > 1 ) {
			mm->set_bb( true, res - 1 );
			mm->set_chi( true, res - 1 );
		}
		if ( res < pose.size() ) {
			mm->set_bb( true, res + 1 );
			mm->set_chi( true, res + 1 );
		}
	}
	
	protocols::simple_moves::MinMoverOP minm( new protocols::simple_moves::MinMover( mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.01, true ) );

	Pose const basis_pose = pose;	
	for ( Size ii = 1; ii <= iterations_; ++ii ) {
		scorefxn->set_weight( fa_rep, scorefxn->get_weight( fa_rep ) + 0.01);

		TR << TR.Blue << "Idealize iteration " << ii << "," << " score " << ( *scorefxn )( pose ) << "." << std::endl;
		TR            << "RMS to starting: " << core::scoring::all_atom_rmsd( pose, basis_pose ) << "." << std::endl; 
		for ( auto const & elem : ref_dofs ) {
			Real const dof_diff = elem.second - start_dofs[ elem.first ];
			Real const target_dof = start_dofs[ elem.first ] + dof_diff * Real( ii ) / Real ( iterations_ );
			
			TR.Trace << "Changing " << ( elem.first.type() == THETA ? "BA" : "BL" ) << " from "
				<< start_dofs[ elem.first ] << " to " << target_dof << std::endl;
			
			pose.conformation().set_dof( elem.first, target_dof );
		}
		
		minm->apply( pose );

		// Rotamer trials
		o2prime_packer.apply( pose );

		
		std::stringstream ss;
		ss << "iter_" << ii << ".pdb";
		pose.dump_scored_pdb( ss.str(), *scorefxn );
	}
	
	
	// Remove added virt
	if ( we_added_the_virt ) {
		pose.conformation().delete_residue_slow( pose.fold_tree().root() );
	}
	
	// Restore original constraints
	pose.constraint_set( cst_set );
	pose.fold_tree( orig_ft );
	
	// Final minimization - optional
	if ( final_minimize_ ) minm->apply( pose );
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
RNAIdealizeMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new RNAIdealizeMover );
}

std::string
RNAIdealizeMoverCreator::keyname() const
{
	return RNAIdealizeMover::class_name();
}

} //protocols
} //farna

