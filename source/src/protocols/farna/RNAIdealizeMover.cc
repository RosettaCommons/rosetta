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
#include <core/pose/util.hh>
#include <core/pose/rna/RNA_SuiteName.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/AtomID.hh>

//#include <protocols/recces/sampler/rna/MC_RNA_KIC_Sampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_KIC_Sampler.hh>

#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMasterMover.hh>

#include <protocols/stepwise/modeler/rna/o2prime/O2PrimePacker.hh>
#include <protocols/stepwise/monte_carlo/mover/TransientCutpointHandler.hh>

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
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

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
	protocols::moves::Mover( RNAIdealizeMover::mover_name() ),
	iterations_( 100 ),
	noise_( false ),
	final_minimize_( false ),
	ang_significance_threshold_( 5 )
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
	ang_significance_threshold_ = tag->getOption< Real >( "ang_significance_threshold", ang_significance_threshold_ );
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

// XRW TEMP std::string
// XRW TEMP RNAIdealizeMover::get_name() const
// XRW TEMP {
// XRW TEMP  return RNAIdealizeMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP RNAIdealizeMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "RNAIdealizeMover";
// XRW TEMP }

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
				+ xyzVector< core::Real >(
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
	Real const bond_tol = 0.02;
	func::FuncOP dist_harm_func( new core::scoring::func::FlatHarmonicFunc( bond_length, bond_length_sd, bond_tol ) );
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

	Real const bond_angle_sd_( radians( 3.0 ) );
	Real const bond_angle = angle_radians(
		ref_pose_.residue( atom_id2.rsd() ).xyz( atom_name2 ),
		ref_pose_.residue( atom_id1.rsd() ).xyz( atom_name1 ),
		ref_pose_.residue( atom_id3.rsd() ).xyz( atom_name3 ) );
	Real const bond_angle_tol( radians( 4.0 ) );
	if ( bond_angle < 0.001 ) TR.Warning << "WHAT THE HELL????????? " << std::endl;

	pose.add_constraint( ConstraintCOP( new AngleConstraint(
		atom_id2, atom_id1, atom_id3,
		//FuncOP( new CircularHarmonicFunc( bond_angle, bond_angle_sd_ ) ),
		FuncOP( new FlatHarmonicFunc( bond_angle, bond_angle_sd_, bond_angle_tol ) ),
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
RNAIdealizeMover::constrain_to_ideal(
	pose::Pose & pose,
	utility::vector1< Size > const & bad_suite_res
) const {

	for ( Size ii = 1; ii <= pose.size(); ++ii )  {
		// If this residue or its successor is in bad suite res, no constraints
		// Update could constrain only PARTICULAR atoms for predecessor/successor situation
		// according to suite torsions
		if ( std::find( bad_suite_res.begin(), bad_suite_res.end(), ii ) != bad_suite_res.end() ) {
			continue;
		}
		if ( std::find( bad_suite_res.begin(), bad_suite_res.end(), ii + 1 ) != bad_suite_res.end() ) {
			continue;
		}
		if ( pose.residue_type( ii ).aa() == core::chemical::aa_vrt ) continue;

		chemical::ResidueType const & residue_type( pose.residue_type( ii ) );

		for ( Size jj = 1; jj <= residue_type.natoms(); ++jj ) {

			AtomID jj_atomid( jj, ii );
			utility::vector1< AtomID > nbrs( pose.conformation().bonded_neighbor_all_res( jj_atomid ) );

			for ( auto const & nbr : nbrs ) {
				if ( nbr.rsd() > pose.size() || nbr.rsd() < 1 ) continue;
				// Also constrain stuff in the atom tree. This way, we don't
				// move all error into atom tree dofs!
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
	using namespace core::chemical::rna;

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

	//scorefxn->set_weight( rna_bond_geometry, 1 );
	scorefxn->set_weight( atom_pair_constraint, 1 );
	//scorefxn->set_weight( angle_constraint, 1 );
	if ( handle_suites_ ) scorefxn->set_weight( dihedral_constraint, 1 );

	// 0. What's the suite situation?
	core::pose::rna::RNA_SuiteName suite_assignment;
	utility::vector1< Size > bad_suite_res;
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		auto assignment = suite_assignment.assign( pose, ii );
		if ( assignment.name == "!!" ) {
			bad_suite_res.push_back( ii );
			TR << "Residue " << ii << " is a suite outlier!" << std::endl;
			auto closest_suite = suite_assignment.closest_suite( pose, ii );
			TR << "Closest suite is " << closest_suite.name << std::endl;

			// Add seven? dihedral constraints
			// Don't bother with circular harmonics yet. Trialing flat harmonic behavior

			if ( handle_suites_ ) {
				AtomID aid1, aid2, aid3, aid4;
				if ( ii > 1 ) {
					pose.conformation().get_torsion_angle_atom_ids( TorsionID( ii - 1, id::BB, DELTA ), aid1, aid2, aid3, aid4 );
					pose.add_constraint( ConstraintCOP( new DihedralConstraint(
						aid1, aid2, aid3, aid4,
						FuncOP( new CircularHarmonicFunc( closest_suite.torsion[1], 0.3 ) ) ) ) );
					pose.conformation().get_torsion_angle_atom_ids( TorsionID( ii - 1, id::BB, EPSILON ), aid1, aid2, aid3, aid4 );
					pose.add_constraint( ConstraintCOP( new DihedralConstraint(
						aid1, aid2, aid3, aid4,
						FuncOP( new CircularHarmonicFunc( closest_suite.torsion[2], 0.3 ) ) ) ) );
					pose.conformation().get_torsion_angle_atom_ids( TorsionID( ii - 1, id::BB, ZETA ), aid1, aid2, aid3, aid4 );
					pose.add_constraint( ConstraintCOP( new DihedralConstraint(
						aid1, aid2, aid3, aid4,
						FuncOP( new CircularHarmonicFunc( closest_suite.torsion[3], 0.3 ) ) ) ) );
				}
				pose.conformation().get_torsion_angle_atom_ids( TorsionID( ii, id::BB, ALPHA ), aid1, aid2, aid3, aid4 );
				pose.add_constraint( ConstraintCOP( new DihedralConstraint(
					aid1, aid2, aid3, aid4,
					FuncOP( new CircularHarmonicFunc( closest_suite.torsion[4], 0.3 ) ) ) ) );
				pose.conformation().get_torsion_angle_atom_ids( TorsionID( ii, id::BB, BETA ), aid1, aid2, aid3, aid4 );
				pose.add_constraint( ConstraintCOP( new DihedralConstraint(
					aid1, aid2, aid3, aid4,
					FuncOP( new CircularHarmonicFunc( closest_suite.torsion[5], 0.3 ) ) ) ) );
				pose.conformation().get_torsion_angle_atom_ids( TorsionID( ii, id::BB, GAMMA ), aid1, aid2, aid3, aid4 );
				pose.add_constraint( ConstraintCOP( new DihedralConstraint(
					aid1, aid2, aid3, aid4,
					FuncOP( new CircularHarmonicFunc( closest_suite.torsion[6], 0.3 ) ) ) ) );
				pose.conformation().get_torsion_angle_atom_ids( TorsionID( ii, id::BB, DELTA ), aid1, aid2, aid3, aid4 );
				pose.add_constraint( ConstraintCOP( new DihedralConstraint(
					aid1, aid2, aid3, aid4,
					FuncOP( new CircularHarmonicFunc( closest_suite.torsion[7], 0.3 ) ) ) ) );

				// Resample

				if ( ! pose.fold_tree().is_cutpoint( ii ) ) {
					//protocols::stepwise::monte_carlo::mover::TransientCutpointHandler tch( ii );
					protocols::stepwise::monte_carlo::mover::TransientCutpointHandler tch( ii-1, ii, true );
					tch.put_in_cutpoints( pose );
				}

				//protocols::recces::sampler::rna::MC_RNA_KIC_Sampler sampler( PoseCOP( new Pose( pose ) ), ii - 1, ii, true );
				protocols::stepwise::sampler::rna::RNA_KIC_Sampler sampler( core::pose::PoseOP( new Pose( pose ) ), ii - 1, ii );
				sampler.init();
				++sampler;
				sampler.apply( pose );

				/*
				if ( ii > 1 && ii < pose.size() ) {
				// be smarter later
				using namespace protocols::stepwise::monte_carlo::mover;
				using namespace protocols::stepwise::monte_carlo::options;
				utility::vector1< Attachment > attachments;
				attachments.emplace_back( ii - 1, BOND_TO_PREVIOUS );
				attachments.emplace_back( ii + 1, BOND_TO_NEXT );
				StepWiseMove stepwise_move( ii, attachments, RESAMPLE_INTERNAL_LOCAL );

				StepWiseMonteCarloOptionsCOP options( new StepWiseMonteCarloOptions );
				StepWiseMasterMover swmm( scorefxn, options );
				swmm.apply( pose, stepwise_move );
				}
				// Add cutpoint and install good suite torsions
				if ( ! pose.fold_tree().is_cutpoint( ii ) ) {
				//protocols::stepwise::monte_carlo::mover::TransientCutpointHandler tch( ii );
				protocols::stepwise::monte_carlo::mover::TransientCutpointHandler tch( ii-1, ii, true );
				tch.put_in_cutpoints( pose );
				}
				if ( ii > 1 ) {
				pose.set_torsion( TorsionID( ii - 1, id::BB, DELTA ), closest_suite.torsion[1] );
				pose.set_torsion( TorsionID( ii - 1, id::BB, EPSILON ), closest_suite.torsion[2] );
				pose.set_torsion( TorsionID( ii - 1, id::BB, ZETA ), closest_suite.torsion[3] );
				}
				pose.set_torsion( TorsionID( ii, id::BB, ALPHA ), closest_suite.torsion[4] );
				pose.set_torsion( TorsionID( ii, id::BB, BETA ), closest_suite.torsion[5] );
				pose.set_torsion( TorsionID( ii, id::BB, GAMMA ), closest_suite.torsion[6] );
				pose.set_torsion( TorsionID( ii, id::BB, DELTA ), closest_suite.torsion[7] );
				*/
			}
		}
	}

	// TODO: don't constrain residues with bad suites nearly as much.....
	constrain_to_ideal( pose, bad_suite_res );

	// Add high-value bounded func coord constraints to every P/CA
	scorefxn->set_weight( coordinate_constraint, 10 );

	Pose const first_basis_pose = pose;
	for ( Size ii = 1; ii <= iterations_; ++ii ) {

		kinematics::MoveMapOP suite_mm( new kinematics::MoveMap );
		suite_mm->set_bb( true );
		suite_mm->set_chi( true );
		//mm->set_bb( false );
		//mm->set_chi( false );
		for ( Size const res : bad_suite_res ) {
			suite_mm->set_bb( true, res );
			suite_mm->set_chi( true, res );

			Size const min_res = res > 4 ? res - 3 : 1;
			Size const max_res = res < pose.size() - 3 ? res + 3 : pose.size();

			for ( Size jj = min_res; jj <= max_res; ++jj ) {
				suite_mm->set_bb( true, jj );
				suite_mm->set_chi( true, jj );
			}
		}

		protocols::simple_moves::MinMoverOP minm( new protocols::simple_moves::MinMover( suite_mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.001, true ) );

		TR << TR.Blue << "Suite-fix iteration " << ii << "," << " score " << ( *scorefxn )( pose ) << "." << std::endl;
		TR            << "RMS to starting: " << core::scoring::all_atom_rmsd( pose, first_basis_pose ) << "." << std::endl;


		minm->apply( pose );

		std::stringstream ss;
		ss << "suite_iter_" << ii << ".pdb";
		pose.dump_scored_pdb( ss.str(), *scorefxn );
	}
	scorefxn->set_weight( dihedral_constraint, 0 );


	// 1. Add virt res
	bool we_added_the_virt = false;
	if ( pose.residue_type( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
		we_added_the_virt = true;
		// attach virt res there
		ResidueOP new_res( ResidueFactory::create_residue( *core::pose::virtual_type_for_pose(pose) ) );
		pose.append_residue_by_jump( *new_res, pose.size() );
		ref_pose_.append_residue_by_jump( *new_res, ref_pose_.size() );
		// make the virt atom the root
		kinematics::FoldTree newF( pose.fold_tree() );
		newF.reorder( pose.size() );
		pose.fold_tree( newF );
		ref_pose_.fold_tree( newF );
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
		//Size const atm_indexP = rsd.has( "P" ) ? rsd.atom_index( "P" ) : rsd.atom_index( "CA" );

		// All coord cst
		for ( Size jj = 1; jj <= pose.residue_type( ii ).nheavyatoms(); ++jj ) {
			pose.add_constraint( ConstraintCOP( new CoordinateConstraint(
				AtomID( jj, ii ),
				AtomID( 1, my_anchor ), rsd.xyz( jj ),
				FuncOP( new FlatHarmonicFunc( 0.0, coord_sdev, coord_tol ) ) ) ) );
		}
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
			//TR << "dis " << ref_pose_.conformation().dof( dis_id ) << " vs " << pose.conformation().dof( dis_id ) << std::endl;
			if ( std::abs( ref_pose_.conformation().dof( dis_id ) - pose.conformation().dof( dis_id ) ) > 0.03 ) {
				TR << "The distance dof defining the position of res " << jj << " atom " << ref_pose_.residue_type( jj ).atom_name( kk ) << " matters." << std::endl;
				TR << ref_pose_.conformation().dof( dis_id ) << " vs " <<  pose.conformation().dof( dis_id ) << std::endl;
				ref_dofs[ dis_id ] = ref_pose_.conformation().dof( dis_id );
				start_dofs[ dis_id ] = pose.conformation().dof( dis_id );
				residues_with_idealizing_dofs.insert( jj );
			}

			kinematics::tree::AtomCOP input_stub_atom2( current_atom->input_stub_atom2() );
			if ( !input_stub_atom2 ) continue;
			if ( input_stub_atom2->is_jump() ) continue;
			if ( input_stub_atom2 == current_atom ) continue;

			// Arbitrary constraints: can't go crazy with sugar torsions
			// or base details, or this gets terrible.
			// TEMP only phos
			//if ( //ref_pose_.residue_type( jj ).atom_name( kk ) != " O2'" &&
			//ref_pose_.residue_type( jj ).atom_name( kk ) != " O3'" &&
			//ref_pose_.residue_type( jj ).atom_name( kk ) != " O5'" &&
			//ref_pose_.residue_type( jj ).atom_name( kk ) != " P  " ) {
			// ref_pose_.residue_type( jj ).atom_name( kk ) != " C3'" ) {
			// continue;
			//}
			if ( ref_pose_.residue_type( jj ).atom_name( kk ) != " C3'" &&
					ref_pose_.residue_type( jj ).atom_name( kk ) != " C4'" &&
					ref_pose_.residue_type( jj ).atom_name( kk ) != " O3'" ) {
				continue;
			}

			auto const ang_id   = DOF_ID( atom_id, THETA );
			// if distance is within 5 degrees we don't care.
			TR << "ang " << jj << ref_pose_.residue_type( jj ).atom_name( kk ) << " " << numeric::conversions::degrees( ref_pose_.conformation().dof( ang_id ) ) << " vs " << numeric::conversions::degrees( pose.conformation().dof( ang_id ) ) << std::endl;
			if ( std::abs( numeric::conversions::degrees( pose.conformation().dof( ang_id ) ) - numeric::conversions::degrees( ref_pose_.conformation().dof( ang_id ) ) ) > ang_significance_threshold_ ) {
				//if ( std::abs( numeric::conversions::degrees( ref_pose_.conformation().dof( ang_id ) ) -
				//     numeric::nearest_angle_degrees( pose.conformation().dof( ang_id ), ref_pose_.conformation().dof( ang_id ) ) ) > 5 ) {
				//if ( std::abs( numeric::conversions::degrees( ref_pose_.conformation().dof( ang_id ) ) -
				//  numeric::nearest_angle_degrees( numeric::conversions::degrees( pose.conformation().dof( ang_id ) ), numeric::conversions::degrees( ref_pose_.conformation().dof( ang_id ) ) ) ) > 5 ) {
				TR << "The angle dof defining the position of res " << jj << " atom " << ref_pose_.residue_type( jj ).atom_name( kk ) << " matters." << std::endl;
				TR << numeric::conversions::degrees( ref_pose_.conformation().dof( ang_id ) ) << " vs " << numeric::conversions::degrees(  pose.conformation().dof( ang_id ) ) << std::endl;
				ref_dofs[ ang_id ] = ref_pose_.conformation().dof( ang_id );
				start_dofs[ ang_id ] = pose.conformation().dof( ang_id );
				residues_with_idealizing_dofs.insert( jj );
			}
		}
	}


	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	mm->set_bb( true );
	mm->set_chi( true );
	//mm->set_bb( false );
	//mm->set_chi( false );
	for ( Size const res : residues_with_idealizing_dofs ) {
		mm->set_bb( true, res );
		mm->set_chi( true, res );

		Size const min_res = res > 4 ? res - 3 : 1;
		Size const max_res = res < pose.size() - 3 ? res + 3 : pose.size();

		for ( Size jj = min_res; jj <= max_res; ++jj ) {
			mm->set_bb( true, jj );
			mm->set_chi( true, jj );
		}
	}

	for ( Size const res : bad_suite_res ) {
		mm->set_bb( true, res - 1 );
		mm->set_bb( true, res );
	}

	protocols::simple_moves::MinMoverOP minm( new protocols::simple_moves::MinMover( mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.001, true ) );

	Pose const basis_pose = pose;
	for ( Size ii = 1; ii <= iterations_; ++ii ) {
		scorefxn->set_weight( fa_rep, scorefxn->get_weight( fa_rep ) + 0.01);

		TR << TR.Blue << "Idealize iteration " << ii << "," << " score " << ( *scorefxn )( pose ) << "." << std::endl;
		TR            << "RMS to starting: " << core::scoring::all_atom_rmsd( pose, basis_pose ) << "." << std::endl;
		for ( auto const & elem : ref_dofs ) {
			// Only interested in halving overall 'error'.
			Real const dof_diff = ( elem.second - start_dofs[ elem.first ] ) / 2.0;
			Real const target_dof = start_dofs[ elem.first ] + dof_diff * Real( ii ) / Real ( iterations_ );

			if ( elem.first.type() == THETA ) {
				TR.Debug << "Changing BA from "
					<< numeric::conversions::degrees( start_dofs[ elem.first ] ) << " to " << numeric::conversions::degrees( target_dof ) << std::endl;
			} else { // D
				TR.Debug << "Changing BL from "
					<< start_dofs[ elem.first ] << " to " << target_dof << std::endl;
			}
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

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP RNAIdealizeMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new RNAIdealizeMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP RNAIdealizeMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return RNAIdealizeMover::mover_name();
// XRW TEMP }

std::string RNAIdealizeMover::get_name() const {
	return mover_name();
}

std::string RNAIdealizeMover::mover_name() {
	return "RNAIdealizeMover";
}

void RNAIdealizeMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "iterations", xsct_positive_integer, "Number of iterations over which to spread idealization" )
		+ XMLSchemaAttribute( "noise", xsct_rosetta_bool, "Salt initial pose with some Gaussian noise" )
		+ XMLSchemaAttribute( "final_minimize", xsct_rosetta_bool, "Minimize the whole pose at once at the end" )
		+ XMLSchemaAttribute( "ang_significance_threshold", xsct_real, "Size of angle deviation to correct (degrees)" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Idealize a pose, moving its bond lengths and angles towards ideal values by repeated relaxation", attlist );
}

std::string RNAIdealizeMoverCreator::keyname() const {
	return RNAIdealizeMover::mover_name();
}

protocols::moves::MoverOP
RNAIdealizeMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RNAIdealizeMover );
}

void RNAIdealizeMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RNAIdealizeMover::provide_xml_schema( xsd );
}


} //protocols
} //farna

