// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license.
// (c) The Rosetta software is developed by the contributing members of the
// (c) Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
// (c) Questions about this can be addressed to University of Washington UW
// (c) TechTransfer, email: license@u.washington.edu.

/// @file   cov_hbs.cc
/// @brief  Sidechain conjugation to acryl amides
/// @author Andy Watkins (amw579@nyu.edu)

// includes
#include <iostream>
#include <fstream>
#include <string>

#include <devel/init.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/ncbb/util.hh>

#include <core/import_pose/import_pose.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/id/TorsionID.hh>
#include <core/id/types.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/func/SumFunc.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Numeric Headers
#include <numeric/conversions.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/RandomTorsionMover.hh>
#include <protocols/simple_moves/a3b_hbs/A3BHbsPatcher.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <numeric/random/random.hh>

using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::a3b_hbs;

using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace cov_creator {
	// pert options
	StringOptionKey const peptide_pdb ( "cov_creator::peptide_pdb" );
	IntegerOptionKey const peptide_res ( "cov_creator::peptide_res" );
	StringOptionKey const protein_pdb ( "cov_creator::protein_pdb" );
	IntegerOptionKey const protein_res ( "cov_creator::protein_res" );
	StringOptionKey const protein_pdb_res ( "cov_creator::protein_pdb_res" );
	StringOptionKey const peptide_pdb_res ( "cov_creator::peptide_pdb_res" );
}

class CovalentPeptidomimeticCreator : public Mover {
	
public:
	
	//default ctor
	CovalentPeptidomimeticCreator(): Mover( "CovalentPeptidomimeticCreator" ){}
	
	//default dtor
	virtual ~CovalentPeptidomimeticCreator(){}
	
	//methods
	
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "CovalentPeptidomimeticCreator"; }
	
};

typedef utility::pointer::shared_ptr< CovalentPeptidomimeticCreator > CovalentPeptidomimeticCreatorOP;
typedef utility::pointer::shared_ptr< CovalentPeptidomimeticCreator const >CovalentPeptidomimeticCreatorCOP;



int main ( int argc, char* argv[] )
{
try {
	//option[ chemical::patch_selectors ].push_back( "CTERM_AMIDATION" );
	
	option.add( cov_creator::peptide_pdb, "peptide pdb").def("");
	option.add( cov_creator::peptide_res, "res id of covalent res on peptide" ).def(5);
	option.add( cov_creator::protein_pdb, "protein pdb" ).def("");
	option.add( cov_creator::protein_res, "res id of covalent res on protein" ).def(1);
	option.add( cov_creator::protein_pdb_res, "protein pdb res" ).def("");
	option.add( cov_creator::peptide_pdb_res, "protein pdb res" ).def("");
	
	
	devel::init(argc, argv);
	
	option[ OptionKeys::out::file::write_pdb_link_records ].value( true );
	
	CovalentPeptidomimeticCreatorOP builder( new CovalentPeptidomimeticCreator() );
	protocols::jd2::JobDistributor::get_instance()->go( builder );
	
} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}

}

void
CovalentPeptidomimeticCreator::apply(
	core::pose::Pose & pose
) {
	using namespace core;
	using namespace utility;
	using namespace scoring;
	using namespace pose;
	using namespace core::chemical;
	using namespace conformation;
	using namespace func;
	using namespace constraints;
	
	using namespace core::id;
	using namespace core::pack;
	using namespace core::pack::task;
	
	ScoreFunctionOP scorefxn = get_score_function();
	
	//Get the residue set we are drawing from.
	core::chemical::ResidueTypeSetCOP residue_set_cap = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	Pose protein;
	Pose peptide;
	Size protein_res = 0;
	Size peptide_res = 0;
	
	ResidueType const & cys = residue_set_cap->name_map( "CYS:S-conjugated" );
	ResidueType const & vdp = residue_set_cap->name_map( "VDP:sidechain_electrophile_conjugated" );
	
	Residue res_cys( cys, true );
	Residue res_vdp( vdp, true );
	
	Size resi_cys = 0;
	Size resi_vdp = 0;
	
	if ( option[ cov_creator::protein_pdb ].user() ||
		 option[ cov_creator::peptide_pdb ].user() ) {
		pose.clear();
		
		import_pose::pose_from_pdb( protein, option[ cov_creator::protein_pdb ].value() );
		import_pose::pose_from_pdb( peptide, option[ cov_creator::peptide_pdb ].value() );

		runtime_assert( protein.conformation().num_chains() == 1 );

		pose = protein;
		
		protein_res = option[ cov_creator::protein_res ].value();
		peptide_res = option[ cov_creator::peptide_res ].value();
		
		pose.append_residue_by_jump( peptide.residue( 1 ), 1 );
		for ( Size i = 2; i <= peptide.total_residue(); ++i ) {
			pose.append_residue_by_bond( peptide.residue( i ), false );
			if ( i == peptide_res ) {
				if ( pose.residue( pose.total_residue() ).type().name3() == "CYS" ) {
					resi_cys = pose.total_residue();
					replace_pose_residue_copying_existing_coordinates(
							pose, pose.total_residue(), cys );
				} else {
					resi_vdp = pose.total_residue();
					replace_pose_residue_copying_existing_coordinates(
							pose, pose.total_residue(), vdp );
				}
			}
		}
		peptide_res = protein.total_residue() + peptide_res;
		
	} else {
		// Single complex provided via -s
		
		if ( option[ cov_creator::protein_res ].user() ||
		 option[ cov_creator::peptide_res ].user() ) {
			protein_res = option[ cov_creator::protein_res ].value();
			peptide_res = option[ cov_creator::peptide_res ].value();
		} else {
			//protein_res = pose.pdb_info()->pdb2pose( option[ cov_creator::protein_pdb_res ].value() );
			//peptide_res = pose.pdb_info()->pdb2pose( option[ cov_creator::peptide_pdb_res ].value() );
		}
	}
	
	kinematics::MoveMapOP pert_mm( new kinematics::MoveMap() );
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( i == protein_res || i == peptide_res ) {
			if ( pose.residue( i ).type().name3() == "CYS" ) {
				resi_cys = i;
				replace_pose_residue_copying_existing_coordinates(
						pose, i, cys );
			} else {
				resi_vdp = i;
				replace_pose_residue_copying_existing_coordinates(
						pose, i, vdp );
			}
			pert_mm->set_bb( i, false );
			pert_mm->set_chi( i, true );
		} else {
			pert_mm->set_bb( i, false );
			pert_mm->set_chi( i, false );
		}
	}
	
	Size vdp_connid = pose.residue( resi_vdp ).type().residue_connection_id_for_atom( pose.residue( resi_vdp ).atom_index( "CZ" ) );
	Size cys_connid = pose.residue( resi_cys ).type().residue_connection_id_for_atom( pose.residue( resi_cys ).atom_index( "SG" ) );
	
	// This seems absurd, but the pose/conformation only gives us const access to the residues.
	// (residues_[] is private)
	// So to set connection partners we have to replace residues...
	
	ResidueOP new_cys = ResidueFactory::create_residue( cys, pose.residue( resi_cys ), pose.conformation() );
	core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( resi_cys ), *new_cys, pose.conformation() );
	new_cys->residue_connection_partner( cys_connid, resi_vdp, vdp_connid );
	pose.conformation().replace_residue( resi_cys, *new_cys, false );
	ResidueOP new_vdp = ResidueFactory::create_residue( vdp, pose.residue( resi_vdp ), pose.conformation() );
	core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( resi_vdp ), *new_vdp, pose.conformation() );
	new_cys->residue_connection_partner( vdp_connid, resi_cys, cys_connid );
	pose.conformation().replace_residue( resi_vdp, *new_vdp, false );
	
	/*cys.residue_connection_partner( cys_connid, resi_vdp, vdp_connid );
	replace_pose_residue_copying_existing_coordinates( pose, resi_cys, cys );
	vdp.residue_connection_partner( cys_connid, resi_cys, cys_connid );
	replace_pose_residue_copying_existing_coordinates( pose, resi_vdp, vdp );
	*/
	//pose.residue( resi_vdp ).residue_connection_partner( vdp_connid, resi_cys, cys_connid);
	//pose.residue( resi_cys ).residue_connection_partner( cys_connid, resi_vdp, vdp_connid);
	pose.conformation().declare_chemical_bond( resi_cys, "SG", resi_vdp, "CZ" );
	
	
	
	using core::pack::task::operation::TaskOperationCOP;
	// create a monte carlo object for the pertubation phase
	moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *scorefxn, 1 ) );
	// create a task factory and task operations
	
	TaskFactoryOP pert_tf(new TaskFactory());
	pert_tf->push_back( operation::TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
	operation::OperateOnCertainResiduesOP pert_res( new operation::OperateOnCertainResidues() );
	utility::vector1< Size > indices;
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( pose.residue( i ).type().name3() == "CYS" ||pose.residue( i ).type().name3() == "VDP" ) {
			indices.push_back( i );
		}
	}
	pert_res->residue_indices( indices );
	pert_tf->push_back( pert_res );
	operation::RestrictToRepackingOP pert_rtrp( new operation::RestrictToRepacking() );
	pert_tf->push_back( pert_rtrp );
	
	// create a rotamer trials mover
	simple_moves::RotamerTrialsMoverOP pert_rt(new simple_moves::EnergyCutRotamerTrialsMover( scorefxn, pert_tf, pert_mc, 0.1 /*energycut*/ ) );
	protocols::simple_moves::MinMoverOP min_mover( new simple_moves::MinMover( pert_mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, true ) );
	
	core::id::AtomID const SG( cys.atom_index( "SG" ), resi_cys );
	core::id::AtomID const CZ( vdp.atom_index( "CZ"  ), resi_vdp );
	core::scoring::func::HarmonicFuncOP harm( new core::scoring::func::HarmonicFunc( 1.81, 0.02 ) );
	core::scoring::constraints::AtomPairConstraintOP bond( new core::scoring::constraints::AtomPairConstraint( SG, CZ, harm ) );
	pose.add_constraint( bond );
	
	if ( scorefxn->get_weight( atom_pair_constraint ) == 0.0 ) {
		scorefxn->set_weight( atom_pair_constraint, 1.0 );
	}
	std::cout << "Beginning cycles of rotamer trials and minimization on conjugated residues " << std::endl;
	for ( Size i = 1; i <= 10; ++i ) {
		//pert_rt->apply( pose );
		std::cout << "After RT " << i << ": score " << ( *scorefxn )( pose ) << std::endl;
		min_mover->apply( pose );
		std::cout << "After min " << i << ": score " << ( *scorefxn )( pose ) << std::endl;
	}
	
	
	/*
	 if ( protein.residue( protein_res ).type().name3() == "CYS" ) {
		// Replace with CYS:S-conjugated
		replace_pose_residue_copying_existing_coordinates(
				pose, protein_res, cys );
		// Now we append the OTHER residue.
		pose.append_residue_by_bond( res_vdp, true, 3, protein_res, 3 );
		resi_cys = protein_res;
		resi_vdp = protein_res + 1;
		
	} else if ( protein.residue( protein_res ).type().name3() == "VDP" ) {
		replace_pose_residue_copying_existing_coordinates(
				pose, protein_res, vdp );
		// Now we append the OTHER residue.
		pose.append_residue_by_bond( res_cys, true, 3, protein_res, 3 );
		resi_vdp = protein_res;
		resi_cys = protein_res + 1;
	}
	
	// Now we have to append and prepend the rest of the peptide.
	// Append residue by bond doesn't put the CONN3 partner near in seq but at end
	// so use that res number
	Size resi = pose.total_residue();//protein_res + 1;
	for ( Size i = peptide_res-1; i >= 1; --i) {
		pose.prepend_polymer_residue_before_seqpos( peptide.residue(i), resi, false );
	}
	
	resi = pose.total_residue();//protein_res + 1;
	for ( Size i = peptide_res+1; i <= peptide.total_residue(); ++i, ++resi ) {
		pose.append_polymer_residue_after_seqpos( peptide.residue(i), resi, false );
	}
	
	core::id::AtomID const C( cys.atom_index( "C" ), resi_cys );
	core::id::AtomID const CA( cys.atom_index( "CA" ), resi_cys );
	core::id::AtomID const CB( cys.atom_index( "CB" ), resi_cys );
	core::id::AtomID const SG( cys.atom_index( "SG" ), resi_cys );
	core::id::AtomID const CZ( vdp.atom_index( "CZ"  ), resi_vdp );
	core::id::AtomID const CE2( vdp.atom_index( "CE2" ), resi_vdp );
	core::id::AtomID const CD( vdp.atom_index( "CD"  ), resi_vdp );
	
	//pose.conformation().set_torsion_angle( C, CA, CB, SG, numeric::conversions::radians(106.5) );
	//pose.conformation().set_torsion_angle( CA, CB, SG, CZ, numeric::conversions::radians(-60.0) );
	//pose.conformation().set_torsion_angle( CB, SG, CZ, CE2, numeric::conversions::radians(180.0) );
	//pose.conformation().set_torsion_angle( SG, CZ, CE2, CD, numeric::conversions::radians(100.5) );
	//pose.conformation().set_bond_angle( CB, SG, CZ, numeric::conversions::radians(90) );
	pose.conformation().set_bond_length( SG, CZ, 1.8 );
	*/
	
	

}
