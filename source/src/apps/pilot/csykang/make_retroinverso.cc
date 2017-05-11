// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: se ...
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/csykang/make_retroinverso.cc
/// @brief Generates a retroinverso version of a pose.
/// @author Christine S. Kang (csykang@uw.edu)

// devel headers
#include <devel/init.hh>

// utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/backtrace.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>

// core headers
#include <core/pose/Pose.hh>
#include <core/pose/util.tmpl.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/types.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/chemical/AA.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResidueLevelTask.hh>

// unit headers
#include <protocols/relax/FastRelax.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>


static THREAD_LOCAL basic::Tracer TR("make_retroinverso");


void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	option.add_relevant( in::file::s );
	//option.add_relevant( in::file::l );

}


/// @brief Given a reference pose, make a newpose where all of the L/D amino acids have been swapped and the sequence runs backwards
/// @details This changes the residue type, but does not fix the position of the atoms.
void
build_inverso_sequence(
	core::pose::PoseCOP refpose,
	core::pose::PoseOP newpose
){
	newpose->clear();
	for ( core::Size i(refpose->total_residue());
			i>=1;
			--i // equivalent to i=i-1
			) {
		core::chemical::ResidueTypeSetCOP standard_residues( refpose->residue_type_set_for_pose( core::chemical::FULL_ATOM_t ) );
		core::chemical::ResidueTypeCOP curtype( standard_residues->get_representative_type_base_name( refpose->residue_type(i).base_name() ) );
		core::conformation::ResidueOP seq_rsd( core::conformation::ResidueFactory::create_residue( *curtype ) );
		core::conformation::ResidueOP flip_rsd( seq_rsd->clone_flipping_chirality( *standard_residues ));

		if ( i==refpose->total_residue() ) { //This is to add the first residue.
			newpose->append_residue_by_jump( *flip_rsd,1 );
		} else {
			newpose->append_residue_by_bond( *flip_rsd,true );
		}
	}
}


/// @brief Mirrors the pose so that the atom positions are correct
/// @details Does not change the residue identities.
void
update_atom_positions(
	core::pose::PoseOP newpose
){
	for ( core::Size r(1),rmax(newpose->total_residue()); // avoids calling the total_residue() repeatedly, instead calls it once and stores the value in "rmax"
			r<=rmax;
			++r
			) {
		for ( core::Size a(1),amax(newpose->residue(r).natoms()); //"newpose->residue(r)" synonymous to "(*newpose).residue(r)"
				a<=amax;
				++a
				) {
			numeric::xyzVector< core::Real > atompos(newpose->xyz(core::id::AtomID(a,r)));
			atompos.x(atompos.x()*(-1.0));
			newpose->set_xyz(core::id::AtomID(a,r),atompos);
		}
	}
}


/// @brief Sets pose conformation for the main chain.
void
set_inverso_conformation_main_chain(
	core::pose::PoseOP refpose,
	core::pose::PoseOP newpose
){
	if ( !( refpose->total_residue() == newpose->total_residue() ) ) utility::exit(__FILE__, __LINE__, "New pose must be the same number of residues as the reference pose!" );
	for ( core::Size r(1),rmax(newpose->total_residue()); // avoids calling the total_residue() repeatedly, instead calls it once and stores the value in "rmax"
			r<=rmax;
			++r
			) {
		newpose->set_phi(r,refpose->psi(rmax-r+1));
		newpose->set_psi(r,refpose->phi(rmax-r+1));
		if ( r<rmax ) newpose->set_omega(r,refpose->omega(rmax-r));
	}
}



/// @brief Sets pose conformation for the side chains.
void
set_inverso_conformation_side_chain(
	core::pose::PoseOP refpose,
	core::pose::PoseOP newpose
){
	if ( !( refpose->total_residue() == newpose->total_residue() ) ) utility::exit(__FILE__, __LINE__, "New pose must be the same number of residues as the reference pose!" );
	for ( core::Size r(1),rmax(newpose->total_residue()); // avoids calling the total_residue() repeatedly, instead calls it once and stores the value in "rmax"
			r<=rmax;
			++r
			) {
		// alternatively, a little bit faster would be to use " core::Size const rref(rmax-r+1); " where rref is the index for refpose
		for ( core::Size c(1),cmax(newpose->residue(r).nchi());
				c<=cmax;
				++c
				) {
			if ( c==1 ) {
				utility::vector1< core::Real > chi_atoms_vector(newpose->residue_type(r).chi_atoms(1));
				// core::Size new_value(newpose->residue_type(r).atom_index("C"));     made the next line a one-liner that includes this line inside of it.
				chi_atoms_vector[1] = newpose->residue_type(r).atom_index("C");

				utility::vector1< core::id::AtomID > atom_id_no(4);
				for ( core::Size a(1);
						a<=4;
						++a
						) {
					runtime_assert_string_msg( atom_id_no.size() == 4 , "Residues only have up to 4 chi." );
					atom_id_no[a] = core::id::AtomID(chi_atoms_vector[a],r);
				}
				newpose->conformation().set_torsion_angle(atom_id_no[1],atom_id_no[2],atom_id_no[3],atom_id_no[4],numeric::conversions::radians(refpose->chi(1,rmax-r+1)));
			} else {
				newpose->set_chi(c,r,refpose->chi(c,rmax-r+1));
			}
		}
	}
}


/// @brief Aligns the newpose and refpose
void
align_conformations(
	core::pose::PoseOP refpose,
	core::pose::PoseOP newpose
){
	if ( !( refpose->total_residue() == newpose->total_residue() ) ) utility::exit(__FILE__, __LINE__, "New pose must be the same number of residues as the reference pose!" );
	core::id::AtomID_Map< core::id::AtomID > aID_map;
	core::pose::initialize_atomid_map(aID_map, *newpose, core::id::BOGUS_ATOM_ID);

	for ( core::Size r(1),rmax(newpose->total_residue());
			r<=rmax;
			++r
			) {
		core::Size const rref(rmax-r+1);
		core::Size const newN(newpose->residue_type(r).atom_index("N"));
		core::Size const refC(refpose->residue_type(r).atom_index("C"));
		aID_map[core::id::AtomID(newN,r)] = core::id::AtomID(refC,rref);
	}
	core::scoring::superimpose_pose(*newpose, *refpose, aID_map);
}


/// @brief Fixes the prolines
void
fix_prolines(
	core::pose::PoseOP newpose
){
	//Creating Rosetta's scorefunction.  This can be specified with flags by the user, and will otherwise default to the score12 scorefunction:
	core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function());
	core::scoring::ScoreFunctionOP sfxn_constrained( core::scoring::get_score_function()); //For final fastrelax; includes constraint linking N- and C-termini.
	sfxn_constrained->set_weight(core::scoring::atom_pair_constraint, 1.0);
	sfxn_constrained->set_weight(core::scoring::dihedral_constraint, 1.0);
	sfxn_constrained->set_weight(core::scoring::angle_constraint, 1.0);

	protocols::relax::FastRelax frlx( sfxn_constrained, 5 ); // will do 5 rounds of fast relaxation.

	core::kinematics::MoveMapOP move( new core::kinematics::MoveMap); // move map tells operator which parts of a structure (or what dof) can move and which must remain fixed
	move->set_bb(false); // Don't allow backbone motion
	move->set_chi(false); // Don't allow side-chain motion
	move->set_jump(false); // Don't bother with "jumps" (relative motion of different molecules relative one another)
	core::select::residue_selector::ResidueIndexSelectorOP ld_pro(new core::select::residue_selector::ResidueIndexSelector);
	for ( core::Size r(1), rmax(newpose->total_residue() );
			r<=rmax;
			++r
			) {
		if ( newpose->residue_type(r).aa() == core::chemical::aa_pro || newpose->residue_type(r).aa() == core::chemical::aa_dpr ) {
			move->set_chi(r, true); // Allow side-chain motion for pro and dpr
			ld_pro->append_index(r);
		}
	}
	core::select::residue_selector::NotResidueSelectorOP not_ld_pro( new core::select::residue_selector::NotResidueSelector( ld_pro) );

	core::pack::task::operation::PreventRepackingRLTOP prevent_repack(new core::pack::task::operation::PreventRepackingRLT);
	core::pack::task::operation::OperateOnResidueSubsetOP operate_res_subset(new core::pack::task::operation::OperateOnResidueSubset( prevent_repack, not_ld_pro ));
	core::pack::task::operation::RestrictToRepackingOP restrict_repack( new core::pack::task::operation::RestrictToRepacking );

	core::pack::task::TaskFactoryOP task_factory( new core::pack::task::TaskFactory );
	task_factory->push_back( operate_res_subset );
	task_factory->push_back( restrict_repack );

	frlx.set_task_factory( task_factory );

	frlx.set_movemap( move );

//	newpose->dump_pdb("temp1.pdb"); //DELETE LATER, only for testing
	frlx.apply( *newpose );
//	newpose->dump_pdb("temp2.pdb"); //DELETE LATER, only for testing
}


//OPT_KEY( Boolean, myboolean_opt ); this is how I would add a custom option

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init( argc, argv );
		register_options();


		//  NEW_OPT( myboolean_opt, brief, default ); when I add a custom option "OPT_KEY ...", NEW_OPT makes the custom option allowed for this app

		if ( ! option [ in::file::s ].user() ) {
			utility_exit_with_message("Please specify -s to specify the input PDB.");
		}

		core::pose::PoseOP refpose(core::import_pose::pose_from_file(option [in::file::s]()[1]));
		core::pose::PoseOP newpose(new core::pose::Pose);

		// Setup retroinverso sequence.

		build_inverso_sequence(refpose,newpose);
		update_atom_positions(newpose);
		set_inverso_conformation_main_chain(refpose,newpose);
		set_inverso_conformation_side_chain(refpose,newpose);
		align_conformations(refpose,newpose);
		fix_prolines(newpose);

		// Output to pdb.

		newpose->dump_pdb("output.pdb");

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
