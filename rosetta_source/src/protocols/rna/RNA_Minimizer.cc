// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief protocols that are specific to RNA_Minimizer
/// @detailed
/// @author Rhiju Das


#include <protocols/rna/RNA_Minimizer.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>

#include <protocols/viewer/viewers.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>


#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>

//Minimizer stuff
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

//Packer stuff
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

#include <numeric/conversions.hh>

// External library headers


//C++ headers
#include <vector>
#include <string>
#include <sstream>
// AUTO-REMOVED #include <fstream>

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif


// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

//Auto Headers
#include <core/chemical/AtomType.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/annotated_sequence.hh>
#include <numeric/xyz.functions.hh>



using namespace core;
using basic::T;

static basic::Tracer TR( "protocols.rna.rna_minimizer" ) ;

namespace protocols {
namespace rna {

RNA_Minimizer::RNA_Minimizer():
	Mover(),
	deriv_check_( false ),
	use_coordinate_constraints_( true ),
	skip_o2star_trials_( false ),
	vary_bond_geometry_( false )
{
	Mover::type("RNA_Minimizer");
	allow_insert_.dimension( 0 );
	if ( basic::options::option[ basic::options::OptionKeys::score::weights ].user() ) {
		scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( basic::options::option[ basic::options::OptionKeys::score::weights ]() );
	} else {
		scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::RNA_HIRES_WTS );
	}
}

/// @details  Apply the RNA full atom minimizer.
///
void RNA_Minimizer::apply( core::pose::Pose & pose	)
{

	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	TR << "Check it! SEQUENCE " << pose.sequence() << std::endl;

	time_t pdb_start_time = time(NULL);

	//	protocols::viewer::add_conformation_viewer( pose.conformation(), "minimizer", 400, 400 );

	if (pose.constraint_set()->has_constraints() ){
		scorefxn_->set_weight( atom_pair_constraint, 1.0 );
		scorefxn_->set_weight( angle_constraint, 1.0 );
	}

	if ( vary_bond_geometry_ ){
		scorefxn_->set_weight( rna_bond_geometry, 1.0 );
	}

	(*scorefxn_)(pose);
	//scorefxn_->show( std::cout, pose );

	/////////////////////////////////////////////////////
	//Need to be careful with coordinate constraints -- remove them later???
	/////////////////////////////////////////////////////
	scoring::constraints::ConstraintSetOP save_pose_constraints = pose.constraint_set()->clone();
	if (use_coordinate_constraints_) core::scoring::constraints::add_coordinate_constraints( pose );

	/////////////////////////////////////////////////////
	// minimzer
	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, deriv_check_, deriv_check_ );
	//	MinimizerOptions options( "dfpmin_armijo_nonmonotone", dummy_tol, use_nblist, deriv_check_, deriv_check_ );
	options.nblist_auto_update( true );

	kinematics::MoveMap mm;

	setup_movemap( mm, pose );

// 	mm.set_bb( false );
// 	mm.set_chi( false );
// 	mm.set_jump( false );
// 	//id::TorsionID torsion_id( 1, id::CHI, 1 );
// 	id::TorsionID torsion_id( 1, id::BB, 3 );
// 	id::DOF_ID dof_id( pose.conformation().dof_id_from_torsion_id( torsion_id) );
// 	mm.set( dof_id, true );

	//Could/should be private data.
	Real const coord_cst_weight = 0.1;
	Size const rounds( option[ basic::options::OptionKeys::rna::minimize_rounds ] );
	Real const fa_rep_final( scorefxn_->get_weight( fa_rep ) );

	for (Size r = 1; r <= rounds; r++ ) {

		ScoreFunctionOP minimize_scorefxn_ = scorefxn_->clone();

		Real const suppress = static_cast<Real>(r)/rounds;
		minimize_scorefxn_->set_weight( fa_rep, fa_rep_final * suppress  );

		//		pose.dump_pdb( "before_o2star_trials.pdb" );
		if (!skip_o2star_trials_) o2star_trials( pose, minimize_scorefxn_ );
		//		pose.dump_pdb( "after_o2star_trials.pdb" );

		//scorefxn_->show( std::cout, pose );

		//Prevent explosions on first minimize.
		if (r == 1 && use_coordinate_constraints_) minimize_scorefxn_->set_weight( coordinate_constraint, coord_cst_weight );
		TR << "Minimizing..." << std::endl;

		minimizer.run( pose, mm, *minimize_scorefxn_, options );

		//		io::pdb::dump_pdb( pose, "minimize_round"+string_of( r)+".pdb" );
		minimize_scorefxn_->show( std::cout, pose );

	}

	time_t pdb_end_time = time(NULL);

	//scorefxn_->show( std::cout, pose );

	pose.constraint_set( save_pose_constraints );
	(*scorefxn_)( pose );

	scorefxn_->show( std::cout, pose );

	TR << "RNA minimizer finished in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;


}

std::string
RNA_Minimizer::get_name() const {
	return "RNA_Minimizer";
}

///////////////////////////////////////////////////////////////////////////////
// Make this its own Mover?
void
RNA_Minimizer::o2star_trials(
  core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP const & packer_scorefxn_,
	bool const do_pack_instead_of_rotamer_trials /* = false */ ) const
{

	//TR << "Repacking 2'-OH ... " << std::endl;

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();

	for (Size i = 1; i <= pose.total_residue(); i++) {
		if ( !pose.residue(i).is_RNA() ) continue;
		task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
		//		task->nonconst_residue_task(i).or_ex4( true );
		task->nonconst_residue_task(i).or_include_current( true );
		// How about bump check?
	}

	TR << "Orienting 2' hydroxyls..." << std::endl;

	if (do_pack_instead_of_rotamer_trials ){
		pack::pack_rotamers( pose, *packer_scorefxn_, task);
	} else {
		pack::rotamer_trials( pose, *packer_scorefxn_, task);
	}
}

void
RNA_Minimizer::setup_movemap( kinematics::MoveMap & mm, pose::Pose & pose ) {

	using namespace core::id;

	Size const nres( pose.total_residue() );

	if ( allow_insert_.size() == 0 ) {
		allow_insert_.dimension( nres );
		for (Size i = 1; i <= nres; i++ ) allow_insert_( i ) = true;
	}

	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( false );

	for  (Size i = 1; i <= nres; i++ )  {
		if ( !allow_insert_( i ) ) continue;
		mm.set_bb(  i, true );
		mm.set_chi( i, true );

		//Chi (glycosidic bond to base)
		//mm.set( pose.conformation().dof_id_from_torsion_id( id::TorsionID( i, id::CHI, 1 ) ), true );
		//2'-hydroxyl
		//		mm.set( pose.conformation().dof_id_from_torsion_id( id::TorsionID( i, id::CHI, 4 ) ), true );

	}

	for (Size n = 1; n <= pose.num_jump(); n++ ){
		Size const jump_pos1( pose.fold_tree().upstream_jump_residue( n ) );
		Size const jump_pos2( pose.fold_tree().downstream_jump_residue( n ) );
		if (allow_insert_( jump_pos1 ) || allow_insert_( jump_pos2 ) ) 	 mm.set_jump( n, true );
	}

	if ( vary_bond_geometry_ ) {
		// Let additional degrees of freedom vary -- but apply constraints to stay near
		// ideal bond lengths and angles!
		pose::Pose pose_reference;
		create_pose_reference( pose, pose_reference );
		vary_bond_geometry( mm, pose, pose_reference );
	}

}


//////////////////////////////////////////////////////////////////////////////////////////////
bool
check_in_bonded_list( core::id::AtomID const & atom_id1,
											core::id::AtomID const & atom_id2,
											utility::vector1< std::pair< core::id::AtomID, core::id::AtomID > > & bonded_atom_list ) {

	for (Size n = 1; n <= bonded_atom_list.size(); n++ ) {
		if( atom_id1 == bonded_atom_list[ n ].first && atom_id2 == bonded_atom_list[ n ].second ) return true;
		if( atom_id2 == bonded_atom_list[ n ].first && atom_id1 == bonded_atom_list[ n ].second ) return true;
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////
bool
check_in_bond_angle_list( core::id::AtomID const & atom_id1,
												 core::id::AtomID const & atom_id2,
												 core::id::AtomID const & atom_id3,
												 utility::vector1< std::pair< core::id::AtomID, std::pair< core::id::AtomID, core::id::AtomID > > > & bond_angle_list ) {

	for (Size n = 1; n <= bond_angle_list.size(); n++ ) {
		if( atom_id1 == bond_angle_list[ n ].first ) {
			if( atom_id2 == bond_angle_list[ n ].second.first && atom_id3 == bond_angle_list[ n ].second.second ) return true;
			if( atom_id3 == bond_angle_list[ n ].second.first && atom_id2 == bond_angle_list[ n ].second.second ) return true;
		}
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
add_bond_constraint( core::id::AtomID const & atom_id1,
										 core::id::AtomID const & atom_id2,
										 utility::vector1< std::pair< core::id::AtomID, core::id::AtomID > > & bonded_atom_list,
										 core::pose::Pose const & pose,
										 core::pose::Pose const & pose_reference,
										 core::scoring::constraints::ConstraintSetOP & cst_set )
{

	using namespace core::scoring;
	using namespace core::scoring::constraints;

	std::string const & atom_name1 = pose.residue( atom_id1.rsd() ).atom_name( atom_id1.atomno() );
	std::string const & atom_name2 = pose.residue( atom_id2.rsd() ).atom_name( atom_id2.atomno() );

	if ( !pose_reference.residue( atom_id1.rsd() ).has( atom_name1 ) ) return;
	if ( !pose_reference.residue( atom_id2.rsd() ).has( atom_name2 ) ) return;

	if ( !check_in_bonded_list( atom_id1, atom_id2,  bonded_atom_list ) ) {
		bonded_atom_list.push_back( std::make_pair( atom_id1, atom_id2 ) );

		Real const bond_length_sd_( 0.05 );

		Real const bond_length = ( pose_reference.residue( atom_id1.rsd() ).xyz( atom_name1 ) -
															 pose_reference.residue( atom_id2.rsd() ).xyz( atom_name2 ) ).length();


		FuncOP dist_harm_func_( new HarmonicFunc( bond_length, bond_length_sd_ ));

		cst_set->add_constraint( new AtomPairConstraint( atom_id1 ,
																										 atom_id2,
																										 dist_harm_func_,
																										 rna_bond_geometry ) )
;
		if ( false ) {
			TR << "PUTTING CONSTRAINT ON DISTANCE: " <<
				atom_id2.rsd() << " " << atom_name1 << "; "  <<
				atom_id1.rsd() << " " << atom_name2 << " "  <<
				bond_length <<
				std::endl;
		}
	}

}


//////////////////////////////////////////////////////////////////////////////////////////////
void
add_bond_angle_constraint( core::id::AtomID const & atom_id1,
													 core::id::AtomID const & atom_id2,
													 core::id::AtomID const & atom_id3,
													 utility::vector1< std::pair < core::id::AtomID, std::pair< core::id::AtomID, core::id::AtomID > > > & bond_angle_list,
													 core::pose::Pose const & pose,
													 core::pose::Pose const & pose_reference,
													 core::scoring::constraints::ConstraintSetOP & cst_set )
{

	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace numeric::conversions;

	if (atom_id2 == atom_id3) return;

	std::string const & atom_name1 = pose.residue( atom_id1.rsd() ).atom_name( atom_id1.atomno() );
	std::string const & atom_name2 = pose.residue( atom_id2.rsd() ).atom_name( atom_id2.atomno() );
	std::string const & atom_name3 = pose.residue( atom_id3.rsd() ).atom_name( atom_id3.atomno() );

	if ( !pose_reference.residue( atom_id1.rsd() ).has( atom_name1 ) ) return;
	if ( !pose_reference.residue( atom_id2.rsd() ).has( atom_name2 ) ) return;
	if ( !pose_reference.residue( atom_id3.rsd() ).has( atom_name3 ) ) return;

	if ( !check_in_bond_angle_list( atom_id1, atom_id2, atom_id3, bond_angle_list ) ) {
		bond_angle_list.push_back( std::make_pair( atom_id1, std::make_pair( atom_id2, atom_id3 ) ) );


		Real const bond_angle_sd_( radians( 5.0 ) );

		Real const bond_angle = angle_radians(
																								 pose_reference.residue( atom_id2.rsd() ).xyz( atom_name2 )  ,
																								 pose_reference.residue( atom_id1.rsd() ).xyz( atom_name1 )  ,
																								 pose_reference.residue( atom_id3.rsd() ).xyz( atom_name3 )
																					 );

		if (bond_angle < 0.001 ) TR << "WHAT THE HELL????????? " << std::endl;

		FuncOP angle_harm_func_( new HarmonicFunc( bond_angle, bond_angle_sd_ ));
		cst_set->add_constraint( new AngleConstraint(
																								 atom_id2 , atom_id1, atom_id3, angle_harm_func_,	rna_bond_geometry ) );

		if ( false ) {
			TR << "PUTTING CONSTRAINT ON ANGLE: " <<
				atom_id2.rsd() << " " << pose_reference.residue( atom_id2.rsd() ).atom_name( atom_id2.atomno() ) << "; "  <<
				atom_id1.rsd() << " " << pose_reference.residue( atom_id1.rsd() ).atom_name( atom_id1.atomno() ) << "; "  <<
				atom_id3.rsd() << " " << pose_reference.residue( atom_id3.rsd() ).atom_name( atom_id3.atomno() ) << " ==> "  << degrees( bond_angle ) << " " << degrees( bond_angle_sd_ ) <<
				std::endl;
		}

	}

}

////////////////////////////////////////////////
bool
check_if_really_connected(
  core::pose::Pose const & pose,
	core::id::AtomID const & atom_id1,
	core::id::AtomID const & atom_id2)
{
	if (  atom_id1.rsd() == atom_id2.rsd() ) return true;


	core::kinematics::tree::Atom const * atom1 ( & pose.atom_tree().atom( atom_id1 ) );
	core::kinematics::tree::Atom const * atom2 ( & pose.atom_tree().atom( atom_id2 ) );

	if ( atom1->parent() == atom2 ) return true;
	if ( atom2->parent() == atom1 ) return true;

	return false;
}

////////////////////////////////////////////////
bool
i_want_this_atom_to_move( conformation::Residue const & residue2, Size const & k )
{

	if (k > residue2.first_sidechain_atom() &&
			k != scoring::rna::first_base_atom_index( residue2 ) ) return false;

	if ( residue2.atom_type( k ).name() == "VIRT" ) {
		//		std::cout << "Is this virtual? " << residue2.atom_name( k ) << std::endl;
		return false;
	}

	return true;

}

bool
i_want_this_atom_to_move( pose::Pose const & pose, core::id::AtomID const & atom_id )
{

	return i_want_this_atom_to_move( pose.residue( atom_id.rsd() ) ,
																	 atom_id.atomno() );

}
//////////////////////////////////////////////////////////////////////////////
// Following has not (yet) been carefully debugged.
void
RNA_Minimizer::vary_bond_geometry(
																	core::kinematics::MoveMap & mm,
																	pose::Pose & pose,
																	pose::Pose const & pose_reference ) {

	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace numeric::conversions;

	ConstraintSetOP cst_set = pose.constraint_set()->clone();
	pose.constraint_set( cst_set );

	Size const nres( pose.total_residue() );

	std::map< AtomID, utility::vector1< AtomID > > lists_of_angle_bonded_atoms;

	for  (Size i = 1; i <= nres; i++ )  {

		if ( !allow_insert_( i ) ) continue;

		conformation::Residue const & residue( pose.residue( i )  );

		for (Size j = 1; j <= residue.natoms(); j++ ) {

			if ( !i_want_this_atom_to_move( residue, j ) ) continue;

			core::kinematics::tree::Atom const * current_atom ( & pose.atom_tree().atom( AtomID(j,i) ) );
			if ( current_atom->is_jump() ) continue;

			///////////////////
			core::kinematics::tree::Atom const * input_stub_atom1( current_atom->input_stub_atom1() );
			if ( !input_stub_atom1 ) continue;
			if ( !i_want_this_atom_to_move( pose, input_stub_atom1->id() ) ) continue;
			mm.set( DOF_ID( AtomID( j, i ), D ), true );
			if ( input_stub_atom1->is_jump() ) continue;

			core::kinematics::tree::Atom const * input_stub_atom2( current_atom->input_stub_atom2() );

			///////////////////
			if ( !input_stub_atom2 ) continue;
			if ( input_stub_atom2 == current_atom ) continue;
			if ( !i_want_this_atom_to_move( pose, input_stub_atom2->id() ) ) continue;
			mm.set( DOF_ID( AtomID( j, i ), THETA ), true );
			if ( input_stub_atom2->is_jump() ) continue;

			///////////////////
			core::kinematics::tree::Atom const * input_stub_atom3( current_atom->input_stub_atom3() );
			if ( !input_stub_atom3 ) continue;
			if ( !i_want_this_atom_to_move( pose, input_stub_atom3->id() ) ) continue;
			if ( input_stub_atom3 == current_atom ) continue;
			mm.set( DOF_ID( AtomID( j, i ), PHI ), true );

		}
	}

	utility::vector1< std::pair< AtomID, AtomID > >  bond_list;
	utility::vector1< std::pair< AtomID, std::pair< AtomID, AtomID > > > bond_angle_list;

 	for  (Size i = 1; i <= nres; i++ )  {

		//Go through all bonds in pose...
		if ( !allow_insert_( i ) ) continue;

		conformation::Residue const & residue( pose.residue( i )  );

		for (Size j = 1; j <= residue.natoms(); j++ ) {

			if ( !i_want_this_atom_to_move( residue, j ) ) continue;

			AtomID atom_id1( j, i );

			utility::vector1< AtomID >  nbrs(	pose.conformation().bonded_neighbor_all_res( atom_id1 ) );

			// Bond lengths.
			for ( Size n = 1; n <= nbrs.size(); n++ ) {

				AtomID const & atom_id2( nbrs[ n ] );

				conformation::Residue const & residue2( pose.residue( atom_id2.rsd() ) ) ;
				Size const & k( atom_id2.atomno() ) ;

				if ( ! check_if_really_connected( pose, atom_id1, atom_id2) ) continue;

				if ( i_want_this_atom_to_move( residue2, k ) )  {
					add_bond_constraint( atom_id1, atom_id2,
															 bond_list,
															 pose, pose_reference, cst_set );
				}

			}

			// Bond angles
			for ( Size m = 1; m <= nbrs.size(); m++ ) {
				AtomID const & atom_id2( nbrs[ m ] );

				conformation::Residue const & residue2( pose.residue( atom_id2.rsd() ) ) ;
				if ( ! check_if_really_connected( pose, atom_id1, atom_id2) ) continue;

				Size const & k( atom_id2.atomno() ) ;

				for ( Size n = 1; n <= nbrs.size(); n++ ) {
					AtomID const & atom_id3( nbrs[ n ] );

					conformation::Residue const & residue3( pose.residue( atom_id3.rsd() ) ) ;
					if ( ! check_if_really_connected( pose, atom_id1, atom_id3) ) continue;

					Size const & q( atom_id3.atomno() ) ;

					if ( i_want_this_atom_to_move( residue2, k ) &&
							 i_want_this_atom_to_move( residue3, q ) )  {
						add_bond_angle_constraint( atom_id1, atom_id2, atom_id3, bond_angle_list,
																			 pose, pose_reference, cst_set );
					}

				}
			}


		}

	}

	pose.constraint_set( cst_set );

}


///////////////////////////////////////////////////
void
apply_ideal_coordinates_for_alternative_pucker( pose::Pose const & pose, pose::Pose & pose_reference )
{

	using namespace core::scoring::rna;
	RNA_TorsionPotential const rna_torsion_potential;

	Real const DELTA_CUTOFF( rna_torsion_potential.delta_cutoff() );

	for (Size n = 1; n <= pose.total_residue(); n++ ) {
		Real const delta = pose.residue( n ).mainchain_torsion( DELTA );
		if ( delta > DELTA_CUTOFF ) {
			//apply_ideal_c2endo_sugar_coords( pose_reference, pose_reference, n );
			apply_ideal_c2endo_sugar_coords( pose_reference, n );
		}
	}

}

///////////////////////////////////////////////////
void
RNA_Minimizer::create_pose_reference(
  pose::Pose const & pose,
	pose::Pose & pose_reference )
{
	using namespace core::chemical;
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	core::pose::make_pose_from_sequence( pose_reference, pose.sequence(),	*rsd_set );
	apply_ideal_coordinates_for_alternative_pucker( pose, pose_reference );
}


///////////////////////////////////////////////////
 void
 RNA_Minimizer::set_score_function( core::scoring::ScoreFunctionOP const & scorefxn ){
	 scorefxn_ = scorefxn->clone();
}

} // namespace rna
} // namespace protocols
