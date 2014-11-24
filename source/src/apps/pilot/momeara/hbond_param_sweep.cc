// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/momeara/hbond_param_sweep.cc
/// @brief Generate geometric conformations of hydrogen bonds

/// @author Matthew O'Meara



// Unit headers
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBondTypeManager.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/LinearPenaltyFunction.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/moves/ReportToDB.hh>
#include <protocols/simple_moves/MinMover.hh>


// Project Headers
#include <core/chemical/types.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/chemical/VariantType.hh>
#include <core/chemical/types.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>


#include <devel/init.hh>

#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>

// C++ Headers
#include <string>
#include <iostream>

using namespace std;
using namespace core;
using namespace chemical;
using namespace conformation;
using namespace kinematics;
using namespace pose;
using namespace scoring;
using namespace scoring::hbonds;
using namespace scoring::constraints;
using namespace protocols::moves;

using pose::PoseOP;
using id::AtomID;
using utility::vector1;


#ifdef DB_SQLITE3
OPT_1GRP_KEY( String, in, feature_schema )
OPT_1GRP_KEY( String, feat_stats, sample_source )
OPT_1GRP_KEY( String, out, feature_database )
#endif // DB_SQLITE3
OPT_1GRP_KEY( Integer, feat_stats, grid_steps )

OPT_1GRP_KEY( Real, toy_relax, constraint_slope )


namespace apps{
namespace pilot{
namespace momeara{


static thread_local basic::Tracer tr( "core.scoring.hbonds.HBondDatabase" );

	void register_options() {
		OPT(out::prefix);

		using namespace basic::options::OptionKeys;
#ifdef DB_SQLITE3
		NEW_OPT(in::feature_schema, "schema file for feature database", "schema.sql" );
		NEW_OPT(feat_stats::sample_source,"Identifier for the sample source", "toy" );
		NEW_OPT(out::feature_database, "sqlite3 database file", "feature_database.sl3" );
#endif //DB_SQLITE3
		NEW_OPT(feat_stats::grid_steps, "Number of points in each dimension of the sample grid", 6);
		NEW_OPT(toy_relax::constraint_slope, "LinearPenality constraint slope parameter", 30);

	}


class HBondConformation {

public:
  HBondConformation(){}


  ~HBondConformation(){}


  HBondConformation( HBondConformation const & src);

  HBondConformation const &
  operator=(HBondConformation const & src);


	/*
		Inputs: hydrogen bond donor, hydrogen bond acceptor, sequence separation
    Output: A minimal pose containing that makes the hydrogen bond
		        More precisely it has sequene 'xA..D', A and D are the acceptor
						and donor residues.  There fold tree has a jump from the acceptor
						residue to the donor hydrogen.

		This is a tool to diagnose the scoring overlab between hydrogen
		bonding and other scoring terms.

		Building the pose and the conformation of the residues is not
		robust:
		  - The default configuration of residues are not energetically
		favorable.  So either configurations need to be checked by hand or
		minimimization should be applied.  To apply minimization will
		probably require adding more residues to fully define torsion
		angles?
		  - Only seq_sep_other has been even partially tested.
	*/

  PoseOP
  make_conformation(
		    HBSeqSep const seq_sep,
		    Residue & don_rsd,
		    Size const hatm,
		    Residue & acc_rsd,
		    Size const & aatm,
		    RT const & hbond_geometry) const {
    ResidueTypeSetCAP RTS( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

		int ss;
		switch(seq_sep){
		case seq_sep_other: ss = 100; break;
		case seq_sep_M4:    ss = 6;   break;
		case seq_sep_M3:    ss = 7;   break;
		case seq_sep_M2:    ss = 8;   break;
		case seq_sep_PM1:   ss = 9;   break;
		case seq_sep_P2:    ss = 11;  break;
		case seq_sep_P3:    ss = 12;  break;
		case seq_sep_P4:    ss = 13;  break;
		default:
			tr << "Cannot build pose: Unrecognized HBSeqSep type '"
				 << HBondTypeManager::name_from_seq_sep_type( seq_sep ) << endl;
		}


    PoseOP pose( new Pose() );
    Residue dummy( RTS->name_map( "VRT" ), true/*no-op param*/ );

		pose->append_residue_by_bond(dummy);

		tr.Debug << "Pose1: " << *pose << endl;
		string const aatm_name(acc_rsd.atom_name(aatm));
		string const hatm_name(don_rsd.atom_name(hatm));
		Size don_pos(0), acc_pos(0);
		switch(seq_sep){
		case seq_sep_other:
			pose->append_residue_by_jump( acc_rsd, 1 );
			acc_pos = 2;
			pose->append_residue_by_jump( don_rsd, 2, aatm_name, hatm_name, true );
			don_pos = 3;
			break;
		case seq_sep_M4:
			pose->append_residue_by_jump( acc_rsd, 1 );
			acc_pos = 2;
			pose->append_residue_by_bond( *dummy.clone(), 2 );
			pose->append_residue_by_bond( *dummy.clone(), 3 );
			pose->append_residue_by_bond( *dummy.clone(), 4 );
			pose->append_residue_by_jump( don_rsd, 2, aatm_name, hatm_name, false );
			don_pos = 6;
			break;
		case seq_sep_M3:
			pose->append_residue_by_jump( acc_rsd, 1 );
			acc_pos = 2;
			pose->append_residue_by_bond( *dummy.clone(), 2 );
			pose->append_residue_by_bond( *dummy.clone(), 3 );
			pose->append_residue_by_jump( don_rsd, 2, aatm_name, hatm_name, false );
			don_pos = 5;
			break;
		case seq_sep_M2:
			pose->append_residue_by_jump( acc_rsd, 1 );
			acc_pos = 2;
			pose->append_residue_by_bond( *dummy.clone(), 2 );
			pose->append_residue_by_jump( don_rsd, 2, aatm_name, hatm_name, false );
			don_pos = 4;
			break;
		case seq_sep_PM1:
			pose->append_residue_by_jump( acc_rsd, 1 );
			acc_pos = 2;
			pose->append_residue_by_jump( don_rsd, 2, aatm_name, hatm_name, false );
			don_pos = 3;
			break;
		case seq_sep_P2:
			pose->append_residue_by_jump( don_rsd, 1 );
			don_pos = 2;
			pose->append_residue_by_bond( *dummy.clone(), 2 );
			pose->append_residue_by_jump( acc_rsd, 2, aatm_name, hatm_name, false );
			acc_pos = 4;
			break;
		case seq_sep_P3:
			pose->append_residue_by_jump( don_rsd, 1 );
			don_pos = 2;
			pose->append_residue_by_bond( *dummy.clone(), 2 );
			pose->append_residue_by_bond( *dummy.clone(), 3 );
			pose->append_residue_by_jump( acc_rsd, 2, aatm_name, hatm_name, false );
			acc_pos = 5;
			break;
		case seq_sep_P4:
			pose->append_residue_by_jump( don_rsd, 1 );
			don_pos = 2;
			pose->append_residue_by_bond( *dummy.clone(), 2 );
			pose->append_residue_by_bond( *dummy.clone(), 3 );
			pose->append_residue_by_bond( *dummy.clone(), 4 );
			pose->append_residue_by_jump( acc_rsd, 2, aatm_name, hatm_name, false );
			acc_pos = 6;
		}
		don_rsd.seqpos(don_pos);
		acc_rsd.seqpos(acc_pos);

		// Try to set reasonable residue conformations
		if( acc_rsd.nchi() > 0 ){
			pose->set_torsion( id::TorsionID( acc_pos, id::CHI, 1), 180);
		}
		if( don_rsd.nchi() > 0 ){
			pose->set_torsion( id::TorsionID( don_pos, id::CHI, 1), 180);
		}
		// this is for serine donor
		pose->set_torsion( id::TorsionID( don_pos, id::CHI, 2), 180);


		AtomID aatm_id(aatm, acc_pos);
		AtomID hatm_id(hatm, don_pos);

		AtomTree const & at = pose->conformation().atom_tree();
		Stub acc_AT_frame( at.atom( aatm_id ).get_stub());
		Stub don_AT_frame( at.atom( hatm_id ).get_stub());

		Stub acc_hbond_frame( build_acceptor_frame(acc_rsd, aatm));
		Stub don_hbond_frame( build_donor_frame(don_rsd, hatm));

		RT acc_AT_to_hbond( acc_AT_frame, acc_hbond_frame );
		RT don_hbond_to_AT( don_hbond_frame, don_AT_frame );

		// Make hbond jump
		//a*(acc_AT_to_hbond*hbond_geometry*don_AT_to_hbond^-1)=
		//(((a*acc_AT_to_hbond)*hbond_geometry)*don_hbond_to_AT)=
		Stub a,b,c;
		acc_AT_to_hbond.make_jump( acc_AT_frame, a );
		hbond_geometry.make_jump( a, b );
		Jump hbond_jump(acc_AT_frame, b);
		pose->set_jump(2, hbond_jump);

		tr.Debug << "acc_AT_frame:" << acc_AT_frame << endl;
		tr.Debug << "acc_AT_to_hbond" << acc_AT_to_hbond << endl;
		tr.Debug << "a:" << a << endl;
		tr.Debug << "hbond_geometry:" << hbond_geometry << endl;
		tr.Debug << "b:" << b << endl;

		return pose;
  }

	Stub
	build_donor_frame(
	  Residue const & rsd,
		Size const & hatm) const
	{
		Vector const & Hxyz( rsd.atom( hatm ).xyz() );
		Vector const & Dxyz( rsd.atom(rsd.atom_base( hatm )).xyz());
		Vector const & DBxyz( rsd.atom(rsd.atom_base(rsd.atom_base(hatm))).xyz());
		return Stub(Hxyz /*center*/, Hxyz, Hxyz - Dxyz, Dxyz - DBxyz);
	}

	Stub
	build_acceptor_frame(
    Residue const & rsd,
    Size const & aatm) const
	{
		Vector const & Axyz( rsd.atom( aatm ).xyz() );
		Vector const & ABxyz( rsd.atom( rsd.atom_base( aatm ) ).xyz() );
		Vector const & AB2xyz( rsd.atom(rsd.abase2(aatm)).xyz());
		Hybridization const & acc_hybrid( rsd.atom_type(aatm).hybridization());

		switch(acc_hybrid){
		case SP2_HYBRID:
			return Stub(Axyz, Axyz, ABxyz, AB2xyz );
			break;
		case SP3_HYBRID:
			return Stub(Axyz, Axyz, AB2xyz, ABxyz );
			break;
		case RING_HYBRID:
			return Stub(Axyz, Axyz, Real(0.5)*(ABxyz+AB2xyz), ABxyz);
			break;
		default:
			tr << "Cannot build acceptor frame: Unrecognized acceptor hybridization '"
				 << HBondTypeManager::name_from_hybridization_type( acc_hybrid ) << endl;
		}
		return Stub();
	}


		// Lock hbond atoms in place quite rigidy
	void
	relax_pose_around_hbond(
	  PoseOP pose,
    ScoreFunctionOP scfxn,
		Size const don_seqpos,
		Size const hatm,
		Size const acc_seqpos,
		Size const aatm,
		Real const well_depth=0,
		Real slope = 400 ) const
	{

		slope = basic::options::option[ basic::options::OptionKeys::toy_relax::constraint_slope ]();

		Residue don_rsd( pose->residue(don_seqpos));
		Residue acc_rsd( pose->residue(acc_seqpos));

		Vector const & Axyz( acc_rsd.atom( aatm ).xyz() );
		Vector const & ABxyz( acc_rsd.atom( acc_rsd.atom_base( aatm ) ).xyz() );
		Vector const & AB2xyz( acc_rsd.atom(acc_rsd.abase2(aatm)).xyz());
		Vector const & Hxyz( don_rsd.atom( hatm ).xyz() );
		Vector const & Dxyz( don_rsd.atom(don_rsd.atom_base( hatm )).xyz());

		AtomID const A(aatm,acc_rsd.seqpos());
		AtomID const AB(acc_rsd.atom_base(aatm), acc_rsd.seqpos());
		AtomID const AB2(acc_rsd.abase2(aatm), acc_rsd.seqpos());
		AtomID const H(hatm, don_rsd.seqpos());
		AtomID const D(don_rsd.atom_base(hatm), don_rsd.seqpos());

		ConstraintSetOP cst_set( pose->constraint_set()->clone() );
		Real half_width = 0; // -> have the penalty begin as soon as it has moved

		// use constraints to make assert that the hbond has not moved
		{ // acceptor -- hydrogen
			Real x_middle = distance(Axyz, Hxyz);
			FuncOP const func( new LinearPenaltyFunction(
												 x_middle, well_depth, half_width, slope));
			tr.Debug << "acceptor -- hydrogen constraint before: " << func->func( Axyz.distance(Hxyz)) << endl;
			cst_set->add_constraint( new AtomPairConstraint( A, H, func ) );
		}
		{ // acceptor -- donor
			Real x_middle = distance(Axyz, Dxyz);
			FuncOP const func( new LinearPenaltyFunction(
                         x_middle, well_depth, half_width, slope));
			tr.Debug << "acceptor -- donor constraint before: " << func->func( Axyz.distance(Dxyz)) << endl;
			cst_set->add_constraint( new AtomPairConstraint( A, D, func ) );
		}
		{ // acceptor base -- hydrogen
			Real x_middle = distance(ABxyz, Hxyz);
			FuncOP const func( new LinearPenaltyFunction(
                         x_middle, well_depth, half_width, slope));
			tr.Debug << "acceptor base -- hydrogen constraint before: " << func->func( ABxyz.distance(Hxyz)) << endl;
			cst_set->add_constraint( new AtomPairConstraint( AB, H, func ) );
		}
		{// acceptor base -- donor
			Real x_middle = distance(ABxyz, Dxyz);
			FuncOP const func( new LinearPenaltyFunction(
                         x_middle, well_depth, half_width, slope));
			tr.Debug << "acceptor base -- donor constraint before: " << func->func( ABxyz.distance(Dxyz)) << endl;
			cst_set->add_constraint( new AtomPairConstraint( AB, D, func ) );
		}
		{	// acceptor base 2 -- hydrogen
			Real x_middle = distance(AB2xyz, Hxyz);
			FuncOP const func( new LinearPenaltyFunction(
                         x_middle, well_depth, half_width, slope));
			tr.Debug << "acceptor base 2 -- hydrogen constraint before: " << func->func( AB2xyz.distance(Hxyz)) << endl;
			cst_set->add_constraint( new AtomPairConstraint( AB2, H, func ) );
		}
		{ // acceptor base 2 -- donor
			Real x_middle = distance(AB2xyz, Dxyz);
			FuncOP const func( new LinearPenaltyFunction(
                         x_middle, well_depth, half_width, slope));
			tr.Debug << "acceptor base 2 -- donor constraint before: " << func->func( AB2xyz.distance(Dxyz)) << endl;
			cst_set->add_constraint( new AtomPairConstraint( AB2, D, func ) );
		}
		// don't optimize for the constraints, just make the jump immobile
		// pose->constraint_set( cst_set );


		Real tolerance( .00001 );

		// assert that the contraints begin satisfied
		cst_set->setup_for_scoring(*pose, *scfxn);
		EnergyMap cst_emap_initial;
		cst_set->residue_pair_energy(acc_rsd, don_rsd, *pose, *scfxn, cst_emap_initial);
		EnergyMap cst_emap_alt_initial;
		pose->constraint_set()->residue_pair_energy(acc_rsd, don_rsd, *pose, *scfxn, cst_emap_alt_initial);
		tr.Debug << "Constraint energy before: " << cst_emap_initial.sum() << endl;
		tr.Debug << "Constraint energy before (alternate): " << cst_emap_alt_initial.sum() << endl;

		MoveMapOP movemap( new MoveMap()); // allow all degrees of freedom
		movemap->set_bb(true);
		movemap->set_chi(true);
		movemap->set_jump(false);
		string min_type("dfpmin");
		bool use_nb_list( true );  //default for Classic Relax
		protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover(movemap, scfxn, min_type, tolerance, use_nb_list ));

		Real score_before( (*scfxn)(*pose));
		tr.Debug << "score before minimizing:" << endl;
		scfxn->show(tr.Debug, *pose);
		Pose initial_pose = *pose; // deep copy

		min_mover->apply( *pose );

		tr.Debug << "all_atom_rmds before and after minimizing:" << all_atom_rmsd(initial_pose, *pose) << endl;
		tr.Debug << "score after minimizing:" << endl;
		scfxn->show(tr.Debug, *pose);
		tr << "Minimizing improved the score from " << score_before << " to " << (*scfxn)(*pose) << endl;

		// Check that the hydrogen bond geometry has not moved very far
		Real threshold = .01;
		Residue new_acc_rsd( pose->residue( acc_rsd.seqpos() ) );
		Residue new_don_rsd( pose->residue( don_rsd.seqpos() ) );

		// assert that the contraints begin satisfied
		cst_set->setup_for_scoring(*pose, *scfxn);
		EnergyMap cst_emap_after;
		cst_set->residue_pair_energy(new_acc_rsd, new_don_rsd, *pose, *scfxn, cst_emap_after);
		tr.Debug << "Constraint energy after: " << cst_emap_after.sum() << endl;
		//		assert(cst_emap_after.sum() < cst_emap_initial.sum() + tolerance );

		Vector const & new_Axyz( new_acc_rsd.atom( aatm ).xyz() );
		Vector const & new_ABxyz(new_acc_rsd.atom(
														 new_acc_rsd.atom_base(aatm)).xyz());
		Vector const & new_AB2xyz( new_acc_rsd.atom(new_acc_rsd.abase2(aatm)).xyz());
		Vector const & new_Hxyz( new_don_rsd.atom(hatm).xyz() );
		Vector const & new_Dxyz( new_don_rsd.atom(
                             new_don_rsd.atom_base(hatm)).xyz());
		{
			Real delta( std::abs(distance(Axyz,Hxyz) - distance(new_Axyz,new_Hxyz)));
			if(delta > threshold){
				stringstream message;
				message << "Relax pose around hbond moved the structure:";
				message << "delta(distance(Axyz,Hxyz))=" << delta;
				utility_exit_with_message( message.str() );
			}
		}
		{
			Real delta( std::abs(distance(Axyz,Dxyz) - distance(new_Axyz,new_Dxyz)));
			if(delta > threshold){
				stringstream message;
				message << "Relax pose around hbond moved the structure:";
				message << "delta(distance(Axyz,Dxyz))=" << delta;
				utility_exit_with_message( message.str() );
			}
		}
		{
			Real delta( std::abs(distance(ABxyz,Hxyz) - distance(new_ABxyz,new_Hxyz)));
			if(delta > threshold){
				stringstream message;
				message << "Relax pose around hbond moved the structure:";
				message << "delta(distance(ABxyz,Hxyz))=" << delta;
				utility_exit_with_message( message.str() );
			}
		}
		{
			Real delta( std::abs(distance(ABxyz,Dxyz) - distance(new_ABxyz,new_Dxyz)));
			if(delta > threshold){
				stringstream message;
				message << "Relax pose around hbond moved the structure:";
				message << "delta(distance(ABxyz,Dxyz))=" << delta;
				utility_exit_with_message( message.str() );
			}
		}
		{
			Real delta( std::abs(distance(AB2xyz,Hxyz) - distance(new_AB2xyz,new_Hxyz)));
			if(delta > threshold){
				stringstream message;
				message << "Relax pose around hbond moved the structure:";
				message << "delta(distance(AB2xyz,Hxyz))=" << delta;
				utility_exit_with_message( message.str() );
			}
		}
		{
			Real delta( std::abs(distance(AB2xyz,Dxyz) - distance(new_AB2xyz,new_Dxyz)));
			if(delta > threshold){
				stringstream message;
				message << "Relax pose around hbond moved the structure:";
				message << "Struture: " << pose->pdb_info()->name() << endl;
				message << "delta(distance(AB2xyz,Dxyz))=" << delta;
				utility_exit_with_message( message.str() );
			}
		}

	}

	void
	example_OH_pose() const
	{
    ResidueTypeSetCAP residue_type_set(
		  ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

		HBSeqSep seq_sep( seq_sep_other );

		Residue acc_res( *ResidueFactory::create_residue( residue_type_set->name_map( "GLY" ) ) );
		Size aatm( acc_res.atom_index( "O" ) ); // backbone oxygen
		Stub acc_frame; // identity frame

		Residue don_res( *ResidueFactory::create_residue( residue_type_set->name_map( "SER" ) ) );
		Size hatm( don_res.atom_index( "HG" ) );  // OH hydrogen
		PoseOP pose;


 		for (Size i=0; i <3; i++){
			for (Size j=0; j <3; j++){
				for (Size k=0; k <3; k++){

					Stub don_frame( Vector(i,j,k)/*center*/,
													Vector(0,0,1)/*a*/,
													Vector(0,1,0)/*b*/,
													Vector(1,0,0)/*c*/);
					RT hbond_geometry(acc_frame, don_frame);
					pose = make_conformation(seq_sep,
																	 don_res, hatm, acc_res, aatm,
																	 hbond_geometry);
					stringstream fname_stringstream;
					fname_stringstream << "/tmp/example_OH_conformation_cartesian_"
														 << i <<"-"<< j <<"-"<< k <<".pdb";
					tr.Debug << "outputing structure '" << fname_stringstream.str() << "'" << endl;
					pose->dump_pdb( fname_stringstream.str() );
				}
			}
		}

	}


	vector1< PoseOP >
	hbond_param_sweep(
		HBSeqSep const seq_sep,
		ResidueOP const & don_rsd,
		Size const hatm,
		ResidueOP const & acc_rsd,
		Size const & aatm,
		string const & pose_name_prefix,
		Size steps
#ifdef DB_SQLITE3
		,
		string const & /*schema_fname*/,
		string const & sample_source,
		string const & database_fname
#endif // DB_SQLITE3
) const {

		ScoreFunctionOP scfxn( get_score_function() );
		scfxn->set_weight( atom_pair_constraint, 1 );


#ifdef DB_SQLITE3
		ReportToDB report_to_db_mover(
			database_fname,
			"hbond_parameter_sweep",
			sample_source,
			scfxn);
#endif // DB_SQLITE3


		Stub acc_frame; // identity frame
		Stub don_frame;
		PoseOP pose;


		Real AHdist_min(MIN_R), AHdist_max(MAX_R); Size AHdist_steps(steps);
		Real cosBAH_min(MIN_xH), cosBAH_max(MAX_xH); Size cosBAH_steps(steps);
		//Real chi_min(MIN_xC), chi_max(MAX_xC); Size chi_steps(steps);

		Real cosAHD_min(MIN_xD), cosAHD_max(MAX_xD); Size cosAHD_steps(steps);

		// Variables to complete the rigid body motions of an hbond
		//		Real AHchi(0), HDchi(0);

		Size nmodels( AHdist_steps*cosBAH_steps*cosAHD_steps );
		tr << "creating " << nmodels << " hydrogen bond models:" << endl;
		tr << "\tAHdist: ["<<AHdist_min<<", " << AHdist_max << "], n_step: " << AHdist_steps << endl;
		tr << "\tcosBAH: ["<<cosBAH_min<<", " << cosBAH_max << "], n_step: " << cosBAH_steps << endl;
		//		tr << "\tchi:    ["<<chi_min<<", " << chi_max << "], n_step: " << chi_steps << endl;
		tr << "\tcosAHD: ["<<cosAHD_min<<", " << cosAHD_max << "], n_step: " << cosAHD_steps << endl;


		vector1< PoseOP > features;
		for( Real AHdist = AHdist_min; AHdist <= AHdist_max;
				 AHdist += (AHdist_max - AHdist_min)/static_cast<Real>(AHdist_steps-1)){
			for( Real cosBAH = cosBAH_min; cosBAH <= cosBAH_max;
					 cosBAH += (cosBAH_max - cosBAH_min)/static_cast<Real>(cosBAH_steps-1)){
				Real chi = 0;
				//				for( Real chi = chi_min; chi <= chi_max;
				//		 chi += (chi_max - chi_min)/static_cast<Real>(chi_steps-1)){
					for( Real cosAHD = cosAHD_min - .0001; cosAHD <= cosAHD_max;
							 cosAHD += (cosAHD_max - cosAHD_min)/static_cast<Real>(cosAHD_steps-1)){

						// In the frame at the acceptor
						Vector Hxyz(cosBAH*AHdist,
												sin(chi)*sqrt(1-cosBAH*cosBAH)*AHdist,
												cos(chi)*sqrt(1-cosBAH*cosBAH)*AHdist);

						Real ADdist, cosBAD;
						Vector Dxyz;
						if ( cosBAH < .5 ) {
							// make AHchi = -180

							// triangle AHD the AHD angle is
							ADdist = sqrt(AHdist*AHdist + 1 +2*AHdist*cosAHD);
							cosBAD = (cosBAH*AHdist + cosBAH*cosAHD + sqrt(1-cosBAH*cosBAH)*sqrt(1-cosAHD*cosAHD))/ADdist;
							//cosBAD = (cosBAH*AHdist + cos(acos(cosBAH)-acos(cosAHD)))/ADdist;

							Dxyz = Vector(cosBAD*ADdist,
														sin(chi)*sqrt(1-cosBAD*cosBAD)*ADdist,
														cos(chi)*sqrt(1-cosBAD*cosBAD)*ADdist);
						} else {
							// make AHchi = 0
							ADdist = sqrt(AHdist*AHdist + 1 +2*AHdist*cosAHD);
							cosBAD = (cosBAH*AHdist + cosBAH*cosAHD - sqrt(1-cosBAH*cosBAH)*sqrt(1-cosAHD*cosAHD))/ADdist;
							Dxyz = Vector(cosBAD*ADdist,
														sin(chi)*sqrt(1-cosBAD*cosBAD)*ADdist,
														cos(chi)*sqrt(1-cosBAD*cosBAD)*ADdist);
						}

						tr.Debug << "hbond geomerty:" << endl;
						tr.Debug << "\tAHdist: " << AHdist << endl;
						tr.Debug << "\tcosBAH: " << cosBAH << endl;
						tr.Debug << "\tchi:    " << chi << endl;
						tr.Debug << "\tcosAHD: " << cosAHD << endl;
						tr.Debug << "\tHxyz:   " << "(" << Hxyz.x() << "," << Hxyz.y() << "," << Hxyz.z() << ")" << endl;
						tr.Debug << "\tAHdist: " << AHdist << endl;
						tr.Debug << "\tcosBAD: " << cosBAD << endl;
						tr.Debug << "\tDxyz: " << Dxyz.x() << "," << Dxyz.y() << "," << Dxyz.z() << ")" << endl;

						Stub don_frame( Hxyz,
														Hxyz,
														Dxyz,
														2*Hxyz);
						RT hbond_geometry(acc_frame, don_frame);

						pose = make_conformation(seq_sep,
																		 *don_rsd, hatm, *acc_rsd, aatm,
																		 hbond_geometry);

						stringstream pose_name;
						pose_name << pose_name_prefix << "_";
						pose_name << don_rsd->type().name() << "_" << hatm << "_";
						pose_name << acc_rsd->type().name() << "_" << aatm << "_";
						pose_name << AHdist << "-" << cosBAH << "-" << chi << "-" << cosAHD;
						PDBInfoOP pdb_info( new PDBInfo( *pose ) );
						assert( pdb_info->nres() == pose->total_residue() );
						pdb_info->name( pose_name.str());
					  pose->pdb_info( pdb_info );
						// features.push_back(pose);

						relax_pose_around_hbond(pose, scfxn,
																		don_rsd->seqpos(), hatm,
																		acc_rsd->seqpos(), aatm);

						pose->dump_pdb( pose_name.str() + ".pdb" );
#ifdef DB_SQLITE3
						report_to_db_mover.apply(*pose);
#endif // DB_SQLITE3

					}
					//				}
			}
		}
		return features;
	}

	ResidueOP
	build_residue(
		ResidueTypeSetCAP residue_type_set,
    std::string const & res_name3,
		bool const uterm_cap = true,
		bool const lterm_cap = true){

		ResidueTypeCOPs possible_types(residue_type_set->name3_map( res_name3 ));
		for ( ResidueTypeCOPs::const_iterator
						type_iter = possible_types.begin(), type_end = possible_types.end();
					type_iter != type_end; ++type_iter ){
			if (uterm_cap != (*type_iter)->has_variant_type( UPPER_TERMINUS )){
				continue;
			}
			if (lterm_cap != (*type_iter)->has_variant_type( LOWER_TERMINUS )){
				continue;
			}

			tr << "Building residue with type '" << (*type_iter)->name() << "'" << endl;
			return new Residue( **type_iter, true );
		}
		stringstream message;
		message << "Unrecognized residue " << res_name3 << " with uterm_cap: '" << uterm_cap << "' lterm_cap: '" << lterm_cap << "'" << endl;
		utility_exit_with_message(message.str());
		return NULL;
	}

	vector1<PoseOP>
	run_example_hbond_sweep(
	  string const & don_res_name,
    string const & hatm_name,
		string const & acc_res_name,
		string const & aatm_name,
		Size steps){

		ResidueTypeSetCAP residue_type_set(
		  ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
		HBSeqSep seq_sep( seq_sep_other );

		ResidueOP acc_res( build_residue(residue_type_set, acc_res_name));
		Size aatm( acc_res->atom_index( aatm_name ) ); // backbone oxygen
		ResidueOP don_res( build_residue(residue_type_set, don_res_name));
		Size hatm( don_res->atom_index( hatm_name ) );  // OH hydrogen



		string const pose_name_prefix( "hbond");
		utility::vector1<PoseOP> features
			(hbond_param_sweep
			 (seq_sep, don_res, hatm, acc_res, aatm, pose_name_prefix, steps
#ifdef DB_SQLITE3
				,
				basic::options::option[ basic::options::OptionKeys::in::feature_schema ](),
				basic::options::option[ basic::options::OptionKeys::feat_stats::sample_source](),
				basic::options::option[ basic::options::OptionKeys::out::feature_database ]()
#endif  // DB_SQLITE3
				));

		stringstream fname_prefix;
		if ( basic::options::option[ basic::options::OptionKeys::out::prefix ].user() ){
			fname_prefix << basic::options::option[ basic::options::OptionKeys::out::prefix ]();
			fname_prefix << "/";
		} else {
			fname_prefix << "";
		}

		if(features.size() > 0){
			tr << "Outputing " << features.size() << " structures like: '"
				 << fname_prefix.str() << features[1]->pdb_info()->name() << ".pdb'" << endl;
		}

		for( utility::vector1<PoseOP>::const_iterator
					 feat_iter = features.begin(), end_iter = features.end();
				 feat_iter != end_iter; ++feat_iter ){
			stringstream fname_ss;
			fname_ss << fname_prefix.str()
							 << (*feat_iter)->pdb_info()->name() << ".pdb";
			(*feat_iter)->dump_pdb( fname_ss.str() );
		}
		return features;
	}

	void
	score_conformations(
	  vector1<PoseOP> features){

		ScoreFunctionOP scfxn( get_score_function() );

		stringstream scores;
		scfxn->show_line_headers(scores);
		scores << endl;
		for( utility::vector1<PoseOP>::const_iterator
					 feat_iter = features.begin(), end_iter = features.end();
				 feat_iter != end_iter; ++feat_iter ){
			(*scfxn)(**feat_iter);
			scores << (*feat_iter)->pdb_info()->name() << " ";
			scfxn->show_line(scores,**feat_iter);
			scores << endl;
		}
		tr << scores.str() << endl;
	}


#ifdef DB_SQLITE3

	void
	report_features_to_db(
		vector1<PoseOP> features,
		string const & /*schema_fname*/,
		string const & sample_source,
		string const & database_fname) const
	{
		ScoreFunctionOP scfxn( get_score_function() );
		ReportToDB report_to_db_mover(
			 database_fname,
			 "hbond_parameter_sweep",
			 sample_source,
			 scfxn);

		for( utility::vector1<PoseOP>::const_iterator
					 feat_iter = features.begin(), end_iter = features.end();
				 feat_iter != end_iter; ++feat_iter ){
			report_to_db_mover.apply(**feat_iter);
		}
	}
#endif // DB_SQLITE3


}; // HBConformation


} // momeara
} // pilot
} // apps


int
main( int argc, char ** argv )
{
	try {
	apps::pilot::momeara::register_options();

	devel::init(argc, argv);

	apps::pilot::momeara::HBondConformation hb_conformation;

	hb_conformation.example_OH_pose();

	Size steps( basic::options::option[ basic::options::OptionKeys::feat_stats::grid_steps]());

	vector1<PoseOP> features
		(hb_conformation.run_example_hbond_sweep
		 ("SER", "HG", "GLY", "O", steps));

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}



	return 0;
}
