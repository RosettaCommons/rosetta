// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_set/RotamerSet_.cc
/// @brief  amino acid rotamer set class implementation
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/pack/rotamer_set/rotamer_building_functions.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/dunbrack/ChiSet.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>

// Project Headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>

#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/kinematics/Stub.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/graph/Graph.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/database/open.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/constants.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// External headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray3D.hh>

// C++ headers
#include <string>
#include <iostream>
#include <fstream>


namespace core {
namespace pack {
namespace rotamer_set {


static THREAD_LOCAL basic::Tracer tt( "core.pack.rotamer_set.rotamer_building_functions", basic::t_info );


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DNA-specific methods
void
read_DNA_rotlib(
	utility::io::izstream & lib_stream,
	utility::vector1< DihedralSet* > & library
)
{
	using namespace std;
	const string delimiters(" ");
	string line;

	while ( !lib_stream.eof() ) {

		DihedralSet *dihedrals = new DihedralSet;

		// Read next entry
		getline( lib_stream, line );

		// Skip delimiters at beginning
		Size last_pos = line.find_first_not_of( delimiters, 0 );
		// Find first "non-delimiter"
		Size pos     = line.find_first_of( delimiters, last_pos );

		while ( string::npos != pos || string::npos != last_pos ) {
			// Found a token, add it to the vector.
			const string next_dih( line.substr( last_pos, pos-last_pos) );
			dihedrals->push_back( strtod( next_dih.c_str(), NULL) );

			// Skip delimiters.  Note the "not_of"
			last_pos = line.find_first_not_of( delimiters, pos );

			// Find next "non-delimiter"
			pos = line.find_first_of( delimiters, last_pos );
		}

		if ( dihedrals->size() < 7 ) {
			delete dihedrals;
		} else {
			library.push_back( dihedrals );
		}
	}
}

void
build_lib_dna_rotamers(
	utility::vector1< DihedralSet* > const & library,
	Size const resid,
	pose::Pose const & pose,
	chemical::ResidueTypeCOP concrete_residue,
	utility::vector1< conformation::ResidueOP > & rotamers
)
{
	using namespace basic;
	using namespace pose;
	using namespace conformation;
	using namespace optimization;
	using namespace id;
	using namespace std;
	using namespace scoring;
	using namespace constraints;

	bool const minimize_rotamers( true );
	Real const distance_cutoff( 2.0 );
	Real min_tolerance( 1e-2 );
	Real const pi( numeric::constants::d::pi );

	Real const std_length( .145 );
	Real const std_angle (  .38 );
	Real const dist_tol  (   .2 ); // A
	Real const angle_tol (   30 ); // degrees
	Real const rad2deg   ( 180/pi );

	Residue const & existing_residue( pose.residue( resid ) );
	debug_assert( existing_residue.is_NA() && concrete_residue->is_NA() );

	bool const simple_way( existing_residue.is_terminus() );

	// a residue of the correct type aligned to the existing backbone
	ResidueOP rot = ResidueFactory::create_residue( *concrete_residue, existing_residue, pose.conformation() );
	rot->set_chi( 1, existing_residue.chi(1) );

	if ( simple_way ) {
		rotamers.push_back( rot );
	} else {
		pose::Pose mini_pose;
		mini_pose.append_residue_by_bond( pose.residue( resid-1 ) );
		mini_pose.append_residue_by_bond( *rot );
		mini_pose.append_residue_by_jump( pose.residue( resid+1 ), 1 );

		Size const seqpos(2);
		Residue const & rsd( mini_pose.residue(seqpos  ) ), next_rsd( mini_pose.residue(seqpos+1) );

		debug_assert( rsd.is_NA() && mini_pose.residue(seqpos-1).is_NA() );

		Vector const P_O3( ( next_rsd.xyz("P") - rsd.xyz("O3'") ) );
		Real const target_distance( P_O3.length() );

		Real const P_bond_angle (angle_of( next_rsd.xyz("O5'"), next_rsd.xyz("P")  ,  rsd.xyz("O3'")  ));
		Real const O3_bond_angle(angle_of( next_rsd.xyz("P")  ,      rsd.xyz("O3'"),  rsd.xyz("C3'")  ));

		// tt << "P-O3 dist.=" << target_distance << ", P_bond_angle=" << P_bond_angle << ", O3_bond_angle=" << O3_bond_angle << std::endl;

		ConstraintSetOP cst_set( new ConstraintSet() );
		func::FuncOP f_target_distance_std_length( new func::HarmonicFunc( target_distance, std_length ) );
		cst_set->add_constraint(
			ConstraintCOP( ConstraintOP( new AtomPairConstraint(
			AtomID(rsd.atom_index("O3'"),seqpos),
			AtomID(next_rsd.atom_index("P"),seqpos+1),
			f_target_distance_std_length
			) ) )
		);

		func::FuncOP f_P_bond_angle_std_angle( new func::HarmonicFunc( P_bond_angle, std_angle ) );
		cst_set->add_constraint(
			ConstraintCOP( ConstraintOP( new AngleConstraint(
			AtomID(next_rsd.atom_index("O5'"),seqpos+1),
			AtomID(next_rsd.atom_index("P"),seqpos+1),
			AtomID(rsd.atom_index("O3'"),seqpos),
			f_P_bond_angle_std_angle
			) ) )
		);

		func::FuncOP f_O3_bond_angle_std_angle( new func::HarmonicFunc( O3_bond_angle, std_angle ) );
		cst_set->add_constraint(
			ConstraintCOP( ConstraintOP( new AngleConstraint(
			AtomID(next_rsd.atom_index("P"),seqpos+1),
			AtomID(next_rsd.atom_index("O3'"),seqpos),
			AtomID(rsd.atom_index("C3'"),seqpos),
			f_O3_bond_angle_std_angle
			) ) )
		);

		// could also add angle constraints
		mini_pose.constraint_set( cst_set );

		ScoreFunction sf;
		sf.set_weight( atom_pair_constraint, 1.0 );
		sf.set_weight( angle_constraint    , 1.0 );

		// setup the options
		MinimizerOptions options( "lbfgs_armijo_nonmonotone" /* or "linmin" */, min_tolerance, true /*use_nblist*/, false /*deriv_check*/, false /*no verbose-deriv-check, is default*/ );

		kinematics::MoveMap mm;
		// allow minimizer to change these five dihedrals
		mm.set( TorsionID(seqpos-1, BB, 6), true );
		mm.set( TorsionID(seqpos  , BB, 1), true );
		mm.set( TorsionID(seqpos  , BB, 2), true );
		mm.set( TorsionID(seqpos  , BB, 3), true );
		mm.set( TorsionID(seqpos  , BB, 4), true );

		for ( utility::vector1< DihedralSet* >::const_iterator lib_iter = library.begin(), end = library.end(); lib_iter != end; ++lib_iter ) {

			DihedralSet const *dihedrals(*lib_iter);

			// Set mainchain dihedral angles (zeta, alpha - epsilon)
			Size const nbb( rsd.n_mainchain_atoms() );
			for ( int i=1; i<=6; ++i ) {
				TorsionID const id ( ( i==1 ) ? ( TorsionID( seqpos-1, BB, nbb ) ) : ( TorsionID( seqpos, BB, i-1 ) ) );
				mini_pose.set_torsion( id, (*dihedrals)[i] );
			}
			// Set chi
			TorsionID const id ( seqpos, CHI, 1 );
			mini_pose.set_torsion( id, (*dihedrals)[7] );

			Vector const rot_P_O3( ( next_rsd.xyz("P") - mini_pose.residue(seqpos).xyz("O3'") ) );

			// consider rotamers only if P-O3' distance is "similar" to that in the original structure.
			if ( fabs( target_distance - rot_P_O3.length() ) <  distance_cutoff ) {

				bool accept_rotamer = true;

				if ( minimize_rotamers ) {
					//      mini_pose.dump_pdb( "before_"+string_of(rotamers.size())+".pdb" );

					//      PROF_START( TEST1 );
					AtomTreeMinimizer().run( mini_pose, mm, sf, options );
					//      PROF_STOP ( TEST1 );

					//      tt << "LIBR.: ";
					//      for ( int i=1; i<=5; ++i ) {
					//       tt << F(10,3,(*dihedrals)[i]);
					//      }

					//      tt << std::endl << "AFTER: ";
					//      for (Size i=1; i <= 5; ++i) {
					//       TorsionID const id ( ( i==1 ) ? ( TorsionID( seqpos-1, BB, nbb ) ) : ( TorsionID( seqpos, BB, i-1 ) ) );
					//       tt << F(10,3, mini_pose.torsion( id ) );
					//      }

					//      tt << std::endl << "DIFF : ";
					for ( int i=1; i<=5; ++i ) {
						TorsionID const id ( ( i==1 ) ? ( TorsionID( seqpos-1, BB, nbb ) ) : ( TorsionID( seqpos, BB, i-1 ) ) );
						Real diff_i = fabs( mini_pose.torsion( id ) - (*dihedrals)[i] );

						if ( diff_i>180 ) {
							diff_i = 360-diff_i;
						}

						if ( diff_i > angle_tol ) {
							accept_rotamer = false;
						}

						//       tt << F(10,3,diff_i);
					}
					//      tt << std::endl << std::endl;
				}

				if ( accept_rotamer ) {
					Residue const &
						cur_rsd( mini_pose.residue( seqpos ) ),
						cur_next_rsd( mini_pose.residue(seqpos+1) );
					Vector const cur_P_O3( ( cur_next_rsd.xyz("P") - cur_rsd.xyz("O3'") ) );
					Real const cur_P_bond_angle ( angle_of( cur_next_rsd.xyz("O5'"), cur_next_rsd.xyz("P")  ,  cur_rsd.xyz("O3'")  ) );
					Real const cur_O3_bond_angle( angle_of( cur_next_rsd.xyz("P")  ,      cur_rsd.xyz("O3'"),  cur_rsd.xyz("C3'")  ) );

					Real const diff_length ( fabs( cur_P_O3.length()  - target_distance ) );

					Real diff_P_angle ( rad2deg * fabs( cur_P_bond_angle  -  P_bond_angle ) );
					Real diff_O3_angle( rad2deg * fabs( cur_O3_bond_angle - O3_bond_angle ) );
					if ( diff_P_angle > 180 ) {
						diff_P_angle = 360 - diff_P_angle;
					}
					if ( diff_O3_angle > 180 ) {
						diff_O3_angle = 360 - diff_O3_angle;
					}

					if ( ( diff_P_angle > angle_tol ) || ( diff_O3_angle > angle_tol ) || ( diff_length > dist_tol ) ) {
						accept_rotamer = false;
					}

					//      if ( accept_rotamer )
					//       tt << "ACCEPTED: Opt. diff: P-O3="  << diff_length << ", P_bond diff="  << diff_P_angle << ", O3_bond diff=" << diff_O3_angle << std::endl;
					//      else
					//       tt << "REJECTED: Opt. diff: P-O3="  << diff_length << ", P_bond diff="  << diff_P_angle << ", O3_bond diff=" << diff_O3_angle << std::endl;

					//      mini_pose.dump_pdb( "after_"+string_of(rotamers.size())+".pdb" );
				}

				if ( accept_rotamer ) {
					ResidueOP new_rot( mini_pose.residue(seqpos).clone() );

					// may have to do still more in current code
					new_rot->seqpos( existing_residue.seqpos() );
					new_rot->chain ( existing_residue.chain () );
					new_rot->copy_residue_connections( existing_residue ); // AARGHH!!
					rotamers.push_back( new_rot );
				}
			}
		}
	}

	tt << "Built " << rotamers.size() << " rotamers for " << resid << " " << concrete_residue->name() << std::endl;
}

// move this helper function
// better yet, stuff inside new rotamer-building class
void
build_random_dna_rotamers(
	Size const resid,
	pose::Pose const & pose,
	chemical::ResidueTypeCOP concrete_residue,
	pack::task::ExtraRotSample const & level,
	utility::vector1< conformation::ResidueOP > & rotamers
)
{
	using namespace pose;
	using namespace conformation;
	using namespace id;
	using namespace pack::task;


	Real const max_rotation( 5.0 ); // degrees
	bool const tether_bond_distance( true );
	bool const tether_bond_angle( false );

	int nrot_to_build(1);
	switch ( level ) {
	case NO_EXTRA_CHI_SAMPLES : // 0
		nrot_to_build = 1;
		break;
	case EX_ONE_STDDEV : // 1
		nrot_to_build = 10;
		break;
	case EX_ONE_HALF_STEP_STDDEV : // 2
		nrot_to_build = 25;
		break;
	case EX_TWO_FULL_STEP_STDDEVS : // 3
		nrot_to_build = 50;
		break;
	case EX_TWO_HALF_STEP_STDDEVS : // 4
		nrot_to_build = 100;
		break;
	case EX_FOUR_HALF_STEP_STDDEVS : // 5
		nrot_to_build = 250;
		break;
	case EX_THREE_THIRD_STEP_STDDEVS : // 6
		nrot_to_build = 500;
		break;
	case EX_SIX_QUARTER_STEP_STDDEVS : // 7
		nrot_to_build = 1000;
		break;
	case ExtraRotSampleCardinality :
	default :
		std::cerr << "Error in RotamerSet_::set_extrachi_samples, invalid ExtraChiSample type" << std::endl;
		utility_exit();
		break;
	}

	Residue const & existing_residue( pose.residue( resid ) );
	debug_assert( existing_residue.is_NA() && concrete_residue->is_NA() );

	if ( existing_residue.is_terminus() ) nrot_to_build = 1;

	//  basic::T( "core.pack.rotamer_set", basic::t_info ) << "Building " << nrot_to_build << " DNA rotamers for " << resid <<
	//  " " << concrete_residue->name() << '\n';

	// a residue of the correct type aligned to the existing backbone
	ResidueOP rot = ResidueFactory::create_residue( *concrete_residue, existing_residue, pose.conformation() );
	rot->set_chi( 1, existing_residue.chi(1) );

	// tt << "Pushing DNA rotamer " << std::endl;
	// for( Size ii = 1 ; ii <= concrete_residue->natoms() ; ++ii) {
	//  Vector coords( rot->xyz(ii) );
	//  tt << "Atom " << concrete_residue->atom_name( ii ) <<
	//   "  " << coords[0] <<
	//   "  " << coords[1] <<
	//   "  " << coords[2] << std::endl;
	// }

	rotamers.push_back( rot );

	if ( nrot_to_build > 1 ) {

		///////////////////////////////////////////////////////////////////
		// new fancy way
		// make a minipose
		pose::Pose mini_pose;
		mini_pose.append_residue_by_bond( pose.residue( resid-1 ) );
		mini_pose.append_residue_by_bond( *rot );
		mini_pose.append_residue_by_jump( pose.residue( resid+1 ), 1 );

		Size const seqpos(2);
		Residue const &
			prev_rsd( mini_pose.residue(seqpos-1) ),
			rsd( mini_pose.residue(seqpos  ) ),
			next_rsd( mini_pose.residue(seqpos+1) );

		debug_assert( rsd.is_DNA() && prev_rsd.is_DNA() );

		utility::vector1< Vector > bonds, atoms;

		atoms.push_back( prev_rsd.xyz("O3'") );
		atoms.push_back(      rsd.xyz(  "P") );
		atoms.push_back(      rsd.xyz("O5'") );
		atoms.push_back(      rsd.xyz("C5'") );
		atoms.push_back(      rsd.xyz("C4'") );

		Vector const r( rsd.xyz("O3'") );

		Vector const n1( ( next_rsd.xyz(  "P") -      rsd.xyz("O3'" ) ).normalized() );
		Vector const n2( ( next_rsd.xyz("O5'") - next_rsd.xyz(  "P" ) ).normalized() );

		Real l1( 0.0 );
		utility::vector1< Real > v1(4), v2(4);

		for ( int i=1; i<= 4; ++i ) {
			Vector const bi( ( atoms[i+1] - atoms[i] ).normalized() );
			Vector const ci( atoms[i+1] );
			v1[i] = n1.dot( bi.cross( r - ci ) );
			v2[i] = n2.dot( bi.cross( r - ci ) );
			l1 += v1[i] * v1[i];
		}

		// want to choose e[i] so that e dot v1 = 0 and e dot v2 = 0
		// eg choose two indices randomly
		// choose random values for those indices
		// now choose the unique values of the other two e[i]'s that will give e.v1= e.v2=0
		//
		// even simpler, orthonormalize v1 and v2, then we start with e' and subtract off dot(e',v1)*v1 and
		// dot(e',v2)*v2 and we are golden.


		// normalize v1, compute dot product with v2
		l1 = std::sqrt( l1 );
		Real v1_dot_v2(0.0);
		for ( int i=1; i<= 4; ++i ) {
			v1[i] /= l1;
			v1_dot_v2 += v1[i] * v2[i];
		}

		// now subtract
		Real l2(0.0);
		for ( int i=1; i<= 4; ++i ) {
			v2[i] -= v1_dot_v2 * v1[i];
			l2 += v2[i] * v2[i];
		}
		// normalize v2
		l2 = std::sqrt(l2);
		for ( int i=1; i<= 4; ++i ) v2[i] /= l2;

		// confirm orthogonality
		Real tmp1(0.0), tmp2(0.0), tmp12(0.0);
		for ( int i=1; i<= 4; ++i ) {
			tmp1  += v1[i] * v1[i];
			tmp2  += v2[i] * v2[i];
			tmp12 += v1[i] * v2[i];
		}
		debug_assert( std::abs(tmp1-1)<1e-3 && std::abs(tmp2-1)<1e-3 && std::abs(tmp12)<1e-3 );


		// now try a bunch of moves:
		for ( int r=1; r<= nrot_to_build-1; ++r ) {


			// choose random starting points
			utility::vector1< Real > e(4);
			Real dot1(0.0), dot2(0.0);
			for ( int i=1; i<= 4; ++i ) {
				e[i] = max_rotation - 2*numeric::random::rg().uniform()*max_rotation;
				dot1 += e[i] * v1[i];
				dot2 += e[i] * v2[i];
			}

			for ( int i=1; i<= 4; ++i ) {
				if ( tether_bond_distance ) e[i] -= dot1 * v1[i];
				if ( tether_bond_angle    ) e[i] -= dot2 * v2[i];
			}

			{ //debug
				dot1 = 0; dot2 = 0;
				for ( int i=1; i<= 4; ++i ) {
					dot1 += e[i] * v1[i];
					dot2 += e[i] * v2[i];
				}
				if ( tether_bond_distance ) debug_assert( std::abs( dot1 ) < 1e-3 );
				if ( tether_bond_angle    ) debug_assert( tether_bond_distance && std::abs( dot2 ) < 1e-3 );
			}


			Size const nbb( rsd.n_mainchain_atoms() );
			for ( int i=1; i<= 4; ++i ) {

				TorsionID const id ( ( i==1 ) ? ( TorsionID( seqpos-1, BB, nbb ) ) : ( TorsionID( seqpos, BB, i-1 ) ) );
				TorsionID const id0( ( i==1 ) ? ( TorsionID(  resid-1, BB, nbb ) ) : ( TorsionID(  resid, BB, i-1 ) ) );
				mini_pose.set_torsion( id, pose.torsion( id0 ) + e[i] );
			}

			if ( true ) { // chi
				TorsionID const id ( seqpos, CHI, 1 );
				TorsionID const id0(  resid, CHI, 1 );
				mini_pose.set_torsion( id, pose.torsion( id0 ) + max_rotation - 2*numeric::random::rg().uniform()*max_rotation );
			}

			rotamers.push_back( mini_pose.residue(seqpos).clone() );

			// this is a little tricky... maybe have to do more initialization here to match the "rotamer-constructor"
			// behavior??
			debug_assert( Size(rot->seqpos()) == resid );
			rotamers[ rotamers.size() ]->seqpos( resid );
			rotamers[ rotamers.size() ]->chain( rot->chain() );
			rotamers[ rotamers.size() ]->copy_residue_connections( existing_residue );

			// check dr:
			//Vector const r0( mini_pose.residue(seqpos).xyz("O4'") - pose.residue(resid).xyz("O4'") );
			//Vector const r ( mini_pose.residue(seqpos).xyz("O3'") - pose.residue(resid).xyz("O3'") );
			//tt << "r: " << tag << F(9,3,r.length()) << F(9,3,r.dot( n1 )) << F(9,3,r.dot( n2 ) ) << std::endl;

			// confirm that phosphate is in the same place
			debug_assert( mini_pose.residue(seqpos  ).xyz("P").distance( pose.residue(resid  ).xyz("P")) < 1e-2 &&
				mini_pose.residue(seqpos+1).xyz("P").distance( pose.residue(resid+1).xyz("P")) < 1e-2 );

		} // r =1,nrot_to_build
	}
}

void
build_dna_rotamers(
	Size const resid,
	pose::Pose const & pose,
	chemical::ResidueTypeCOP concrete_residue,
	pack::task::PackerTask const & task,
	utility::vector1< conformation::ResidueOP > & rotamers
)
{
	using conformation::ResidueOP;

	utility::vector1< ResidueOP > new_rotamers;

	pack::task::ExtraRotSample const level
		( task.residue_task( resid ).extrachi_sample_level( true /*buried*/, 1 /*chi*/, *concrete_residue ) );

	// need to improve this logic:
	if ( level == task::EX_SIX_QUARTER_STEP_STDDEVS ) {

		utility::io::izstream lib_stream;
		basic::database::open( lib_stream, "rotamer/dna/VQ-DNA-64.rotlib" );

		// Read rotamers from rotamer library
		utility::vector1< DihedralSet* > library;
		read_DNA_rotlib( lib_stream, library );

		build_lib_dna_rotamers( library, resid, pose, concrete_residue, new_rotamers );

		// Close stream and free library memory
		lib_stream.close();
		for ( utility::vector1< DihedralSet* >::iterator lib_iter = library.begin(), end = library.end(); lib_iter!=end; ++lib_iter ) {
			delete( *lib_iter );
		}

	} else {
		build_random_dna_rotamers( resid, pose, concrete_residue, level, new_rotamers );
	}

	// check for compatibility
	for ( Size ii=1; ii<= new_rotamers.size(); ++ii ) {
		ResidueOP rot( new_rotamers[ii] );

		//    for( Size i = 1 ; i <= rot->natoms() ; ++i ) {
		//      tt << "precompat " << rot->atom_name(i) <<
		//          " x " << rot->xyz(i)[0] <<
		//          " y " << rot->xyz(i)[1] <<
		//          " z " << rot->xyz(i)[2] << std::endl;
		//    }


		if ( task.rotamer_couplings_exist() ) {
			bool ok( true );
			RotamerCouplings const & couplings( * (task.rotamer_couplings() ) );
			int const other_resid( couplings[ resid ].first );
			if ( other_resid && !task.pack_residue( other_resid ) ) {
				conformation::ResidueMatcher const & matcher( *(couplings[resid].second ) );
				if ( !matcher( *rot, pose.residue( other_resid ) ) ) {
					ok = false;
				}
			}
			if ( ok ) {
				rotamers.push_back( rot );
			}
		} else {
			//tt << "added dna rotamer: " << existing_residue.seqpos() << ' ' << existing_residue.name() << ' '<<
			// concrete_residue->name() << std::endl;
			rotamers.push_back( rot );
		} // do rotamer couplings exist?
	}
}


/// @brief Internal, recursive function to implement the external interface of build_rotamers_from_rotamer_bins(), below
void
build_rotamers_from_rotamer_bins( conformation::Residue const & residue,
	utility::vector1< conformation::ResidueOP > & rotamers,
	uint current_chi_index,
	utility::vector1< uint > & current_bin_indices )
{
	using namespace std;

	Size const n_chis( residue.type().nchi() );

	if ( current_chi_index <= n_chis ) {
		RotamerBins const bins_for_current_chi( residue.chi_rotamers( current_chi_index ) );
		Size const n_bins_for_current_chi( bins_for_current_chi.size() );
		if ( n_bins_for_current_chi != 0 ) {
			for ( uint i( 1 ); i <= n_bins_for_current_chi; ++i ) {
				current_bin_indices.at( current_chi_index ) = i;
				build_rotamers_from_rotamer_bins( residue, rotamers, current_chi_index + 1, current_bin_indices );  // Recurse.
			}
		} else {  // Leave the bin index at 0 and drop to next level.
			build_rotamers_from_rotamer_bins( residue, rotamers, current_chi_index + 1, current_bin_indices );  // Recurse
		}
		// If the current chi index is greater than the number of chis, every index has been changed, and it is time to
		// cease recursing and make a rotamer from the current state of the bin indices.
	} else {
		// tt.Trace << "Current bin index for rotamer generation: " << current_bin_indices << endl;
		rotamers.push_back( residue.create_rotamer() );  // Clone current residue.
		// Set the torsions of the new rotamer to values indexed by current_bin_indices.
		tt.Debug << "Selecting ";
		for ( uint i = 1; i <= n_chis; ++i ) {
			RotamerBins bins( residue.chi_rotamers( i ) );  // Get bins for ith torsion angle.
			uint const bin_index( current_bin_indices.at( i ) );  // Get the appropriate bin for this rotamer by index.
			if ( bin_index != 0 ) {  // If index is 0, there are no bins: skip setting this torsion.
				RotamerBin bin = bins[ bin_index ];
				Angle const torsion( bin.first );  // first is the setting; second is the SD.
				tt.Debug << "bin #" << current_bin_indices.at(i) <<
					" from torsion angle chi" << i << " (" << torsion << ")  ";
				rotamers[ rotamers.size() ]->set_chi( i, torsion );  // Set the rotamer's torsion angle.
			}
		}
		tt.Debug << endl;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Make a rotamer (Residue) for every combination of torsion angle in the rotamer bins listed in the params file for the
// given Residue.
/// @param   <residue>: the Residue for which rotamers are to be made
/// @param   <rotamers>: a list of rotamers (Residues) to which this method will add new rotamers
/// @details Uses recursion to make a nested loop n levels deep, where n is the number of side-chain torsions ("chis").
/// @author  Labonte
/// @note    This currently assumes that all torsions are independent.  For now, this is the default method of handling
/// carbohydrate rotamers.  Ideally, we will find an alternative method for handling carbohydrate rotamer libraries....
void build_rotamers_from_rotamer_bins(conformation::Residue const & residue,
	utility::vector1<conformation::ResidueOP> & rotamers ) {

	Size const n_chis( residue.type().nchi() );
	tt.Debug << "Creating list of indices; ";
	utility::vector1< uint > current_bin_indices( n_chis, 0 );  // Initialize with zeroes.
	tt.Debug << "nesting " << n_chis << " levels deep." << std::endl;

	build_rotamers_from_rotamer_bins( residue, rotamers, 1, current_bin_indices );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
debug_dump_rotamers(
	utility::vector1< conformation::ResidueOP > & rotamers
)
{
	using conformation::ResidueOP;

	for ( Size ii=1; ii<= rotamers.size(); ++ii ) {
		ResidueOP rot( rotamers[ii] );

		tt << "Rotamer at position " << rot->seqpos() << std::endl;
		tt << "  of type " << rot->type().name() << std::endl;

		for ( Size jj=1; jj <= rot->natoms() ; ++jj ) {
			conformation::Atom const & atm = rot->atom( jj );
			tt << "Atom " << rot->atom_name( jj ) << " at " << atm.xyz()[0] << "  " << atm.xyz()[1] << "  " << atm.xyz()[2] << std::endl;
		}
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Water-specific methods
conformation::ResidueOP
create_oriented_water_rotamer(
	chemical::ResidueType const & h2o_type,
	Vector const & xyz_atom1,
	Vector const & xyz_O,
	Vector const & xyz_atom2,
	std::string const & name1,
	std::string const & name2,
	conformation::Residue const & tp5 // for the approx geometry of the lone pairs
)
{


	Vector const tp5_O( tp5.xyz("O") ), tp5_atom1( tp5.xyz( name1 ) ), tp5_atom2( tp5.xyz( name2 ) );

	Vector const xyz_midpoint( xyz_O + ( xyz_atom1 - xyz_O ).normalized() + ( xyz_atom2 - xyz_O ).normalized() );
	Vector const tp5_midpoint( tp5_O + ( tp5_atom1 - tp5_O ).normalized() + ( tp5_atom2 - tp5_O ).normalized() );

	kinematics::Stub const xyz_stub( xyz_O, xyz_midpoint, xyz_atom1 );
	kinematics::Stub const tp5_stub( tp5_O, tp5_midpoint, tp5_atom1 );

	conformation::ResidueOP rot( conformation::ResidueFactory::create_residue( h2o_type ) );
	debug_assert( rot->natoms() == 3 );
	rot->set_xyz(  "O", xyz_O );
	rot->set_xyz( "H1", xyz_stub.local2global( tp5_stub.global2local( tp5.xyz("H1") ) ) );
	rot->set_xyz( "H2", xyz_stub.local2global( tp5_stub.global2local( tp5.xyz("H2") ) ) );
	return rot;
}

void
build_fixed_O_water_rotamers_independent(
	Size const seqpos,
	chemical::ResidueType const & h2o_type,
	pose::Pose const & pose,
	utility::graph::GraphCOP packer_neighbor_graph,
	utility::vector1< conformation::ResidueOP > & new_rotamers
)
{
	using conformation::Residue;
	using conformation::ResidueOP;
	using numeric::conversions::degrees;

	static Real const heavy_atom_cutoff2( 3.3 *  3.3 );
	Real const min_theta(  90.0 );
	Real const max_theta( 180.0 );


	debug_assert( h2o_type.natoms() == 3 ); // only works for this guy now

	Residue const & existing_rsd( pose.residue(seqpos) );
	Vector const & xyz_O( existing_rsd.nbr_atom_xyz() );

	chemical::ResidueTypeSetCOP residue_set( h2o_type.residue_type_set() );

	// logic: look for nearby donors/acceptors
	// look for acceptors within 3.2
	// look for donors with:
	//  (a) donor-base distance < 3.2 (could be more subtle by requiring good hydrogen geometry if that H isnt moving)
	//  (b) ...

	utility::vector1< Vector > acceptors, donors;

	for ( utility::graph::Graph::EdgeListConstIter
			ir  = packer_neighbor_graph->get_node( seqpos )->const_edge_list_begin(),
			ire = packer_neighbor_graph->get_node( seqpos )->const_edge_list_end();
			ir != ire; ++ir ) {
		Size const i( (*ir)->get_other_ind( seqpos ) );
		debug_assert( i != seqpos );
		Residue const & rsd( pose.residue(i) );

		// donors
		for ( chemical::AtomIndices::const_iterator
				hnum  = rsd.Hpos_polar().begin(),
				hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
			Vector const & xyz( rsd.xyz( rsd.atom_base( *hnum ) ) );
			if ( xyz_O.distance_squared( xyz ) < heavy_atom_cutoff2 ) {
				donors.push_back( xyz );
			}
		}

		// acceptors
		for ( chemical::AtomIndices::const_iterator
				anum  = rsd.accpt_pos().begin(),
				anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {
			Vector const & xyz( rsd.xyz( *anum ) );
			if ( xyz_O.distance_squared( xyz ) < heavy_atom_cutoff2 ) {
				acceptors.push_back( xyz );
			}
		}
	} // i=1,nres


	ResidueOP tp5( conformation::ResidueFactory::create_residue( residue_set->name_map("TP5") ) );

	for ( Size i=1; i<= donors.size(); ++i ) {
		Vector const & xyz1( donors[i] );

		// donor-donor pairs ////////////////////////////////
		for ( Size j=i+1; j<= donors.size(); ++j ) {
			Vector const & xyz2( donors[j] );

			Real const theta( degrees( angle_of( xyz1, xyz_O, xyz2 ) ) );
			if ( theta > min_theta && theta < max_theta ) {

				new_rotamers.push_back( create_oriented_water_rotamer( h2o_type, xyz1, xyz_O, xyz2, "EP1", "EP2", *tp5 ) );

			}
		}

		// donor-acceptor pairs /////////////////////////////
		for ( Size j=1; j<= acceptors.size(); ++j ) {
			Vector const & xyz2( acceptors[j] );

			Real const theta( degrees( angle_of( xyz1, xyz_O, xyz2 ) ) );
			if ( theta > min_theta && theta < max_theta ) {

				new_rotamers.push_back( create_oriented_water_rotamer( h2o_type, xyz1, xyz_O, xyz2, "EP1", "H1", *tp5 ) );
				new_rotamers.push_back( create_oriented_water_rotamer( h2o_type, xyz1, xyz_O, xyz2, "EP1", "H2", *tp5 ) );

			}
		}
	}

	// acceptor-acceptor pairs ////////////////////////////
	for ( Size i=1; i<= acceptors.size(); ++i ) {
		Vector const & xyz1( acceptors[i] );
		for ( Size j=i+1; j<= acceptors.size(); ++j ) {
			Vector const & xyz2( acceptors[j] );
			Real const theta( degrees( angle_of( xyz1, xyz_O, xyz2 ) ) );

			if ( theta > min_theta && theta < max_theta ) {
				new_rotamers.push_back( create_oriented_water_rotamer( h2o_type, xyz1, xyz_O, xyz2, "H1", "H2", *tp5 ) );
			}

		}
	}


	// update sequence position and chain information
	for ( Size i=1; i<= new_rotamers.size(); ++i ) {
		new_rotamers[i]->seqpos( seqpos );
		new_rotamers[i]->chain( existing_rsd.chain() );
	}


	// now loop over pairs
	tt << "build_water_rotamers seqpos= " << seqpos << " found " << donors.size() << " nearby donors, " <<
		acceptors.size() << " nearby acceptors. Built " << new_rotamers.size() << " rotamers." << std::endl;


}

Vector
build_optimal_water_O_on_donor(
	Vector const & hxyz,
	Vector const & dxyz
)
{

	Real const distance( 2.75 );

	return ( dxyz + distance * ( hxyz - dxyz ).normalized() );
}

void
build_optimal_water_Os_on_acceptor(
	Vector const & a_xyz,
	Vector b1_xyz, // local copy
	Vector const & b2_xyz,
	chemical::Hybridization const & hybrid,
	utility::vector1< Vector > & O_list
)
{
	using numeric::conversions::radians;
	using namespace chemical;

	Real const distance( 2.75 ); // acceptor--O distance

	Real theta(0.0), step_size(0.0);
	utility::vector1< Real > phi_list, phi_steps;
	phi_steps.push_back( -4 );
	phi_steps.push_back( -2 );
	phi_steps.push_back( -1 );
	phi_steps.push_back(  0 );
	phi_steps.push_back(  1 );
	phi_steps.push_back(  2 );
	phi_steps.push_back(  4 );

	switch( hybrid ) {
	case SP2_HYBRID :
		theta = 180.0 - 120.0;
		step_size = 15.0;
		phi_list.push_back(   0.0 );
		phi_list.push_back( 180.0 );
		break;
	case SP3_HYBRID :
		theta = 180.0 - 109.0;
		step_size = 10.0;
		phi_list.push_back( 120.0 );
		phi_list.push_back( 240.0 );
		break;
	case RING_HYBRID :
		b1_xyz = 0.5 * ( b1_xyz + b2_xyz );
		theta = 0.0;
		phi_steps.clear();
		phi_steps.push_back( 0.0 );
		phi_list.push_back( 0.0 ); // doesnt matter
		step_size = 0.0; // doesnt matter
		break;
	default :
		tt << "Bad hybridization type for acceptor " << hybrid << '\n';
		utility_exit();
	}

	kinematics::Stub stub( a_xyz, b1_xyz, b2_xyz );
	for ( Size i=1; i<= phi_list.size(); ++i ) {
		for ( Size j=1; j<= phi_steps.size(); ++j ) {
			O_list.push_back( stub.spherical( radians( phi_list[i] + step_size * phi_steps[j]), radians( theta ), distance));
		}
	}
}

void
build_donor_donor_waters(
	conformation::Residue const & rsd1,
	Size const hatm1,
	conformation::Residue const & rsd2,
	Size const hatm2,
	chemical::ResidueType const & h2o_type,
	conformation::Residue const & tp5,
	Size const nstep,
	utility::vector1< conformation::ResidueOP > & new_waters
)
{
	using namespace scoring::hbonds;
	HBondDatabaseCOP hb_database( HBondDatabase::get_database() );
	HBondOptions hbondopts; // default ctor, reads from command line

	Real const o12_dis2_cutoff( 25.0 ); // wild guess
	Real const hbond_energy_threshold( -0.10 );

	Vector const & hatm1_xyz( rsd1.xyz( hatm1 ) );
	Vector const & datm1_xyz( rsd1.xyz( rsd1.atom_base( hatm1 ) ) );
	Vector const & hatm2_xyz( rsd2.xyz( hatm2 ) );
	Vector const & datm2_xyz( rsd2.xyz( rsd2.atom_base( hatm2 ) ) );

	Vector const o1( build_optimal_water_O_on_donor( hatm1_xyz, datm1_xyz ) );
	Vector const o2( build_optimal_water_O_on_donor( hatm2_xyz, datm2_xyz ) );

	Real const o12_dis2( o1.distance_squared( o2 ) );
	if ( o1.distance_squared( o2 ) > o12_dis2_cutoff ) return;

	debug_assert( nstep >= 2 && nstep%2 == 0 );
	for ( Size step=0; step<= nstep; ++step ) {

		// calculate an intermediate water position
		Vector const xyz_O( o1 + ( Real(step)/nstep) * ( o2 - o1 ) );

		conformation::ResidueOP
			rot( create_oriented_water_rotamer( h2o_type, hatm1_xyz, xyz_O, hatm2_xyz, "EP1", "EP2", tp5 ) );

		/// evaluate the energy of the two hbonds
		Real energy1(0.0), energy2(0.0);
		HBEvalTuple hbe1(
			get_hb_don_chem_type(rsd1.atom_base(hatm1),rsd1),
			hbacc_H2O,
			seq_sep_other);
		hb_energy_deriv( *hb_database, hbondopts, hbe1, datm1_xyz, hatm1_xyz,
			rot->xyz("O"), rot->xyz("H1"), rot->xyz("H2"), energy1 );

		//if ( energy1 > hbond_energy_threshold ) return;

		HBEvalTuple hbe2(
			get_hb_don_chem_type(rsd2.atom_base(hatm2),rsd2),
			hbacc_H2O,
			seq_sep_other);
		hb_energy_deriv( *hb_database, hbondopts, hbe2, datm2_xyz, hatm2_xyz,
			rot->xyz("O"), rot->xyz("H1"), rot->xyz("H2"), energy2 );

		if ( energy1 > hbond_energy_threshold ) return;
		if ( energy2 > hbond_energy_threshold ) return;

		tt << "hbenergies: dd step= " << step << " nstep= " << nstep <<
			" E1= " << ObjexxFCL::format::F(9,3,energy1) << " E2= " << ObjexxFCL::format::F(9,3,energy2) <<
			" o12_distance= " << ObjexxFCL::format::F(9,3,std::sqrt( o12_dis2 )) << '\n';


		// accept this rotamer
		new_waters.push_back( rot );
	}
}

void
build_donor_acceptor_waters(
	conformation::Residue const & rsd1,
	Size const hatm1,
	conformation::Residue const & rsd2,
	Size const aatm2,
	chemical::ResidueType const & h2o_type,
	conformation::Residue const & tp5,
	Size const nstep,
	utility::vector1< conformation::ResidueOP > & new_waters
)
{
	using namespace scoring::hbonds;
	HBondDatabaseCOP hb_database( HBondDatabase::get_database());
	HBondOptions hbondopts; // default ctor, reads from command line

	Real const o12_dis2_cutoff( 25.0 ); // wild guess
	Real const hbond_energy_threshold( -0.10 );

	Vector const & hatm1_xyz( rsd1.xyz( hatm1 ) );
	Vector const & datm1_xyz( rsd1.xyz( rsd1.atom_base( hatm1 ) ) );
	Vector const & aatm2_xyz( rsd2.xyz( aatm2 ) );
	Vector const & aatm2_base_xyz( rsd2.xyz( rsd2.atom_base( aatm2 ) ) );
	Vector const & aatm2_base2_xyz( rsd2.xyz( rsd2.abase2( aatm2 ) ) );

	Vector const o1( build_optimal_water_O_on_donor( hatm1_xyz, datm1_xyz ) );
	utility::vector1< Vector > o2_pos_list;
	chemical::Hybridization const & hybrid( rsd2.atom_type(aatm2).hybridization() );
	build_optimal_water_Os_on_acceptor( aatm2_xyz, aatm2_base_xyz, aatm2_base2_xyz, hybrid, o2_pos_list );

	for ( Size jj=1; jj<= o2_pos_list.size(); ++jj ) {
		Vector const & o2( o2_pos_list[jj] );

		Real const o12_dis2( o1.distance_squared( o2 ) );
		if ( o1.distance_squared( o2 ) > o12_dis2_cutoff ) continue;

		debug_assert( nstep>= 2 && nstep%2 == 0 );

		for ( Size step=0; step<= nstep; ++step ) {

			// calculate an intermediate water position
			Vector const xyz_O( o1 + ( Real(step)/nstep) * ( o2 - o1 ) );

			conformation::ResidueOP rot
				( create_oriented_water_rotamer( h2o_type, hatm1_xyz, xyz_O, aatm2_xyz, "EP1", "H1", tp5 ) );

			/// evaluate the energy of the two hbonds

			// donor to water
			Real energy1(0.0);
			HBEvalTuple const hbe_type( get_hb_don_chem_type( rsd1.atom_base( hatm1 ), rsd1 ), hbacc_H2O, seq_sep_other);
			hb_energy_deriv( *hb_database, hbondopts, hbe_type, datm1_xyz, hatm1_xyz,
				rot->xyz("O"), rot->xyz("H1"), rot->xyz("H2"), energy1 );

			//if ( energy1 > hbond_energy_threshold ) continue;

			Real energy2(0.0);
			{ // water to acceptor
				using namespace chemical;
				HBEvalTuple const hbe_type
					( hbdon_H2O, get_hb_acc_chem_type( aatm2, rsd2 ) , seq_sep_other);
				hb_energy_deriv( *hb_database, hbondopts, hbe_type, rot->xyz("O"), rot->xyz("H1"),
					aatm2_xyz, aatm2_base_xyz, aatm2_base2_xyz, energy2 );
			}

			if ( energy1 > hbond_energy_threshold ) continue;
			if ( energy2 > hbond_energy_threshold ) continue;

			tt << "hbenergies: da step= " << step << " nstep= " << nstep <<
				" E1= " << ObjexxFCL::format::F(9,3,energy1) << " E2= " << ObjexxFCL::format::F(9,3,energy2) <<
				" o12_distance= " << ObjexxFCL::format::F(9,3,std::sqrt( o12_dis2 )) << '\n';


			// accept this rotamer
			new_waters.push_back( rot );
			new_waters.push_back( create_oriented_water_rotamer( h2o_type, hatm1_xyz, xyz_O, aatm2_xyz,
				"EP1", "H2", tp5 ) );
			//    {
			//     using namespace ObjexxFCL;
			//     using namespace ObjexxFCL::format;

			//     std::string const filename( "donor_acceptor_test_"+lead_zero_string_of(counter,4)+"_"+
			//                   lead_zero_string_of(step,4)+".pdb" );

			//     tt << "hbenergy_da: step " << I(4,step) << F(9,3,energy1) << F(9,3,energy2) <<
			//      F(9,3,std::sqrt(o12_dis2)) << ' '<< filename << '\n';

			//     pose::Pose pose;
			//     pose.append_residue_by_jump( rsd1, 1 );
			//     pose.append_residue_by_jump( rsd2, 1 );
			//     pose.append_residue_by_jump( *rot, 1 );
			//     pose.dump_pdb( filename );
			//    }

		}
	}
}

void
build_acceptor_acceptor_waters(
	conformation::Residue const & rsd1,
	Size const aatm1,
	conformation::Residue const & rsd2,
	Size const aatm2,
	chemical::ResidueType const & h2o_type,
	conformation::Residue const & tp5,
	Size const nstep,
	utility::vector1< conformation::ResidueOP > & new_waters
)
{
	using namespace scoring::hbonds;
	HBondDatabaseCOP hb_database( HBondDatabase::get_database());
	HBondOptions hbondopts; // default ctor, reads from command line

	Real const o12_dis2_cutoff( 25.0 ); // wild guess
	Real const hbond_energy_threshold( -0.10 );

	Vector const & aatm1_xyz( rsd1.xyz( aatm1 ) );
	Vector const & aatm1_base_xyz( rsd1.xyz( rsd1.atom_base( aatm1 ) ) );
	Vector const & aatm1_base2_xyz( rsd1.xyz( rsd1.abase2( aatm1 ) ) );

	Vector const & aatm2_xyz( rsd2.xyz( aatm2 ) );
	Vector const & aatm2_base_xyz( rsd2.xyz( rsd2.atom_base( aatm2 ) ) );
	Vector const & aatm2_base2_xyz( rsd2.xyz( rsd2.abase2( aatm2 ) ) );

	utility::vector1< Vector > o1_pos_list;
	chemical::Hybridization const & hybrid1( rsd1.atom_type(aatm1).hybridization() );
	build_optimal_water_Os_on_acceptor( aatm1_xyz, aatm1_base_xyz, aatm1_base2_xyz, hybrid1, o1_pos_list );

	utility::vector1< Vector > o2_pos_list;
	chemical::Hybridization const & hybrid2( rsd2.atom_type(aatm2).hybridization() );
	build_optimal_water_Os_on_acceptor( aatm2_xyz, aatm2_base_xyz, aatm2_base2_xyz, hybrid2, o2_pos_list );

	for ( Size ii=1; ii<= o1_pos_list.size(); ++ii ) {
		Vector const & o1( o1_pos_list[ii] );

		for ( Size jj=1; jj<= o2_pos_list.size(); ++jj ) {
			Vector const & o2( o2_pos_list[jj] );

			Real const o12_dis2( o1.distance_squared( o2 ) );
			if ( o1.distance_squared( o2 ) > o12_dis2_cutoff ) continue;

			debug_assert( nstep>= 2 && nstep%2 == 0 );

			for ( Size step=0; step<= nstep; ++step ) {

				// calculate an intermediate water position
				Vector const xyz_O( o1 + ( Real(step)/nstep) * ( o2 - o1 ) );


				conformation::ResidueOP rot
					( create_oriented_water_rotamer( h2o_type, aatm1_xyz, xyz_O, aatm2_xyz, "H1", "H2", tp5 ) );

				/// evaluate the energy of the two hbonds
				Real energy1(0.0), energy2( 0.0 );

				{ // water to acceptor1
					using namespace chemical;
					HBEvalTuple const hbe_type( hbdon_H2O, get_hb_acc_chem_type( aatm1, rsd1 ), seq_sep_other);
					hb_energy_deriv( *hb_database, hbondopts, hbe_type, rot->xyz("O"), rot->xyz("H1"),
						aatm1_xyz, aatm1_base_xyz, aatm1_base2_xyz, energy1 );
				}

				{ // water to acceptor
					using namespace chemical;
					HBEvalTuple const hbe_type
						( hbdon_H2O, get_hb_acc_chem_type( aatm2, rsd2 ), seq_sep_other);
					hb_energy_deriv( *hb_database, hbondopts, hbe_type, rot->xyz("O"), rot->xyz("H2"),
						aatm2_xyz, aatm2_base_xyz, aatm2_base2_xyz, energy2 );
				}

				if ( energy1 > hbond_energy_threshold ) continue;
				if ( energy2 > hbond_energy_threshold ) continue;

				tt << "hbenergies: aa step= " << step << " nstep= " << nstep <<
					" E1= " << ObjexxFCL::format::F(9,3,energy1) << " E2= " << ObjexxFCL::format::F(9,3,energy2) <<
					" o12_distance= " << ObjexxFCL::format::F(9,3,std::sqrt( o12_dis2 )) << '\n';


				// accept this rotamer
				new_waters.push_back( rot );
			} // step=0,nstep
		}
	}
}

void
build_moving_O_bridge_waters(
	conformation::Residue const & rsd1,
	Size const anchor_atom,
	conformation::Residue const & rsd2,
	chemical::ResidueType const & h2o_type,
	conformation::Residue const & tp5,
	Size const nstep,
	//  utility::vector1< id::AtomID > & atom_ids,
	utility::vector1< conformation::ResidueOP > & new_rotamers
)
{
	Real const dis2_cutoff( 36.0 );

	utility::vector1< id::AtomID > atom_ids;

	tt << "build_moving_O_bridge_waters " << rsd1.seqpos() << ' ' << rsd1.name() << ' ' << rsd2.seqpos() <<
		' ' << rsd2.name() << '\n';

	// get info about the anchor atom
	bool const anchor_is_donor( rsd1.atom_type( anchor_atom ).is_polar_hydrogen() );
	Vector const & anchor_xyz( rsd1.xyz( anchor_atom ) );

	// look for donors or acceptors that could be making an interaction with the water
	//
	// donors
	for ( chemical::AtomIndices::const_iterator
			hnum  = rsd2.Hpos_polar().begin(),
			hnume = rsd2.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const hatm( *hnum );

		if ( rsd2.xyz( hatm ).distance_squared( anchor_xyz ) > dis2_cutoff ) continue;

		Size const old_size( new_rotamers.size() );

		if ( anchor_is_donor ) {
			build_donor_donor_waters( rsd1, anchor_atom, rsd2, hatm, h2o_type, tp5, nstep, new_rotamers );
		} else {
			build_donor_acceptor_waters( rsd2, hatm, rsd1, anchor_atom, h2o_type, tp5, nstep, new_rotamers );
		}

		if ( new_rotamers.size() > old_size ) {
			tt << "built " << new_rotamers.size() - old_size << " rotamers between " <<
				rsd1.seqpos() << ' ' << rsd1.name() << ' ' << rsd1.atom_name( anchor_atom ) << " and " <<
				rsd2.seqpos() << ' ' << rsd2.name() << ' ' << rsd2.atom_name( hatm ) << '\n';
			atom_ids.push_back( id::AtomID( hatm, rsd2.seqpos() ) );
		}
	}

	// acceptors
	for ( chemical::AtomIndices::const_iterator
			anum  = rsd2.accpt_pos().begin(),
			anume = rsd2.accpt_pos().end(); anum != anume; ++anum ) {
		Size const aatm( *anum );

		if ( rsd2.xyz( aatm ).distance_squared( anchor_xyz ) > dis2_cutoff ) continue;

		Size const old_size( new_rotamers.size() );

		if ( anchor_is_donor ) {
			build_donor_acceptor_waters( rsd1, anchor_atom, rsd2, aatm, h2o_type, tp5, nstep, new_rotamers );
		} else {
			build_acceptor_acceptor_waters( rsd1, anchor_atom, rsd2, aatm, h2o_type, tp5, nstep, new_rotamers );
		}

		if ( new_rotamers.size() > old_size ) {
			tt << "built " << new_rotamers.size() - old_size << " rotamers between " <<
				rsd1.seqpos() << ' ' << rsd1.name() << ' ' << rsd1.atom_name( anchor_atom ) << " and " <<
				rsd2.seqpos() << ' ' << rsd2.name() << ' ' << rsd2.atom_name( aatm ) << '\n';
			atom_ids.push_back( id::AtomID( aatm, rsd2.seqpos() ) );
		}
	}
}

void
build_moving_O_water_rotamers_dependent(
	RotamerSets const & rotsets,
	//Size const seqpos_water,
	WaterAnchorInfo const & water_info,
	chemical::ResidueType const & h2o_type,
	pack::task::PackerTask const & task,
	pose::Pose const & pose,
	utility::graph::GraphCOP packer_neighbor_graph,
	utility::vector1< conformation::ResidueOP > & new_rotamers
)
{
	using conformation::Residue;
	using conformation::ResidueOP;
	using numeric::conversions::degrees;

	debug_assert( h2o_type.natoms() == 3 ); // only works for this guy now

	Size nstep( water_info.nstep() );
	if ( nstep < 2 ) nstep = 2;
	if ( nstep%2 == 1 ) ++nstep;

	// who are we "attached to" ?
	//

	Size const i( water_info.anchor_residue() );

	tt << "build_moving_O_water_rotamers_dependent: anchor_pos=  " << i << '\n';

	// build a tp5 water for geometry calculations below
	ResidueOP tp5( conformation::ResidueFactory::create_residue( h2o_type.residue_type_set()->name_map("TP5") ) );

	// build a list of residue/rotamers at this dna position to loop over
	utility::vector1< conformation::ResidueCOP > rsd1_list;
	if ( task.pack_residue(i) ) {
		// get the rotamerset
		RotamerSetCOP rotset( rotsets.rotamer_set_for_residue( i ) );
		for ( Size ii=1; ii<= rotset->num_rotamers(); ++ii ) {
			rsd1_list.push_back( rotset->rotamer(ii) );
		}
	} else {
		rsd1_list.push_back( pose.residue(i).get_self_ptr() );
	}

	// here we assume that the two residues bridged by water must be neighbors...
	// this may not be completely true: D-O + A-O could be bigger than 5.5, but prob not much
	//
	// may want to refine this

	for ( utility::graph::Graph::EdgeListConstIter
			jr  = packer_neighbor_graph->get_node( i )->const_edge_list_begin(),
			jre = packer_neighbor_graph->get_node( i )->const_edge_list_end();
			jr != jre; ++jr ) {
		Size const j( (*jr)->get_other_ind( i ) );
		debug_assert( i != j );

		//if ( !pose.residue(j).is_protein() ) continue; // first just protein-DNA bridges
		if ( pose.residue(j).name() == "TP3" ) continue;

		// fixed-fixed interactions were already hit inside the first pass of building
		if ( !task.pack_residue( i ) && !task.pack_residue( j ) ) continue;

		utility::vector1< conformation::ResidueCOP > rsd2_list;
		if ( task.pack_residue(j) ) {
			// get the rotamerset
			RotamerSetCOP rotset( rotsets.rotamer_set_for_residue( j ) );
			for ( Size jj=1; jj<= rotset->num_rotamers(); ++jj ) {
				rsd2_list.push_back( rotset->rotamer(jj) );
			}
		} else {
			rsd2_list.push_back( pose.residue(j).get_self_ptr() );
		}


		/// now loop over residue pairs

		/// now build the waters
		for ( Size ii=1; ii<= rsd1_list.size(); ++ii ) {
			Residue const & rsd1( *rsd1_list[ii] );
			if ( ! water_info.attaches_to_residue_type( rsd1.type() ) ) continue;
			Size const anchor_atom( water_info.anchor_atom( rsd1.type() ) );

			for ( Size jj=1; jj<= rsd2_list.size(); ++jj ) {
				build_moving_O_bridge_waters( rsd1, anchor_atom, *( rsd2_list[jj] ), h2o_type, *tp5, nstep, new_rotamers );
			}
		}
	} // nbrs of anchor_pos
}

void
build_moving_O_water_rotamers_independent(
	WaterAnchorInfo const & water_info,
	chemical::ResidueType const & h2o_type,
	pack::task::PackerTask const & task,
	pose::Pose const & pose,
	utility::graph::GraphCOP packer_neighbor_graph,
	utility::vector1< conformation::ResidueOP > & new_rotamers
)
{
	using conformation::Residue;
	using conformation::ResidueOP;
	using numeric::conversions::degrees;

	debug_assert( h2o_type.natoms() == 3 ); // only works for this guy now

	Size nstep( water_info.nstep() );
	if ( nstep < 2 ) nstep = 2;
	if ( nstep%2 == 1 ) ++nstep;

	Size const i( water_info.anchor_residue() );

	tt << "build_moving_O_water_rotamers_independent: anchor_pos=  " << i << '\n';

	if ( task.pack_residue(i) ) return; // have to wait

	Residue const & rsd1( pose.residue(i) );
	if ( ! water_info.attaches_to_residue_type( rsd1.type() ) ) return;
	Size const anchor_atom( water_info.anchor_atom( rsd1.type() ) );

	// build a tp5 water for geometry calculations below
	ResidueOP tp5( conformation::ResidueFactory::create_residue( h2o_type.residue_type_set()->name_map("TP5") ) );

	// here we assume that the two residues bridged by water must be neighbors...
	// this may not be completely true: D-O + A-O could be bigger than 5.5, but prob not much
	//
	// may want to refine this

	for ( utility::graph::Graph::EdgeListConstIter
			jr  = packer_neighbor_graph->get_node( i )->const_edge_list_begin(),
			jre = packer_neighbor_graph->get_node( i )->const_edge_list_end();
			jr != jre; ++jr ) {
		Size const j( (*jr)->get_other_ind( i ) );
		debug_assert( i != j );

		//if ( !pose.residue(j).is_protein() ) continue; // first just protein-DNA bridges
		if ( pose.residue(j).name() == "TP3" ) continue;

		if ( task.pack_residue( j ) ) continue;

		build_moving_O_bridge_waters( rsd1, anchor_atom, pose.residue(j), h2o_type, *tp5, nstep, new_rotamers );

	} // nbrs of anchor_pos
}

void
build_independent_water_rotamers(
	Size const seqpos_water,
	chemical::ResidueType const & h2o_type,
	pack::task::PackerTask const & task,
	pose::Pose const & pose,
	utility::graph::GraphCOP packer_neighbor_graph,
	utility::vector1< conformation::ResidueOP > & new_rotamers
)
{
	//using core::pose::datacache::CacheableDataType::WATER_PACKING_INFO;

	{
		tt << "water " << seqpos_water << " nbrs:";
		for ( utility::graph::Graph::EdgeListConstIter
				jr  = packer_neighbor_graph->get_node( seqpos_water )->const_edge_list_begin(),
				jre = packer_neighbor_graph->get_node( seqpos_water )->const_edge_list_end();
				jr != jre; ++jr ) {
			tt << ' ' << (*jr)->get_other_ind( seqpos_water );
		}
		tt << '\n';
	}


	if ( pose.data().has( core::pose::datacache::CacheableDataType::WATER_PACKING_INFO ) ) {
		WaterPackingInfo const & water_info
			( static_cast< WaterPackingInfo const & >( pose.data().get( core::pose::datacache::CacheableDataType::WATER_PACKING_INFO ) ) );

		build_moving_O_water_rotamers_independent( water_info[ seqpos_water], h2o_type, task, pose,
			packer_neighbor_graph, new_rotamers );
	} else {

		build_fixed_O_water_rotamers_independent( seqpos_water, h2o_type, pose, packer_neighbor_graph, new_rotamers );

	}
}

void
build_dependent_water_rotamers(
	RotamerSets const & rotsets,
	Size const seqpos_water,
	chemical::ResidueType const & h2o_type,
	pack::task::PackerTask const & task,
	pose::Pose const & pose,
	utility::graph::GraphCOP packer_neighbor_graph,
	utility::vector1< conformation::ResidueOP > & new_rotamers
)
{
	//using core::pose::datacache::CacheableDataType::WATER_PACKING_INFO;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::WATER_PACKING_INFO ) ) {
		WaterPackingInfo const & water_info
			( static_cast< WaterPackingInfo const & >( pose.data().get( core::pose::datacache::CacheableDataType::WATER_PACKING_INFO ) ) );

		build_moving_O_water_rotamers_dependent( rotsets, /*seqpos_water,*/ water_info[ seqpos_water], h2o_type, task, pose,
			packer_neighbor_graph, new_rotamers );
	} else {

		// doesnt exist yet
		//build_fixed_O_water_rotamers_dependent( rotsets, seqpos_water, h2o_type, pose, packer_neighbor_graph,
		// new_rotamers );

	}


}


} // rotamer_set
} // pack
} // core
