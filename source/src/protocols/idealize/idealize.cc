// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/idealize/idealize.cc
/// @brief
/// @author

// Unit Headers
#include <protocols/idealize/idealize.hh>


// // Rosetta Headers
#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>

#include <core/pose/PDBInfo.hh>


#include <basic/basic.hh>
#include <basic/Tracer.hh> // tracer output

#include <core/pose/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <ObjexxFCL/string.functions.hh>

// Numeric headers
#include <numeric/random/random.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

#include <core/chemical/VariantType.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>
namespace protocols {
namespace idealize {

using namespace core;
using namespace ObjexxFCL;

static THREAD_LOCAL basic::Tracer TR( "protocols.idealize" );

/// helper
void
dihedral_distance(
	pose::Pose const & pose1,
	pose::Pose const & pose2,
	utility::vector1< bool > const & use_pos,
	Real & avg_bb_angle_dev,
	Real & max_bb_angle_dev,
	Real & avg_chi_angle_dev,
	Real & max_chi_angle_dev
) {
	using basic::subtract_degree_angles;

	avg_bb_angle_dev = 0.0;
	max_bb_angle_dev = 0.0;

	avg_chi_angle_dev = 0.0;
	max_chi_angle_dev = 0.0;

	Size nchi_dihedrals(0), nbb_dihedrals(0);

	for ( Size pos = 1; pos <= pose1.total_residue(); ++pos ) {
		if ( ! use_pos[ pos ] ) continue;

		conformation::Residue const & rsd1( pose1.residue( pos ) );
		conformation::Residue const & rsd2( pose2.residue( pos ) );

		if ( !rsd1.is_polymer() ) continue;

		Size const nbb ( rsd1.n_mainchain_atoms() );
		Size const nchi( rsd1.nchi() );
		runtime_assert( rsd2.is_polymer() && rsd2.nchi() ==  nchi && rsd2.n_mainchain_atoms() == nbb );

		// first the bb dev's
		for ( Size i=1; i<= nbb; ++i ) {
			if ( ( i ==     1 && rsd1.is_lower_terminus() ) ||
					( i >= nbb-1 && rsd1.is_upper_terminus() ) ) {
				continue;
			}
			Real const dev( std::abs( subtract_degree_angles( rsd1.mainchain_torsion( i ), rsd2.mainchain_torsion(i) ) ) );
			//if ( dev > 0.01 ) std::cout << "bbdev: " << pos << ' ' << pose1.residue(pos).name1() << ' ' << i << ' ' <<
			//          dev << ' ' << rsd1.mainchain_torsion( i ) << ' ' <<  rsd2.mainchain_torsion(i) << std::endl;
			avg_bb_angle_dev += dev;
			++nbb_dihedrals;
			max_bb_angle_dev = std::max( max_bb_angle_dev, dev );
		}

		for ( Size i=1; i<= nchi; ++i ) {
			Real const dev( std::abs( subtract_degree_angles( rsd1.chi( i ), rsd2.chi( i ) ) ) );
			avg_chi_angle_dev += dev;
			++nchi_dihedrals;
			max_chi_angle_dev = std::max( max_chi_angle_dev, dev );
		}
	}

	avg_bb_angle_dev  /=  nbb_dihedrals;
	avg_chi_angle_dev /= nchi_dihedrals;
	TR.flush();
}

// positions within window of the idealized positions will move during minimization
void
basic_idealize(
	pose::Pose & pose,
	utility::vector1< Size > pos_list, // local copy
	scoring::ScoreFunction const & scorefxn,
	bool const fast,
	bool const chainbreaks,
	bool const cis_omega
) {
	using namespace optimization;
	using namespace optimization::symmetry;
	using namespace id;
	using namespace ObjexxFCL::format;
	using scoring::all_atom_rmsd;

	Size const window_width( 3 ); // window:  from pos-window_width to pos+window_width


	pose::Pose const start_pose( pose );
	Size const nres ( pose.total_residue() );

	// keep chainbreaks if they exist
	if ( chainbreaks ) {

		// squared distance at which bond is considered discontinuous
		Real const chain_break_cutoff = { 4.0 };

		// find chain breaks to add cut points
		bool new_cutpoint = false;
		pose::PDBInfoCOP pdbinfo = pose.pdb_info();
		kinematics::FoldTree f( pose.fold_tree() );
		for ( Size i = 1; i < nres; ++i ) {
			if ( f.is_cutpoint(i) ) continue;
			bool chain_break = false;
			Size j = i+1;
			if ( pdbinfo->number(i)+1 != pdbinfo->number(j) ) {
				TR.Info << "non-sequential at res nums " << i << '-' << j << std::endl;
				TR.Info << "non-sequential pdb res nums " << pdbinfo->number(i) << pdbinfo->chain(i) <<
					'-' << pdbinfo->number(j) << pdbinfo->chain(j) << std::endl;
				chain_break = true;
			} else {
				conformation::Residue const & rsd = pose.residue(i);
				conformation::Residue const & next_rsd = pose.residue(j);
				if ( rsd.is_polymer() && next_rsd.is_polymer() ) {
					Real dist_squared = rsd.atom( rsd.upper_connect_atom() ).xyz().distance_squared(next_rsd.atom( next_rsd.lower_connect_atom() ).xyz());
					if ( dist_squared > chain_break_cutoff ) {
						TR.Info << "chain break at res nums: " << i << '-' << j << ' ' << std::sqrt(dist_squared) << std::endl;
						TR.Info << "chain break pdb res nums: " << pdbinfo->number(i) << pdbinfo->chain(i) <<
							'-' << pdbinfo->number(j) << pdbinfo->chain(j) << std::endl;
						chain_break = true;
					} else if ( dist_squared < 0.1 ) {
						TR.Info << "zero length bond at res nums: " << i << '-' << j << std::endl;
						TR.Info << "zero length bond pdb res nums: " << pdbinfo->number(i) << pdbinfo->chain(i) <<
							'-' << pdbinfo->number(j) << pdbinfo->chain(j) << std::endl;
						chain_break = true;
					}
				}
			}
			if ( chain_break ) {
				// add cutpoint
				f.new_jump( i, j, i );
				pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, i );
				pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, i+1 );
				new_cutpoint = true;
				TR.Info << "added cutpoint at: " << i << std::endl;
			}
		}
		if ( new_cutpoint ) pose.fold_tree( f );
	}

	Size const njump( pose.num_jump() );

	// setup the minimizer options
	MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.001, true /*use_nblist*/, false /*deriv_check*/ );
	//MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.001, true /*use_nblist*/, true /*deriv_check*/, true );
	kinematics::MoveMap final_mm;

	bool const lastjumpmin (
		pose.residue( nres ).aa() == core::chemical::aa_vrt &&
		pose.fold_tree().upstream_jump_residue( njump ) == int(nres)
	);

	TR.Info << "lastjumpmin: " << lastjumpmin << std::endl;
	utility::vector1< bool > idealized( nres, false );

	// get symmetry info
	core::conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		TR.Info << "setting up symmetric idealize " << std::endl;
		core::conformation::symmetry::SymmetricConformation & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info();
	}

	while ( !pos_list.empty() ) {
		Size const seqpos( pos_list[ static_cast< Size >( numeric::random::rg().uniform() * pos_list.size() + 1 ) ] );
		pos_list.erase( std::find( pos_list.begin(), pos_list.end(), seqpos ) );

		// idealize the mainchain + sidechain
		//pose.dump_pdb( "pre_idl_"+right_string_of(seqpos,4,'0')+".pdb" );
		if ( seqpos > (Size)pose.conformation().size() ) { continue; }
		conformation::idealize_position( seqpos, pose.conformation() );
		//pose.dump_pdb( "post_idl_"+right_string_of(seqpos,4,'0')+".pdb" );
		idealized[ seqpos ] = true;

		if ( cis_omega && pose.residue(seqpos).is_protein() && pose.residue(seqpos).aa() != chemical::aa_pro ) {
			core::Real omega_i = fmod( pose.omega(seqpos), 360.0); if ( omega_i <0 ) omega_i+=360;
			if ( omega_i<90 || omega_i>270 )  {
				TR << "Fix cis omega at position " << seqpos << std::endl;
				pose.set_omega(seqpos, 180.0);
			}
		}

		// setup the window of positions to minimize, also records flexible positions for the final minimize
		utility::vector1< Size > window;
		for ( Size i=seqpos; i >= seqpos-window_width; --i ) {
			window.push_back( i );
			if ( i == 1 || pose.residue(i).is_lower_terminus() ) break;
		}
		for ( Size i=seqpos; i <= seqpos+window_width; ++i ) {
			window.push_back( i );
			if ( i == nres || pose.residue(i).is_upper_terminus() ) break;
		}


		kinematics::MoveMap local_mm;
		local_mm.set_chi( seqpos, true );
		final_mm.set_chi( seqpos, true );
		for ( Size ii=1; ii<= window.size(); ++ii ) {
			Size const i( window[ ii ] );
			local_mm.set_bb( i, true );
			final_mm.set_bb( i, true );
			// disallow proline PHI
			if ( pose.residue(i).aa() == chemical::aa_pro ) local_mm.set( TorsionID( phi_torsion, BB, i ), false );
			if ( pose.residue(i).aa() == chemical::aa_pro ) final_mm.set( TorsionID( phi_torsion, BB, i ), false );
		}

		// if jumpmin
		if ( lastjumpmin ) local_mm.set_jump( pose.num_jump(), true );

		// special case for symmetry
		//    - make mm symmetric
		//    - allow symmjumps to minimize
		if ( symm_info ) {
			local_mm.set_jump( true );
			core::pose::symmetry::make_symmetric_movemap( pose, local_mm );
		}

		// dont minimize or calculate stats after each idealization in fast mode
		if ( fast ) {
			TR.Info << "forced ideal geometry on seqpos " << seqpos << std::endl;
			continue;
		}

		Real max_bb_angle_dev, avg_bb_angle_dev, max_chi_angle_dev, avg_chi_angle_dev;

		dihedral_distance( pose, start_pose, idealized, avg_bb_angle_dev, max_bb_angle_dev,
			avg_chi_angle_dev, max_chi_angle_dev );

		TR.Info << "premin:  (pos,rmsd,avg-bb,max-bb,avg-chi,max-chi,score) " <<
			I( 4, seqpos ) << ' ' << pose.residue(seqpos).name1() << F( 9, 3, all_atom_rmsd( pose, start_pose ) ) <<
			F(9,3,avg_bb_angle_dev) << F(9,3,max_bb_angle_dev) <<
			F(9,3,avg_chi_angle_dev) << F(9,3,max_chi_angle_dev) <<
			F(12,3,scorefxn( pose ) ) << std::endl;

		if ( symm_info ) {
			SymAtomTreeMinimizer().run( pose, local_mm, scorefxn, options );
		} else {
			AtomTreeMinimizer().run( pose, local_mm, scorefxn, options );
		}
		//pose.dump_pdb( "post_min_"+right_string_of(seqpos,4,'0')+".pdb" );

		dihedral_distance( pose, start_pose, idealized, avg_bb_angle_dev, max_bb_angle_dev,
			avg_chi_angle_dev, max_chi_angle_dev );

		TR.Info << "postmin: (pos,rmsd,avg-bb,max-bb,avg-chi,max-chi,score) " <<
			I(4,seqpos) << ' ' << pose.residue(seqpos).name1() << F(9,3,all_atom_rmsd(pose,start_pose)) <<
			F(9,3,avg_bb_angle_dev) << F(9,3,max_bb_angle_dev) <<
			F(9,3,avg_chi_angle_dev) << F(9,3,max_chi_angle_dev) <<
			F(12,3,scorefxn( pose ) ) << std::endl;
		scorefxn.show( TR, pose );
	}

	if ( lastjumpmin ) final_mm.set_jump( pose.num_jump(), true );

	// special case for symmetry
	//    - make mm symmetric
	//    - allow symmjumps to minimize
	if ( symm_info ) {
		final_mm.set_jump( true );
		core::pose::symmetry::make_symmetric_movemap( pose, final_mm );
	}

	// final minimize
	Real max_bb_angle_dev, avg_bb_angle_dev, max_chi_angle_dev, avg_chi_angle_dev;
	dihedral_distance( pose, start_pose, idealized, avg_bb_angle_dev, max_bb_angle_dev,
		avg_chi_angle_dev, max_chi_angle_dev );

	TR.Info << "pre-finalmin: (pos,rmsd,avg-bb,max-bb,avg-chi,max-chi,score) " <<
		F(9,3,all_atom_rmsd(pose,start_pose)) <<
		F(9,3,avg_bb_angle_dev) << F(9,3,max_bb_angle_dev) <<
		F(9,3,avg_chi_angle_dev) << F(9,3,max_chi_angle_dev) <<
		F(12,3,scorefxn( pose ) ) << std::endl;

	if ( symm_info ) {
		SymAtomTreeMinimizer().run( pose, final_mm, scorefxn, options );
	} else {
		AtomTreeMinimizer().run( pose, final_mm, scorefxn, options );
	}

	dihedral_distance( pose, start_pose, idealized, avg_bb_angle_dev, max_bb_angle_dev,
		avg_chi_angle_dev, max_chi_angle_dev );

	TR.Info << "post-finalmin: (pos,rmsd,avg-bb,max-bb,avg-chi,max-chi,score) " <<
		F(9,3,all_atom_rmsd(pose,start_pose)) <<
		F(9,3,avg_bb_angle_dev) << F(9,3,max_bb_angle_dev) <<
		F(9,3,avg_chi_angle_dev) << F(9,3,max_chi_angle_dev) <<
		F(12,3,scorefxn( pose ) ) << std::endl;

	scorefxn.show( TR.Info, pose );
	TR.Info << std::endl;

	TR.flush();
} // basic_idealize

} // namespace idealize
} // namespace protocols
