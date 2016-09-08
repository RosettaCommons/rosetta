// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/relax/Ramady.cc
/// @brief Header for the Rana energy repair code, Ramady
/// @author Mike Tyka


#include <protocols/relax/Ramady.hh>

#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/random/random.hh>

//Auto Headers
#include <core/kinematics/Jump.hh>

////////////////////////////////////////////////////////////////////////////////////////////////////

static THREAD_LOCAL basic::Tracer TR( "protocols.relax.Ramady" );

namespace protocols {
namespace relax {
////////////////////////////////////////////////////////////////////////////////////////////////////


void add_coordinate_constraints_to_pose( core::pose::Pose & pose, const core::pose::Pose &constraint_target_pose,  protocols::loops::Loops &exclude_regions ){
	using namespace core;
	using namespace conformation;
	using namespace pose;
	using namespace scoring;
	using namespace constraints;
	using namespace id;
	using namespace kinematics;
	using namespace moves;

	core::Size nnonvrt_cst_target = constraint_target_pose.size();
	core::Size nnonvrt_pose = pose.size();

	while ( pose.residue( nnonvrt_pose ).aa() == core::chemical::aa_vrt ) { nnonvrt_pose--; }
	while ( constraint_target_pose.residue( nnonvrt_cst_target ).aa() == core::chemical::aa_vrt ) { nnonvrt_cst_target--; }

	protocols::loops::Loops coordconstraint_segments;
	coordconstraint_segments = exclude_regions.invert( nnonvrt_cst_target );

	//TR << coordconstraint_segments << std::endl;

	if ( nnonvrt_pose != nnonvrt_cst_target ) {
		std::cerr << "ERROR coord constraint pose length mismatch with input pose: " << nnonvrt_cst_target << " vs. " << nnonvrt_pose << std::endl;
		utility_exit();
	}

	if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
		pose.append_residue_by_jump
			( *ResidueFactory::create_residue( pose.residue(1).residue_type_set()->name_map( "VRT" ) ),
			pose.size()/2 );
	}


	Size nres = pose.size();
	Real const coord_sdev( 0.5 );
	for ( Size i = 1; i<= (Size)nres; ++i ) {
		if ( i==(Size)pose.fold_tree().root() ) continue;
		if ( coordconstraint_segments.is_loop_residue( i ) ) {
			Residue const & nat_i_rsd( pose.residue(i) );
			for ( Size ii = 1; ii<= nat_i_rsd.last_backbone_atom(); ++ii ) {
				func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, coord_sdev ) );
				pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new CoordinateConstraint(
					AtomID(ii,i), AtomID(1,nres), nat_i_rsd.xyz( ii ),
					fx ) ) ) );
			}
		}
	}


}


bool rama_list_pred( const std::pair < core::Size, core::Real > &left, const std::pair < core::Size, core::Real > &right )
{
	return left.second > right.second;
}

void fix_worst_bad_ramas( core::pose::Pose & original_pose, core::Size how_many, core::Real skip_prob, core::Real limit_RMS, core::Real limit_rama  ){

	using namespace core;
	using namespace id;
	using namespace optimization;
	using namespace protocols::moves;
	using namespace core::scoring;
	using namespace core::pose;
	using namespace conformation;
	using namespace kinematics;

	//const core::Real limit_RMS = 0.5;
	//const core::Real limit_rama = 2.0;
	//const core::Real limit_rama_min = 2.0;
	const core::Real limit_rama_min = limit_rama;


	// Original RAFT set
	//-140  153 180 0.135 B
	// -72  145 180 0.155 B
	//-122  117 180 0.073 B
	// -82  -14 180 0.122 A
	// -61  -41 180 0.497 A
	//  57   39 180 0.018 L

	core::Real post_fix_rama_e=0;
	core::Real change_rms=0 ;
	core::Real ok_phi[] =  { -140, -72, -122, -82, -61, 57 } ;
	core::Real ok_psi[] =  {  153, 145,  117, -14, -41, 39 } ;
	core::pose::Pose pose = original_pose;
	core::pose::Pose temp_pose = pose;
	protocols::loops::Loops exclude_regions;
	add_coordinate_constraints_to_pose( pose, original_pose, exclude_regions ) ;
	core::scoring::ScoreFunction rama_scorefxn;
	rama_scorefxn.set_weight( coordinate_constraint, 0.2 );
	rama_scorefxn.set_weight( rama, 1.0 );
	rama_scorefxn.set_weight( omega, 0.2 );
	Energies & energies( pose.energies() );
	rama_scorefxn(pose); //apply score

	MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.02, true /*use_nblist*/, false /*deriv_check*/ );
	kinematics::MoveMap final_mm;
	final_mm.set_bb(true);


	std::vector < std::pair < core::Size, core::Real > > rama_list;

	for ( Size j=1; j<= pose.size(); ++j ) {
		EnergyMap & emap( energies.onebody_energies( j ) );
		if (  emap[ rama ] > limit_rama_min ) {
			rama_list.push_back( std::make_pair( j, emap[ rama ] ) );
		}
	}

	if ( rama_list.size() == 0 ) return;
	// sort the rama_list

	std::sort(rama_list.begin(), rama_list.end(), rama_list_pred);

	for ( auto & g : rama_list ) {
		TR << "RAMALIST: " << g.first << "  " << g.second << std::endl;
	}

	for ( Size g=0; g< how_many; g++ ) {

		if ( g >= rama_list.size() ) break;

		if ( numeric::random::uniform() < skip_prob ) continue;

		core::Size i = rama_list[g].first;
		EnergyMap & emap( energies.onebody_energies( i ) );
		core::Real rama_e =  emap[ rama ];

		TR << "RAMALIST: " << rama_list[g].first << "  " << rama_list[g].second << " RAMA: " << rama_e << " RES: " << i << std::endl;
		// save the angles
		for ( core::Size ir = 1; ir < pose.size(); ir ++ ) {
			temp_pose.set_phi(   ir, pose.phi(   ir ) );
			temp_pose.set_psi(   ir, pose.psi(   ir ) );
			temp_pose.set_omega( ir, pose.omega( ir ) );
		}
		std::vector < core::Size > used_angles;

		while ( used_angles.size() < 6 ) {


			// pick a random reasonable phi/spi pair

			//pose.dump_pdb("ramarep_" + utility::to_string( i ) + "_" + "pre.pdb" );

			const core::Size ok_angles = 6;
			core::Real curphi = pose.phi(i);
			core::Real curpsi = pose.psi(i);
			core::Real newphi = ok_phi[0];
			core::Real newpsi = ok_psi[0];
			core::Real bestdist = 1000000;
			core::Size bestindex = 0;
			TR << "S:" << curphi << " " << curpsi << "  " << std::endl;
			for ( core::Size a = 0; a < ok_angles; a++ ) {
				core::Real diffphi =  ok_phi[a]  - curphi; while ( diffphi > 180 ) diffphi-=360.0;  while ( diffphi < -180 ) diffphi += 360.0;
				core::Real diffpsi =  ok_psi[a]  - curpsi; while ( diffpsi > 180 ) diffpsi-=360.0;  while ( diffpsi < -180 ) diffpsi += 360.0;
				core::Real dist = sqrt( diffphi*diffphi + diffpsi * diffpsi );
				TR << "T:" << ok_phi[a] << "  " << ok_psi[a] << "  " << dist << std::endl;
				if ( (dist < bestdist) || (a == 0) ) {
					if ( std::find(used_angles.begin(), used_angles.end(), a)!=used_angles.end() ) continue;
					newphi = ok_phi[a];
					newpsi = ok_psi[a];
					bestdist = dist;
					bestindex = a;
				}
			}
			used_angles.push_back( bestindex );
			TR << "B:" << newphi << "  " << newpsi << "  " << bestdist << "BEST: " << bestindex << std::endl;

			pose.set_phi( i, newphi );
			pose.set_psi( i, newpsi );
			//pose.dump_pdb("ramarep_" + utility::to_string( i ) + "_" + "mid.pdb" );

			kinematics::MoveMap local_mm;
			local_mm.set_bb(true);
			for ( int ii=-4; ii< 3; ii ++ ) {
				int res = i + ii;

				if ( res < 1 ) continue;
				if ( res > (int)pose.size() ) continue;
				local_mm.set_bb( res, true );
				local_mm.set( TorsionID( omega_torsion, BB, res), false );
				// disallow proline PHI
				if ( pose.residue(res).aa() == chemical::aa_pro ) local_mm.set( TorsionID( phi_torsion, BB, res), false );
			}

			AtomTreeMinimizer().run( pose, local_mm, rama_scorefxn, options );
			AtomTreeMinimizer().run( pose, final_mm, rama_scorefxn, options );

			// Now check score
			rama_scorefxn(pose); //apply score

			post_fix_rama_e =  emap[ rama ];
			change_rms = core::scoring::CA_rmsd( temp_pose, pose );

			TR << "RMS:" << change_rms << std::endl;


			// break out if structure is good
			if ( ( post_fix_rama_e < limit_rama ) &&
					( change_rms < limit_RMS ) ) {
				break;
			}

			// otherwise restore and continue
			for ( core::Size ir = 1; ir < pose.size(); ir ++ ) {
				pose.set_phi(   ir, temp_pose.phi(   ir ) );
				pose.set_psi(   ir, temp_pose.psi(   ir ) );
				pose.set_omega( ir, temp_pose.omega( ir ) );
			}
		}


		TR << "RAMADY " << rama_e << "  " << post_fix_rama_e << "  " << change_rms << std::endl;


		//pose.dump_pdb("ramarep_" + utility::to_string( i ) + "_" + "post.pdb" );
	}


	// Final RamaCheck
	rama_scorefxn(pose); //apply score
	for ( Size i=1; i<= pose.size(); ++i ) {
		if ( !pose.residue_type(i).is_protein() ) continue;
		EnergyMap & emap( energies.onebody_energies( i ) );
		TR << "CHECK: " << i << "  "
			<< pose.phi(i) << "  "
			<< pose.psi(i) << "  "
			<< emap[ rama ] << std::endl;
	}

	// save the angles
	for ( core::Size ir = 1; ir < original_pose.size(); ir ++ ) {
		if ( !pose.residue_type(ir).is_protein() ) continue;
		original_pose.set_phi(   ir, pose.phi(   ir ) );
		original_pose.set_psi(   ir, pose.psi(   ir ) );
		original_pose.set_omega( ir, pose.omega( ir ) );
	}

}


}

}
