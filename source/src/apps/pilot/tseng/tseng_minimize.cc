// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Will Sheffler
/// @brief

// 8/26/2008 by Paul Tseng
// This code tries to "learn" the local minima region of the Rosetta energy landscape,
// starting from the input structure.  In particular, it checks for
// canyon-like landscape, which may be useful for global optimization.
// This is done by starting at the input structure (e.g., experimentally
// determined native structure) and
// perturbing this pocore::Size and then do local optimization to see if we get a
// different local min with nearly equal energy and is not too far away
// (i.e., another pocore::Size on the canyon floor).  If yes, then repeat with this
// local min.
// By varying the direction of perturbation and the size of perturbation,
// we hope to "learn" the shape of the canyon floor and the steepness of the
// walls around it.


// libRosetta headers
// AUTO-REMOVED #include <protocols/simple_moves/ScoreMover.hh>

// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>

// AUTO-REMOVED #include <protocols/jobdist/standard_mains.hh>

#include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/moves/NullMover.hh>
// AUTO-REMOVED #include <protocols/simple_moves/MinMover.hh>
// AUTO-REMOVED #include <protocols/simple_moves/ProlineFixMover.hh>
// AUTO-REMOVED #include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/BackboneMover.hh>
// AUTO-REMOVED #include <protocols/relax_protocols.hh>
#include <core/scoring/dssp/Dssp.hh>

// AUTO-REMOVED #include <core/io/silent/silent.fwd.hh>
// AUTO-REMOVED #include <core/io/silent/ProteinSilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/excn/Exceptions.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/Tracer.hh>
#include <devel/init.hh>

#include <utility/io/izstream.hh>

#include <numeric/random/random.hh>
#include <utility/io/ozstream.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <basic/options/option.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jumping/StrandPairing.hh>




using namespace core;
using namespace basic::options;
using namespace ObjexxFCL;

static thread_local basic::Tracer TR( "protocols.moves.ScoreMover" );

namespace score_app { BooleanOptionKey linmin( "score_app:linmin" ); }

void test( std::string fname ) {

	using namespace protocols;
	using namespace protocols::jobdist;
	using namespace protocols::moves;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::file;
	using namespace core;
	using namespace optimization;
	using namespace scoring;
	using namespace io::pdb;

	ScoreFunctionOP sfxnOP( get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS ) );
	ScoreFunction & sfxn( *sfxnOP );

	AtomTreeMinimizer minimizer;

	core::kinematics::MoveMap movemap;
	movemap.set_bb( true );
	movemap.set_chi( true );

	core::pose::Pose pose;
	core::import_pose::pose_from_pdb(pose,fname);

	Real eps,meps,maxpert,lmdel,pert,diff,score;
	core::Size ntrial,trial,nlm,k,chino;
    // NM is the maximum number of local min to be stored
	core::Size const NM( 1000 );
	FArray1D<Real> sflm( NM );
    // N is the number of residues of the protein
	core::Size const N( pose.total_residue() );
	FArray2D<Real> lm_phi( N, NM );
	FArray2D<Real> lm_psi( N, NM );
	FArray2D<Real> lm_chi( N, NM );
	FArray2D<Real> lm_omega( N, NM );
	FArray1D<Real> x1_phi( N );
	FArray1D<Real> x1_psi( N );
	FArray1D<Real> x1_chi( N );
	FArray1D<Real> x1_omega( N );

    // Prcore::Size out secondary structure info, suggested by James Thompson
	for ( core::Size i = 1; i <= N; ++i )
		std::cout << pose.secstruct(i);
	std::cout << std::endl;
	core::scoring::dssp::Dssp dssp_obj( pose );
	dssp_obj.insert_ss_into_pose( pose );
	for ( core::Size i = 1; i <= N; ++i )
		std::cout << pose.secstruct(i);
	std::cout << std::endl;

    // Lines suggested by James Thompson (
	Size nmoves_( 20 );
	Real  m_Temperature = 0.8;
	core::kinematics::MoveMapOP movemap_ = new core::kinematics::MoveMap();
	simple_moves::SmallMoverOP small_mover( new simple_moves::SmallMover( movemap_,
                                                                         m_Temperature, nmoves_ ) );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 2.0 );
	small_mover->angle_max( 'L', 3.0 );

    // Run local optimization to get to a local min

	std::cout << "Energy of input structure: " << sfxn(pose) <<  " #residues: " << N <<std::endl;
	minimizer.run( pose, movemap, sfxn, MinimizerOptions( "dfpmin_armijo_nonmonotone", 0.001, true ) );
	std::cout << "Energy after local minimization:  " << sfxn(pose) << std::endl;

	ntrial=100;
    // A trial means perturbing phi,psi,chi,omega in a structure and do local optimization.
    // maxpert = maximum amount of pertubation on phi, psi, chi, omega per trial (e.g., 5 degrees)
	maxpert=5.;

    // Parameter for checking if a local min is on the canyon floor.
	eps=0.5;		// eps must be +ve and not too large.
	meps=0.01;

	chino=1;
	trial=0;
	nlm=1;
	for ( core::Size i = 1; i <= N; ++i ) {
		lm_phi(i,nlm)=pose.phi(i);
		lm_psi(i,nlm)=pose.psi(i);
        if ( chino <= pose.residue( i ).nchi() ){
			lm_chi(i,nlm)=pose.chi(chino,i);
		}
		lm_omega(i,nlm)=pose.omega(i);
	}
	sflm(1)=sfxn(pose);

L1:
	if (trial < ntrial) {
		++trial;

        // Perturb current local min x0 to x1 and then minimize to get a "new" local min

		for ( core::Size i = 1; i <= N; ++i ) {
			x1_phi(i)=lm_phi(i,nlm)+maxpert*(numeric::random::rg().uniform()*2.-1.);
			x1_psi(i)=lm_psi(i,nlm)+maxpert*(numeric::random::rg().uniform()*2.-1.);
			x1_chi(i)=lm_chi(i,nlm)+maxpert*(numeric::random::rg().uniform()*2.-1.);		//perturb chi by multiples of 60 deg?
			x1_omega(i)=lm_omega(i,nlm);  //  +maxpert*(numeric::random::rg().uniform()*2.-1.);
			pose.set_phi(i, x1_phi(i) );
			pose.set_psi(i, x1_psi(i) );
            if ( chino <= pose.residue( i ).nchi() ){
				pose.set_chi(chino, i, x1_chi(i) );
			}
			pose.set_omega(i, x1_omega(i) );
		}

		std::cout << "trial: " << trial << " #local min: " << nlm << " energy after perturbation: " << sfxn(pose) << std::endl;
		minimizer.run( pose, movemap, sfxn, MinimizerOptions( "dfpmin_armijo_nonmonotone", 0.001, true ) );

        // Check if "new" local min is near x0 and its energy is not too high (so it's in the valley floor)
		lmdel=0.0;
		pert=0.0;
		for ( core::Size i = 1; i <= N; ++i ) {
            lmdel += (pose.phi(i)-lm_phi(i,nlm))*(pose.phi(i)-lm_phi(i,nlm));
            lmdel += (pose.psi(i)-lm_psi(i,nlm))*(pose.psi(i)-lm_psi(i,nlm));
            lmdel += (pose.chi(chino,i)-lm_chi(i,nlm))*(pose.chi(chino,i)-lm_chi(i,nlm));
            lmdel += (pose.omega(i)-lm_omega(i,nlm))*(pose.omega(i)-lm_omega(i,nlm));
            pert += (x1_phi(i)-lm_phi(i,nlm))*(x1_phi(i)-lm_phi(i,nlm));
            pert += (x1_psi(i)-lm_psi(i,nlm))*(x1_psi(i)-lm_psi(i,nlm));
            pert += (x1_chi(i)-lm_chi(i,nlm))*(x1_chi(i)-lm_chi(i,nlm));
            pert += (x1_omega(i)-lm_omega(i,nlm))*(x1_omega(i)-lm_omega(i,nlm));
		}
		lmdel=std::sqrt(lmdel);
		pert=std::sqrt(pert);
		score=sfxn(pose);

		std::cout << " perturb level: " << pert << " energy after min: " << sfxn(pose) << " change in local min: " << lmdel <<  std::endl;


		if (lmdel<(1+eps)*pert && score<sflm(nlm)+eps*lmdel) {

            // Check if local min has already been found
			k=1;
			Real mindiff_phi=1.0E30f; Real mindiff_psi=1.0E30f; Real mindiff_chi=1.0E30f; Real mindiff_omega=1.0E30f;
			Real maxdiff_phi=0.0; Real maxdiff_psi=0.0; Real maxdiff_chi=0.0; Real maxdiff_omega=0.0;
        L2:
			if ( k <= nlm ) {
				Real diff_phi=0.0; Real diff_psi=0.0; Real diff_chi=0.0; Real diff_omega=0.0;
				for ( core::Size i = 1; i <= N; ++i ) {
                    diff_phi += (pose.phi(i)-lm_phi(i,k))*(pose.phi(i)-lm_phi(i,k));
                    diff_psi += (pose.psi(i)-lm_psi(i,k))*(pose.psi(i)-lm_psi(i,k));
                    diff_chi += (pose.chi(chino,i)-lm_chi(i,k))*(pose.chi(chino,i)-lm_chi(i,k));
                    diff_omega += (pose.omega(i)-lm_omega(i,k))*(pose.omega(i)-lm_omega(i,k));
				}
				diff=std::sqrt(diff_phi+diff_psi+diff_chi+diff_omega);
				mindiff_phi=std::min(mindiff_phi,std::sqrt(diff_phi));
				mindiff_psi=std::min(mindiff_psi,std::sqrt(diff_psi));
				mindiff_chi=std::min(mindiff_chi,std::sqrt(diff_chi));
				mindiff_omega=std::min(mindiff_omega,std::sqrt(diff_omega));
				maxdiff_phi=std::max(maxdiff_phi,std::sqrt(diff_phi));
				maxdiff_psi=std::max(maxdiff_psi,std::sqrt(diff_psi));
				maxdiff_chi=std::max(maxdiff_chi,std::sqrt(diff_chi));
				maxdiff_omega=std::max(maxdiff_omega,std::sqrt(diff_omega));
				if ( diff>meps ) {
					++k;
					goto L2;
				}
			}

            std::cout << " k= " << k << " nlm= " << nlm << std::endl;
            std::cout << " mindiff_phi= " << mindiff_phi << " mindiff_psi= " << mindiff_psi << " mindiff_chi= " << mindiff_chi << " mindiff_omega= " << mindiff_omega << std::endl;
            std::cout << " maxdiff_phi= " << maxdiff_phi << " maxdiff_psi= " << maxdiff_psi << " maxdiff_chi= " << maxdiff_chi << " maxdiff_omega= " << maxdiff_omega << std::endl;


			if ( k > nlm ) {
      			++nlm;
				for ( core::Size i = 1; i <= N; ++i ) {
					lm_phi(i,nlm)=pose.phi(i);
					lm_psi(i,nlm)=pose.psi(i);
                    if ( chino <= pose.residue( i ).nchi() ){
						lm_chi(i,nlm)=pose.chi(chino,i);
					}
					lm_omega(i,nlm)=pose.omega(i);
				}
      			sflm(nlm)=score;
				std::cout << " #canyon floor local min increases to " << nlm << "    Energy(new local min): " << sflm(nlm) << std::endl;
			}
		}
		goto L1;
	}

	std::cout << " #canyon floor local min found: " << nlm << std::endl;


    // write local min to a file, as suggested by Will Sheffler:
	utility::io::ozstream out( "lm_phi.out" );
	for( core::Size k = 1; k <= nlm; k++) {
  		for( core::Size i = 1; i <= N; i++) {
   			out << lm_phi(i,k) << ' ';
		}
  		out << std::endl;
	}
	out.close();


}

int
main( int argc, char * argv [] )
{
    try {
        devel::init( argc, argv );

        using namespace basic::options;
        using namespace basic::options::OptionKeys;
        using namespace utility;

        // test_io();

        // test_sasa_dots();

        // core::pose::Pose native_pose;

        if( option[ in::file::s ].user() ) {
            vector1<file::FileName> files( option[ in::file::s ]() );
            for( size_t i = 1; i <= files.size(); ++i ) {
                test( files[i] );
            }
        } else if( option[ in::file::l ].user() ) {
            vector1<file::FileName> files( option[ in::file::l ]() );
            for( size_t i = 1; i <= files.size(); ++i ) {
                utility::io::izstream list( files[i] );
                std::string fname;
                while( list >> fname ) {
                    // std::cerr << "'" << fname << "'" << std::endl;
                    test( fname );
                }
            }
        }


        return 0;
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cout << "caught exception " << e.msg() << std::endl;
        return -1;
    }
}
