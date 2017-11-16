// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/CombinePoseMover.cc
/// @brief
/// @author Hahnbeom Park

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/simple_moves/CombinePoseMover.hh>

#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <core/scoring/rms_util.hh>
#include <utility/string_util.hh>
#include <numeric/random/random.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace simple_moves {

using namespace core;

static basic::Tracer TR("protocols.simple_moves.CombinePoseMover");

CombinePoseMover::CombinePoseMover( core::scoring::ScoreFunctionCOP sfxn,
	pose::Pose const &pose_ref )
{
	set_default();
	pose_ref_ = pose_ref;
	sfxn_ = sfxn->clone();
	pose_tag_ = "combine";

	// sfxn is not used currently; will be used as minimization added
}

CombinePoseMover::~CombinePoseMover()= default;

void
CombinePoseMover::set_default(){
	max_struct_ = 1;
	max_struct_try_ = 5;
	vdwcut_ = 10.0;
	minfrac_crossover_ = 0.1;
	maxfrac_crossover_ = 0.4;
	rmsdcut_ = 2.0;

	cartesian_crossover_ = false;
	do_minimize_ = false;

	nonideal_ = false;
	sfxn0_ = scoring::ScoreFunctionFactory::create_score_function( "score0" );
	store_silents_ = false;

}

void
CombinePoseMover::apply( pose::Pose &pose2 )
{

	protocols::moves::MoverOP tocen = protocols::moves::MoverOP(
		new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );

	// Fraction should be like below of course
	runtime_assert( maxfrac_crossover_ <= 1.0 );
	runtime_assert( minfrac_crossover_ <= maxfrac_crossover_ );
	runtime_assert( minfrac_crossover_ >= 0.0 );
	runtime_assert( pose_ref_.size() == pose2.size() );

	sampled_structures_.resize( 0 );

	Size const nres( pose_ref().size() );
	Size models_failed( 0 );
	Size models_built( 0 );

	// Clear combined residue indices
	clear_combine_history();

	Real rmsd0, minrmsd, maxrmsd;
	// superimpose
	rmsd0 = scoring::calpha_superimpose_pose( pose2, pose_ref() );
	minrmsd = 0.5*rmsd0;
	maxrmsd = rmsdcut_; // initial
	TR << "Initial RMSD between two poses, set min/max rmsd as: ";
	TR << std::setw(8) << rmsd0 << " " << std::setw(8) << minrmsd << std::endl;

	if ( rmsd0 < 0.5 ) {
		TR << "Too close to each other; return without attempts." << std::endl;
		return;
	}

	pose::Pose pose0( pose_ref() );
	pose::Pose pose_min( pose_ref() );

	Size i_struct( 0 );
	Size ntry( 0 );
	Size n_seg_search( 0 );

	while ( true ) {

		// to prevent from running forever...
		//if( ntry > max_struct()*3 )  0.0; // just send anything

		// Get fraction of crossover
		Real const f_crossover
			= numeric::random::rg().uniform()*(maxfrac_crossover_ - minfrac_crossover_) + minfrac_crossover_;

		// Get Starting/Ending seqpos
		Size n_crossover = (Size)(nres*f_crossover);

		Size res1( numeric::random::rg().uniform()*( nres - n_crossover ) );
		Size res2( res1 + n_crossover - 1 );

		// make sure
		Size const start( res1 >= 1 ? res1 : 1 );
		Size const end( res2 <= nres ? res2 : nres );

		// scan through to get deviation part, filter out if the region is too similar to each other
		core::Real seg_meand( 0.0 );
		for ( core::Size ires = start; ires <= end; ++ires ) {
			Vector crd1 = pose0.residue( ires ).xyz( " CA " );
			Vector crd2 = pose2.residue( ires ).xyz( " CA " );
			seg_meand += crd1.distance(crd2);
		}
		seg_meand /= n_crossover;
		n_seg_search++;
		if ( seg_meand < (minrmsd*nres)/n_crossover && n_seg_search < 10 ) {
			//TR.Debug << "segment dist " << seg_meand << ", too similar: " << res1 << " " << res2 << ", skip" << std::endl;
			continue;
		} else {
			n_seg_search = 0;
			TR.Debug << "segment dist/nsearch " << seg_meand << " / " << n_seg_search;
			TR.Debug << " : " << start << " " << end << ", go!!" << std::endl;
		}

		i_struct++;
		ntry++;

		// Replace residues
		pose::Pose newpose( pose_ref() );

		std::vector< Size > combineres;
		combineres.push_back( start ); combineres.push_back( end );
		append_combine_history( combineres );

		for ( Size ires = start; ires <= end; ++ ires ) {
			if ( cartesian_crossover_ ) {
				// Copy cartesian coordinate
				for ( Size j = 1; j <= newpose.residue(ires).natoms(); ++j ) {
					Vector xyz( pose2.residue( ires ).xyz( j ) );
					core::id::AtomID id( j, ires );
					newpose.set_xyz( id, xyz );
				}

			} else {
				// Copy internal coordinates
				newpose.set_phi  ( ires, pose2.phi( ires ) );
				newpose.set_psi  ( ires, pose2.psi( ires ) );
				newpose.set_omega( ires, pose2.omega( ires ) );
				for ( Size ichi = 1; ichi <= newpose.residue_type( ires ).nchi(); ++ichi ) {
					newpose.set_chi( ichi, ires, pose2.chi( ichi, ires ) );
				}
			}
		}

		// Should we do this only in centroid?
		//if( !newpose.is_centroid() ) tocen->run( newpose );

		// Check if there is very severe clash
		// convert into centroid before scoring
		pose::Pose newpose_cen( newpose );
		if ( !newpose.is_centroid() ) tocen->apply( newpose_cen );

		Real score = sfxn0_->score( newpose_cen );
		Real rmsd = scoring::CA_rmsd( newpose, pose0 );
		Real rmsd2 = scoring::CA_rmsd( newpose, pose2 );

		if ( rmsd < minrmsd || rmsd2 < minrmsd ) {
			TR.Debug << "Reject, too similar to either structure, rmsd to 1/2/cut: ";
			TR.Debug << std::setw(8) << rmsd << " " << std::setw(8) << rmsd2;
			TR.Debug << " " << std::setw(8) << maxrmsd << std::endl;
			continue;
		}

		// Use dynamic vdw cut depending on the size of segment
		Real const vdwcut( n_crossover*2.0 );

		if ( score > vdwcut || rmsd > maxrmsd ) {
			models_failed++;
			TR.Debug << "i/max/start/end: " << i_struct << " / " << max_struct_try() << " " << start << " " << end;
			TR.Debug << ", failed: vdW ";
			TR.Debug << score << " ? " << vdwcut;
			TR.Debug << ", rmsd: " << rmsd << " ? " << maxrmsd << " " << std::endl;
			continue;
		} else {
			TR << "i/max/start/end: " << i_struct << " / " << max_struct_try() << " " << start << " " << end;
			TR << ", success(vdW/cut): ";
			TR << score << " < " << vdwcut << " ";
			TR << ", rmsd: " << rmsd << std::endl;
		}

		// Short fast minimization
		if ( do_minimize_ ) {
			optimization::MinimizerOptions minopt( "lbfgs_armijo_nonmonotone", 0.02, true, false, false );
			minopt.max_iter( 50 );
			optimization::AtomTreeMinimizer minimizer;
			kinematics::MoveMap mm_loc;
			mm_loc.set_jump( true ); mm_loc.set_bb( true ); mm_loc.set_chi( true );
			minimizer.run( newpose, mm_loc, *sfxn_, minopt );
		}

		// Store
		if ( store_silents_ ) {
			core::io::silent::SilentFileOptions opts;
			io::silent::SilentStructOP new_struct = nonideal_ ?
				io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary", opts) :
				io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( opts );

			new_struct->fill_struct( newpose );    // make the silent struct from the copy pose
			// Set donorhistory
			std::string donorhistory = new_struct->get_comment("donorhistory");
			donorhistory = donorhistory + utility::to_string( pose_tag() )
				+ "_" + utility::to_string( start  )
				+ "_" + utility::to_string( end );

			new_struct->erase_comment( "donorhistory" );
			new_struct->add_comment( "donorhistory", donorhistory );
			sampled_structures_.push_back( new_struct );
		}

		models_built++;

		if ( models_built >= max_struct() ) break;

		// if not enough structures sampled, increase rmsdcut
		if ( i_struct%max_struct_try() == 0 ) {
			maxrmsd = 1.0;
			minrmsd = 0.0;
		}
	}

	// Return minimium energy pose
	pose2 = pose_min;
}

} // namespace simple_moves
} // namespace protocol
