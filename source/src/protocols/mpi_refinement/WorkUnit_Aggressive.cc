// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/mpi_refinement/WorkUnit_Sampler.cc
/// @brief
/// @author Hahnbeom Park: Generalized as a "Sampler" from "Loop Hasher"

//#include <protocols/mpi_refinement/WorkUnit_Sampler.hh>
#include <protocols/mpi_refinement/WorkUnit_Aggressive.hh>
#include <protocols/mpi_refinement/util.hh>
#include <protocols/wum/WorkUnitBase.hh>
#include <protocols/wum/SilentStructStore.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/rms_util.hh>
//#include <protocols/wum/SilentStructStore.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <protocols/relax/AtomCoordinateCstMover.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
//
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/CombinePoseMover.hh>
#include <protocols/normalmode/NormalModeRelaxMover.hh>
#include <protocols/normalmode/NormalModeRelaxMover.fwd.hh>

#include <core/pose/selection.hh>
#include <protocols/simple_moves/BackboneMover.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <numeric/random/random.hh>

//Auto Headers
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/excn/Exceptions.hh>
#include <utility> //for std::pair

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

namespace protocols {
namespace mpi_refinement {

static basic::Tracer TR("WorkUnit_Sampler.AggressiveType");

////////////////////////////////////////////
//////// WorkUnit Combine
WorkUnit_CombinePose::WorkUnit_CombinePose( core::Size nstruct,
	bool const cartesian
)
{
	set_defaults();
	cartesian? set_cartesian( 1 ) : set_cartesian( 0 );
	set_nstruct( nstruct );
}

void
WorkUnit_CombinePose::set_defaults(){}

void
WorkUnit_CombinePose::init_from_cmd( const core::Size )
{}

void
WorkUnit_CombinePose::run()
{
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( decoys().size() == 0 ) {
		TR << "Empty WorkUnit ! Cannot execute run() " << std::endl;
		return;
	}

	bool const nonideal( option[ OptionKeys::lh::bss]() );
	bool const cartesian( get_cartesian() );
	core::Size const nstruct( get_nstruct() );

	// Assert if mixing struct has been already sent
	runtime_assert( decoys().size() >= 2 );

	SilentStructCOP start_struct = decoys().get_struct(0);

	core::pose::Pose pose;
	core::pose::Pose pose2;
	decoys().get_pose( 0, pose );
	decoys().get_pose( 1, pose2 );
	std::string tag1( decoys().get_struct(0)->decoy_tag() );
	std::string tag2( decoys().get_struct(1)->decoy_tag() );

	// clear the store of structures
	decoys().clear();

	std::string combtype = cartesian? " cartesian": " torsion";
	TR << "Executing WorkUnit_CombinePose_Mover on " << combtype << "." << std::endl;

	core::Size starttime = time(nullptr);

	core::scoring::ScoreFunctionCOP sfxn_loc = get_energy( "talaris2013_cart" );

	protocols::simple_moves::CombinePoseMover comb( sfxn_loc, pose );
	comb.set_max_struct( nstruct );
	comb.set_max_try( nstruct*2 );
	comb.set_minfrac_crossover( minfrac_ );
	comb.set_maxfrac_crossover( maxfrac_ );
	comb.set_nonideal( nonideal );
	comb.set_cartesian( cartesian );
	comb.set_rmsdcut( 3.0 ); // don't deviate to much
	comb.set_store_silents( true );
	//comb.set_minimize( true ); // use rerelax instead
	comb.apply( pose2 );

	std::vector< SilentStructOP > decoys_out = comb.return_silent();
	TR << "combpose: How many? " << decoys_out.size() << "/" << nstruct << std::endl;

	// this happens sometimes when parents are too close... then just return with startings
	if ( decoys_out.size() < nstruct ) {
		for ( core::Size istruct = 1; istruct <= nstruct - decoys_out.size(); ++istruct ) {
			decoys().store().push_back( start_struct->clone() );
		}

	}

	// Add info here if you want
	// Note that this should go to batchrelax so no minimization/repack will be called here
	for ( core::Size istruct = 0; istruct < decoys_out.size(); ++ istruct ) {
		SilentStructOP ss = decoys_out[istruct];
		std::string tag = "comb_" + tag1 + "_" + tag2 + "_" + ObjexxFCL::string_of( istruct );
		store_to_decoys( start_struct, ss, tag );
	}

	core::Size endtime = time(nullptr);
	TR.Debug << "Build " << decoys().size() << " structures in ";
	TR.Debug << endtime - starttime << " s " << std::endl;

}

////////////////////////////////////////////
//////// WorkUnit NormalMode
WorkUnit_NormalMode::WorkUnit_NormalMode( core::Size const nmodes,
	core::Size const nmtype,
	core::Size const relaxtype,
	core::Real const maxscale )
{
	// nmtype     1: CartCen 2:TorsCen 3:CartFull  4:TorsFull
	// relaxtype  1: Cartmin 2:Torsmin 3:CartExpol 4:TorsExpol

	set_defaults();
	set_nmodes( nmodes );
	set_nmtype( nmtype ); //
	set_relaxtype( relaxtype );
	set_maxscale( maxscale );
}

void
WorkUnit_NormalMode::set_defaults(){}

void
WorkUnit_NormalMode::init_from_cmd( const core::Size )
{}

// The logic it takes options might look complicated because there are so many options to decide;
// # modes to use, solve normal mode in cart or torsion space, relax on cart or torsion space,
// scale to perturb, how to relax after extrapolating, relax on centroid or fullatom...
void
WorkUnit_NormalMode::run()
{
	using namespace protocols::normalmode;
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( decoys().size() == 0 ) {
		TR << "Empty WorkUnit ! Cannot execute run() " << std::endl;
		return;
	}

	// Get struct
	SilentStructCOP start_struct = decoys().get_struct(0);
	core::pose::Pose pose;
	decoys().get_pose( 0, pose );

	// clear the store of structures
	decoys().clear();

	TR << "Executing WorkUnit_NormalMode_Mover" << std::endl;

	core::Size starttime = time(nullptr);

	// 1. Scales/mode setup
	// Use maxmode + mixture of them by scale given as input
	core::Size const nmode( get_nmodes() ); // Pure modes
	core::Size const nmix( 0 );  // mixed modes
	core::Size const ncomb( nmode + nmix );
	utility::vector1< core::Size > modes; // Use top [nmode] modes for both pure/mixing
	for ( core::Size i = 1; i <= nmode; ++i ) modes.push_back( i );

	// distance
	core::Real const dist( option[ lh::NMdist ]() );

	utility::vector1< core::Real > scales; // line-search scales
	core::Real const s = get_maxscale();
	scales.push_back( -s ); scales.push_back( s );

	core::Size const nstruct( scales.size()*ncomb ); // 6 scales * (5 pure + 5 mixed modes)

	// Mode weights
	utility::vector1< utility::vector1< core::Real > > modescales( ncomb );
	// First fill in by pure eig-vectors
	utility::vector1< core::Real > nullv( nmode, 0.0 );
	core::Size i_comb( 1 );
	for ( ; i_comb <= nmode; ++i_comb ) {
		utility::vector1< core::Real > scalev( nullv ); scalev[i_comb] = 1.0;
		modescales[i_comb] = scalev;
	}

	// Add up remainings by mixing
	for ( ; i_comb <= ncomb; ++i_comb ) {
		utility::vector1< core::Real > scalev( nmode, 0.0 );
		for ( core::Size i = 1; i <= nmode; ++i ) scalev[i] = numeric::random::rg().uniform();
		modescales[i_comb] = scalev;
	}

	// 2.Set relaxtype
	// relaxtype  1: Cartmin 2:Torsmin 3:CartExpol 4:TorsExpol
	std::string relaxmode;
	if ( get_relaxtype() <= 2 ) {
		relaxmode = "min";
	} else if ( get_relaxtype() <= 4 ) {
		relaxmode = "extrapolate_and_relax";
	} else {
		TR << "Unknown relaxmode: " << get_relaxtype() << "! using default relaxtype = min." << std::endl;
		relaxmode = "min";
	}

	bool cartmin( false );
	if ( get_relaxtype() == 1 || get_relaxtype() == 3 ) cartmin = true;

	// 3. Set normalmode type
	// nmtype   1: CartCen 2:TorsCen 3:CartFull  4:TorsFull
	bool cartnm( false );
	bool iscen( false );
	if ( get_nmtype() <= 2 ) iscen = true;
	if ( get_nmtype() == 1 || get_nmtype() == 3 ) cartnm = true;

	core::scoring::ScoreFunctionOP sfxn_loc = ( iscen )?
		core::scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart" ):
		core::scoring::ScoreFunctionFactory::create_score_function( "talaris2013" );

	core::scoring::ScoreFunctionOP sfxn_pack
		= core::scoring::ScoreFunctionFactory::create_score_function( "soft_rep" );

	if ( cartmin && sfxn_loc->get_weight( core::scoring::cart_bonded ) == 0.0 ) {
		sfxn_loc->set_weight( core::scoring::cart_bonded, 0.5 ); //make sure!
		sfxn_loc->set_weight( core::scoring::pro_close, 0.0 ); //make sure!
	}

	if ( sfxn_loc->get_weight( core::scoring::elec_dens_fast ) > 0.0 ) {
		TR << "Sampling with elec_dens_fast : " << sfxn_loc->get_weight( core::scoring::elec_dens_fast ) << std::endl;
	}

	// change pose level if necessary
	protocols::moves::MoverOP tofa
		( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD ) );
	protocols::moves::MoverOP tocen
		( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );

	if ( iscen && !pose.is_centroid() ) {
		tocen->apply( pose );
	} else if ( !iscen && pose.is_centroid() ) {
		tofa->apply( pose );
	}

	// finally get normal mode instance
	core::kinematics::MoveMapOP mm = get_movemap( pose, "pack", false );
	/*
	NormalModeRelaxMoverOP NM = ( cartnm )?
	get_NMmover( pose, sfxn_loc, mm, dist, relaxmode, true ) :  // return CNM
	get_NMmover( pose, sfxn_loc, mm, dist, relaxmode, false );  // return TNM
	*/
	NormalModeRelaxMoverOP NM( new NormalModeRelaxMover( sfxn_loc, cartnm, mm, relaxmode, dist ) );

	core::optimization::MinimizerOptionsOP minoption
		( new core::optimization::MinimizerOptions( "lbfgs_armijo_nonmonotone",
		0.001, true, false, false ) );
	minoption->max_iter( 20 );

	NM->set_cartesian_minimize( cartmin );
	NM->set_minoption( minoption );

	core::optimization::AtomTreeMinimizer minimizer;
	protocols::relax::AtomCoordinateCstMover coord_cst_mover;
	coord_cst_mover.cst_sd( 1.0 );

	// 4. Run! (will generate ncomb*6 scales)
	for ( i_comb = 1; i_comb <= ncomb; ++i_comb ) {
		utility::vector1< core::Real > &modescale = modescales[i_comb];

		NM->set_mode( modes, modescale );

		for ( core::Size i = 1; i <= scales.size(); ++i ) {
			core::pose::Pose pose_tmp( pose );
			core::Real const scale = scales[i];

			NM->set_extrapolate_scale( scale );
			NM->apply_on_pose( pose_tmp );

			// Always store in full atom!
			if ( pose_tmp.is_centroid() ) {
				tofa->apply( pose_tmp );
				// at least repack before storing them....
				repack( pose_tmp, sfxn_pack );

				// this is to bring better discrimination by GOAP
				// minimizing can slow down overall procedure
				// also think about moving too much; this is TorsionMinimizer
				// minimize with coordinate cst if you really want to

				if ( option[ lh::minimize_after_nmsearch ]() ) {
					coord_cst_mover.apply( pose );
					minimizer.run( pose_tmp, *mm, *sfxn_pack, *minoption );
				}
			}

			//ss->add_energy( "NMmode", i_comb );
			//ss->add_energy( "NMscale", scales[i] );
			store_to_decoys( start_struct, pose_tmp );
		}

		// Just to make sure
		if ( decoys().store().size() >= nstruct ) break;
	}

	core::Size endtime = time(nullptr);
	TR.Debug << "Build " << decoys().size() << " structures in ";
	TR.Debug << endtime - starttime << " s " << std::endl;

}

////////////////////////////////////////////
//////// WorkUnit RamaPerturber
WorkUnit_RamaPerturber::WorkUnit_RamaPerturber( core::Size const nsteps,
	core::Size const res1,
	core::Size const res2,
	core::Real const kT
)
{
	set_defaults();
	set_nsteps( nsteps );
	set_res1( res1 );
	set_res2( res2 );
	set_kT( kT );
}

void
WorkUnit_RamaPerturber::set_defaults(){}

void
WorkUnit_RamaPerturber::init_from_cmd( const core::Size )
{
}

void
WorkUnit_RamaPerturber::run()
{
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( decoys().size() == 0 ) {
		TR << "Empty WorkUnit ! Cannot execute run() " << std::endl;
		return;
	}

	SilentStructCOP start_struct = decoys().get_struct(0);

	core::pose::Pose pose0;
	decoys().get_pose( 0, pose0 );
	runtime_assert( pose0.is_fullatom() );
	// Assign secstruc by running dssp
	core::scoring::dssp::Dssp dssp( pose0 );
	dssp.insert_ss_into_pose( pose0 );

	protocols::moves::MoverOP tofa
		( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD ) );
	protocols::moves::MoverOP tocen
		( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );

	// clear the store of structures
	decoys().clear();

	// get number of runs
	core::Size starttime = time(nullptr);

	utility::vector1< core::Size > loopres;
	loopres.resize( 0 );
	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_jump(false ); mm->set_bb( false ); mm->set_chi( false );
	mm->set( core::id::THETA, false ); mm->set( core::id::D, false );

	for ( core::Size i = get_res1(); i <= get_res2(); ++i ) {
		loopres.push_back( i );
		mm->set_bb( i, true );
		mm->set_chi( i, true );
	}

	if ( loopres.size() == 0 ) {
		TR << "Empty loop region! nothing to execute." << std::endl;
		return;
	}

	// Setup perturber: bigger pert on Loops; the purpose is to preserve SecStruct
	core::Real rama_kbT( 99.0 );
	simple_moves::ShearMover bbsampler( mm, rama_kbT, 1 ); //perturb once every time
	bbsampler.angle_max('H',   5.0);
	bbsampler.angle_max('E',  15.0);
	bbsampler.angle_max('L', 180.0);

	// Score setup: Cen, downweight hbond & vdw
	core::scoring::ScoreFunctionOP sfxn_cen
		= core::scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart" );
	sfxn_cen->set_weight( core::scoring::hbond_sr_bb, 0.2 );
	sfxn_cen->set_weight( core::scoring::hbond_lr_bb, 0.2 );
	sfxn_cen->set_weight( core::scoring::vdw, 0.05 ); //original 1.0, max clash is then 2~5

	core::scoring::ScoreFunctionOP sfxn_famin
		= core::scoring::ScoreFunctionFactory::create_score_function( "talaris2013" );
	sfxn_famin->set_weight( core::scoring::cart_bonded, 0.5 );

	// Run MC
	core::pose::Pose pose_work( pose0 );
	tocen->apply( pose_work );
	protocols::moves::MonteCarlo mc( pose_work, *sfxn_cen, get_kT() );
	//for( core::Size iiter = 1; iiter <= get_nsteps(); ++ iiter ){
	for ( core::Size iiter = 1; iiter <= 50; ++ iiter ) {
		bbsampler.apply( pose_work );
		mc.boltzmann( pose_work );
	}

	// Once done, minimize with original vdw
	core::optimization::MinimizerOptions minoption( "lbfgs_armijo_nonmonotone",
		0.001, true, false, false );
	core::optimization::AtomTreeMinimizer minimizer;
	minoption.max_iter( 50 );

	sfxn_cen->set_weight( core::scoring::vdw, 1.0 ); //original 1.0, max clash is the
	minimizer.run( pose_work, *mm, *sfxn_cen, minoption );

	// Finally relax in full atom
	tofa->apply( pose_work );
	// pose/loopres/sfxn/nonideal(true)/ramp(true)/eff(false)/envdist(0.0)
	ramp_minpack_loop2( pose_work, loopres, sfxn_famin, false, true, false, 6.0 );

	superimpose_to_ref( pose0, pose_work );

	store_to_decoys( start_struct, pose_work );

	core::Size endtime = time(nullptr);
	TR.Debug << "Build " << decoys().size() << " structures in " << endtime - starttime << " s " << std::endl;

} // class WorkUnit_RamaPerturber

}
}

