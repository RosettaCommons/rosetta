// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/mpi_refinement/WorkUnit_Sampler.cc
/// @brief
/// @author Mike Tyka
/// @author Hahnbeom Park: Generalized as a "Sampler" from "Loop Hasher"

#include <protocols/mpi_refinement/WorkUnit_Sampler.hh>
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
#include <core/scoring/rms_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

//#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <numeric/random/random.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

namespace protocols {
namespace mpi_refinement {

static basic::Tracer TR("WorkUnit_Sampler");

////////////////////////////////////////////
/////// Parent
core::kinematics::MoveMapOP
WorkUnit_Sampler::get_movemap( core::pose::Pose const &pose,
	std::string const mode,
	bool const nonideal ) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	if ( mode.compare("pack") == 0 ) {
		mm->set_jump( true ); mm->set_bb( true ); mm->set_chi( true );

	} else if ( mode.compare("full") == 0 ) {
		mm->set_jump( true ); mm->set_bb( true ); mm->set_chi( true );
		if ( nonideal ) {
			mm->set( core::id::PHI, true ); mm->set( core::id::THETA, true ); mm->set( core::id::D, true );
		}

	} else if ( mode.compare("looponly") == 0 ) {
		std::string loopstr( option[ lh::loop_string ]() );
		std::string segstr( option[ lh::seg_string ]() );
		TR << "Movemap setting on loop/segs: " << loopstr << " / " << segstr << std::endl;
		utility::vector1< core::Size > loopres = loopstring_to_loopvector( loopstr );
		utility::vector1< core::Size > segres  = loopstring_to_loopvector( segstr );

		add_movemap_from_loopres( *mm, pose, loopres, nonideal );
		add_movemap_from_loopres( *mm, pose,  segres, nonideal );
	}

	return mm;
}

void
WorkUnit_Sampler::init_from_cmd( const core::Size ){}

void
WorkUnit_Sampler::store_to_decoys( core::io::silent::SilentStructCOP start_struct,
	core::pose::Pose const pose,
	std::string const additional_tag ){
	// exported them into silentstruct
	core::io::silent::SilentStructOP ss =
		core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary");
	ss->fill_struct( pose );
	ss->add_energy( "samplemethod", start_struct->get_energy("samplemethod") );
	ss->set_decoy_tag( start_struct->decoy_tag() + additional_tag );
	decoys().store().push_back( ss );
}

void
WorkUnit_Sampler::store_to_decoys( core::io::silent::SilentStructCOP start_struct,
	core::io::silent::SilentStructOP ss,
	std::string const additional_tag )
{
	ss->add_energy( "samplemethod", start_struct->get_energy("samplemethod") );
	ss->set_decoy_tag( start_struct->decoy_tag() + additional_tag );
	decoys().store().push_back( ss );
}

void
WorkUnit_Sampler::repack( core::pose::Pose &pose,
	core::scoring::ScoreFunctionOP sfxn )
{
	core::kinematics::MoveMapOP mmpack = get_movemap( pose, "pack", false );
	protocols::simple_moves::PackRotamersMoverOP packer = setup_packer( pose, *mmpack, sfxn );
	packer->apply( pose );
}

// This will only work in Cartesian space because
// local loophash doesn't care about ideal geometry
void
WorkUnit_Sampler::ramp_minpack_loop2( core::pose::Pose &pose,
	utility::vector1< core::Size > const loopres,
	core::scoring::ScoreFunctionCOP sfxn,
	bool const nonideal,
	bool const ramp,
	bool const efficient,
	core::Real dist_cut
)
{
	// Expand loopres into its neighbors
	utility::vector1< core::Size > touched_residue = get_touched_res( pose, loopres );

	// Get movemap
	core::kinematics::MoveMap mm;
	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( true );
	bool env_min = true;
	if ( dist_cut < 3.0 ) env_min = false;

	TR.Debug << "loopres/touched_res: " << loopres.size() << " " << touched_residue.size() << std::endl;
	for ( core::Size ires = 1; ires <= loopres.size(); ++ires ) {
		core::Size const resno( loopres[ ires ] );
		mm.set_bb( resno, true );
		mm.set_chi( resno, true );

		if ( nonideal ) {
			for ( core::Size j=1; j<=pose.residue_type(resno).natoms(); ++j ) {
				mm.set( core::id::DOF_ID(core::id::AtomID(j,resno), core::id::THETA ), true );
				mm.set( core::id::DOF_ID(core::id::AtomID(j,resno), core::id::D ), true );
			}
		}
	}

	for ( core::Size ires = 1; ires <= touched_residue.size(); ++ires ) {
		core::Size const resno( touched_residue[ ires ] );
		if ( env_min ) mm.set_bb( resno, true );
		mm.set_chi( resno, true );

		if ( nonideal ) {
			for ( core::Size j=1; j<=pose.residue_type(resno).natoms(); ++j ) {
				mm.set( core::id::DOF_ID(core::id::AtomID(j,resno), core::id::THETA ), true );
				mm.set( core::id::DOF_ID(core::id::AtomID(j,resno), core::id::D ), true );
			}
		}
	}

	protocols::simple_moves::PackRotamersMoverOP packer = setup_packer( pose, mm, sfxn );
	core::optimization::CartesianMinimizer minimizer1;
	core::optimization::AtomTreeMinimizer minimizer2;
	core::optimization::MinimizerOptions minoption( "lbfgs_armijo_nonmonotone",
		0.00001, true, false, false );

	if ( ramp ) {
		// Ramp schedule
		float w_ramp[] = { 0.02, 0.25, 0.55, 1.0 };
		int max_iter[] = { 50, 50, 100, 200 };

		// Pack & Min
		core::scoring::ScoreFunctionOP sfxn_loc = sfxn->clone();
		if ( nonideal && (*sfxn_loc)[ core::scoring::cart_bonded ] < 1.e-3 ) {
			sfxn_loc->set_weight( core::scoring::cart_bonded, 0.8 );
		}

		for ( int i = 0; i< 4; ++i ) {
			packer->apply( pose );
			minoption.max_iter( (core::Size)(max_iter[i]) );
			sfxn_loc->set_weight( core::scoring::fa_rep,
				(*sfxn)[ core::scoring::fa_rep ] * (core::Real)(w_ramp[i]) );

			if ( nonideal ) {
				minimizer1.run( pose, mm, *sfxn_loc, minoption );
			} else {
				minimizer2.run( pose, mm, *sfxn_loc, minoption );
			}

		}
	} else {
		packer->apply( pose );
		if ( efficient ) {
			minoption.max_iter( 30 );
		} else {
			minoption.max_iter( 200 );
		}

		if ( nonideal ) {
			minimizer1.run( pose, mm, *sfxn, minoption );
		} else {
			minimizer2.run( pose, mm, *sfxn, minoption );
		}
	}
}

void
WorkUnit_Sampler::superimpose_to_ref( core::pose::Pose const &pose_ref,
	core::pose::Pose &pose_work,
	utility::vector1< core::Size > exclude_res ) const
{

	std::map< core::id::AtomID, core::id::AtomID > atom_map;

	for ( core::Size ires = 1; ires <= pose_ref.total_residue(); ++ires ) {
		if ( !pose_ref.residue(ires).has( " CA " ) || exclude_res.contains( ires ) ) continue;
		core::Size iatm = pose_ref.residue(ires).atom_index(" CA ");
		core::id::AtomID atomid1( iatm, ires );
		core::id::AtomID atomid2( iatm, ires );
		atom_map[ atomid1 ] = atomid2;
	}

	core::scoring::superimpose_pose( pose_work, pose_ref, atom_map );
}

core::scoring::ScoreFunctionOP
WorkUnit_Sampler::get_energy( std::string const sfxn_name,
	bool const softpack,
	core::Real const weight_coord_cst ) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( sfxn_name.compare("scorefacts_cart") == 0 ) {
		// hacky way for LJ Hbond...
		option[ corrections::score::lj_hbond_hdis ].value( 2.4 );
		option[ corrections::score::lj_hbond_OH_donor_dis ].value( 3.2 );
	}

	core::scoring::ScoreFunctionOP sfxn =
		( sfxn_name.compare("cen_cart") == 0 )?
		core::scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart" )
		: ( sfxn_name.compare("cen_loop") == 0 )?
		core::scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart" )
		: ( sfxn_name.compare("talaris2013_cart") == 0 )?
		core::scoring::ScoreFunctionFactory::create_score_function( "talaris2013" )
		: ( sfxn_name.compare("talaris2013_cart_softrep") == 0 )?
		core::scoring::ScoreFunctionFactory::create_score_function( "soft_rep" )
		: ( sfxn_name.compare("scorefacts_cart") == 0 )?
		core::scoring::ScoreFunctionFactory::create_score_function( "scorefacts_cart" )
		:
		core::scoring::ScoreFunctionFactory::create_score_function( "talaris2013" );

	// is this important?
	if ( sfxn_name.compare("cen_loop") == 0 ) {
		sfxn->set_weight( core::scoring::cbeta_smooth, 3.0 );
		sfxn->set_weight( core::scoring::cen_env_smooth, 3.0 );
	}

	if ( sfxn_name.compare( "talaris2013_cart" ) == 0 ||
			sfxn_name.compare( "talaris2013_cart_softrep" ) == 0 ) {
		sfxn->set_weight( core::scoring::cart_bonded, 0.5 );
	}
	sfxn->set_weight( core::scoring::pro_close, 0.0 );

	if ( softpack ) {
		sfxn->set_weight( core::scoring::fa_rep, 0.004 );
		sfxn->set_weight( core::scoring::fa_atr, 0.0  );
	}

	if ( weight_coord_cst > 0.0 ) {
		sfxn->set_weight( core::scoring::coordinate_constraint, 1.0 );
	}

	return sfxn;
}

void
WorkUnit_Sampler::revert_facts_params() const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// hacky way for LJ Hbond...
	// any better way?
	option[ corrections::score::lj_hbond_hdis ].value( 1.75 );
	option[ corrections::score::lj_hbond_OH_donor_dis ].value( 2.6 );
}

}
}

