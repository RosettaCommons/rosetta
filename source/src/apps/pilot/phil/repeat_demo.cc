// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/phil/repeat_demo.cc
/// @brief  Example of how to set up a symmetric repeat protein system where you fold through the polypeptide chain
/// @author Phil Bradley (pbradley@fredhutch.org)

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/variant_util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/relax/FastRelax.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

#include <ObjexxFCL/format.hh>

static basic::Tracer TR( "apps.pilot.phil.repeat_demo" );

// Lazy:
using namespace basic::options;
using namespace core;
using namespace pose;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using ObjexxFCL::format::F;

using namespace std;

// OPTIONS
namespace my_options {
IntegerOptionKey nrepeat("my:nrepeat");
}

void
add_my_options()
{
	option.add( my_options::nrepeat, "nrepeat" );
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// fills in a SymmetryInfo object with the necessary info
///
void
setup_repeat_symminfo(
	Size const repeatlen,
	Size const nrepeat,
	Size const base_repeat,
	conformation::symmetry::SymmetryInfo & symminfo
)
{

	Size const base_offset( (base_repeat-1)*repeatlen );

	for ( Size i=1; i<= nrepeat; ++i ) {
		if ( i == base_repeat ) continue;
		Size const offset( (i-1)*repeatlen );
		for ( Size j=1; j<= repeatlen; ++j ) {
			symminfo.add_bb_clone ( base_offset+j, offset+j );
			symminfo.add_chi_clone( base_offset+j, offset+j );
		}
	}

	symminfo.num_virtuals( 1 ); // the one at the end...
	symminfo.set_use_symmetry( true );
	symminfo.set_flat_score_multiply( nrepeat*repeatlen+1, 1 );
	symminfo.torsion_changes_move_other_monomers( true ); // note -- this signals that we are folding between repeats

	/// repeats prior to base_repeat have score_multiply ==> 0
	for ( Size i=1; i<base_repeat; ++i ) {
		Size const offset( (i-1)*repeatlen );
		for ( Size j=1; j<= repeatlen; ++j ) symminfo.set_score_multiply( offset+j, 0 );
	}

	symminfo.update_score_multiply_factor();

}

///////////////////////////////////////////////////////////////////////////////
/// sets up a repeat pose, starting from a non-symmetric pdb with nres=repeatlen*nrepeat
///

void
setup_repeat_pose(
	Size const repeatlen,
	Size const nrepeat,
	Size const base_repeat, // repeat number of the independent repeat (aka "monomer")
	Pose & pose
)
{
	runtime_assert( !pose::symmetry::is_symmetric( pose ) ); // not to begin with...
	runtime_assert( nrepeat * repeatlen == pose.size() );

	runtime_assert( base_repeat > 1 ); // why? well, the base repeat should get the right context info from nbring repeats
	// but note that with base_repeat>1 we probably can't use linmem_ig and there are other places in the code that
	// assume that monomer 1 is the independent monomer. These should gradually be fixed. Let me (PB) know if you run into
	// trouble.

	Size const nres_protein( nrepeat * repeatlen );
	remove_upper_terminus_type_from_pose_residue( pose, pose.size() );
	ResidueOP vrtrsd
		( conformation::ResidueFactory::create_residue( *core::pose::get_restype_for_pose( pose, "VRTBB" ) ) );
	pose.append_residue_by_bond( *vrtrsd, true ); // since is polymer...
	pose.conformation().insert_chain_ending( nres_protein );

	kinematics::FoldTree f( pose.size() );
	f.reorder( pose.size() );
	pose.fold_tree( f );


	conformation::symmetry::SymmetryInfo symminfo;

	setup_repeat_symminfo( repeatlen, nrepeat, base_repeat, symminfo );

	/// now make symmetric
	pose::symmetry::make_symmetric_pose( pose, symminfo );


	runtime_assert( pose::symmetry::is_symmetric( pose ) );
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
rebuild_test()
{
	// Unfortunately linmem_ig assumes the independent monomer is the first one (I think)
	runtime_assert( !option[ OptionKeys::packing::linmem_ig ].user() );

	// this is a bummer: there's some dumb code that doesn't understand that dofs associated with an atom for one residue
	// can correspond to torsion angles for another. I can fix this if it's a big problem (PB).
	runtime_assert( option[ OptionKeys::optimization::old_sym_min ] );

	// this doofy command line flag is the super-secret signal to rotamer building and scorefxn generation that
	// we are using symmetry. This should probably be changed. In the meantime you can always add
	// "-symmetry_definition stoopid" to the command line. If you are mixing symmetry and normal stuff you can
	// fiddle with this global variable. I know this is the most horrible hack in the world, but for once I can
	// claim innocence.... although I should fix it rather than complaining!
	runtime_assert( option[ OptionKeys::symmetry::symmetry_definition ].user() );

	string const filename( start_file() );
	Pose pose;

	import_pose::pose_from_file( pose, filename , core::import_pose::PDB_file);

	pose.dump_pdb( "start.pdb" );
	Pose const pdb_pose( pose );


	Size const nrepeat( option[ my_options::nrepeat ] );

	Size const base_repeat( 3 );
	runtime_assert( base_repeat < nrepeat );

	runtime_assert( pose.size()%nrepeat == 0 );
	Size const repeatlen( pose.size() / nrepeat );

	// this is the first helper function:
	setup_repeat_pose( repeatlen, nrepeat, base_repeat, pose );

	Size const base_offset( (base_repeat-1)*repeatlen );
	TR.Trace << "symmetrize backbone pose" << endl;

	for ( Size i=1; i<= repeatlen; ++i ) {
		Size const pos( i+base_offset );
		pose.set_phi  ( pos, pose.phi  (pos) );
		pose.set_psi  ( pos, pose.psi  (pos) );
		pose.set_omega( pos, pose.omega(pos) );
	}
	pose.dump_pdb( "after_bb_symm.pdb" );

	TR.Trace << "symmetrize sidechains pose" << endl;

	for ( Size i=1; i<= repeatlen; ++i ) {
		Size const pos( i+base_offset );
		ResidueOP oldrsd( pose.residue(pos).clone() );
		pose.replace_residue( pos, *oldrsd, false );
	}

	pose.dump_pdb( "after_chi_symm.pdb" );

	// compute rmsd to pdb pose
	Real const rmsd( CA_rmsd( pose, pdb_pose, 1, nrepeat * repeatlen ) );
	TR.Trace << "rmsd to pdb_pose: " << F(9,3,rmsd)<< endl;

	{ // try minimizing, do deriv_check
		scoring::symmetry::SymmetricScoreFunctionOP fa_scorefxn( new scoring::symmetry::SymmetricScoreFunction() );
		fa_scorefxn->set_weight( fa_atr, 1.0 );
		fa_scorefxn->set_weight( fa_rep, 1.0 );
		//fa_scorefxn->set_weight( fa_elec, 1.0 );

		kinematics::MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
		movemap->set_bb ( true );
		movemap->set_chi( true );
		movemap->set_bb ( pose.size(), false );
		movemap->set_chi( pose.size(), false );

		bool const use_nblist( true ), deriv_check( true ), deriv_check_verbose( false );
		protocols::simple_moves::symmetry::SymMinMoverOP min_mover
			( new protocols::simple_moves::symmetry::SymMinMover(movemap, fa_scorefxn, "dfpmin_armijo_nonmonotone", 1e-2,
			use_nblist, deriv_check, deriv_check_verbose ) );

		min_mover->apply( pose );

		pose.dump_pdb("after_min.pdb");
	}

	// now try a fastrelax
	ScoreFunctionOP fa_scorefxn( get_score_function() );
	runtime_assert( pose::symmetry::is_symmetric( pose ) );
	runtime_assert( pose::symmetry::is_symmetric( *fa_scorefxn ) );
	{
		protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );

		kinematics::MoveMapOP movemap = kinematics::MoveMapOP( new kinematics::MoveMap );
		movemap->set_bb ( true );
		movemap->set_chi( true );
		movemap->set_bb ( pose.size(), false );
		movemap->set_chi( pose.size(), false );
		fastrelax.set_movemap( movemap );
		TR.Trace << "start relaxing" << endl;
		fastrelax.apply( pose );
		TR.Trace << "finished relaxing" << endl;
		Real const rmsd( CA_rmsd( pose, pdb_pose, 1, nrepeat * repeatlen ) );
		TR.Trace << "rmsd to pdb_pose after relax: " << F(9,3,rmsd)<< endl;
		pose.dump_pdb( "after_relax.pdb" );
	}

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	using namespace basic::options;

	try {
		add_my_options();

		devel::init(argc, argv);

		rebuild_test();

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
