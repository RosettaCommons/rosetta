// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief protocols that are specific to relax
/// @details
/// @author Mike Tyka, Monica Berrondo

//#include <protocols/jobdist/Jobs.hh>

#include <protocols/relax/cst_util.hh>
#include <protocols/relax/RelaxProtocolBase.hh>

#include <protocols/relax/AtomCoordinateCstMover.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/rms_util.hh>
//#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/pose/util.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//*only for debug structures
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

#include <basic/Tracer.hh>

#include <core/conformation/Residue.hh>
//#include <protocols/jobdist/Jobs.hh>
//#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.relax.ClassicRelax" );

using namespace core;
using io::pdb::dump_pdb;
////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace relax {


RelaxProtocolBase::RelaxProtocolBase( core::scoring::ScoreFunctionOP score_in ) :
	parent( "RelaxProtocol" ),
	min_type_("dfpmin_armijo_nonmonotone"),
	scorefxn_( score_in ),
	task_factory_(/* NULL */)
{
	set_defaults();
}

RelaxProtocolBase::RelaxProtocolBase( std::string const & movername ) :
	parent( movername ),
	min_type_("dfpmin_armijo_nonmonotone"),
	scorefxn_( /* 0 */ ),
	task_factory_(/* NULL */)
{
	set_defaults();
}

RelaxProtocolBase::RelaxProtocolBase( std::string const & movername, core::scoring::ScoreFunctionOP score_in ) :
	parent( movername ),
	min_type_("dfpmin_armijo_nonmonotone"),
	scorefxn_(score_in ),
	task_factory_(/* NULL */)
{
	set_defaults();
}


RelaxProtocolBase::~RelaxProtocolBase() {}

void RelaxProtocolBase::initialize_movemap(
	core::pose::Pose const & pose,
	core::kinematics::MoveMap & movemap
)
{
	using namespace core::id;

	if ( fix_omega_ ) {
		for ( Size i=1; i<=pose.total_residue(); ++i ) {
			movemap.set( TorsionID(i, BB, 3),false );
		}
	}

	if ( minimize_bond_lengths_ ) {
		// 0 Default  all bondlengths
		// 1          backbone only
		// 2          sidechain only
		// 3          CA only (Ca-C,Ca-N and Ca-Cb)

		if ( minimize_bondlength_subset_ == 0 ) {
			movemap.set( core::id::D, true );
		} else if ( minimize_bondlength_subset_ == 1 ) {
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::chemical::AtomIndices const & ii_mainchain_atoms( pose.residue(ii).mainchain_atoms() );
				for ( Size jj = 1; jj <= ii_mainchain_atoms.size(); ++jj ) {
					//if ( jj == 1 ) {
					// if ( ii > 1 && pose.residue(ii).is_bonded( ii-1 ) && !pose.residue(ii).has_variant_type("CUTPOINT_UPPER")) {
					//  movemap.set( DOF_ID( AtomID( ii_mainchain_atoms[ jj ], ii ), core::id::D ), true );
					// }
					//} else {
					movemap.set( DOF_ID( AtomID( ii_mainchain_atoms[ jj ], ii ), core::id::D ), true );
					//}
				}
			}
		} else if ( minimize_bondlength_subset_ == 2 ) {
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::conformation::Residue const &res_i = pose.residue(ii);
				for ( Size jj = 1; jj <= res_i.natoms(); ++jj ) {
					if ( res_i.atom_is_backbone(jj) ) continue;
					movemap.set( DOF_ID( AtomID( jj, ii ), core::id::D ), true );
				}
			}
		} else if ( minimize_bondlength_subset_ == 3 ) {
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::conformation::Residue const &res_i = pose.residue(ii);
				if ( res_i.type().has( " C  ") ) {
					movemap.set( DOF_ID( AtomID( res_i.atom_index(" C  "), ii ), core::id::D ), true );
				}
				if ( res_i.type().has( " CA ") ) {
					movemap.set( DOF_ID( AtomID( res_i.atom_index(" CA "), ii ), core::id::D ), true );
				}
				if ( res_i.type().has( " CB ") ) {
					movemap.set( DOF_ID( AtomID( res_i.atom_index(" CB "), ii ), core::id::D ), true );
				}
			}
		}
	}

	if ( minimize_bond_angles_ ) {


		// 0 Default  all bondangles
		// 1          backbone only
		// 2          sidechain only
		// 3          tau only
		// 4          Ca-Cb only
		if ( minimize_bondangle_subset_ == 0 ) {
			movemap.set( core::id::THETA, true );
		} else if ( minimize_bondangle_subset_ == 1 ) {
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::chemical::AtomIndices const & ii_mainchain_atoms( pose.residue(ii).mainchain_atoms() );
				for ( Size jj = 1; jj <= ii_mainchain_atoms.size(); ++jj ) {
					//if ( jj == 1 || jj == 2 ) {  //fpd  add jj==2
					// if ( ii > 1 && pose.residue(ii).is_bonded( ii-1 ) && !pose.residue(ii).has_variant_type("CUTPOINT_UPPER")) {
					//  movemap.set( DOF_ID( AtomID( ii_mainchain_atoms[ jj ], ii ), core::id::THETA ), true );
					// }
					//} else {
					movemap.set( DOF_ID( AtomID( ii_mainchain_atoms[ jj ], ii ), core::id::THETA ), true );
					//}
				}
			}
		} else if ( minimize_bondangle_subset_ == 2 ) {
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::conformation::Residue const &res_i = pose.residue(ii);
				for ( Size jj = 1; jj <= res_i.natoms(); ++jj ) {
					if ( res_i.atom_is_backbone(jj) ) continue;
					movemap.set( DOF_ID( AtomID( jj, ii ), core::id::THETA ), true );
				}
			}
		} else if ( minimize_bondangle_subset_ == 3 ) {
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::conformation::Residue const &res_i = pose.residue(ii);
				if ( res_i.type().has( " C  ") ) {
					movemap.set( DOF_ID( AtomID( res_i.atom_index(" C  "), ii ), core::id::THETA ), true );
				}
			}
		} else if ( minimize_bondangle_subset_ == 4 ) {
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::conformation::Residue const &res_i = pose.residue(ii);
				if ( res_i.type().has( " CB ") ) {
					movemap.set( DOF_ID( AtomID( res_i.atom_index(" CB "), ii ), core::id::THETA ), true );
				}
			}
		}
	}

}

void RelaxProtocolBase::set_defaults(){
	set_default_minimization_settings();
	set_default_coordinate_settings();
	set_default_movemap();

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	dry_run_ = basic::options::option[ OptionKeys::run::dry_run ]();
}

void RelaxProtocolBase::set_default_minimization_settings(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	fix_omega_ = option[ OptionKeys::relax::fix_omega]();
	minimize_bond_lengths_ = option[ OptionKeys::relax::minimize_bond_lengths ]();
	minimize_bond_angles_ = option[ OptionKeys::relax::minimize_bond_angles ]();
	minimize_bondlength_subset_ = option[ OptionKeys::relax::minimize_bondlength_subset ]();
	minimize_bondangle_subset_  = option[ OptionKeys::relax::minimize_bondangle_subset ]();

	//fpd extras
	cartesian_ = option[ OptionKeys::relax::cartesian ]();
	if ( option[ OptionKeys::relax::min_type ].user() ) {
		min_type_ = option[ OptionKeys::relax::min_type ]();
	} else if ( cartesian_ || minimize_bond_lengths_ || minimize_bond_angles_ ) {
		min_type_ = "lbfgs_armijo_nonmonotone";  // default is different for cartesian/nonideal minimization
	}

	//fpd
	max_iter_ = 0;  // if nonzero, override MinimizerOption default
}

void RelaxProtocolBase::set_default_coordinate_settings(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	constrain_relax_to_native_coords_ = option[ OptionKeys::relax::constrain_relax_to_native_coords ]();
	constrain_relax_to_start_coords_ = option[ OptionKeys::relax::constrain_relax_to_start_coords ]();
	constrain_coords_ = constrain_relax_to_native_coords_ || constrain_relax_to_start_coords_;
	coord_constrain_sidechains_ = option[ OptionKeys::relax::coord_constrain_sidechains ]();
	if ( coord_constrain_sidechains_ && ! constrain_coords_ ) {
		utility_exit_with_message("Option -relax:coord_constrain_sidechains also requires either -relax:constrain_relax_to_start_coords or -relax:constrain_relax_to_native_coords");
	}
	explicit_ramp_constraints_ = option[ OptionKeys::relax::ramp_constraints ].user();
	ramp_down_constraints_ = option[ OptionKeys::relax::ramp_constraints ]();
	constrain_relax_segments_ = option[ OptionKeys::relax::constrain_relax_segments ].user();
	limit_aroma_chi2_ = option[ OptionKeys::relax::limit_aroma_chi2 ]();
}

void RelaxProtocolBase::set_default_movemap(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	movemap_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap() );

	if ( option[ OptionKeys::in::file::movemap ].user() ) {

		//Allow user settings to be applied before movemap is read in. But, by default, use the movemap file.
		if ( option[ OptionKeys::relax::jump_move ].user() ) {
			movemap_->set_jump( option[ OptionKeys::relax::jump_move ]() );
		}
		if ( option[ OptionKeys::relax::bb_move ].user() ) {
			movemap_->set_bb( option[ OptionKeys::relax::bb_move ]() );
		}
		if ( option[ OptionKeys::relax::chi_move ].user() ) {
			movemap_->set_chi( option[ OptionKeys::relax::chi_move ]() );
		}

		movemap_->init_from_file(option[ OptionKeys::in::file::movemap ]() );

	} else {
		movemap_->set_jump( option[ OptionKeys::relax::jump_move ]() );
		movemap_->set_bb( option[ OptionKeys::relax::bb_move ]() );
		movemap_->set_chi( option[ OptionKeys::relax::chi_move ]() );
	}
}

void RelaxProtocolBase::register_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static bool runonce( false );
	if ( runonce ) return;

	runonce = true;

	option.add_relevant( OptionKeys::relax::fix_omega );
	option.add_relevant( OptionKeys::relax::minimize_bond_lengths );
	option.add_relevant( OptionKeys::relax::minimize_bond_angles );
	option.add_relevant( OptionKeys::relax::minimize_bondlength_subset );
	option.add_relevant( OptionKeys::relax::minimize_bondlength_subset );

}


//mjo TODO: please move this with the rest of the disulfide code

////////////////////////////////////////////////////////////////////////////////////////////////////////
// apply_disulfides
//  this is really just "repack disulfides", where a repack is called using only the FA disulfide score
//
//  Anyway, the function as is doesn't do much of anything right now that the subsequent repack steps
//  wouldn't also do. It's in here more as a placeholder, as I'm planning to rewrite it into a minimize
//  with constraints type disulfide minimizer
//
// -rvernon
////////////////////////////////////////////////////////////////////////////////////////////////////////
void RelaxProtocolBase::apply_disulfides( core::pose::Pose & pose ){
	using namespace moves;
	using namespace scoring;
	using namespace core::pose::datacache;


	if ( ( basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].user() ) ||
			( basic::options::option[ basic::options::OptionKeys::in::detect_disulf ].user() ) ) {

		std::string weight_set("score12_justdisulfides");
		core::scoring::ScoreFunctionOP disulf_score_only(scoring::ScoreFunctionFactory::create_score_function( weight_set));

		protocols::simple_moves::PackRotamersMoverOP full_repack;
		core::pack::task::PackerTaskOP task;
		task = pack::task::TaskFactory::create_packer_task( pose );

		utility::vector1<bool> allow_repack(pose.total_residue(), false);

		for ( Size i = 1; i<= pose.total_residue() ; ++i ) {
			allow_repack[i] = movemap_->get_chi(i);
		}

		if ( basic::options::option[ basic::options::OptionKeys::relax::chi_move].user() ) {
			bool const repack = basic::options::option[ basic::options::OptionKeys::relax::chi_move]();
			allow_repack.assign( pose.total_residue(), repack);
		}

		task->initialize_from_command_line().restrict_to_repacking().restrict_to_residues(allow_repack);
		task->or_include_current( true );
		if ( basic::options::option[ basic::options::OptionKeys::symmetry::symmetry_definition ].user() )  {
			full_repack = protocols::simple_moves::PackRotamersMoverOP( new simple_moves::symmetry::SymPackRotamersMover( disulf_score_only, task ) );
		} else {
			full_repack = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover( disulf_score_only, task ) );
		}

		( *disulf_score_only )( pose );

		full_repack->apply(pose);
	}
}


void RelaxProtocolBase::set_up_constraints( core::pose::Pose &pose, core::kinematics::MoveMap& local_movemap ) {
	using namespace conformation;
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose::datacache;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace protocols::moves;
	using namespace core::scoring;

	if ( constrain_coords_ ) {
		//fpd  Make ramping on by default if one of -constrain_relax_* is specified
		//fpd  Let it be overridden if '-ramp_constraints false' is specified
		if ( !option[ OptionKeys::relax::ramp_constraints ].user() ) {
			ramp_down_constraints_ = true;
		}

		protocols::relax::AtomCoordinateCstMover coord_cst_mover;
		if ( constrain_relax_to_native_coords_ ) {
			if ( get_native_pose() ) {
				coord_cst_mover.set_refstruct( get_native_pose() );
			} else {
				utility_exit_with_message("Native pose needed for OptionKeys::relax::constrain_relax_to_native_coords");
			}
		}

		if ( constrain_relax_segments_ ) {
			coord_cst_mover.set_loop_segments( protocols::loops::LoopsCOP( protocols::loops::LoopsOP( new protocols::loops::Loops(  option[ OptionKeys::relax::constrain_relax_segments ]() ) ) ) );
		}

		//if -relax::coord_cst_width is given, the code instead uses _bounded_ constraints on
		//heavy atom positions, of the specified width
		coord_cst_mover.cst_sidechain( coord_constrain_sidechains_ );
		coord_cst_mover.cst_sd( option[ OptionKeys::relax::coord_cst_stdev ] );
		coord_cst_mover.bounded( option[ OptionKeys::relax::coord_cst_width ].user() );
		coord_cst_mover.cst_width( option[ OptionKeys::relax::coord_cst_width ]() );
		coord_cst_mover.ambiguous_hnq( option[ OptionKeys::packing::flip_HNQ ]() );

		// Add virtual root if one doesn't exist
		if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
			core::pose::addVirtualResAsRoot(pose);
		}

		// Actually apply the constraints
		coord_cst_mover.apply( pose );

		if ( get_scorefxn()->get_weight( coordinate_constraint ) == 0 ) {
			get_scorefxn()->set_weight( coordinate_constraint, 0.5 );
		}

		// Should this be only if we've added the virtual root?
		local_movemap.set_jump( pose.num_jump(), true );
	} // if constrain_coords_

	// Support for RosettaScripts
	if ( cst_files_.size() > 0 ) {
		// To preserve? let's just turn off
		//core::scoring::constraints::ConstraintSetOP
		// save_pose_constraint_set = pose.constraint_set()->clone();

		for ( Size i_cst = 1; i_cst <= cst_files_.size(); ++i_cst ) {
			std::string const filename = cst_files( i_cst );
			ConstraintSetOP user_csts
				= ConstraintIO::get_instance()->read_constraints_new( filename,
				ConstraintSetOP( new ConstraintSet ), pose );
			pose.constraint_set( user_csts );
		}
	} // if constrain_user_defined_

	if ( option[ OptionKeys::relax::sc_cst_maxdist ].user() && option[ OptionKeys::relax::sc_cst_maxdist ]() > 0 ) {
		// derive a set of side-chain restraints
		Real const upper_dist_cutoff( option[ OptionKeys::relax::sc_cst_maxdist ]() );
		derive_sc_sc_restraints( pose, upper_dist_cutoff );

		if ( get_scorefxn()->get_weight( atom_pair_constraint ) == 0 ) {
			get_scorefxn()->set_weight( atom_pair_constraint, 2.0 );
		}
	}

} // setup_up_constraints


void RelaxProtocolBase::output_debug_structure( core::pose::Pose & pose, std::string prefix ) {
	using namespace core::io::silent;
	using namespace basic::options;
	if ( option[ basic::options::OptionKeys::out::file::silent ].user() ) {
		std::string silent_file="_"+prefix;
		std::string output_file_name= protocols::jd2::current_output_name();
		if ( output_file_name != "" ) {
			silent_file = output_file_name + silent_file;
		} else silent_file = "bla"+silent_file;

		SilentFileData sfd;
		//filename might have been changed -- e.g., to also have an MPI rank in there

		//  ProteinSilentStruct pss;
		io::silent::SilentStructOP pss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
		pss->fill_struct( pose, get_current_tag() );

		sfd.write_silent_struct( *pss, silent_file, false /* bWriteScoresOnly */ );
	} // if option[ out::file::silent ].user()
}

} // namespace relax
} // namespace protocols
