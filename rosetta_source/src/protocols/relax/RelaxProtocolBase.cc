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
/// @detailed
/// @author Mike Tyka, Monica Berrondo

//#include <protocols/jobdist/Jobs.hh>

#include <protocols/relax/cst_util.hh>
#include <protocols/relax/RelaxProtocolBase.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/io/pdb/pose_io.hh>

#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/rms_util.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>
#include <protocols/moves/PackRotamersMover.hh>
#include <protocols/moves/symmetry/SymPackRotamersMover.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

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

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

#include <basic/Tracer.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <utility/io/mpistream.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.relax.ClassicRelax");

using namespace core;
using io::pdb::dump_pdb;
////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace relax {



RelaxProtocolBase::RelaxProtocolBase( core::scoring::ScoreFunctionOP score_in ) :
	parent( "RelaxProtocol" ),
	min_type_("dfpmin_armijo_nonmonotone"),
	scorefxn_( score_in )
{
	set_defaults();
}

RelaxProtocolBase::RelaxProtocolBase( std::string const & movername ) :
	parent( movername ),
	min_type_("dfpmin_armijo_nonmonotone"),
	scorefxn_( 0 )
{
	set_defaults();
}

RelaxProtocolBase::RelaxProtocolBase( std::string const & movername, core::scoring::ScoreFunctionOP score_in ) :
	parent( movername ),
	min_type_("dfpmin_armijo_nonmonotone"),
	scorefxn_(score_in )
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
	if ( minimize_bond_lengths_ ) {
		movemap.set( core::id::D, true );
	} else if ( minimize_mainchain_bond_lengths_ ) {
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			core::chemical::AtomIndices const & ii_mainchain_atoms( pose.residue(ii).mainchain_atoms() );
			for ( Size jj = 1; jj <= ii_mainchain_atoms.size(); ++jj ) {
				if ( jj == 1 ) {
					if ( ii > 1 && pose.residue(ii).is_bonded( ii-1 ) && !pose.residue(ii).has_variant_type("CUTPOINT_UPPER")) {
						movemap.set( DOF_ID( AtomID( ii_mainchain_atoms[ jj ], ii ), core::id::D ), true );
					}
				} else {
					movemap.set( DOF_ID( AtomID( ii_mainchain_atoms[ jj ], ii ), core::id::D ), true );
				}
			}
		}
	}

	if ( minimize_bond_angles_ ) {
		movemap.set( core::id::THETA, true );
	} else if ( minimize_mainchain_bond_angles_ ) {
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			core::chemical::AtomIndices const & ii_mainchain_atoms( pose.residue(ii).mainchain_atoms() );
			for ( Size jj = 1; jj <= ii_mainchain_atoms.size(); ++jj ) {
				if ( jj == 1 ) {
					if ( ii > 1 && pose.residue(ii).is_bonded( ii-1 ) && !pose.residue(ii).has_variant_type("CUTPOINT_UPPER")) {
						movemap.set( DOF_ID( AtomID( ii_mainchain_atoms[ jj ], ii ), core::id::THETA ), true );
					}
				} else {
					movemap.set( DOF_ID( AtomID( ii_mainchain_atoms[ jj ], ii ), core::id::THETA ), true );
				}
			}
		}
	}

}

void RelaxProtocolBase::set_scorefxn( core::scoring::ScoreFunctionOP score ) {
	scorefxn_ = score;
}

core::scoring::ScoreFunctionOP RelaxProtocolBase::get_scorefxn() {
	return scorefxn_;
}

const core::scoring::ScoreFunctionCOP RelaxProtocolBase::get_scorefxn() const {
	return scorefxn_;
}

void RelaxProtocolBase::set_defaults(){
	set_default_minimization_settings();
	set_default_coordinate_settings();
	set_default_movemap();
	task_factory_ = 0;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	dry_run_ = basic::options::option[ OptionKeys::run::dry_run ]();
}

void RelaxProtocolBase::set_default_minimization_settings(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	minimize_bond_lengths_ = option[ OptionKeys::relax::minimize_bond_lengths ]();
	minimize_bond_angles_ = option[ OptionKeys::relax::minimize_bond_angles ]();
	minimize_mainchain_bond_lengths_ = option[ OptionKeys::relax::minimize_mainchain_bond_lengths ]();
	minimize_mainchain_bond_angles_  = option[ OptionKeys::relax::minimize_mainchain_bond_angles ]();

	//fpd extras
	cartesian_ = option[ OptionKeys::relax::cartesian ]();
	if ( option[ OptionKeys::relax::min_type ].user() )
		min_type_ = option[ OptionKeys::relax::min_type ]();
	else if (cartesian_)
		min_type_ = "lbfgs_armijo_nonmonotone";  // default is different for cartesian minimization
}

void RelaxProtocolBase::set_default_coordinate_settings(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	constrain_relax_to_native_coords_ = option[ OptionKeys::relax::constrain_relax_to_native_coords ]();
	constrain_relax_to_start_coords_ = option[ OptionKeys::relax::constrain_relax_to_start_coords ]();
	constrain_coords_ = constrain_relax_to_native_coords_ || constrain_relax_to_start_coords_;
	ramp_down_constraints_ = option[ OptionKeys::relax::ramp_constraints ]();
	constrain_relax_segments_ = option[ OptionKeys::relax::constrain_relax_segments ].user();
	limit_aroma_chi2_ = option[ OptionKeys::relax::limit_aroma_chi2 ]();
}

void RelaxProtocolBase::set_default_movemap(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	movemap_ = new core::kinematics::MoveMap();
	movemap_->set_jump( option[ OptionKeys::relax::jump_move ]() );
	movemap_->set_bb( option[ OptionKeys::relax::bb_move ]() );
	movemap_->set_chi( option[ OptionKeys::relax::chi_move ]() );
}

void RelaxProtocolBase::set_movemap( core::kinematics::MoveMapOP movemap ) {
	movemap_ = movemap;
}

void RelaxProtocolBase::set_min_type( std::string min_type ) {
	min_type_ = min_type;
}

void RelaxProtocolBase::set_task_factory( core::pack::task::TaskFactoryOP taskf ) {
	task_factory_ = taskf;
}

void RelaxProtocolBase::minimize_bond_lengths( bool setting ) { minimize_bond_lengths_ = setting; }
void RelaxProtocolBase::minimize_bond_angles( bool setting ) { minimize_bond_angles_ = setting; }
void RelaxProtocolBase::minimize_mainchain_bond_lengths( bool setting ) { minimize_mainchain_bond_lengths_ = setting; }
void RelaxProtocolBase::minimize_mainchain_bond_angles( bool setting ) { minimize_mainchain_bond_angles_ = setting; }

bool RelaxProtocolBase::minimize_bond_lengths() const { return minimize_bond_lengths_;}
bool RelaxProtocolBase::minimize_bond_angles() const { return minimize_bond_angles_;}
bool RelaxProtocolBase::minimize_mainchain_bond_lengths() const { return minimize_mainchain_bond_lengths_;}
bool RelaxProtocolBase::minimize_mainchain_bond_angles() const { return minimize_mainchain_bond_angles_;}

void RelaxProtocolBase::register_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static bool runonce( false );
	if ( runonce ) return;

	runonce = true;

	option.add_relevant( OptionKeys::relax::minimize_bond_lengths );
	option.add_relevant( OptionKeys::relax::minimize_bond_angles );
	option.add_relevant( OptionKeys::relax::minimize_mainchain_bond_lengths );
	option.add_relevant( OptionKeys::relax::minimize_mainchain_bond_angles );

}

core::kinematics::MoveMapOP
RelaxProtocolBase::get_movemap()
{
	return movemap_;
}

core::pack::task::TaskFactoryOP const &
RelaxProtocolBase::get_task_factory() const
{
	return task_factory_;
}

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

		moves::PackRotamersMoverOP full_repack;
		core::pack::task::PackerTaskOP task;
		task = pack::task::TaskFactory::create_packer_task( pose );

		utility::vector1<bool> allow_repack(pose.total_residue(), false);

		for ( Size i = 1; i<= pose.total_residue() ; ++i ) {
			allow_repack[i] = movemap_->get_chi(i);
		}

		if (basic::options::option[ basic::options::OptionKeys::relax::chi_move].user() ){
			bool const repack = basic::options::option[ basic::options::OptionKeys::relax::chi_move]();
			allow_repack.assign( pose.total_residue(), repack);
		}

		task->initialize_from_command_line().restrict_to_repacking().restrict_to_residues(allow_repack);
		task->or_include_current( true );
		if ( basic::options::option[ basic::options::OptionKeys::symmetry::symmetry_definition ].user() )  {
			full_repack = new moves::symmetry::SymPackRotamersMover( disulf_score_only, task );
		} else {
			full_repack = new moves::PackRotamersMover( disulf_score_only, task );
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
	using namespace core::id;
	using namespace protocols::moves;
	using namespace core::scoring;

	core::pose::Pose constraint_target_pose = pose;

	//fpd  Make ramping on by default if one of -constrain_relax_* is specified
	//fpd  Let it be overridden if '-ramp_constraints false' is specified
	if (constrain_coords_ && !option[ OptionKeys::relax::ramp_constraints ].user() )
		ramp_down_constraints_ = true;

	if( constrain_relax_to_native_coords_ ){
		if ( get_native_pose() ) {
			constraint_target_pose = *get_native_pose();
			if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
				/// Align the native pose to the input pose to avoid rotation/translation based
				///  errors.
				///fpd  (Only if not rooted on a VRT to avoid problems with density/symmetry)
				core::scoring::calpha_superimpose_pose( constraint_target_pose, pose );
			}
		} else {
			std::cerr << "Native pose needed for OptionKeys::relax::constrain_relax_to_native_coords" << std::endl;
			utility_exit();
		}
	}

	protocols::loops::Loops coordconstraint_segments_;
	if( constrain_relax_segments_ ){
		coordconstraint_segments_.read_loop_file(  option[ OptionKeys::relax::constrain_relax_segments ]() );
	}else{
		coordconstraint_segments_.add_loop( protocols::loops::Loop( 1, pose.total_residue(), 1 ) );
	}

	if ( constrain_coords_ ) {
		core::Size nnonvrt_cst_target = constraint_target_pose.total_residue();
		core::Size nnonvrt_pose = pose.total_residue();

		while ( pose.residue( nnonvrt_pose ).aa() == core::chemical::aa_vrt ) { nnonvrt_pose--; }
		while ( constraint_target_pose.residue( nnonvrt_cst_target ).aa() == core::chemical::aa_vrt ) { nnonvrt_cst_target--; }

		//fpd This isn't strictly necessary, the cst target only needs to define all constrained residues
		//fpd Probably still a good sanity check
		if ( nnonvrt_pose != nnonvrt_cst_target ) {
			std::cerr << "ERROR coord constraint pose length mismatch with input pose: " << nnonvrt_cst_target << " vs. " << nnonvrt_pose << std::endl;
			utility_exit();
		}

		if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
			pose.append_residue_by_jump
				( *ResidueFactory::create_residue( pose.residue(1).residue_type_set().name_map( "VRT" ) ),
					pose.total_residue()/2 );
		}

		//fpd if -relax::coord_cst_width is given, the code instead uses _bounded_ constraints on
		//fpd    CA positions, of the specified width
		if (!option[ OptionKeys::relax::coord_cst_width ].user() ) {
			Size nres = pose.total_residue();
			Real const coord_sdev( option[ OptionKeys::relax::coord_cst_stdev ] );
				// default is 0.5 (from idealize) -- maybe too small
			for ( Size i = 1; i<= nres; ++i ) {
				if ( (Size)i==(Size)pose.fold_tree().root() ) continue;
				if( coordconstraint_segments_.is_loop_residue( i ) ) {
					Residue const & nat_i_rsd( constraint_target_pose.residue(i) );
					for ( Size ii = 1; ii<= nat_i_rsd.last_backbone_atom(); ++ii ) {
						pose.add_constraint( new CoordinateConstraint(
							AtomID(ii,i), AtomID(1,nres), nat_i_rsd.xyz( ii ),
							new HarmonicFunc( 0.0, coord_sdev ) ) );
					}
				}
			}
		} else {
			Size nres = pose.total_residue();
			Real const cst_width( option[ OptionKeys::relax::coord_cst_width ]() );
			Real const coord_sdev( option[ OptionKeys::relax::coord_cst_stdev ] );
			for ( Size i = 1; i<= nres - 1; ++i ) {
				if( coordconstraint_segments_.is_loop_residue( i ) ) {
					Residue const & nat_i_rsd( constraint_target_pose.residue(i) );
					for ( Size ii = 1; ii<= nat_i_rsd.last_backbone_atom(); ++ii ) {
						pose.add_constraint( new CoordinateConstraint(
							AtomID(ii,i), AtomID(1,nres), nat_i_rsd.xyz( ii ),
							new BoundFunc( 0, cst_width, coord_sdev, "xyz" )) );
					}
				}
			}
		}
		if ( get_scorefxn()->get_weight( coordinate_constraint ) == 0 ) {
			get_scorefxn()->set_weight( coordinate_constraint, 0.5 );
		}

		local_movemap.set_jump( pose.num_jump(), true );
	} // if constrain_coords_


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
		if ( get_current_job() && get_current_job()->output_file_name() != "" ) {
			silent_file = get_current_job()->output_file_name()+silent_file;
		} else silent_file = "bla"+silent_file;

		SilentFileData sfd;
		//filename might have been changed -- e.g., to also have an MPI rank in there

		//		ProteinSilentStruct pss;
		io::silent::SilentStructOP pss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
		pss->fill_struct( pose, get_current_tag() );

		sfd.write_silent_struct( *pss, silent_file, false /* bWriteScoresOnly */ );
	} // if option[ out::file::silent ].user()
}

} // namespace relax
} // namespace protocols
