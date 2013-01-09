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
		// 0 Default  all bondlengths
		// 1          backbone only
		// 2          sidechain only
		// 3          CA only (Ca-C,Ca-N and Ca-Cb)

		if (minimize_bondlength_subset_ == 0) {
			movemap.set( core::id::D, true );
		} else if ( minimize_bondlength_subset_ == 1) {
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::chemical::AtomIndices const & ii_mainchain_atoms( pose.residue(ii).mainchain_atoms() );
				for ( Size jj = 1; jj <= ii_mainchain_atoms.size(); ++jj ) {
					//if ( jj == 1 ) {
					//	if ( ii > 1 && pose.residue(ii).is_bonded( ii-1 ) && !pose.residue(ii).has_variant_type("CUTPOINT_UPPER")) {
					//		movemap.set( DOF_ID( AtomID( ii_mainchain_atoms[ jj ], ii ), core::id::D ), true );
					//	}
					//} else {
						movemap.set( DOF_ID( AtomID( ii_mainchain_atoms[ jj ], ii ), core::id::D ), true );
					//}
				}
			}
		} else if ( minimize_bondlength_subset_ == 2) {
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::conformation::Residue const &res_i = pose.residue(ii);
				for ( Size jj = 1; jj <= res_i.natoms(); ++jj ) {
					if (res_i.atom_is_backbone(jj)) continue;
					movemap.set( DOF_ID( AtomID( jj, ii ), core::id::D ), true );
				}
			}
		} else if ( minimize_bondlength_subset_ == 3) {
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::conformation::Residue const &res_i = pose.residue(ii);
				if (res_i.type().has_atom_name( " C  "))
					movemap.set( DOF_ID( AtomID( res_i.atom_index(" C  "), ii ), core::id::D ), true );
				if (res_i.type().has_atom_name( " CA "))
					movemap.set( DOF_ID( AtomID( res_i.atom_index(" CA "), ii ), core::id::D ), true );
				if (res_i.type().has_atom_name( " CB "))
					movemap.set( DOF_ID( AtomID( res_i.atom_index(" CB "), ii ), core::id::D ), true );
			}
		}
	}

	if ( minimize_bond_angles_ ) {


		// 0 Default  all bondangles
		// 1          backbone only
		// 2          sidechain only
		// 3          tau only
		// 4          Ca-Cb only
		if (minimize_bondangle_subset_ == 0) {
			movemap.set( core::id::THETA, true );
		} else if ( minimize_bondangle_subset_ == 1) {
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::chemical::AtomIndices const & ii_mainchain_atoms( pose.residue(ii).mainchain_atoms() );
				for ( Size jj = 1; jj <= ii_mainchain_atoms.size(); ++jj ) {
					//if ( jj == 1 || jj == 2 ) {  //fpd  add jj==2
					//	if ( ii > 1 && pose.residue(ii).is_bonded( ii-1 ) && !pose.residue(ii).has_variant_type("CUTPOINT_UPPER")) {
					//		movemap.set( DOF_ID( AtomID( ii_mainchain_atoms[ jj ], ii ), core::id::THETA ), true );
					//	}
					//} else {
						movemap.set( DOF_ID( AtomID( ii_mainchain_atoms[ jj ], ii ), core::id::THETA ), true );
					//}
				}
			}
		} else if ( minimize_bondangle_subset_ == 2) {
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::conformation::Residue const &res_i = pose.residue(ii);
				for ( Size jj = 1; jj <= res_i.natoms(); ++jj ) {
					if (res_i.atom_is_backbone(jj)) continue;
					movemap.set( DOF_ID( AtomID( jj, ii ), core::id::THETA ), true );
				}
			}
		} else if ( minimize_bondangle_subset_ == 3) {
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::conformation::Residue const &res_i = pose.residue(ii);
				if (res_i.type().has_atom_name( " C  "))
					movemap.set( DOF_ID( AtomID( res_i.atom_index(" C  "), ii ), core::id::THETA ), true );
			}
		} else if ( minimize_bondangle_subset_ == 4) {
			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				core::conformation::Residue const &res_i = pose.residue(ii);
				if (res_i.type().has_atom_name( " CB "))
					movemap.set( DOF_ID( AtomID( res_i.atom_index(" CB "), ii ), core::id::THETA ), true );
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
	minimize_bondlength_subset_ = option[ OptionKeys::relax::minimize_bondlength_subset ]();
	minimize_bondangle_subset_  = option[ OptionKeys::relax::minimize_bondangle_subset ]();

	//fpd extras
	cartesian_ = option[ OptionKeys::relax::cartesian ]();
	if ( option[ OptionKeys::relax::min_type ].user() )
		min_type_ = option[ OptionKeys::relax::min_type ]();
	else if (cartesian_ || minimize_bond_lengths_ || minimize_bond_angles_)
		min_type_ = "lbfgs_armijo_nonmonotone";  // default is different for cartesian/fbal minimization
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
	movemap_ = new core::kinematics::MoveMap();
	movemap_->set_jump( option[ OptionKeys::relax::jump_move ]() );
	movemap_->set_bb( option[ OptionKeys::relax::bb_move ]() );
	movemap_->set_chi( option[ OptionKeys::relax::chi_move ]() );
    if (option[ OptionKeys::in::file::movemap ].user()) {
        movemap_->init_from_file(option[ OptionKeys::in::file::movemap ]() );
    }
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
void RelaxProtocolBase::minimize_bondangle_subset( int setting ) { minimize_bondangle_subset_ = setting; }
void RelaxProtocolBase::minimize_bondlength_subset( int setting ) { minimize_bondlength_subset_ = setting; }

bool RelaxProtocolBase::minimize_bond_lengths() const { return minimize_bond_lengths_;}
bool RelaxProtocolBase::minimize_bond_angles() const { return minimize_bond_angles_;}
int RelaxProtocolBase::minimize_bondlength_subset() const { return minimize_bondlength_subset_;}
int RelaxProtocolBase::minimize_bondangle_subset() const { return minimize_bondangle_subset_;}

void RelaxProtocolBase::register_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static bool runonce( false );
	if ( runonce ) return;

	runonce = true;

	option.add_relevant( OptionKeys::relax::minimize_bond_lengths );
	option.add_relevant( OptionKeys::relax::minimize_bond_angles );
	option.add_relevant( OptionKeys::relax::minimize_bondlength_subset );
	option.add_relevant( OptionKeys::relax::minimize_bondlength_subset );

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

		if (basic::options::option[ basic::options::OptionKeys::relax::chi_move].user() ){
			bool const repack = basic::options::option[ basic::options::OptionKeys::relax::chi_move]();
			allow_repack.assign( pose.total_residue(), repack);
		}

		task->initialize_from_command_line().restrict_to_repacking().restrict_to_residues(allow_repack);
		task->or_include_current( true );
		if ( basic::options::option[ basic::options::OptionKeys::symmetry::symmetry_definition ].user() )  {
			full_repack = new simple_moves::symmetry::SymPackRotamersMover( disulf_score_only, task );
		} else {
			full_repack = new protocols::simple_moves::PackRotamersMover( disulf_score_only, task );
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
	core::id::SequenceMapping seq_map; // A mapping of pose -> constraint_target_pose numbering

	//fpd  Make ramping on by default if one of -constrain_relax_* is specified
	//fpd  Let it be overridden if '-ramp_constraints false' is specified
	if (constrain_coords_ && !option[ OptionKeys::relax::ramp_constraints ].user() ) {
		ramp_down_constraints_ = true;
	}

	if( constrain_relax_to_native_coords_ ){
		if ( get_native_pose() ) {
			constraint_target_pose = *get_native_pose();
		} else {
			utility_exit_with_message("Native pose needed for OptionKeys::relax::constrain_relax_to_native_coords");
		}
		// TODO: Allow for input of an alignment on commandline
		if (  pose.total_residue() == constraint_target_pose.total_residue() &&
				( !coord_constrain_sidechains_ || pose.sequence() != constraint_target_pose.sequence() ) ) {
			// We match in size and (for sidechains) sequence - we're looking at the traditional 1:1 mapping.
			seq_map = core::id::SequenceMapping::identity( pose.total_residue() );
		} else {
			// Try to match on a PDB-identity basis, or a sequence alignment basis if that fails.
			TR << "Length " << (coord_constrain_sidechains_?"and/or identities ":"") <<
					"of input structure and native don't match - aligning on PDB identity or sequence." << std::endl;
			seq_map = core::pose::sequence_map_from_pdbinfo( pose, constraint_target_pose );
		}
		// Align the native pose to the input pose to avoid rotation/translation based
		//  errors.
		//fpd  (Only if not already rooted on a VRT to avoid problems with density/symmetry)
		if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
			core::id::SequenceMapping rev_seq_map( seq_map ); // constraint_target_pose -> pose mapping
			rev_seq_map.reverse();
			core::sequence::calpha_superimpose_with_mapping(constraint_target_pose, pose, rev_seq_map);
		}
	} else {
		// Aligning to input - mapping is 1:1
		seq_map = core::id::SequenceMapping::identity( pose.total_residue() );
	}

	protocols::loops::Loops coordconstraint_segments_;
	if( constrain_relax_segments_ ){
		coordconstraint_segments_ = protocols::loops::Loops(  option[ OptionKeys::relax::constrain_relax_segments ]() );
	}else{
		coordconstraint_segments_.add_loop( protocols::loops::Loop( 1, pose.total_residue(), 1 ) );
	}

	if ( constrain_coords_ ) {
		//if -relax::coord_cst_width is given, the code instead uses _bounded_ constraints on
		//heavy atom positions, of the specified width
		Real const coord_sdev( option[ OptionKeys::relax::coord_cst_stdev ] );
    bool bounded( option[ OptionKeys::relax::coord_cst_width ].user() );
		Real const cst_width( option[ OptionKeys::relax::coord_cst_width ]() );

		// Add virtual root
		if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
			core::pose::addVirtualResAsRoot(pose);
		}

		core::Size nres = pose.total_residue();

		for ( Size i = 1; i<= nres; ++i ) {
			if ( pose.fold_tree().is_root(i) ) continue; // Skip root virtual atom.

			if ( coordconstraint_segments_.is_loop_residue( i ) ) {
				Size j(seq_map[i]);
				if( j == 0 ) continue;
				assert( j <= constraint_target_pose.total_residue() ); // Should be, if map was set up properly.

				Residue const & pose_i_rsd( pose.residue(i) );
				Residue const & targ_j_rsd( constraint_target_pose.residue(j) );
				core::Size last_atom( pose_i_rsd.last_backbone_atom() );
				core::Size last_targ_atom( targ_j_rsd.last_backbone_atom() );
				bool use_atom_names(false);
				if ( coord_constrain_sidechains_ ) {
					last_atom = pose_i_rsd.nheavyatoms();
					last_targ_atom = targ_j_rsd.nheavyatoms();
					use_atom_names = pose_i_rsd.name() != targ_j_rsd.name(); // Don't bother with lookup if they're the same residue type.
				}
				if ( !use_atom_names && last_atom != last_targ_atom ) {
					TR.Warning << "Warning: Coordinate constraint reference residue has different number of " << (coord_constrain_sidechains_?"heavy":"backbone") << " atoms: ref. "
						<< targ_j_rsd.name() << " (res " << j << ") versus  " << pose_i_rsd.name() << " (res " << i << "). - skipping." << std::endl;
					continue;
				}
				for ( Size ii = 1; ii<= last_atom; ++ii ) {
					Size jj(ii);
					if ( use_atom_names ) {
						std::string atomname( pose_i_rsd.atom_name(ii) );
						if ( ! targ_j_rsd.has(atomname) ) {
							TR.Debug << "Skip adding coordinate constraints for atom " << atomname << " of residue " << i << " (" << pose_i_rsd.name() <<
								") - not found in residue " << j << " (" << targ_j_rsd.name() << ") of reference structure." << std::endl;
							continue;
						}
						jj = targ_j_rsd.atom_index( atomname );
					}
					core::scoring::constraints::FuncOP function;
					if( bounded ) {
						function = new BoundFunc( 0, cst_width, coord_sdev, "xyz" );
					} else {
						function = new HarmonicFunc( 0.0, coord_sdev );
					}
					pose.add_constraint( new CoordinateConstraint(
						AtomID(ii,i), AtomID(1,nres), targ_j_rsd.xyz( jj ), function ) );
				} // for atom
			} // if(loop)
		} // for residue

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
		std::string output_file_name= protocols::jd2::current_output_name();
		if ( output_file_name != "" ) {
			silent_file = output_file_name + silent_file;
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
