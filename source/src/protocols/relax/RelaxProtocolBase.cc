// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_protocols
/// @brief protocols that are specific to relax
/// @details
/// @author Mike Tyka, Monica Berrondo

//#include <protocols/jobdist/Jobs.hh>

#include <protocols/relax/cst_util.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/util.hh>

#include <protocols/relax/AtomCoordinateCstMover.hh>

#include <core/kinematics/MoveMap.hh>



#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pose/util.hh>
#include <core/io/pdb/build_pose_as_is.hh>

#include <core/select/movemap/MoveMapFactory.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/palette/PackerPalette.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>

#include <protocols/loops/Loops.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//*only for debug structures
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// ObjexxFCL Headers

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

#include <basic/Tracer.hh>

#include <core/conformation/Residue.hh>
#include <protocols/jd2/util.hh>

#include <utility>
#include <utility/vector1.hh>

#include <basic/citation_manager/CitationCollectionBase.hh> // AUTO IWYU For CitationCollectionList

static basic::Tracer TR( "protocols.relax.ClassicRelax" );

using namespace core;
////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace relax {


RelaxProtocolBase::RelaxProtocolBase( core::scoring::ScoreFunctionOP score_in ) :
	parent( "RelaxProtocol" ),
	min_type_("lbfgs_armijo_nonmonotone"),
	scorefxn_(std::move( score_in )),
	task_factory_(/* NULL */)
{
	set_defaults();
}

RelaxProtocolBase::RelaxProtocolBase( std::string const & movername ) :
	parent( movername ),
	min_type_("lbfgs_armijo_nonmonotone"),
	scorefxn_( /* 0 */ ),
	task_factory_(/* NULL */)
{
	set_defaults();
}

RelaxProtocolBase::RelaxProtocolBase( std::string const & movername, core::scoring::ScoreFunctionOP score_in ) :
	parent( movername ),
	min_type_("lbfgs_armijo_nonmonotone"),
	scorefxn_(std::move(score_in )),
	task_factory_(/* NULL */)
{
	set_defaults();
}


RelaxProtocolBase::~RelaxProtocolBase() = default;

void RelaxProtocolBase::initialize_movemap(
	core::pose::Pose const & pose,
	core::kinematics::MoveMap & movemap
)
{
	using namespace core::id;

	// If the user has provided a MoveMapFactory, then ignore all of the other
	// behaviors that are set for this instance.
	if ( movemap_factory_ ) {
		movemap_factory_->edit_movemap_given_pose( pose, movemap );
		if ( fix_omega_ || minimize_bond_lengths_ || minimize_bond_angles_ ) {
			TR.Warning << "Because a MoveMapFactory is set for the RelaxProtocolBase, the" <<
				" movemap-altering behaviors controlled by [" << ( fix_omega_ ? " fix_omega_" : "" ) <<
				( minimize_bond_lengths_ ? " minimize_bond_lengths_" : "" ) <<
				( minimize_bond_angles_ ? " minimize_bond_angles_" : "" ) << "] flag(s), which has (have) been set to true," <<
				" are being ignored. These behaviors must be set through the MoveMapFactory if you want"
				" to enable them." << std::endl;
		}
	} else {

		if ( fix_omega_ ) {
			for ( core::Size i=1; i<=pose.size(); ++i ) {
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
				for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
					core::chemical::AtomIndices const & ii_mainchain_atoms( pose.residue(ii).mainchain_atoms() );
					for ( core::Size jj = 1; jj <= ii_mainchain_atoms.size(); ++jj ) {
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
				for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
					core::conformation::Residue const &res_i = pose.residue(ii);
					for ( core::Size jj = 1; jj <= res_i.natoms(); ++jj ) {
						if ( res_i.atom_is_backbone(jj) ) continue;
						movemap.set( DOF_ID( AtomID( jj, ii ), core::id::D ), true );
					}
				}
			} else if ( minimize_bondlength_subset_ == 3 ) {
				for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
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
				for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
					core::chemical::AtomIndices const & ii_mainchain_atoms( pose.residue(ii).mainchain_atoms() );
					for ( core::Size jj = 1; jj <= ii_mainchain_atoms.size(); ++jj ) {
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
				for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
					core::conformation::Residue const &res_i = pose.residue(ii);
					for ( core::Size jj = 1; jj <= res_i.natoms(); ++jj ) {
						if ( res_i.atom_is_backbone(jj) ) continue;
						movemap.set( DOF_ID( AtomID( jj, ii ), core::id::THETA ), true );
					}
				}
			} else if ( minimize_bondangle_subset_ == 3 ) {
				for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
					core::conformation::Residue const &res_i = pose.residue(ii);
					if ( res_i.type().has( " C  ") ) {
						movemap.set( DOF_ID( AtomID( res_i.atom_index(" C  "), ii ), core::id::THETA ), true );
					}
				}
			} else if ( minimize_bondangle_subset_ == 4 ) {
				for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
					core::conformation::Residue const &res_i = pose.residue(ii);
					if ( res_i.type().has( " CB ") ) {
						movemap.set( DOF_ID( AtomID( res_i.atom_index(" CB "), ii ), core::id::THETA ), true );
					}
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
	ramp_down_constraints_ = option[ OptionKeys::relax::ramp_constraints ]();
	if ( constrain_coords_ ) {
		if ( !option[ OptionKeys::relax::ramp_constraints ].user() ) {
			ramp_down_constraints_ = true;
		}
	}
	constrain_relax_segments_ = option[ OptionKeys::relax::constrain_relax_segments ].user();
	limit_aroma_chi2_ = option[ OptionKeys::relax::limit_aroma_chi2 ]();
}

void RelaxProtocolBase::set_default_movemap(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	movemap_ = utility::pointer::make_shared< core::kinematics::MoveMap >();

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

core::kinematics::MoveMapCOP
RelaxProtocolBase::get_movemap() const { return movemap_; }

core::kinematics::MoveMapOP
RelaxProtocolBase::get_movemap() { return movemap_; }

const core::scoring::ScoreFunctionCOP
RelaxProtocolBase::get_scorefxn() const { return scorefxn_; }

core::pack::task::TaskFactoryOP const &
RelaxProtocolBase::get_task_factory() const { return task_factory_; }

void RelaxProtocolBase::set_movemap( core::kinematics::MoveMapOP movemap )
{
	movemap_ = movemap;
}

void RelaxProtocolBase::set_movemap_factory( core::select::movemap::MoveMapFactoryOP mm_factory )
{
	movemap_factory_ = mm_factory;
}

void RelaxProtocolBase::set_scorefxn( core::scoring::ScoreFunctionOP scorefxn ) { scorefxn_ = scorefxn; }
void RelaxProtocolBase::set_task_factory( core::pack::task::TaskFactoryOP task_factory ) { task_factory_ = task_factory; }

void RelaxProtocolBase::cartesian( bool newval ) { cartesian_ = newval; }
void RelaxProtocolBase::min_type( std::string min_type ) { min_type_ = min_type; }
void RelaxProtocolBase::max_iter( core::Size max_iter ) { max_iter_ = max_iter; }
void RelaxProtocolBase::dry_run( bool setting ) { dry_run_ = setting; }

void RelaxProtocolBase::constrain_relax_to_native_coords( bool constrain_relax_to_native_coords ) { constrain_relax_to_native_coords_ = constrain_relax_to_native_coords; }
void RelaxProtocolBase::constrain_relax_to_start_coords(  bool constrain_relax_to_start_coords ) { constrain_relax_to_start_coords_ = constrain_relax_to_start_coords; }
void RelaxProtocolBase::constrain_coords( bool constrain_coords ) { constrain_coords_ = constrain_coords; }
void RelaxProtocolBase::coord_constrain_sidechains( bool coord_constrain_sidechains ) {
	coord_constrain_sidechains_ = coord_constrain_sidechains;
}
void RelaxProtocolBase::constrain_relax_segments( bool constrain_relax_segments ) { constrain_relax_segments_ = constrain_relax_segments; }
void RelaxProtocolBase::ramp_down_constraints( bool ramp_down_constraints ) {
	ramp_down_constraints_ =  ramp_down_constraints;
}

void RelaxProtocolBase::minimize_bond_lengths( bool minimize_bond_lengths ) { minimize_bond_lengths_ = minimize_bond_lengths; }
void RelaxProtocolBase::minimize_bond_angles( bool minimize_bond_angles ) { minimize_bond_angles_ = minimize_bond_angles; }

/// @brief Provide the citation.
void
RelaxProtocolBase::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {

	citations.add( scorefxn_ );

	if ( task_factory_ != nullptr ) {
		if ( task_factory_->has_packer_palette() ) {
			debug_assert( task_factory_->packer_palette() != nullptr ); //Should be true.
			citations.add( task_factory_->packer_palette() );
		}
		for ( auto const & entry: *task_factory_ ) {
			citations.add( entry );
		}
	}
}

core::scoring::ScoreFunctionOP RelaxProtocolBase::get_scorefxn() { return scorefxn_; }

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
			( basic::options::option[ basic::options::OptionKeys::in::detect_disulf ]() ) ) {  //changed to make sure it doesn't flag "-detect_disulf false" as true

		std::string weight_set("score12_justdisulfides");
		core::scoring::ScoreFunctionOP disulf_score_only(scoring::ScoreFunctionFactory::create_score_function( weight_set));

		protocols::minimization_packing::PackRotamersMoverOP full_repack;
		core::pack::task::PackerTaskOP task;
		task = pack::task::TaskFactory::create_packer_task( pose );

		utility::vector1<bool> allow_repack(pose.size(), false);

		for ( core::Size i = 1; i<= pose.size() ; ++i ) {
			allow_repack[i] = movemap_->get_chi(i);
		}

		if ( basic::options::option[ basic::options::OptionKeys::relax::chi_move].user() ) {
			bool const repack = basic::options::option[ basic::options::OptionKeys::relax::chi_move]();
			allow_repack.assign( pose.size(), repack);
		}

		task->initialize_from_command_line().restrict_to_repacking().restrict_to_residues(allow_repack);
		task->or_include_current( true );
		full_repack = utility::pointer::make_shared< protocols::minimization_packing::PackRotamersMover >( disulf_score_only, task );

		( *disulf_score_only )( pose );

		full_repack->apply(pose);
	}
}


void RelaxProtocolBase::set_up_constraints(
	core::pose::Pose & pose,
	core::kinematics::MoveMap & local_movemap
) {
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
		//  if ( !option[ OptionKeys::relax::ramp_constraints ].user() ) {
		//   ramp_down_constraints_ = true;
		//  }

		protocols::relax::AtomCoordinateCstMover coord_cst_mover;
		if ( constrain_relax_to_native_coords_ ) {
			std::string const & native_pdb_fname = basic::options::option[ basic::options::OptionKeys::in::file::native ]();
			if ( native_pdb_fname.empty() ) {
				utility_exit_with_message("Native pose needs to be specified with -in:file:native for -relax::constrain_relax_to_native_coords");
			}
			core::pose::PoseOP ref_pose( new core::pose::Pose() );
			// We need to load this as a fullatom structure (the default)
			// Relax is intrinsically a fullatom protocol, so loading it as centroid doesn't make sense.
			core::io::pdb::build_pose_from_pdb_as_is( *ref_pose, native_pdb_fname );
			coord_cst_mover.set_refstruct( ref_pose );
		}

		if ( constrain_relax_segments_ ) {
			coord_cst_mover.set_loop_segments( utility::pointer::make_shared< protocols::loops::Loops >(  option[ OptionKeys::relax::constrain_relax_segments ]() ) );
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

		for ( core::Size i_cst = 1; i_cst <= cst_files_.size(); ++i_cst ) {
			std::string const filename = cst_files( i_cst );
			ConstraintSetOP user_csts
				= ConstraintIO::get_instance()->read_constraints_new( filename,
				utility::pointer::make_shared< ConstraintSet >(), pose );
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

	//cyclic

	if ( option[ OptionKeys::relax::cyclic_peptide ]() ) {
		cyclize_pose( pose );
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

		SilentFileOptions opts;
		SilentFileData sfd( opts );
		//filename might have been changed -- e.g., to also have an MPI rank in there

		//  ProteinSilentStruct pss;
		io::silent::SilentStructOP pss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( opts );
		pss->fill_struct( pose, get_current_tag() );

		sfd.write_silent_struct( *pss, silent_file, false /* bWriteScoresOnly */ );
	} // if option[ out::file::silent ].user()
}

} // namespace relax
} // namespace protocols
