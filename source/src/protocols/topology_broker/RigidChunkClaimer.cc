// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @detailed responsibilities:
/// @author Oliver Lange

// Unit Headers
#include <protocols/topology_broker/RigidChunkClaimer.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/claims/JumpClaim.hh>
#include <protocols/topology_broker/claims/CutClaim.hh>
#include <protocols/topology_broker/claims/BBClaim.hh>
#include <protocols/topology_broker/SequenceNumberResolver.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/id/DOF_ID.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <protocols/comparative_modeling/util.hh>
// AUTO-REMOVED #include <core/kinematics/Exceptions.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.tmpl.hh>
#include <protocols/loops/LoopsFileIO.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/comparative_modeling/ThreadingJob.hh>
#include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/chemical/VariantType.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/xyz.functions.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <vector>
// AUTO-REMOVED #include <iterator>

#include <protocols/jd2/Job.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>



static thread_local basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

using namespace core;

//helper function
protocols::loops::Loops generate_rigid_from_alignment( pose::Pose query_pose, core::sequence::SequenceAlignment const& align, Size min_loop_size ) {
	using core::Size;
	using namespace basic::options;

	loops::LoopsOP loops = comparative_modeling::loops_from_alignment( query_pose.total_residue(), align, min_loop_size );

	// this is now done in select_parts()
 	// randomly grow loops by N residues (4 is a good amount)
	return loops->invert( query_pose.total_residue() );
}

// having AdjacentJumps causes problems in fix_internal_coords_of_siblings.
// can have more than two child on C(=0) for instance
RigidChunkClaimer::RigidChunkClaimer()
	: bExclusive_( true ),
		bAllowAdjacentJumps_( false ),
		bUseInputPose_( true ),
		bRigidInRelax_( false ) {}

RigidChunkClaimer::RigidChunkClaimer( pose::Pose const& input_pose, loops::Loops rigid ) :
	input_pose_( input_pose ),
	centroid_input_pose_( input_pose ),
	rigid_core_( rigid ),
	bExclusive_( true ),
	bAllowAdjacentJumps_( false ),
	bUseInputPose_( true ),
	bRigidInRelax_( false )
{
	if ( centroid_input_pose_.total_residue() && centroid_input_pose_.is_fullatom() )
		core::util::switch_to_residue_type_set(	centroid_input_pose_, chemical::CENTROID );
}

void RigidChunkClaimer::set_defaults() {
	Parent::set_defaults();
	bUseInputPose_ = true;
	bUseThreadingJobLoops_ = false;
	min_loop_size_ = 0;
	using namespace basic::options;
	random_grow_loops_by_ = option[ OptionKeys::loops::random_grow_loops_by ]();
}

void RigidChunkClaimer::receive_message( ClaimerMessage& cm ) {
	if ( typeid( cm ) == typeid( CM_SuggestFixResidue ) ) {
		CM_SuggestFixResidue& msg = dynamic_cast< CM_SuggestFixResidue& >( cm );

		//find good residue
		if ( !current_rigid_core_.size() ) return;
		msg.good_fix_pos_ =	current_rigid_core_.begin()->start() + numeric::nint( 0.5 * current_rigid_core_.begin()->size() ) - 1;
		msg.received_by_me( this );
	}
}

bool RigidChunkClaimer::read_tag( std::string tag, std::istream& is )
{
	loops::PoseNumberedLoopFileReader reader;
	reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );

	if ( tag == "pdb" || tag == "PDB" || tag == "pdb:" || tag == "PDB_FILE" ) {
		std::string file;
		is >> file;
		core::import_pose::pose_from_pdb( input_pose_,
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ),
			file );
		runtime_assert( input_pose_.is_fullatom() );
	} else if ( tag == "REGION" ) {
		loops::SerializedLoopList loops = reader.read_pose_numbered_loops_file( is, type(), false /*no strict checking */ );
		rigid_core_ = loops::Loops( loops );
	} else if ( tag == "region_file" || tag == "REGION_FILE" ) {
		std::string file;
		is >> file;
		std::ifstream infile( file.c_str() );

		if (!infile.good()) {
			utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + file + "'" );
		}
		loops::SerializedLoopList loops = reader.read_pose_numbered_loops_file( infile, file, false /*no strict checking */ );
		rigid_core_ = loops::Loops( loops ); // <==
	} else if ( tag == "loop_file" || tag == "LOOP_FILE" ) {
		std::string file;
		is >> file;
		std::ifstream infile( file.c_str() );

		if (!infile.good()) {
			utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + file + "'" );
		}
		loops::SerializedLoopList loops = reader.read_pose_numbered_loops_file( infile, file, false /*no strict checking */ );
		protocols::loops::Loops loop_defs = loops::Loops( loops ); // <==
		loop_defs = loop_defs.invert( input_pose_.total_residue() );
		tr << "Rigid core: " << input_pose_.total_residue() << std::endl << loop_defs << std::endl;
		rigid_core_ = loop_defs;
	} else if ( tag == "NO_USE_INPUT_POSE" ) {
		bUseInputPose_ = false;
	} else if ( tag == "USE_THREADING_LOOPS" ) {
		bUseThreadingJobLoops_ = true;
	} else if ( tag == "MIN_LOOP_SIZE" ) {
		is >> min_loop_size_;
	} else if ( tag == "KEEP_FLEXIBLE" ) {
		runtime_assert( false );
		//right now not supported
		// have problems with JumpClaims being reissued ...
		// in this line: current_jumps_.push_back( foreign_claim.clone() ); //o
		// jumps get added to our own list ... in finalize_claims they will be send back to the Broker, who will call allow_claim again for this jump
		bExclusive_ = false;
	} else if ( tag == "RANDOM_GROW_LOOP_BY" ) {
		is >>  random_grow_loops_by_;
	} else if ( tag == "RIGID_IN_RELAX" ) {
		bRigidInRelax_ = true;
	}
	else return Parent::read_tag( tag, is );
	return true;
}

void RigidChunkClaimer::init_after_reading() {
	tr.Debug << type() << " initialized with input_pdb: " << input_pose_.sequence() << " and regions " << rigid_core_ << std::endl;
	centroid_input_pose_=input_pose_;
	if ( centroid_input_pose_.total_residue() && centroid_input_pose_.is_fullatom() ) core::util::switch_to_residue_type_set(	centroid_input_pose_, chemical::CENTROID );
}

void RigidChunkClaimer::select_parts() {
	current_rigid_core_.clear();
	int attempts = 0;
	if ( rigid_core_.size() == 0 ) return;//nothing to select
	while ( current_rigid_core_.size() == 0 && attempts < 50 ) {
		++attempts;
		for ( loops::Loops::const_iterator it = rigid_core_.begin(), eit = rigid_core_.end();
					it != eit; ++it ) {
			if ( numeric::random::rg().uniform() >= it->skip_rate() )  {
				current_rigid_core_.push_back( *it );
			}
		}
	}
	if ( current_rigid_core_.size() == 0 ) {
		current_rigid_core_ = rigid_core_;
	}
	if ( random_grow_loops_by_ > 0 ) {
		core::Size nres( current_rigid_core_[ current_rigid_core_.size() ].stop() + 200 ); //it doesn't matter for this where exactly nres is.
		loops::Loops loops( current_rigid_core_.invert( nres ) );
		loops.grow_all_loops( nres, random_grow_loops_by_ );
		tr.Info << "Enlarged loops: " << std::endl;
		tr.Info << loops << std::endl;
		current_rigid_core_ = loops.invert( nres );
	}
}

///@detail generate exclusive backbone claims for each residue of the rigid-chunk
/// jumps are not exclusive and are added later in final_claims --- want to reuse other jumps if possible
void RigidChunkClaimer::generate_claims( claims::DofClaims& new_claims ) {
	using namespace loops;
	tr.Trace << "rigid chunk -- generate claim " << std::endl;

	//stochastically select rigid_core ( if skip-rate is set >0, otherwise all loops are selected )
	tr.Trace << "region selected: " << current_rigid_core_ << std::endl;

	//	new_claims.push_back( new CutBiasClaim( *this ) ); we don't need this claim type --- always call manipulate_cut_bias
	for ( Loops::const_iterator loop_it = current_rigid_core_.begin(); loop_it != current_rigid_core_.end(); ++loop_it ) {
		for ( Size pos = loop_it->start(); pos <= loop_it->stop(); ++pos ) {
			new_claims.push_back( new claims::BBClaim( this,
																								 std::make_pair( label(), pos ),
																								 claims::DofClaim::EXCLUSIVE ) );
		}
	}
}

void RigidChunkClaimer::new_decoy() {
	select_parts(); //set members current_XXX
	current_jump_calculator_ = new JumpCalculator( current_rigid_core_, bAllowAdjacentJumps_ );
	current_jumps_.clear();
}

void RigidChunkClaimer::new_decoy( core::pose::Pose const& pose ) {
	//should we read input structures dofs?
	tr.Debug << "New decoy:" << std::endl;
	if ( bUseInputPose_ ) {
		input_pose_ = pose;
		centroid_input_pose_=input_pose_;
		if ( centroid_input_pose_.total_residue() && centroid_input_pose_.is_fullatom() ) core::util::switch_to_residue_type_set(	centroid_input_pose_, chemical::CENTROID );

		// use loops from ThreadingJob ???
		if ( bUseThreadingJobLoops_ ) {
			using namespace protocols::jd2;
			protocols::comparative_modeling::ThreadingJobCOP job = dynamic_cast< protocols::comparative_modeling::ThreadingJob const*  >( JobDistributor::get_instance()->current_job()->inner_job().get() );
			if ( job ) {
				tr.Debug << "------------------found ThreadingJob ... get loops " << std::endl;
				rigid_core_ = generate_rigid_from_alignment( input_pose_, job->alignment(), min_loop_size_ );
			}
		} //bUseThreading
		tr.Debug << "RigidChunk defined for " << rigid_core_ << std::endl;
	}
	new_decoy();
}

bool RigidChunkClaimer::allow_claim( claims::DofClaim const& foreign_claim ) {
	if ( foreign_claim.owner() == this ) return true; // always allow your own claims!

	// check foreign claim

	claims::BBClaimCOP bb_ptr( dynamic_cast< const claims::BBClaim* >( &foreign_claim ) );

	if ( bb_ptr && current_rigid_core_.is_loop_residue( bb_ptr->global_position() ) ) {
		if ( bExclusive_ ) { 	// if we want exclusive claim this is not acceptable
			return false;
		} else {
			// allow only the weakest claim. We want to initialize ourselves... don't know if we need to be so restrictive!
			if ( !(foreign_claim.right() == claims::DofClaim::CAN_INIT) ) return false;
		}
	} // DofClaim::BB

	claims::JumpClaimCOP jump_ptr( dynamic_cast< const claims::JumpClaim* >( &foreign_claim ) );

	if ( jump_ptr ) {
		runtime_assert( current_jump_calculator_ );
		if ( current_jump_calculator_->irrelevant_jump( jump_ptr->global_pos1(), jump_ptr->global_pos2() ) ) {
			return true;
		} else if ( !current_jump_calculator_->good_jump( jump_ptr->global_pos1(), jump_ptr->global_pos2() ) ) {
			return false;
		} else if ( bExclusive_ ) { // but probably a good jump --- since it has a physical reason.
			//reclaim the claim
			current_jumps_.push_back( new claims::JumpClaim( this, jump_ptr->local_pos1(), jump_ptr->local_pos2(), claims::DofClaim::EXCLUSIVE ) );
			return false;
		} else {
			current_jumps_.push_back( foreign_claim.clone() ); //ok - remember this jump, since it connects rigid1 and rigid2
		}
	} // DofClaim::JUMP

	claims::CutClaimCOP cut_ptr( dynamic_cast< const claims::CutClaim* >( &foreign_claim ) );

	if ( cut_ptr ) {
		for ( loops::Loops::const_iterator region = current_rigid_core_.begin(); region != current_rigid_core_.end(); ++region ) {

			//TODO: ensure that the label setting code is correctly functioning in this claimer

			claims::LocalPosition cut_position = cut_ptr->get_position();
			Size absolute_cut_position = broker().sequence_number_resolver().find_global_pose_number( cut_position );

			if ( absolute_cut_position >= region->start() &&
					 absolute_cut_position < region->stop() ) // cut claim can be at the chunk end
				return false; // no cuts within our rigid-core boundaries
		}
	}
	return true;
}

bool RigidChunkClaimer::accept_declined_claim( claims::DofClaim const& was_declined ) {
	tr.Warning << "[WARNING] RigidChunkClaimer couldn't get " << was_declined << std::endl;
	return false; // no tolerance here --- don't accept any of these
}

void RigidChunkClaimer::finalize_claims( claims::DofClaims& new_claims ) {
	claims::DofClaims extra_jumps;
	current_jump_calculator_->generate_rigidity_jumps( this, extra_jumps, label() );

	std::copy( extra_jumps.begin(), extra_jumps.end(), std::back_inserter( current_jumps_ ) );
	std::copy( current_jumps_.begin(), current_jumps_.end(), std::back_inserter( new_claims ) );
}

void fix_internal_coords_of_siblings( pose::Pose& pose, pose::Pose const& ref_pose, id::AtomID atom, id::AtomID ref_atom ) {
	runtime_assert( atom.rsd() >= 1 && atom.rsd() <= pose.total_residue() );
	runtime_assert( pose.conformation().atom_tree().has( atom ) );
	runtime_assert( ref_pose.conformation().atom_tree().has( ref_atom ) );

	bool has_par1( pose.conformation().atom_tree().atom( atom ).parent() );
	bool ref_has_par1( ref_pose.conformation().atom_tree().atom( ref_atom ).parent() );

	//folding direction matters for the angle we have to set...hence find the parent atoms and get the angle
	core::id::AtomID par1O;
	core::id::AtomID ref_par1O;
	if ( has_par1 && ref_has_par1 ) {
		par1O=pose.conformation().atom_tree().atom( atom ).parent()->id();
		std::string const & aname(pose.residue(par1O.rsd()).atom_name(par1O.atomno()));
		ref_par1O=core::id::AtomID( ref_pose.residue( par1O.rsd() ).atom_index( aname ), par1O.rsd() );
	}	else {
		tr.Warning << "cannot fix internal coords of " << atom << " in RigidChunk because 1st parent is missing " << std::endl;
		return;
	}
	bool has_par2( pose.conformation().atom_tree().atom( par1O ).parent() );
	bool ref_has_par2( ref_pose.conformation().atom_tree().atom( ref_par1O ).parent() );
	core::id::AtomID par2O;
	core::id::AtomID ref_par2O;
 	if ( has_par2 && ref_has_par2 ) {
		par2O=pose.conformation().atom_tree().atom( par1O ).parent()->id();
		std::string const & aname(pose.residue(par2O.rsd()).atom_name(par2O.atomno()));
		ref_par2O=core::id::AtomID( ref_pose.residue( par2O.rsd() ).atom_index( aname ), par2O.rsd() );
	} else {
		tr.Warning << "cannot fix internal coords of " << atom << " in RigidChunk because 2nd parent is missing " << std::endl;
		return;
	}
	runtime_assert( ref_pose.conformation().atom_tree().has( ref_par1O ) );
	runtime_assert( ref_pose.conformation().atom_tree().has( ref_par2O ) );
	runtime_assert( pose.conformation().atom_tree().has( par1O ) );
	runtime_assert( pose.conformation().atom_tree().has( par2O ) );

	core::Real angle( numeric::angle_radians( ref_pose.xyz( ref_atom ), ref_pose.xyz( ref_par1O ), ref_pose.xyz( ref_par2O ) ) );
	tr.Trace << "ref angle direct: " << angle << std::endl;
	pose.conformation().set_bond_angle(  par2O, par1O, atom, angle );

	id::DOF_ID torsion_offset_dof( atom, id::PHI );
	id::DOF_ID ref_torsion_offset_dof( ref_atom, id::PHI );
	core::Real value( ref_pose.conformation().atom_tree().dof( ref_torsion_offset_dof ) );
	pose.conformation().set_dof( torsion_offset_dof, value );
}

void fix_mainchain_connect( pose::Pose& pose, pose::Pose const& ref_pose, core::Size upper_residue ) {
	core::conformation::Residue const & prev_rsd( ref_pose.residue( upper_residue-1 ) );
	core::conformation::Residue const &      rsd( ref_pose.residue( upper_residue ) );
	core::Size const nbb_prev( prev_rsd.n_mainchain_atoms() );
	core::id::AtomID bbM1   ( prev_rsd.mainchain_atom( nbb_prev-2 ),  upper_residue-1 );
	core::id::AtomID bb0    ( prev_rsd.mainchain_atom( nbb_prev-1 ),  upper_residue-1 );
	core::id::AtomID bb1    ( prev_rsd.mainchain_atom( nbb_prev   ),  upper_residue-1 );
	core::id::AtomID bb2    (      rsd.mainchain_atom(        1   ),  upper_residue   );
	core::id::AtomID bb3    (      rsd.mainchain_atom(        2   ),  upper_residue   );
	core::id::AtomID bb4    (      rsd.mainchain_atom(        3   ),  upper_residue   );

	core::conformation::Residue const & ref_resi = ref_pose.residue( upper_residue );
	tr.Trace << "mainchain torsion: ref: " << ref_resi.mainchain_torsion( 1 ) << " atom-tree: "
					 << ref_pose.conformation().torsion_angle( bb1, bb2, bb3, bb4 ) << std::endl;

	core::conformation::Residue const & resi = pose.residue( upper_residue );
	tr.Trace << "mainchain torsion (before): conf: " << resi.mainchain_torsion( 1 ) << " atom-tree: "
					 << pose.conformation().torsion_angle( bb1, bb2, bb3, bb4 ) << std::endl;

	pose.conformation().set_bond_length( bb1, bb2, ref_pose.conformation().bond_length( bb1, bb2 ) );
	pose.conformation().set_bond_angle ( bb0, bb1, bb2, ref_pose.conformation().bond_angle( bb0, bb1, bb2 ) );
	pose.conformation().set_bond_angle ( bb1, bb2, bb3, ref_pose.conformation().bond_angle( bb1, bb2, bb3 ) );
	pose.conformation().set_torsion_angle( bbM1, bb0, bb1, bb2, ref_pose.conformation().torsion_angle( bbM1, bb0, bb1, bb2 ) );
	pose.conformation().set_torsion_angle( bb0, bb1, bb2, bb3, ref_pose.conformation().torsion_angle( bb0, bb1, bb2, bb3 ) );
	pose.conformation().set_torsion_angle( bb1, bb2, bb3, bb4, ref_pose.conformation().torsion_angle( bb1, bb2, bb3, bb4 ) );

	core::conformation::Residue const & new_resi = pose.residue( upper_residue ); //this should trigger update of coords and torsions
	tr.Trace << "mainchain torsion (after): conf: " << new_resi.mainchain_torsion( 1 ) << " atom-tree: "
					 << pose.conformation().torsion_angle( bb1, bb2, bb3, bb4 ) << std::endl;

	if ( prev_rsd.has( "O" ) ) {
		core::id::AtomID ref_atomO( prev_rsd.atom_index( "O" ), upper_residue-1 );
		core::id::AtomID atomO( pose.residue_type( upper_residue-1 ).atom_index( "O" ), upper_residue-1 );
		fix_internal_coords_of_siblings( pose, ref_pose, atomO, ref_atomO );
	}
	if ( rsd.has( "H" ) ) {
		core::id::AtomID ref_atomH( rsd.atom_index( "H" ), upper_residue );
		core::id::AtomID atomH( new_resi.atom_index( "H" ), upper_residue );
		runtime_assert( new_resi.has( "H" ) );
		fix_internal_coords_of_siblings( pose, ref_pose, atomH, ref_atomH );
	}

	if ( tr.Trace.visible() ) {
		bool ideal1( core::pose::is_ideal_position( upper_residue, ref_pose ) );
		if ( ideal1 && !core::pose::is_ideal_position( upper_residue, pose ) ) {
			tr.Warning << " pose in RigidChunkClaimer is not ideal at position " << upper_residue << " although template pose was ideal there " << std::endl;
		}

		bool ideal2( core::pose::is_ideal_position( upper_residue-1, ref_pose ) );
		if ( ideal2 && !core::pose::is_ideal_position( upper_residue-1, pose ) ) {
			tr.Warning << " pose in RigidChunkClaimer is not ideal at position " << upper_residue-1 << " although template pose was ideal there " << std::endl;
		}
	}
}

void copy_internal_coords( pose::Pose& pose, pose::Pose const& ref_pose, loops::Loops core ) {
	///fpd if there are post modifications to pose (not in ref_pose), we can't just copy ref_pose->pose
	///fpd    instead ... make xyz copy in rigid regions
	for ( loops::Loops::const_iterator region = core.begin(); region != core.end(); ++region ) {
		for (Size i=region->start(); i<=region->stop(); ++i) {
			core::conformation::Residue const &rsd_i = ref_pose.residue(i);
			pose.replace_residue ( i , rsd_i , false );
		}
	}

	if ( tr.Trace.visible() ) {
		tr.Trace << pose.fold_tree() << std::endl;
		tr.Trace << ref_pose.fold_tree() << std::endl;
	}

	///fpd fix connections
	///fpd this requires that the input pose have one flanking residue on each side of each region
	for ( loops::Loops::const_iterator region = core.begin(); region != core.end(); ++region ) {
		Size loop_start = region->start();
		Size loop_stop  = region->stop();

		bool lower_connect = ( loop_start > 1
													 && !pose.residue(loop_start).is_lower_terminus()
													 && !pose.fold_tree().is_cutpoint( loop_start-1 ) );
		bool upper_connect = ( loop_stop < pose.total_residue()
													 && !pose.residue(loop_stop).is_upper_terminus()
													 && !pose.fold_tree().is_cutpoint( loop_stop ) );

		if ( lower_connect ) {
			tr.Trace << "fixing lower connection for " << loop_start << std::endl;
			fix_mainchain_connect( pose, ref_pose, loop_start );
		} else {
			tr.Trace << "NOT fixing lower connection for " << loop_start << std::endl;
		}

		if ( upper_connect ) {
			tr.Trace << "fixing upper connection for " << loop_stop << std::endl;
			fix_mainchain_connect( pose, ref_pose, loop_stop+1 );
		} else {
			tr.Trace << "NOT fixing upper connection for " << loop_stop << std::endl;
		}
	}
}

void RigidChunkClaimer::adjust_relax_movemap(  core::kinematics::MoveMap& mm ) const {
	if ( bRigidInRelax_ ) {
		rigid_core_.switch_movemap( mm, id::BB, false );
	}
}

void RigidChunkClaimer::initialize_dofs( core::pose::Pose& pose, claims::DofClaims const& /* init_claims*/, claims::DofClaims& /*failed_init_claims*/ ) {
	//need to copy coords and jumps --- if chunks were idealized no problem .... but non-idealized stuff ?
	//also take care of fullatom vs centroid...
	//	tr.Warning << "[WARNING] *** use input structure of RigidChunkClaimer --- NEEDS TO BE IDEALIZED !!! *** \n";
	// to do this without idealized :

	// get rigid-body reorientation for first residue... this must be applied to all residues in the chunk,
	// and then just set from the coordinates.

	//in a sense this is not necessary since we asked for exclusive rights for the stuff in RIGID ...
	/*	for ( DofClaims::const_iterator claim=init_claims.begin(); claim!=init_claims.end(); ++claim ) {
		if ( (*claim)->owner() != this ) continue;

		}*/
	//	superimpose_chain (pose, input_pose_, rigid_core_ );


	//need to have same number of residues for fold-tree transfer...
	// would be nice to drop this restriction but for now, fill up with missing density...

	//fpd runtime_assert( pose.total_residue() == centroid_input_pose_.total_residue() );
	//fpd   really, we just have to make sure that #residues in the input pose > index of last rigid chunk
	//fpd   (strictly greater-than since we have to have a flanking res on each side of each region)
	//fpd   we still need missing dens in the gaps (but not at c term now!)
	core::Size lastChunk=1;
	for ( loops::Loops::const_iterator it = current_rigid_core_.begin(); it!=current_rigid_core_.end(); ++it )
		lastChunk = std::max( lastChunk , it->stop() );
	runtime_assert ( lastChunk <= centroid_input_pose_.total_residue() );

	bool missing_density( false );
	//sanity check: no missing density in backbon in any of the rigid_core residues?
	for ( loops::Loops::const_iterator it = current_rigid_core_.begin(); it!=current_rigid_core_.end(); ++it ) {
		for ( Size pos = it->start(); pos <=it->stop(); ++pos ) {
			// Do we really have Sidechains ?
			// check this my making sure that no SC atom is more than 20A (?) away from CA
			numeric::xyzVector< core::Real> ca_pos = input_pose_.residue( pos ).atom("CA").xyz();
			numeric::xyzVector< core::Real> n_pos = input_pose_.residue( pos ).atom("N").xyz();
			numeric::xyzVector< core::Real> o_pos = input_pose_.residue( pos ).atom("O").xyz();
			if ( ( n_pos - ca_pos).length() > 20 || ( ( n_pos - o_pos ).length() > 20 ) ) {
				tr.Error << "missing backbone in rigid-chunk at " << pos << std::endl;
				missing_density = true;
			}
		}
	}
	if ( missing_density ) throw utility::excn::EXCN_BadInput( " missing density in backbone of rigid-chunk... check your LOOP definition");
	//centroid_input_pose_.fold_tree( pose.fold_tree() );
	runtime_assert( !pose.is_fullatom()  );
	runtime_assert( !centroid_input_pose_.is_fullatom() );
	copy_internal_coords( pose, centroid_input_pose_, current_rigid_core_ );
}

///@brief multiply your bias to this -- if its zero don't change that, i.e., multiply only
void RigidChunkClaimer::manipulate_cut_bias( utility::vector1< core::Real >& cut_bias ) {
	current_rigid_core_.transfer_to_residue_vector( cut_bias, 0.0 );
}


void RigidChunkClaimer::switch_to_fullatom( core::pose::Pose& pose , utility::vector1< bool > bNeedToRepack ) const {
	loops::Loops const& region( current_rigid_core_ );

	// copy sidechain torsions from input pose
	tr.Debug << "copy side chains for residues with * / missing density residues with - ";
	for ( loops::Loops::const_iterator it = region.begin(); it!=region.end(); ++it ) {
		for ( Size pos = it->start(); pos <=it->stop(); ++pos ) {
			bNeedToRepack[ pos ] = false; //in principle our residues don't need a repack since we have a side-chains for them.
			// Do we really have Sidechains ?
			// check this my making sure that no SC atom is more than 20A (?) away from CA
			numeric::xyzVector< core::Real> ca_pos = input_pose_.residue( pos ).atom("CA").xyz();
			for ( Size j = 1; j<=input_pose_.residue( pos ).natoms(); ++j ) {
				if ( ( ca_pos - input_pose_.residue( pos ).atom( j ).xyz()).length() > 20 ) {
					tr.Debug << "-" << pos << " ";
					bNeedToRepack[ pos ] = true;
					break; //one bad atom is enough
				}
			}
			//copy sidechains only for non-loop regions
			if ( !bNeedToRepack[ pos ] ) {
				tr.Debug <<  "*" << pos << " ";
				bool const lower_cut ( pose.residue( pos ).has_variant_type( chemical::CUTPOINT_LOWER ) );
				bool const upper_cut ( pose.residue( pos ).has_variant_type( chemical::CUTPOINT_UPPER ) );
				//bool const disulf ( pose.residue( pos ).has_variant_type( chemical::DISULFIDE ) );
				pose.replace_residue( pos, input_pose_.residue( pos ), true /*orient backbone*/ );
				if ( lower_cut ) core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, pos );
				if ( upper_cut ) core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, pos );
			}
		} //for all rigid residues
	}
	tr.Debug << " that have not moved from start-structure" << std::endl;
} //switch to fullatom


// ================================================================================
// ======================== JumpCalculator ========================================
// ================================================================================
RigidChunkClaimer::JumpCalculator::JumpCalculator( loops::Loops const& rigid, bool bAllowAdjacentJumps )
	: rigid_ ( rigid ),
		visited_( rigid.size(), 0 ),
		new_nr_( 1 ),
		bAllowAdjacentJumps_( bAllowAdjacentJumps ) {}

bool is_not_neighbor_to_rigid( loops::Loops const& rigid, Size pos1, Size pos2 ) {
	Size up1 = rigid.loop_index_of_residue( pos1-1 );
	Size in1 = rigid.loop_index_of_residue( pos1 );
	Size down1 = rigid.loop_index_of_residue( pos1+1 );

	Size down2 = rigid.loop_index_of_residue( pos2+1 );
	Size in2 = rigid.loop_index_of_residue( pos2 );
	Size up2 = rigid.loop_index_of_residue( pos2-1 );
	if ( !in1 && ( up1 && down1 ) ) return false;
	if ( !in2 && ( up2 && down2 ) ) return false;
	return true;
}

bool connects_rigid_regions( loops::Loops const& rigid, Size pos1, Size pos2 ) {
	//TODO: this is probably easier with label checks...
	Size rigid1 = rigid.loop_index_of_residue( pos1 );
	Size rigid2 = rigid.loop_index_of_residue( pos2 );
	return rigid1 && rigid2;
}

bool RigidChunkClaimer::JumpCalculator::irrelevant_jump( Size global_start, Size global_end ) {
	//TODO make better use of local positions
	if( tr.Trace.visible() ) {
		tr.Trace << "Irrelevant_jump check for " << global_start << "->" << global_end << std::endl;
		tr.Trace << "connects_rigid: " << connects_rigid_regions( rigid_,	global_start, global_end ) << std::endl;
		tr.Trace << "is not_neighbor_to_rigid: " << is_not_neighbor_to_rigid( rigid_, global_start, global_end ) << std::endl;
		tr.Trace << "bAllowAdjacent: " << bAllowAdjacentJumps_ << std::endl;
	}

	if ( !connects_rigid_regions( rigid_, global_start, global_end ) ) {
		return bAllowAdjacentJumps_ || is_not_neighbor_to_rigid( rigid_, global_start, global_end );//jump doesn't connect two rigid regions -- irrelevant
	}
	return false; //either connects rigid regions or is neighbor to rigid region
}

///@brief check if this (relevant) jump is compatible with our rigid-structure
///   not on the same continuous stretch of rigid residues  ( we don't allow cuts within rigid stretches )
///   not connecting already conntected rigid stretches
/// if it connects two unconnected rigid stretches ---> its a good jump we'll keep it,
/// *** --> update visited_ ***
bool RigidChunkClaimer::JumpCalculator::good_jump( core::Size global_start, core::Size global_end ) {
	//TODO make better use of local positions

	Size up_loop( rigid_.loop_index_of_residue( global_start ) );
	Size down_loop( rigid_.loop_index_of_residue( global_end ) );

	//we arrive here only if jump is not irrelevant...
	//if we don't allow adjacent jump that means this jump connects rigid regions or is a bad neighbor
	if ( !bAllowAdjacentJumps_ && !connects_rigid_regions( rigid_, global_start, global_end ) ) return false;

	runtime_assert( up_loop && down_loop ); //since this would be irrelevant --- already checked.

	//don't allow jumps within same loop -- that will make a nasty cut within rigid region
	if ( up_loop == down_loop ) return false;

	// at this point rigid1 and rigid2 refer to rigid regions but not the same --- this might be useful if we need to connect rigid regions

	runtime_assert( visited_.size() >= up_loop );
	runtime_assert( visited_.size() >= down_loop );
	// jump touches unvisited regions or visited but yet disconnected regions
	if ( !visited_[ up_loop ] || !visited_[ down_loop ] || ( visited_[ up_loop ] != visited_[ down_loop ] ) ) {

		// decide upon 3 cases: both nodes unvisited, 1 node visited, both nodes visited
		// case0 : both new tag with new jump_nr
		Size visit_nr = new_nr_++;
		// case1 : both visited--> replace all higher visit_nr by lower visit_nr
		if ( visited_[ up_loop ] && visited_[ down_loop ] ) {
			Size old_visit_nr = visited_[ down_loop ]; //arbitrary choice
			visit_nr = visited_[ up_loop ];
			for ( Size i=1; i<=visited_.size(); i++ ) {
				if ( visited_[ i ] == old_visit_nr ) visited_[ i ] = visit_nr;
			}
		} else if ( visited_[ up_loop ] || visited_[ down_loop ]) {
			// case2: one already visited the other is still zero and thus neutral to addition
			visit_nr = visited_[ up_loop ] + visited_[ down_loop ];
		} // case3: none visited
		visited_[ up_loop ] = visit_nr;
		visited_[ down_loop ] = visit_nr;
		return true;
	} // jump between different regions
	return false;
}

///@detail generate a list of Jumps (Size tupels) that fix the remaining part of the chunk
void
RigidChunkClaimer::JumpCalculator::generate_rigidity_jumps( RigidChunkClaimer* parent_claimer, claims::DofClaims& extra_jumps, std::string label ) {
	if( visited_.size() == 0 ){ // No rigid chunks ??
		return;
	}

	//now we have a connection pattern based on the jumps already present.
	//take a visited region and make it the root-reg
	Size root_reg = 0;
	for ( Size region = 1; region <= visited_.size(); region++ ) {
		if ( visited_[ region ] ) {
			root_reg = region;
			break;
		}
	}

	// if no rigid regions are yet connected, define one arbitrarily as the root-reg
	if ( root_reg == 0 ) {
		root_reg = 1;
		runtime_assert( visited_.size() > 0 );
		visited_[ root_reg ] = 1;
	}

	loops::Loops::LoopList rigid_loops = rigid_.loops(); // loops in sequence that correspond to the regions

	//	take middle of this loop piece. ... there might be better ways to make the extra jumps...
	Size const anchor( static_cast< Size >( 0.5*(rigid_loops[ root_reg ].stop()
				- rigid_loops[ root_reg ].start()) ) + rigid_loops[ root_reg ].start() );

	for ( Size region = 1; region <= visited_.size(); region++ ) {
		Size old_visited = visited_[ region ];
		if ( visited_[ region ] != visited_[ root_reg ] ) {
			Size target_pos (	rigid_loops[ region ].start()
							 + static_cast< Size >( 0.5*( rigid_loops[ region ].stop()-rigid_loops[ region ].start() ) ) );

			extra_jumps.push_back( new claims::JumpClaim( parent_claimer,
																										std::make_pair( label, anchor),
																										std::make_pair( label, target_pos),
																										claims::DofClaim::EXCLUSIVE ) );
			visited_[ region ] = visited_[ root_reg ];

			if ( old_visited ) { //if we connected a cluster make sure to update all its nodes
				for ( Size i=1; i<=visited_.size(); i++ ) {
					if ( visited_[ i ] == old_visited ) visited_[ i ] = visited_[ root_reg ];
				}
			}
		}
	} // for region
} //generate_rigidity_jumps

} //topology_broker
} //protocols
