// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/BuildManager.hh
/// @brief a container for managing BuildInstructions
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/build/BuildManager.hh>

// package headers
#include <protocols/forge/build/Interval.hh>
#include <protocols/forge/build/BuildInstruction.hh>
#include <protocols/forge/methods/util.hh>
#include <protocols/forge/methods/fold_tree_functions.hh>

// project headers

#include <core/kinematics/MoveMap.hh>
//testing
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/util.hh>

#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/option.hh>


#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/id/SequenceMapping.hh>
#include <basic/Tracer.hh>
#ifdef WIN32
// apparently this is required for a Visual Studio build.
#include <core/conformation/Residue.hh>
#endif

// utility headers
#include <utility/exit.hh>

// C++ headers
#include <algorithm>

#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace build {


// Tracer instance for this file
// Named after the original location of this code
static THREAD_LOCAL basic::Tracer TR( "protocols.forge.build.BuildManager" );


/// @brief default constructor
BuildManager::BuildManager() :
	Super(),
	modify_was_successful_( false )
{}


/// @brief copy constructor
BuildManager::BuildManager( BuildManager const & rval ) :
	Super( rval ),
	original2modified_( rval.original2modified_ ),
	seqmap_( rval.seqmap_.get() ? new SequenceMapping( *rval.seqmap_ ) : 0 ),
	modify_was_successful_( rval.modify_was_successful_ )
{
	// add instructions in the exact same order
	for ( BIOPConstIterator i = rval.begin(), ie = rval.end(); i != ie; ++i ) {
		add( *i );
	}

	// rebuild the dependencies
	reconstruct_dependencies( rval.instruction_dependencies_ );
}


/// @brief default destructor
BuildManager::~BuildManager() {}


/// @brief copy assignment
BuildManager & BuildManager::operator =( BuildManager const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );

		clear(); // clear all instructions

		// add the instructions in the same order
		for ( BIOPConstIterator i = rval.begin(), ie = rval.end(); i != ie; ++i ) {
			add( *i );
		}

		// rebuild the dependencies
		reconstruct_dependencies( rval.instruction_dependencies_ );

		original2modified_ = rval.original2modified_;
		seqmap_ = SequenceMappingOP( rval.seqmap_.get() ? new SequenceMapping( *rval.seqmap_ ) : 0 );
		modify_was_successful_ = rval.modify_was_successful_;
	}
	return *this;
}


/// @brief clone this object
BuildManagerOP BuildManager::clone() const {
	return BuildManagerOP( new BuildManager( *this ) );
}


/// @brief create a new instance of this type of object
BuildManagerOP BuildManager::create() const {
	return BuildManagerOP( new BuildManager() );
}


/// @brief reset all accounting info (intervals, positions, etc) to initial
///  state
void BuildManager::reset_accounting() {
	for ( BIOPIterator i = instructions_.begin(), ie = instructions_.end(); i != ie; ++i ) {
		(**i).reset_accounting();
	}
	original2modified_.clear();
	seqmap_.reset(); // to NULL
	modify_was_successful_ = false;
}


/// @brief add an instruction directly (no copy)
void BuildManager::add( BuildInstructionOP bi ) {
	// A different option is to clone the instruction, which is how it was
	// originally done prior to the advent of the dependency tracking.  Storing
	// the pointer directly made the initial code for letting the user
	// construct the dependency graph (BuildManager::create_directed_dependency)
	// somewhat easier.  Consider changing this and how the graph is
	// tracked/constructed if need be.
	instructions_.push_back( bi );
}


/// @brief clear all instructions
void BuildManager::clear() {
	instructions_.clear();
	instruction_dependencies_.clear();
}


/// @brief create a directed dependency: instruction 'u' must complete
///  before instruction 'v' can complete, i.e. 'v' depends on 'u'
void BuildManager::create_directed_dependency( BuildInstructionOP u, BuildInstructionOP v ) {
	// find the index of the instructions in the instructions_ array
	BIOPIterator u_iter = find_instruction( u );
	BIOPIterator v_iter = find_instruction( v );
	runtime_assert( u_iter != instructions_.end() );
	runtime_assert( v_iter != instructions_.end() );

	// link the two instructions, safe even if two instructions already linked
	v->add_dependency_to( BuildInstructionAP( u ) );

	// add to the edge list
	instruction_dependencies_.push_back(
		std::make_pair(
		std::distance( instructions_.begin(), u_iter ) + 1,
		std::distance( instructions_.begin(), v_iter ) + 1
		)
	);
}


/// @brief clear all dependencies
/// @return number of dependencies dropped
BuildManager::Size BuildManager::clear_dependencies() {
	Size const n = instruction_dependencies_.size();

	// clear for each instruction
	for ( BIOPIterator i = instructions_.begin(), ie = instructions_.end(); i != ie; ++i ) {
		(**i).clear_dependencies();
	}

	// clear the list
	instruction_dependencies_.clear();

	return n;
}


/// @brief modify the pose using the instructions in this container
/// @param[in,out] pose the Pose to modify
/// @return a map translating original residue -> modified residue for
///  positions that existed within both the original Pose and modified Pose
BuildManager::Original2Modified BuildManager::modify( Pose & pose ) {
	using namespace protocols::forge::build::BuildInstructionState;
	using protocols::forge::methods::closed_range;

	//testing
	using core::chemical::ResidueTypeCOPs;
	using core::kinematics::Edge;
	using core::kinematics::FoldTree;
	using namespace core::conformation;
	using core::conformation::get_anchor_atomno;

	// sanity check
	if ( !compatibility_check() ) {
		// full stop
		utility_exit_with_message( "ERROR: BuildManager::modify() : BuildInstructions incompatible!" );
	}
	//special case for two chain build -- danger, primitive at the moment, only
	//works with non-N,C term extension.  Try placing in SegRebuld
	if ( basic::options::option[basic::options::OptionKeys::remodel::two_chain_tree].user() ) {
		Size second_start = basic::options::option[basic::options::OptionKeys::remodel::two_chain_tree];
		Size nres( pose.size());

		//FoldTree f(pose.fold_tree());
		FoldTree f;
		//make cutpoint
		f.add_edge(1, second_start-1, Edge::PEPTIDE);
		f.add_edge(second_start, nres, Edge::PEPTIDE);
		f.add_edge(second_start-1,second_start,1);//jump across the cut
		//make cutpoint
		//f.cut_edge(second_start-1);
		f.reorder(nres);
		pose.fold_tree(f);
		//using add_variant_type_to_pose_residue;
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, second_start-1);
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, second_start);
		pose.conformation().declare_chemical_bond( second_start-1, pose.residue( second_start-1 ).atom_name( pose.residue( second_start-1 ).upper_connect_atom() ),
																							 second_start, pose.residue( second_start ).atom_name( pose.residue( second_start ).lower_connect_atom() ) );


		protocols::loops::Loops chain_def_loops;
		chain_def_loops.add_loop(protocols::loops::Loop(second_start-1,second_start, second_start-1));
		//chain_def_loops.add_loop(protocols::loops::Loop(second_start,pose.size()));

		//add virtual residue, as star foldtree requires it
		if ( pose.residue( pose.size()).aa() != core::chemical::aa_vrt ) {
			pose.append_residue_by_jump(*core::conformation::ResidueFactory::create_residue( pose.residue(1).residue_type_set()->name_map("VRT")), second_start);
		}

		//update foldtree to new foldtree
		f = pose.fold_tree();
		TR << "before star" << f << std::endl;

		//update nres to include the new residue
		nres =  pose.size();

		f.reorder(nres);
		pose.fold_tree(f);
		protocols::forge::methods::make_star_foldtree(pose, chain_def_loops);
		TR << "after star" << pose.fold_tree() << std::endl;
	}


	// each call to the Manager's modify() is for a fresh Pose, so reset
	// the accounting information
	reset_accounting();

	// store old nres for constructing SequenceMapping
	Size const old_nres = pose.size();

	// gather residues of original Pose for mapping purposes
	Positions old_r = closed_range( Size( 1 ), pose.size() ); // old residues

	// attach all instructions as length observers and set flag
	// so detach does not occur after modify
	for ( BIOPIterator i = instructions_.begin(), ie = instructions_.end(); i != ie; ++i ) {
		(**i).attach_to( pose );
		(**i).detach_after_modify( false );
	}

	// Now perform all instructions; note that modify() will automatically
	// re-attach here but that's fine.
	Size instructions_remaining = instructions_.size();
	Size instructions_remaining_prior_loop = instructions_remaining; // safeguard against cycles in the dependency topology
	BIOPIterator iii = instructions_.begin(), iiie = instructions_.end();

	while ( instructions_remaining > 0 ) {

		// reset iterator if at the end
		if ( iii == iiie ) {
			iii = instructions_.begin();

			// For now we do the following loop counter check in lieu of actually
			// checking the topology of the dependency graph for cycles.  If the
			// number of instructions remaining in the current loop is equal to the
			// number of instructions remaining in the prior loop, there is a problem.
			if ( instructions_remaining == instructions_remaining_prior_loop ) {
				TR.Fatal << "FATAL: Infinite loop will occur -- check the topology of the dependency graph and make sure there are no cycles." << std::endl;
				utility_exit();
			}

			instructions_remaining_prior_loop = instructions_remaining;
		}

		BuildInstruction & instruction = (**iii);

		// if instruction has not run yet, call modify()
		if ( !instruction.modify_was_successful() ) {
			instruction.modify( pose );
			TR << pose.fold_tree() << std::endl;

			// if modify finished successfully, decrement the counter
			if ( instruction.modify_was_successful() ) {
				--instructions_remaining;
			}
		}

		// increment iterator
		++iii;
	}

	// Manually unhook length observers.  Also reset detach_after_modify flag
	// to true as a safety, since we are not necessarily storing clones of the
	// instructions and we don't want to be inadvertantly causing weird behavior
	// downstream if the instructions are somehow reused.
	for ( BIOPIterator i = instructions_.begin(), ie = instructions_.end(); i != ie; ++i ) {
		(**i).detach_from();
		(**i).detach_after_modify( true ); // safety
	}

	// gather residues of newly modified Pose for mapping purposes
	Positions new_r = closed_range( Size( 1 ), pose.size() ); // new residues

	// remove original deleted positions from original residues
	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = **i;
		Positions const odp = bi.original_deleted_positions();

		for ( Positions::const_iterator j = odp.begin(), je = odp.end(); j != je; ++j ) {
			old_r.erase( *j );
		}
	}

	// remove newly created positions from modified residues
	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = **i;
		Positions const np = bi.new_positions();

		for ( Positions::const_iterator j = np.begin(), je = np.end(); j != je; ++j ) {
			new_r.erase( *j );
		}
	}

	runtime_assert( old_r.size() == new_r.size() );

	// create old -> new map for unmodified regions
	original2modified_.clear();
	for ( Positions::const_iterator o = old_r.begin(), n = new_r.begin(), oe = old_r.end(), ne = new_r.end();
			o != oe && n != ne; ++o, ++n ) {
		original2modified_[ *o ] = *n;
	}
	runtime_assert( original2modified_.size() == old_r.size() && original2modified_.size() == new_r.size() ); // paranoia

	// flip switch so position mapping functions are turned on
	modify_was_successful_ = true;

	/* ALL CODE AFTER THIS POINT MUST OCCUR AFTER modify_was_successful_ IS SET TRUE */

	// create the new SequenceMapping consisting of oldnew() plus
	// old2new_region_endpoints() information
	seqmap_ = SequenceMappingOP( new SequenceMapping() );
	for ( Size r = 1; r <= old_nres; ++r ) {
		Original2Modified::const_iterator i = original2modified_.find( r );
		if ( i != original2modified_.end() ) {
			seqmap_->push_back( i->second );
		} else {
			seqmap_->push_back( 0 );
		}
	}

	Original2Modified o2n_re = original2modified_interval_endpoints();
	for ( Original2Modified::const_iterator i = o2n_re.begin(), ie = o2n_re.end(); i != ie; ++i ) {
		(*seqmap_)[ i->first ] = i->second;
	}

	return original2modified_;
}


/// @brief a dry run of modify() with an all-ala helical Pose of the given length
/// @param[in] nres The length of the dummy structure to use.
/// @return The final length of the modified Pose.
/// @remarks Use this to do a fake run of modify() if you need any
///  position or mapping information prior to actually calling modify().
BuildManager::Size BuildManager::dummy_modify( Size const nres ) {
	using core::pose::make_pose_from_sequence;

	runtime_assert( !instructions_.empty() );

	// fake a poly-ala helical dummy pose
	Pose dummy_pose;
	core::pose::make_pose_from_sequence(
		dummy_pose,
		String( nres, 'A' ),
		( **instructions_.begin() ).residue_type_set()
	);

	for ( core::Size i = 1, ie = dummy_pose.size(); i <= ie; ++i ) {
		dummy_pose.set_secstruct( i, 'H' );
	}

	for ( Size i = 1, ie = dummy_pose.size(); i <= ie; ++i ) {
		dummy_pose.set_phi( i, -60.0 );
		dummy_pose.set_psi( i, -45.0 );
		dummy_pose.set_omega( i, 180.0 );
	}

	modify( dummy_pose );

	return dummy_pose.size();
}


/// @brief check if instruction regions are compatible with each other
/// @return true if regions compatible, false if regions incompatible
bool BuildManager::compatibility_check() const {
	for ( BIOPConstIterator i = begin(), ie = --end(); i != ie; ++i ) {
		for ( BIOPConstIterator j = i + 1, je = end(); j != je; ++j ) {
			if ( !(**i).compatible_with( **j ) ) {
				return false;
			}
		} // foreach j
	} // foreach i

	return true;
}


/// @brief return the combined movemap from all instructions in this manager
BuildManager::MoveMap BuildManager::movemap() const {
	return *movemap_as_OP();
}

BuildManager::MoveMapOP BuildManager::movemap_as_OP() const {
	MoveMapOP combined_mm( new MoveMap() );

	if ( !modify_was_successful_ ) {
		return combined_mm;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		combined_mm->import( (**i).movemap() );
	}

	return combined_mm;
}


/// @brief SequenceMapping consistent with the original -> modified mapping from the
///  most recent modify() call
/// @return valid Sequence mapping if modify() was called; otherwise returns
///  NULL
/// @remarks This mapping contains the same information as original2modified()
///  combined with original2modified_interval_endpoints().
BuildManager::SequenceMappingCOP BuildManager::sequence_mapping() const {
	return seqmap_;
}


/// @brief return a map translating original residue -> modified residue for
///  positions that existed within both the original Pose and modified Pose
/// @return map; empty if modify() has not yet been called
BuildManager::Original2Modified const & BuildManager::original2modified() const {
	return original2modified_;
}


/// @brief return a map translating original intervals to modified intervals
/// @remarks modified intervals with no equivalent original intervals (e.g. cases
///  such as insertions) will not appear in this map
/// @return map; empty if modify() has not yet been called
BuildManager::Interval2Interval BuildManager::original2modified_intervals() const {
	Interval2Interval o2m;

	// if modify() not called, the mapping is unknown
	if ( !modify_was_successful_ ) {
		return o2m;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = **i;

		if ( bi.original_interval_valid() ) {
			o2m[ bi.original_interval() ] = bi.interval();
		}
	}

	return o2m;
}


/// @brief return a map translating original interval endpoints to modified
///  interval endpoints
/// @remarks modified intervals with no equivalent original interval (e.g. cases
///  such as insertions) will not appear in this map
/// @return map; empty if modify() has not yet been called
BuildManager::Original2Modified BuildManager::original2modified_interval_endpoints() const {
	Original2Modified o2m;

	// if modify() not called, the mapping is unknown
	if ( !modify_was_successful_ ) {
		return o2m;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = **i;

		if ( bi.original_interval_valid() ) {
			o2m[ bi.original_interval().left ] = bi.interval().left;
			o2m[ bi.original_interval().right ] = bi.interval().right;
		}
	}

	return o2m;
}


/// @brief return a map translating modified intervals to original intervals
/// @remarks modified intervals with no equivalent original interval (e.g. cases
///  such as insertions) will not appear in this map
/// @return map; empty if modify() has not yet been called
BuildManager::Interval2Interval BuildManager::modified2original_intervals() const {
	Interval2Interval m2o;

	// if modify() not called, the mapping is unknown
	if ( !modify_was_successful_ ) {
		return m2o;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = **i;

		if ( bi.original_interval_valid() ) {
			m2o[ bi.interval() ] = bi.original_interval();
		}
	}

	return m2o;
}


/// @brief return a map translating modified interval endpoints to original
///  interval endpoints
/// @remarks modified intervals with no equivalent original interval (e.g. cases
///  such as insertions) will not appear in this map
BuildManager::Modified2Original BuildManager::modified2original_interval_endpoints() const {
	Modified2Original m2o;

	// if modify() not called, the mapping is unknown
	if ( !modify_was_successful_ ) {
		return m2o;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = **i;

		if ( bi.original_interval_valid() ) {
			m2o[ bi.interval().left ] = bi.original_interval().left;
			m2o[ bi.interval().right ] = bi.original_interval().right;
		}
	}

	return m2o;
}


/// @brief return all modified intervals
/// @remarks Since this encompasses everything this is typically not useful
///  except for overall tracking purposes.
/// @return If modify() has not been called will return an empty set.
std::set< Interval > BuildManager::intervals() const {
	std::set< Interval > intervals;

	if ( !modify_was_successful_ ) {
		return intervals;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = (**i);
		intervals.insert( bi.interval() );
	}

	return intervals;
}


/// @brief return modified intervals that have no equivalent original interval
///  in their BuildInstructions (original_interval_valid() = false)
/// @remarks This is for cases such as insertions where there
///  is no equivalent original region.
/// @return If modify() has not been called will return an empty set.
std::set< Interval > BuildManager::intervals_without_valid_original_equivalents() const {
	std::set< Interval > intervals;

	if ( !modify_was_successful_ ) {
		return intervals;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = (**i);

		// only return those modified intervals whose original interval is
		// invalid
		if ( !bi.original_interval_valid() ) {
			intervals.insert( bi.interval() );
		}
	}

	return intervals;
}


/// @brief return all intervals containing positions that were pre-existing
///  in the original Pose prior to calling modify()
/// @return If modify() has not been called will return an empty set.
std::set< Interval > BuildManager::intervals_containing_preexisting_positions() const {
	std::set< Interval > intervals;

	if ( !modify_was_successful_ ) {
		return intervals;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = (**i);
		if ( !bi.preexisting_positions().empty() ) {
			intervals.insert( bi.interval() );
		}
	}

	return intervals;
}


/// @brief return all intervals containing positions that are "new" and did
///  not exist in the original Pose
/// @return If modify() has not been called will return an empty set.
std::set< Interval > BuildManager::intervals_containing_new_positions() const {
	std::set< Interval > intervals;

	if ( !modify_was_successful_ ) {
		return intervals;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = (**i);
		if ( !bi.new_positions().empty() ) {
			intervals.insert( bi.interval() );
		}
	}

	return intervals;
}


/// @brief return all intervals containing positions with defined conformation
/// @return If modify() has not been called will return an empty set.
std::set< Interval > BuildManager::intervals_containing_defined_positions() const {
	std::set< Interval > intervals;

	if ( !modify_was_successful_ ) {
		return intervals;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = (**i);
		if ( !bi.defined_positions().empty() ) {
			intervals.insert( bi.interval() );
		}
	}

	return intervals;
}


/// @brief return all intervals containing positions with undefined conformation
/// @remarks typically used to define intervals appropriate for loop modeling
/// @return If modify() has not been called will return an empty set.
std::set< Interval > BuildManager::intervals_containing_undefined_positions() const {
	std::set< Interval > intervals;

	if ( !modify_was_successful_ ) {
		return intervals;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = (**i);
		if ( !bi.undefined_positions().empty() ) {
			intervals.insert( bi.interval() );
		}
	}

	return intervals;
}


/// @brief return all original intervals containing positions that will
///  be kept by the BuildInstructions
/// @remarks returns valid data even without calling modify()
std::set< Interval > BuildManager::original_intervals_containing_kept_positions() const {
	std::set< Interval > intervals;

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = (**i);
		if ( !bi.original_kept_positions().empty() ) {
			intervals.insert( bi.original_interval() );
		}
	}

	return intervals;
}


/// @brief return all original intervals containing positions that will
///  be deleted by the BuildInstructions
/// @remarks returns valid data even without calling modify()
std::set< Interval > BuildManager::original_intervals_containing_deleted_positions() const {
	std::set< Interval > intervals;

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = (**i);
		if ( !bi.original_deleted_positions().empty() ) {
			intervals.insert( bi.original_interval() );
		}
	}

	return intervals;
}


/// @brief return all positions within the modified intervals
/// @remarks Since this encompasses everything this is typically not useful
///  except for overall tracking purposes.
/// @return If modify() has not been called will return an empty set.
BuildManager::Positions BuildManager::positions() const {
	using protocols::forge::methods::insert_closed_range;

	Positions p;

	if ( !modify_was_successful_ ) {
		return p;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		Interval const ival = (**i).interval();
		insert_closed_range( ival.left, ival.right, p );
	}

	return p;
}


/// @brief return the set of positions within the new regions that were
///  pre-existing in the original Pose prior to calling modify()
/// @return If modify() has not been called will return an empty set.
BuildManager::Positions BuildManager::preexisting_positions() const {
	Positions preexisting;

	if ( !modify_was_successful_ ) {
		return preexisting;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		Positions const bi_pp = (**i).preexisting_positions();
		preexisting.insert( bi_pp.begin(), bi_pp.end() );
	}

	return preexisting;
}


/// @brief return a copy of the set of positions that are "new" and did
///  not exist in the original Pose.
/// @return If modify() has not been called will return an empty set.
BuildManager::Positions BuildManager::new_positions() const {
	Positions newp;

	if ( !modify_was_successful_ ) {
		return newp;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		Positions const bi_np = (**i).new_positions();
		newp.insert( bi_np.begin(), bi_np.end() );
	}

	return newp;
}


/// @brief return a copy of the set of positions within the newly modified
///  regions that have a defined conformation.  E.g. existing or copied residues.
/// @return If modify() has not been called will return an empty set.
BuildManager::Positions BuildManager::defined_positions() const {
	Positions defined;

	if ( !modify_was_successful_ ) {
		return defined;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		Positions const bi_dp = (**i).defined_positions();
		defined.insert( bi_dp.begin(), bi_dp.end() );
	}

	return defined;
}


/// @brief return a copy of the set of positions within the newly modified
///  regions that have an undefined conformation.  E.g. newly created residues.
/// @return If modify() has not been called will return an empty set.
BuildManager::Positions BuildManager::undefined_positions() const {
	Positions undefined;

	if ( !modify_was_successful_ ) {
		return undefined;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		Positions const bi_up = (**i).undefined_positions();
		undefined.insert( bi_up.begin(), bi_up.end() );
	}

	return undefined;
}


/// @brief the positions representing the union of all intervals containing
///  positions with undefined conformation
/// @remarks Useful as a reference for defining neighborhoods around
///  loop modeled regions.
/// @return If modify() has not been called will return an empty set.
BuildManager::Positions BuildManager::union_of_intervals_containing_undefined_positions() const {
	using protocols::forge::methods::insert_closed_range;

	Positions p;

	if ( !modify_was_successful_ ) {
		return p;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = (**i);
		Interval const ival = bi.interval();

		if ( !bi.undefined_positions().empty() ) {
			insert_closed_range( ival.left, ival.right, p );
		}
	}

	return p;
}


/// @brief return the set of positions within the original intervals that
///  will be kept by the BuildInstructions
/// @remarks returns valid data even without calling modify()
BuildManager::Positions BuildManager::original_kept_positions() const {
	Positions kept;

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		Positions const bi_kp = (**i).original_kept_positions();
		kept.insert( bi_kp.begin(), bi_kp.end() );
	}

	return kept;
}


/// @brief return set of positions within the original intervals that will
///  be deleted by the BuildInstructions
/// @remarks returns valid data even without calling modify()
BuildManager::Positions BuildManager::original_deleted_positions() const {
	Positions deleted;

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		Positions const bi_dp = (**i).original_deleted_positions();
		deleted.insert( bi_dp.begin(), bi_dp.end() );
	}

	return deleted;
}


/// @brief return a map from modified intervals to the set of pre-existing
///  positions inside them
/// @return If modify() has not been called will return an empty map.
BuildManager::Interval2Positions BuildManager::modified_i2p_preexisting() const {
	Interval2Positions i2p;

	if ( !modify_was_successful_ ) {
		return i2p;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = **i;
		i2p.insert( std::make_pair( bi.interval(), bi.preexisting_positions() ) );
	}

	return i2p;
}


/// @brief return a map from modified intervals to the set of "new"
///  positions inside them that were not present in the original Pose
/// @return If modify() has not been called will return an empty map.
BuildManager::Interval2Positions BuildManager::modified_i2p_new() const {
	Interval2Positions i2p;

	if ( !modify_was_successful_ ) {
		return i2p;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = **i;
		i2p.insert( std::make_pair( bi.interval(), bi.new_positions() ) );
	}

	return i2p;
}


/// @brief return a map from modified intervals to the set of
///  positions inside them that have defined conformation
/// @return If modify() has not been called will return an empty map.
BuildManager::Interval2Positions BuildManager::modified_i2p_defined() const {
	Interval2Positions i2p;

	if ( !modify_was_successful_ ) {
		return i2p;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = **i;
		i2p.insert( std::make_pair( bi.interval(), bi.defined_positions() ) );
	}

	return i2p;
}


/// @brief return a map from modified intervals to the set of
///  positions inside them that have undefined conformation
/// @return If modify() has not been called will return an empty map.
BuildManager::Interval2Positions BuildManager::modified_i2p_undefined() const {
	Interval2Positions i2p;

	if ( !modify_was_successful_ ) {
		return i2p;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = **i;
		i2p.insert( std::make_pair( bi.interval(), bi.undefined_positions() ) );
	}

	return i2p;
}


/// @brief return a map from modified intervals to their individual movemaps
/// @return If modify() has not been called will return an empty map.
BuildManager::Interval2MoveMap BuildManager::modified_interval2movemap() const {
	Interval2MoveMap i2m;

	if ( !modify_was_successful_ ) {
		return i2m;
	}

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = **i;
		i2m.insert( std::make_pair( bi.interval(), bi.movemap() ) );
	}

	return i2m;
}


/// @brief return a map from original intervals to the set of positions
///  inside them that will be kept by the BuildInstructions
/// @remarks returns valid data even without calling modify()
BuildManager::Interval2Positions BuildManager::original_i2p_kept() const {
	Interval2Positions i2p;

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = **i;
		i2p.insert( std::make_pair( bi.original_interval(), bi.original_kept_positions() ) );
	}

	return i2p;
}


/// @brief return a map from original intervals to the set of positions
///  inside them that will be deleted by the BuildInstructions
/// @remarks returns valid data even without calling modify()
BuildManager::Interval2Positions BuildManager::original_i2p_deleted() const {
	Interval2Positions i2p;

	for ( BIOPConstIterator i = begin(), ie = end(); i != ie; ++i ) {
		BuildInstruction const & bi = **i;
		i2p.insert( std::make_pair( bi.original_interval(), bi.original_deleted_positions() ) );
	}

	return i2p;
}


/// @brief find the given instruction
/// @return iterator pointing to the BuildInstructionOP if found, otherwise
///  the 'end' iterator
BuildManager::BIOPIterator BuildManager::find_instruction( BuildInstructionCOP u ) {
	return std::find( instructions_.begin(), instructions_.end(), u );
}


/// @brief find the given instruction
/// @return iterator pointing to the BuildInstructionCOP if found, otherwise
///  the 'end' iterator
BuildManager::BIOPConstIterator BuildManager::find_instruction( BuildInstructionCOP u ) const {
	return std::find( instructions_.begin(), instructions_.end(), u );
}


/// @brief find the edge specifying the given dependency
/// @return iterator pointing to the DependencyEdge if found, otherwise the
///  'end' iterator
BuildManager::DependencyEdges::iterator BuildManager::find_dependency( BuildInstructionCOP u, BuildInstructionCOP v ) {
	for ( DependencyEdges::iterator i = instruction_dependencies_.begin(), ie = instruction_dependencies_.end(); i != ie; ++i ) {
		if ( instructions_[ i->first ] == u && instructions_[ i->second ] == v ) {
			return i;
		}
	}

	return instruction_dependencies_.end();
}


/// @brief find the given dependency
/// @return const iterator pointing to the DependencyEdge if found, otherwise the
///  'end' const iterator
BuildManager::DependencyEdges::const_iterator BuildManager::find_dependency( BuildInstructionCOP u, BuildInstructionCOP v ) const {
	for ( DependencyEdges::const_iterator i = instruction_dependencies_.begin(), ie = instruction_dependencies_.end(); i != ie; ++i ) {
		if ( instructions_[ i->first ] == u && instructions_[ i->second ] == v ) {
			return i;
		}
	}

	return instruction_dependencies_.end();
}


/// @brief clear the current dependency list and reconstruct the dependencies
///  using the given list
void BuildManager::reconstruct_dependencies( DependencyEdges const & dependency_list ) {
	clear_dependencies();

	// add the new dependencies
	for ( DependencyEdges::const_iterator i = dependency_list.begin(), ie = dependency_list.end(); i != ie; ++i ) {
		create_directed_dependency( instructions_[ i->first ], instructions_[ i->second ] );
	}
}


} // namespace build
} // namespace forge
} // namespace protocols
