// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/forge/build/BuildManager.hh
/// @brief a container for managing BuildInstructions
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_BuildManager_hh
#define INCLUDED_protocols_forge_build_BuildManager_hh

// unit headers
#include <protocols/forge/build/BuildManager.fwd.hh>

// package headers
#include <core/types.hh>
#include <protocols/forge/build/Interval.fwd.hh>
#include <protocols/forge/build/BuildInstruction.fwd.hh>

#if defined(WIN_PYROSETTA) || defined(WIN32)
#include <protocols/forge/build/BuildInstruction.hh>
#endif

// project headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>
#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief a container for managing BuildInstructions
/// @note Compatibility checks wrt dependencies currently do not exist.
///  It remains to be seen how to handle this.
class BuildManager : public utility::pointer::ReferenceCount {


private: // typedefs


	typedef utility::pointer::ReferenceCount Super;


public: // typedefs


	typedef core::Size Size;
	typedef core::kinematics::MoveMap MoveMap;
	typedef core::kinematics::MoveMapOP MoveMapOP;
	typedef core::pose::Pose Pose;

	typedef core::id::SequenceMapping SequenceMapping;
	typedef core::id::SequenceMappingOP SequenceMappingOP;
	typedef core::id::SequenceMappingCOP SequenceMappingCOP;

	typedef utility::vector1< BuildInstructionOP > BuildInstructionOPs;
	typedef BuildInstructionOPs::iterator BIOPIterator;
	typedef BuildInstructionOPs::const_iterator BIOPConstIterator;
	typedef std::pair< Size, Size > DependencyEdge;
	typedef utility::vector1< DependencyEdge > DependencyEdges;
	typedef std::set< Size > Positions;
	typedef std::map< Size, Size > Original2Modified;
	typedef std::map< Size, Size > Modified2Original;
	typedef std::map< Interval, Positions > Interval2Positions;
	typedef std::map< Interval, Interval > Interval2Interval;
	typedef std::map< Interval, MoveMap > Interval2MoveMap;

	typedef std::string String;


public: // construct/destruct


	/// @brief default constructor
	BuildManager();


	/// @brief copy constructor
	BuildManager( BuildManager const & rval );


	/// @brief default destructor
	virtual
	~BuildManager();


public: // assignment


	/// @brief copy assignment
	BuildManager & operator =( BuildManager const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	BuildManagerOP clone() const;


	/// @brief create a new instance of this type of object
	virtual
	BuildManagerOP create() const;


public: // mutators


	/// @brief reset all accounting info (intervals, positions, etc) to initial
	///  state
	void reset_accounting();


public: // instruction interface


	/// @brief add an instruction directly (no copy)
	void add( BuildInstructionOP bi );


	/// @brief clear all instructions
	void clear();


	/// @brief current number of instructions
	inline
	Size size() const {
		return instructions_.size();
	}


	/// @brief no instructions?
	inline
	bool empty() const {
		return instructions_.empty();
	}


	/// @brief const iterator pointing to the first instruction
	inline
	BIOPConstIterator begin() const {
		return instructions_.begin();
	}


	/// @brief const iterator pointing just beyond the last instruction
	inline
	BIOPConstIterator end() const {
		return instructions_.end();
	}


public: // instruction dependency interface


	/// @brief create a directed dependency: instruction 'u' must complete
	///  before instruction 'v' can complete, i.e. 'v' depends on 'u'
	void create_directed_dependency( BuildInstructionOP u, BuildInstructionOP v );


	/// @brief the number of dependencies currently defined (i.e. # of edges in
	///  the dependency graph)
	inline
	Size n_dependencies() const {
		return instruction_dependencies_.size();
	}


	/// @brief have dependencies been defined?
	inline
	bool dependencies_exist() const {
		return !instruction_dependencies_.empty();
	}


	/// @brief clear all dependencies
	/// @return number of dependencies dropped
	Size clear_dependencies();


public: // pose modification interface


	/// @brief modify the pose using the instructions in this container
	/// @param[in,out] pose the Pose to modify
	/// @return a map translating original residue -> modified residue for
	///  positions that existed within both the original Pose and modified Pose
	Original2Modified modify( Pose & pose );


	/// @brief a dry run of modify() with an all-ala helical Pose of the given
	/// length
	/// @param[in] nres The length of the dummy structure to use.
	/// @return The final length of the modified Pose.
	/// @remarks Use this to do a fake run of modify() if you need any
	///  position or mapping information prior to actually calling modify().
	Size dummy_modify( Size const nres );


	/// @brief check if instruction regions are compatible with each other
	/// @return true if regions compatible, false if regions incompatible
	bool compatibility_check() const;


public: // movemap


	/// @brief return the combined movemap from all instructions in this manager
	/// @return If modify() has not been called will return an empty MoveMap.
	MoveMap movemap() const;

	MoveMapOP movemap_as_OP() const;


public: // mapping


	/// @brief SequenceMapping consistent with the original -> modified mapping from the
	///  most recent modify() call
	/// @return valid Sequence mapping if modify() was called; otherwise returns
	///  NULL
	/// @remarks This mapping contains the same information as original2modified()
	///  combined with original2modified_interval_endpoints().
	SequenceMappingCOP sequence_mapping() const;


	/// @brief return a map translating original residue -> modified residue for
	///  positions that existed within both the original Pose and modified Pose
	/// @return map; empty if modify() has not yet been called
	Original2Modified const & original2modified() const;


	/// @brief return a map translating original intervals to modified intervals
	/// @remarks modified intervals with no equivalent original interval (e.g. cases
	///  such as insertions) will not appear in this map
	/// @return map; empty if modify() has not yet been called
	Interval2Interval original2modified_intervals() const;


	/// @brief return a map translating original interval endpoints to modified
	///  interval endpoints
	/// @remarks modified intervals with no equivalent original interval (e.g. cases
	///  such as insertions) will not appear in this map
	/// @return map; empty if modify() has not yet been called
	Original2Modified original2modified_interval_endpoints() const;


	/// @brief return a map translating modified intervals to original intervals
	/// @remarks modified intervals with no equivalent original interval (e.g. cases
	///  such as insertions) will not appear in this map
	/// @return map; empty if modify() has not yet been called
	Interval2Interval modified2original_intervals() const;


	/// @brief return a map translating modified interval endpoints to original
	///  interval endpoints
	/// @remarks modified intervals with no equivalent original interval (e.g. cases
	///  such as insertions) will not appear in this map
	/// @return map; empty if modify() has not yet been called
	Modified2Original modified2original_interval_endpoints() const;


public: // intervals in the newly modified Pose


	/// @brief return all modified intervals
	/// @remarks Since this encompasses everything this is typically not useful
	///  except for overall tracking purposes.
	/// @return If modify() has not been called will return an empty set.
	std::set< Interval > intervals() const;


	/// @brief return modified intervals that have no equivalent original interval
	///  in their BuildInstructions (original_interval_valid() = false)
	/// @remarks This is for cases such as insertions where there
	///  is no equivalent original region.
	/// @return If modify() has not been called will return an empty set.
	std::set< Interval > intervals_without_valid_original_equivalents() const;


	/// @brief return all intervals containing positions that were pre-existing
	///  in the original Pose prior to calling modify()
	/// @return If modify() has not been called will return an empty set.
	std::set< Interval > intervals_containing_preexisting_positions() const;


	/// @brief return all intervals containing positions that are "new" and did
	///  not exist in the original Pose
	/// @return If modify() has not been called will return an empty set.
	std::set< Interval > intervals_containing_new_positions() const;


	/// @brief return all intervals containing positions with defined conformation
	/// @return If modify() has not been called will return an empty set.
	std::set< Interval > intervals_containing_defined_positions() const;


	/// @brief return all intervals containing positions with undefined conformation
	/// @remarks typically used to define intervals appropriate for loop modeling
	/// @return If modify() has not been called will return an empty set.
	std::set< Interval > intervals_containing_undefined_positions() const;


public: // original intervals in the original Pose


	/// @brief return all original intervals containing positions that will
	///  be kept by the BuildInstructions
	/// @remarks returns valid data even without calling modify()
	std::set< Interval > original_intervals_containing_kept_positions() const;


	/// @brief return all original intervals containing positions that will
	///  be deleted by the BuildInstructions
	/// @remarks returns valid data even without calling modify()
	std::set< Interval > original_intervals_containing_deleted_positions() const;


public: // positions in the newly modified Pose


	/// @brief return all positions within the modified intervals
	/// @remarks Since this encompasses everything this is typically not useful
	///  except for overall tracking purposes.
	/// @return If modify() has not been called will return an empty set.
	Positions positions() const;


	/// @brief return the set of positions within the new regions that were
	///  pre-existing in the original Pose prior to calling modify()
	/// @return If modify() has not been called will return an empty set.
	Positions preexisting_positions() const;


	/// @brief return a copy of the set of positions that are "new" and did
	///  not exist in the original Pose.
	/// @return If modify() has not been called will return an empty set.
	Positions new_positions() const;


	/// @brief return a copy of the set of positions within the newly modified
	///  regions that have a defined conformation.  E.g. existing or copied residues.
	/// @return If modify() has not been called will return an empty set.
	Positions defined_positions() const;


	/// @brief return a copy of the set of positions within the newly modified
	///  regions that have an undefined conformation.  E.g. newly created residues.
	/// @return If modify() has not been called will return an empty set.
	Positions undefined_positions() const;


	/// @brief the positions representing the union of all intervals containing
	///  positions with undefined conformation
	/// @remarks Useful as a reference for defining neighborhoods around
	///  loop modeled regions.
	/// @return If modify() has not been called will return an empty set.
	Positions union_of_intervals_containing_undefined_positions() const;


public: // positions in the original unmodified Pose


	/// @brief return the set of positions within the original intervals that
	///  will be kept by the BuildInstructions
	/// @remarks returns valid data even without calling modify()
	Positions original_kept_positions() const;


	/// @brief return set of positions within the original intervals that will
	///  be deleted by the BuildInstructions
	/// @remarks returns valid data even without calling modify()
	Positions original_deleted_positions() const;


public: // map modified intervals to position information


	/// @brief return a map from modified intervals to the set of pre-existing
	///  positions inside them
	/// @return If modify() has not been called will return an empty map.
	Interval2Positions modified_i2p_preexisting() const;


	/// @brief return a map from modified intervals to the set of "new"
	///  positions inside them that were not present in the original Pose
	/// @return If modify() has not been called will return an empty map.
	Interval2Positions modified_i2p_new() const;


	/// @brief return a map from modified intervals to the set of
	///  positions inside them that have defined conformation
	/// @return If modify() has not been called will return an empty map.
	Interval2Positions modified_i2p_defined() const;


	/// @brief return a map from modified intervals to the set of
	///  positions inside them that have undefined conformation
	/// @return If modify() has not been called will return an empty map.
	Interval2Positions modified_i2p_undefined() const;


	/// @brief return a map from modified intervals to their individual movemaps
	/// @return If modify() has not been called will return an empty map.
	Interval2MoveMap modified_interval2movemap() const;


public: // map original intervals to original position information


	/// @brief return a map from original intervals to the set of positions
	///  inside them that will be kept by the BuildInstructions
	/// @remarks returns valid data even without calling modify()
	Interval2Positions original_i2p_kept() const;


	/// @brief return a map from original intervals to the set of positions
	///  inside them that will be deleted by the BuildInstructions
	/// @remarks returns valid data even without calling modify()
	Interval2Positions original_i2p_deleted() const;


private: // instruction handling


	/// @brief find the given instruction
	/// @return iterator pointing to the BuildInstructionOP if found, otherwise
	///  the 'end' iterator
	BIOPIterator find_instruction( BuildInstructionCOP u );


	/// @brief find the given instruction
	/// @return iterator pointing to the BuildInstructionCOP if found, otherwise
	///  the 'end' iterator
	BIOPConstIterator find_instruction( BuildInstructionCOP u ) const;


private: // dependency handling


	/// @brief find the given dependency
	/// @return iterator pointing to the DependencyEdge if found, otherwise the
	///  'end' iterator
	DependencyEdges::iterator find_dependency( BuildInstructionCOP u, BuildInstructionCOP v );


	/// @brief find the given dependency
	/// @return const iterator pointing to the DependencyEdge if found, otherwise the
	///  'end' const iterator
	DependencyEdges::const_iterator find_dependency( BuildInstructionCOP u, BuildInstructionCOP v ) const;


	/// @brief clear the current dependency list and reconstruct the dependencies
	///  using the given list
	void reconstruct_dependencies( DependencyEdges const & dependency_list );


private: // data


	/// @brief the list of BuildInstructions to apply
	BuildInstructionOPs instructions_;


	/// @brief list recording BuildInstruction inter-dependencies by index into
	///  instructions_ array
	/// @remarks For each pair "first,second" in the list, the instruction with
	///  index "first" must complete before the instruction with index "second"
	///  I.e. "second" depends on "first".
	DependencyEdges instruction_dependencies_;


	/// @brief map translation old residue -> new residue for
	///  non-modified regions only
	/// @remarks only filled in after modify() called
	Original2Modified original2modified_;


	/// @brief SequenceMapping consistent with the old -> new mapping from the
	///  most recent modify() call
	SequenceMappingOP seqmap_;


	/// @brief indicates modify() has been called and succeeded
	bool modify_was_successful_;


};


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_BuildManager_HH */
