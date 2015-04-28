// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/components/VarLengthBuild.fwd.hh
/// @brief  Component that performs a simplified version of a protocol for
///         variable length remodeling of protein backbone segments.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_components_VarLengthBuild_hh
#define INCLUDED_protocols_forge_components_VarLengthBuild_hh

// unit headers
#include <protocols/forge/components/VarLengthBuild.fwd.hh>

// package headers
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/Interval.fwd.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.fwd.hh>
#include <protocols/forge/remodel/RemodelData.hh>
#ifdef WIN32
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh> // WIN32 INCLUDE
#endif

// project headers
#include <core/fragment/ConstantLengthFragSet.fwd.hh>
#include <core/fragment/OrderedFragSet.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/Mover.hh>

// utility headers

// C++ headers
#include <string>

#include <core/fragment/FrameList.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace components {


/// @brief Component that performs a protocol for user-specified variable length
///  remodeling of protein backbone segments.
/// @details This bootstrap implementation performs centroid level loop building
///  using 9,3,1-mer fragments and running a LoopMover with the Remodel 'remodel_cen'
///  score function.  The default loop mover is the forge RemodelLoopMover.
///  It does not yet handle extensions or continuous replacements.
///  This class is under heavy development, expect many changes to both API and
///  internals.
class VarLengthBuild : public protocols::moves::Mover {


private: // typedefs


	typedef protocols::moves::Mover Super;


public: // typedefs

	typedef std::string String;

	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::fragment::ConstantLengthFragSetOP ConstantLengthFragSetOP;
	typedef core::fragment::FrameList FrameList;
	typedef core::fragment::OrderedFragSetOP OrderedFragSetOP;
	typedef core::kinematics::MoveMap MoveMap;
	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef protocols::loops::Loops Loops;
	typedef protocols::forge::remodel::RemodelData RemodelData;
	typedef protocols::moves::MoverOP MoverOP;

	typedef protocols::forge::build::BuildManager BuildManager;
	typedef protocols::forge::build::Interval Interval;
	typedef protocols::forge::remodel::RemodelConstraintGeneratorOP RemodelConstraintGeneratorOP;

	typedef BuildManager::BIOPConstIterator BIOPConstIterator;
	typedef BuildManager::Original2Modified Original2Modified;
	typedef BuildManager::Positions Positions;


public: // construct/destruct


	/// @brief default constructor
	VarLengthBuild();


	/// @brief BuildManager constructor
	VarLengthBuild( BuildManager const & manager );

	/// @brief BuildManager + RemodelData constructor
	VarLengthBuild( BuildManager const & manager , RemodelData const & remodel_data); //copy remodeldata


	/// @brief copy constructor
	VarLengthBuild( VarLengthBuild const & rval );


	/// @brief default destructor
	~VarLengthBuild();


private: // assignment


	/// @brief disallow copy assignment
	/// @remarks Mover base class prevents this from working...
	VarLengthBuild & operator =( VarLengthBuild const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	MoverOP clone() const;


	/// @brief create a new instance of this type of object
	virtual
	MoverOP fresh_instance() const;


public: // accessors


	/// @brief build manager
	inline
	BuildManager const & manager() const {
		return manager_;
	}


	/// @brief at the end of apply(), reset the Pose to the original Pose if
	///  mover was not successful?  (default true)
	inline
	bool recover_original_on_failure() const {
		return recover_original_on_failure_;
	}


	/// @brief the string id of the loop mover to use
	///  (default "RemodelLoopMover")
	/// @return "RemodelLoopMover" for the forge RemodelLoopMover, otherwise a
	///  string recognized by create_loop_mover() in the "LoopMoverFactory".
	inline
	String const & loop_mover_str() const {
		return loop_mover_str_;
	}


	/// @brief cache fragments after picking? (default true)
	/// @details If true, will cache fragments and reuse them upon each
	///  apply() call instead of repicking every time.
	inline
	bool cache_fragments() const {
		return cache_fragments_;
	}


	/// @brief option dictating whether to keep VallLibrary in memory or clear
	///  it under certain circumstances after picking fragments.
	///  (default KEEP_IN_MEMORY)
	inline
	VLB_VallMemoryUsage::Enum vall_memory_usage() const {
		return vall_memory_usage_;
	}


	/// @brief the number of fragments to pick at each position
	///  (default 200)
	inline
	Size num_fragpick() const {
		return num_fragpick_;
	}


	/// @brief also use fragments spanning the entire length of a loop?
	///  (default false)
	inline
	bool use_fullmer() const {
		return use_fullmer_;
	}


	/// @brief full sequence string corresponding to original input pose
	///  used to pick sequence biased fragments; if empty, sequence bias
	///  is not used when picking fragments
	inline
	String const & original_sequence() {
		return original_sequence_;
	}


	/// @brief full sequence string with length corresponding to the *new*
	///  modified pose used to pick secondary structure biased fragments.
	/// @remarks IMPORTANT: This is an override!  If this string is filled,
	///  it will be used as the string to pick secondary structure biased
	///  fragments without modification.  All secondary structure settings
	///  that might be taken from the original Pose or from the BuildInstructions
	///  *will be ignored*.  The length of this string must be equal to the
	///  length of the *NEW* modified pose, otherwise the protocol will stop
	///  with error -- you can use BuildManager::dummy_modify() to help
	///  figure things out.
	inline
	String const & new_secondary_structure_override() const {
		return new_secondary_structure_override_;
	}


	/// @brief full amino acid string with length corresponding to the *new*
	///  modified pose used to pick sequence biased fragments.
	/// @remarks IMPORTANT: This is an override!  If this string is filled,
	///  it will be used as the string to pick sequence biased
	///  fragments without modification.  All sequence settings
	///  that might be taken from the original Pose or from the BuildInstructions
	///  *will be ignored*.  The length of this string must be equal to the
	///  length of the *NEW* modified pose, otherwise the protocol will stop
	///  with error -- you can use BuildManager::dummy_modify() to help
	///  figure things out.
	inline
	String const & new_sequence_override() const {
		return new_sequence_override_;
	}


	/// @brief return the highest linear chainbreak score a chainbreak can have
	///  and still be considered closed
	/// @remarks default 0.07
	inline
	Real max_linear_chainbreak() const {
		return max_linear_chainbreak_;
	}


	/// @brief Flag to turn on restart mode, in which VLB assumes that the Pose
	///  fed to it during apply() has already been modified by the manager.
	///  (default False)
	/// @remarks In restart mode, VLB only runs the manager's dummy modify
	///  capability during apply() to get the mapping info.
	inline
	bool restart_mode() const {
		return restart_mode_;
	}


public: // mutators


	/// @brief set ScoreFunction used during build
	void scorefunction( ScoreFunction const & sfx );


	/// @brief set ScoreFunction used during build
	void scorefunction( ScoreFunctionOP const & sfx );


	/// @brief set build manager; also clears any cached fragments
	void manager( BuildManager const & manager );


	/// @brief at the end of apply(), reset the Pose to the original Pose if
	///  mover was not successful?
	inline
	void recover_original_on_failure( bool const flag ) {
		recover_original_on_failure_ = flag;
	}


	/// @brief set the loop mover to use via string
	/// @details use "RemodelLoopMover" for the forge RemodelLoopMover,
	///  otherwise set it to a string recognized by
	///  create_loop_mover() in the "LoopMoverFactory".
	inline
	void loop_mover_str( String const & str ) {
		loop_mover_str_ = str;
	}


	/// @brief cache fragments after picking?
	/// @details If true, will cache fragments and reuse them upon each
	///  apply() call instead of repicking every time.
	inline
	void cache_fragments( bool const flag ) {
		cache_fragments_ = flag;
	}


	/// @brief option dictating whether to keep VallLibrary in memory or clear
	///  it under certain circumstances after picking fragments.
	inline
	void vall_memory_usage( VLB_VallMemoryUsage::Enum const level ) {
		vall_memory_usage_ = level;
	}


	/// @brief the number of fragments to pick at each position (default 200)
	inline
	void num_fragpick( Size const num ){
		num_fragpick_ = num;
	}


	/// @brief also use fragments spanning the entire length of a loop?
	inline
	void use_fullmer( bool const flag ) {
		use_fullmer_ = flag;
	}


	/// @brief full sequence string corresponding to original input pose
	///  used to pick sequence biased fragments; if empty, sequence bias
	///  is not used when picking fragments
	inline
	void original_sequence( String const & seq ) {
		original_sequence_ = seq;
	}


	/// @brief full sequence string with length corresponding to the *new*
	///  modified pose used to pick secondary structure biased fragments.
	/// @param[in] str String with length equals to the *new* modified pose.
	///  See remarks for help on how to determine this.  String is allowed
	///  to be empty, in which case it will clear the setting.
	/// @remarks IMPORTANT: This is an override!  If this string is filled,
	///  it will be used as the string to pick secondary structure biased
	///  fragments without modification.  All secondary structure settings
	///  that might be taken from the original Pose or from the BuildInstructions
	///  *will be ignored*.  The length of this string must be equal to the
	///  length of the *NEW* modified pose, otherwise the protocol will stop
	///  with error -- you can use BuildManager::dummy_modify() to help
	///  figure things out.
	inline
	void new_secondary_structure_override( String const & str ) {
		new_secondary_structure_override_ = str;
	}


	/// @brief full amino acid string with length corresponding to the *new*
	///  modified pose used to pick sequence biased fragments.
	/// @param[in] str String with length equals to the *new* modified pose.
	///  See remarks for help on how to determine this.  String is allowed
	///  to be empty, in which case it will clear the setting.
	/// @remarks IMPORTANT: This is an override!  If this string is filled,
	///  it will be used as the string to pick sequence biased
	///  fragments without modification.  All sequence settings
	///  that might be taken from the original Pose or from the BuildInstructions
	///  *will be ignored*.  The length of this string must be equal to the
	///  length of the *NEW* modified pose, otherwise the protocol will stop
	///  with error -- you can use BuildManager::dummy_modify() to help
	///  figure things out.
	inline
	void new_sequence_override( String const & str ) {
		new_sequence_override_ = str;
	}


	/// @brief set the highest linear chainbreak score a chainbreak can have
	///  and still be considered closed
	inline
	void max_linear_chainbreak( Real const tol ) {
		max_linear_chainbreak_ = tol;
	}

	inline
	void loop_mover_fold_tree_constant( bool const flag ){
		loop_mover_fold_tree_constant_ = flag;
	}

	/// @brief Flag to turn on restart mode, in which VLB assumes that the Pose
	///  fed to it during apply() has already been modified by the manager.
	///  (default False)
	/// @remarks In restart mode, VLB only runs the manager's dummy modify
	///  capability during apply() to get the mapping info.
	inline
	void restart_mode( bool const flag ) {
		restart_mode_ = flag;
	}

	inline
	void ignore_cmdline_enzdes_cstfile( bool const flag ) {
		ignore_cmdline_enzdes_cstfile_ = flag;
	}


public: //constraint / setup mover management management

	void
	clear_rcgs();

	void
	add_rcg( RemodelConstraintGeneratorOP rcg );

	void
	clear_setup_movers();

	void
	add_setup_mover( moves::MoverOP mover_in );

	void
	clear_user_provided_movers();

	void
	add_user_provided_mover( moves::MoverOP mover_in );


public: // fragment management


	/// @brief clear any currently cached fragments
	void clear_fragments();

	/// @brief set abego definition for fragments
	void set_abego( utility::vector1< String > const & abego ) {
		abego_ = abego;
	}


public: // main operations


	/// @brief run protocol on given Pose
	/// @return if procedure successful, return Pose with modifications and a
	///  sealed fold tree, otherwise return Pose with modifications and the
	///  in-progress cut fold tree
	/// @remarks Before invoking this function it's best to make sure
	///  the secondary structure in the Pose is marked via the method
	///  that you would prefer, e.g. by Dssp (protocols::jumping::Dssp),
	///  by the old Rosetta++ binning method (core::pose::set_ss_from_phipsi)
	///  or by external method such as reading in a file.
	virtual
	void apply( Pose & pose );

	virtual std::string get_name() const;

protected: // main operations


	/// @brief run centroid level protocol on given Pose
	/// @return true if regions modeled within tolerances, false otherwise
	virtual
	bool centroid_build(
		Pose & pose
	);


	/// @brief return the appropriate loop mover
	/// @param[in] loops The loops to model.
	/// @param[in] false_mm Enforce False settings in this MoveMap.  Currently
	///  only useful with the RemodelLoopMover.
	virtual
	MoverOP loop_mover_instance(
		loops::LoopsOP const loops,
		MoveMap const & false_mm
	);


protected: // fragment management


	/// @brief pick fragments of size full, 9, 3, 1
	/// @param[in] complete_ss The complete secondary structure string, typically from a Pose.
	/// @param[in] complete_aa The complete amino acid string, typically from a Pose;
	///            can be empty.  If empty, sequence bias is not used to pick fragments.
	/// @param[in] complete_abego The complete abego string, typically from setter, set_abego
	/// @param[in] interval The interval [left, right] to pick fragments from; Pose
	///  numbering (i.e. 1-based indexing).
	/// @param[in] n_frags The number of fragments to pick per position.
	void pick_all_fragments(
		String const & complete_ss,
		String const & complete_aa,
		utility::vector1< String > const & complete_abego,
		Interval const & interval,
		Size const n_frags
	);


	/// @brief pick fragments of a given length, padding when necessary
	/// @param[in] complete_ss The complete secondary structure string, typically from a Pose.
	/// @param[in] complete_aa The complete amino acid string, typically from a Pose;
	///            can be empty.  If empty, sequence bias is not used to pick fragments.
	/// @param[in] complete_abego The complete abego string, typically from a setter, set_abego
	/// @param[in] interval The interval [left, right] to pick fragments from; Pose
	///  numbering (i.e. 1-based indexing).
	/// @param[in] frag_length The desired length of the fragments
	/// @param[in] n_frags The number of fragments to pick per position.
	FrameList pick_fragments(
		String const & complete_ss,
		String const & complete_aa,
		utility::vector1< String > const & complete_abego,
		Interval const & interval,
		Size const frag_length,
		Size const n_frags
	);


protected:  //constraint management


	void setup_remodel_constraints( Pose & pose );


	void remove_remodel_constraints( Pose & pose );


private: // data


	/// @brief manages desired build instructions
	BuildManager manager_;


	/// @brief the ScoreFunction used during build
	ScoreFunctionOP sfx_;

	/// @remodelData, used when constraints are defined through Remodel
	protocols::forge::remodel::RemodelData remodel_data_;


	/// @brief at the end of apply(), reset the Pose to the original Pose if
	///  mover was not successful?  (default true)
	bool recover_original_on_failure_;


	/// @brief the string ID of the LoopMover used during build
	///  (default "RemodelLoopMover" )
	/// @details use "RemodelLoopMover" for the forge RemodelLoopMover,
	///  otherwise set it to a string recognized by
	///  create_loop_mover() in the "LoopMoverFactory".
	String loop_mover_str_;


	/// @brief cache fragments after picking? (default true)
	/// @details If true, will cache fragments and reuse them upon each
	///  apply() call instead of repicking every time.
	bool cache_fragments_;


	/// @brief option dictating whether to keep VallLibrary in memory or clear
	///  it under certain circumstances after picking fragments.
	///  (default KEEP_IN_MEMORY)
	VLB_VallMemoryUsage::Enum vall_memory_usage_;


	/// @brief the number of fragments to pick at each position (default 200)
	Size num_fragpick_;


	/// @brief also use fragments spanning the entire length of a loop?
	///  (default false)
	bool use_fullmer_;


	/// @brief full sequence string corresponding to original input pose
	///  used to pick sequence biased fragments; if empty, sequence bias
	///  is not used when picking fragments
	String original_sequence_;


	/// @brief full sequence string with length corresponding to the *new*
	///  modified pose used to pick secondary structure biased fragments.
	/// @remarks IMPORTANT: This is an override!  If this string is filled,
	///  it will be used as the string to pick secondary structure biased
	///  fragments without modification.  All secondary structure settings
	///  that might be taken from the original Pose or from the BuildInstructions
	///  *will be ignored*.  The length of this string must be equal to the
	///  length of the *NEW* modified pose, otherwise the protocol will stop
	///  with error -- you can use BuildManager::dummy_modify() to help
	///  figure things out.
	String new_secondary_structure_override_;


	/// @brief full amino acid string with length corresponding to the *new*
	///  modified pose used to pick sequence biased fragments.
	/// @remarks IMPORTANT: This is an override!  If this string is filled,
	///  it will be used as the string to pick sequence biased
	///  fragments without modification.  All sequence settings
	///  that might be taken from the original Pose or from the BuildInstructions
	///  *will be ignored*.  The length of this string must be equal to the
	///  length of the *NEW* modified pose, otherwise the protocol will stop
	///  with error -- you can use BuildManager::dummy_modify() to help
	///  figure things out.
	String new_sequence_override_;


	/// @brief the highest linear chainbreak score a chainbreak can have
	///  and still be considered closed
	Real max_linear_chainbreak_;


	/// @brief internal flag indicating if fragments for all intervals have been
	///  picked
	bool fragments_picked_;


	/// @brief full-mer to use during fragment insertion
	OrderedFragSetOP fragfull_;


	/// @brief 9-mers to use during fragment insertion
	ConstantLengthFragSetOP frag9_;


	/// @brief 3-mers to use during fragment insertion
	ConstantLengthFragSetOP frag3_;


	/// @brief 1-mers to use during fragment insertion
	ConstantLengthFragSetOP frag1_;

	/// @brief abego for picking fragments
	utility::vector1< String > abego_;

	/// @brief collection of RCGs to manage remodel constraints
	utility::vector1< RemodelConstraintGeneratorOP > rcgs_;

	/// @brief collection of movers that get a chance to modify
	/// the pose after the new lenght has been setup but before
	/// remodeling starts
	utility::vector1< moves::MoverOP > setup_movers_;

	/// @brief collection of movers that mess with the pose
	/// during RemodelLoopMover apply where n is determined by the variable below
	/// gets piped straight into RemodelLoopMover, only has an effect
	/// if this is the loop mover used
	utility::vector1< moves::MoverOP > user_provided_movers_;

	/// @brief determines how often the above movers get called
	/// in the loop mover, also only has an effect if the
	/// RemodelLoopMover is used (default 3)
	Size user_provided_movers_apply_cycle_;

	/// @brief tell RemodelLoopMover to keep the foldtree untouched
	/// also only has an effect if RemodelLoopmover is used
	/// note: will only tell loop mover is this is set to true
	bool loop_mover_fold_tree_constant_;

	/// @brief Flag to turn on restart mode, in which VLB assumes that the Pose
	///  fed to it during apply() has already been modified by the manager.
	///  (default False)
	/// @remarks In restart mode, VLB only runs the manager's dummy modify
	///  capability during apply() to get the mapping info.
	bool restart_mode_;


	/// @brief flo jan 11 reading an enzdes cstfile in the middle of
	/// remodel messes up a remodel run that was started from enzdes
	bool ignore_cmdline_enzdes_cstfile_;

	Size repeat_tail_length_;


};


} // components
} // forge
} // protocols


#endif /* INCLUDED_protocols_forge_components_VarLengthBuild_HH */
