// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/AnchoredDesign/AnchorMoversData.hh
/// @brief header for AnchorMoversData
/// @author Steven Lewis

#ifndef INCLUDED_protocols_anchored_design_AnchorMoversData_hh
#define INCLUDED_protocols_anchored_design_AnchorMoversData_hh

// Unit Headers
#include <protocols/anchored_design/AnchorMoversData.fwd.hh>
#include <protocols/anchored_design/Anchor.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh> //ultimately, non OPed loop objects are contained
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/fragment/FragSet.fwd.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1_bool.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/keys/Key3Tuple.hh> //container holds Loop, movemap, movemap (OPs)

#include <utility/vector1.hh>


// C++ Headers

namespace protocols {
namespace anchored_design {

/// @brief This data class wraps all the data needed for the AnchoredDesign movers.
/// @details This data class keeps the anchor, mobile loops, and movemaps associated with the loops all together.
///It generates its own movemaps.
class AnchorMoversData : public utility::pointer::ReferenceCount
{
public:
	//////////////////////////////////typedefs//////////////////////////////////////////////
	//these are internal data structures
	typedef utility::keys::Key3Tuple< protocols::loops::Loop, core::kinematics::MoveMapOP, core::kinematics::MoveMapOP> Loop_mm_tuple;
	typedef utility::vector1< Loop_mm_tuple > Loop_mm_tuples;

	//////////////////mechanicals (ctor, dtor, etc)/////////////////////////////////////
	/// @brief empty constructor is empty - don't use it unless you intend to set alllllllll the data manually...
	AnchorMoversData();

	/// @brief constructor takes an anchor and loop object and sets up internals reasonably; options boolean is optionally optional.  If you use this constructor, you will need to later manually give fragments to the AnchorMoversData object or use one of its fragments functions to determine what type of fragments to use.
	AnchorMoversData(
		protocols::anchored_design::AnchorCOP anchor,
		protocols::loops::Loops const & loops,
		bool const options=false
	);

	/// @brief constructor takes a pose to generate its own internals (also using the option system)
	AnchorMoversData( core::pose::Pose const & pose );

	/// @brief virtual dtors make c++ happy
	virtual ~AnchorMoversData();

	/// @brief copy ctor
	AnchorMoversData( AnchorMoversData const & rhs );

	/// @brief assignment operator
	AnchorMoversData & operator=( AnchorMoversData const & rhs );

	AnchorMoversDataOP clone() const;

	/// @brief randomly reset loop cutpoints.  Useful only when starting structure is well-closed.  Best for MPI-style runs
	void pick_new_cutpoints( bool reset_always );

	/// @brief set up constraints for anchor_noise_constraints_mode
	void anchor_noise_constraints_setup( core::pose::Pose & pose );

	///////////////////set functions/////////////////////////////////////////////////////////
	/// @brief set fragments object
	void set_frags( core::fragment::FragSetOP );
	/// @brief brew its own design frags (LLLLL secondary structure)
	core::fragment::FragSetCOP autogenerate_design_frags();
	/// @brief brew its own sequence-specific frags
	core::fragment::FragSetCOP autogenerate_constseq_frags(std::string const & seq);
	/// @brief figure out what type of fragments to use, and use them; pose needed to make packertask to detect if design is occuring
	void autogenerate_frags( core::pose::Pose const & pose );
	/// @brief set packertask factory
	void set_task_factory( core::pack::task::TaskFactoryOP );
	/// @brief set fullatom scorefunction
	void set_fullatom_scorefunction( core::scoring::ScoreFunctionOP );
	/// @brief set centroid scorefunction
	void set_centroid_scorefunction( core::scoring::ScoreFunctionOP );
	/// @brief set up kinematics' loops and anchors
	void set_loops_and_anchor( protocols::anchored_design::AnchorCOP anchor, protocols::loops::Loops loops);

	/////////////////loops and movemap functions//////////////////////////////////////////////
	/// @brief access anchored loop
	protocols::loops::Loop const & anchored_loop() const; // {return anchor_loop_index_;}

	/// @brief access for movemap that covers all loops; fullatom
	core::kinematics::MoveMapOP movemap_fa_all() const;// {return movemap_fa_all_;}
	/// @brief access for movemap that covers all loops; centroid
	core::kinematics::MoveMapOP movemap_cen_all() const;// {return movemap_cen_all_;}

	/// @brief number of loops/mms
	inline core::Size num_loops() const { return loops_and_fa_mms_.size(); } //WHY NOT CEN HERE?
	/// @brief accessor for a loop
	protocols::loops::Loop const & loop( core::Size i ) const;// { return loops_and_fa_mms_[ i ].key1(); }
	/// @brief accessor for omega-variable movemap (most movers); fullatom phase (may allow anchor movement w/constraints)
	core::kinematics::MoveMapOP movemap_fa( core::Size i ) const;// { return loops_and_fa_mms_[ i ].key2(); }
	/// @brief accessor for omega-fixed movemap (appropriate for CCD movers); fullatom phase (may allow anchor movement w/constraints)
	core::kinematics::MoveMapOP movemap_fa_omegafixed( core::Size i ) const;// { return loops_and_fa_mms_[ i ].key3(); }
	/// @brief accessor for omega-variable movemap (most movers); centroid phase (no anchor movement)
	core::kinematics::MoveMapOP movemap_cen( core::Size i ) const;// { return loops_and_cen_mms_[ i ].key2(); }
	/// @brief accessor for omega-fixed movemap (appropriate for CCD movers); centroid phase
	core::kinematics::MoveMapOP movemap_cen_omegafixed( core::Size i ) const;// { return loops_and_cen_mms_[ i ].key3(); }
	/// @brief accessor for loops object
	protocols::loops::Loops const & loops() const;


	/////////////////////////////anchor functions//////////////////////////////////////////////
	/// @brief access for anchor start
	core::Size anchor_start() const;
	/// @brief access for anchor end
	core::Size anchor_end() const;

	////////////////////////scfxn, packertask, fragments accessors//////////////////////////////
	/// @brief access fragments object
	core::fragment::FragSetCOP get_frags() const;// { return fragset_; }
	/// @brief access packertask factory
	core::pack::task::TaskFactoryCOP get_task_factory() const;// { return task_factory_; }
	/// @brief access packertask factory
	core::pack::task::TaskFactoryCOP get_late_factory() const;// { return task_factory_; }
	/// @brief access fullatom scorefunction
	core::scoring::ScoreFunctionOP get_fullatom_scorefunction() const;// { return fullatom_scorefunction_; }
	/// @brief access centroid scorefunction
	core::scoring::ScoreFunctionOP get_centroid_scorefunction() const;// { return centroid_scorefunction_; }
	/// @brief access centroid scorefunction for minimization
	core::scoring::ScoreFunctionOP get_centroid_scorefunction_min() const;// { return centroid_scorefunction_min_; }
	/// @brief return packertask from factory
	core::pack::task::PackerTaskOP get_task(core::pose::Pose const & pose) const;


	////////////////option system replacement/////////////////////////////
	/// @brief dye position used in dye modeling publication
	void set_akash_dyepos(core::Size const akash_dyepos);
	/// @brief used for unbound mode
	void set_unbound_mode(bool unbound_mode);
	/// @brief used to test anchoring via constraints
	void set_anchor_via_constraints(bool anchor_via_constraints);
	/// @brief VDW weight in centroid scorefunction
	void set_VDW_weight(core::Real VDW_weight);
	/// @brief chainbreak weight in fullatom scorefunction
	void set_chainbreak_weight(core::Real chainbreak_weight);
	/// @brief allow anchor to repack
	void set_allow_anchor_repack(bool allow_anchor_repack);
	/// @brief resfile for design
	void set_resfile_1(std::string const & resfile_1);
	/// @brief later-stage resfile if desired
	void set_resfile_2(std::string const & resfile_2);
	// /// @brief whether to automatically initialize from the options system; defaults to true
	//void set_autoinitialize(bool const autoinitialize);
	/// @brief copy of cmdline option loop_file
	//void set_loop_file(std::string const & loop_file); //This setter in un-provided because it cannot be used; the variable is only used in the option-system-dependent ctor; it's only a class variable to fit the rules
	/// @brief copy of cmdline option frag3
	void set_frag3(std::string const & frag3);
	/// @brief do not use fragments?
	void set_no_frags(bool no_frags);
	/// @brief special anchor_noise_constraints_mode
	void set_anchor_noise_constraints_mode(bool anchor_noise_constraints_mode);
	/// @brief special super_secret_fixed_interface_mode
	void set_super_secret_fixed_interface_mode(bool super_secret_fixed_interface_mode);

	/// @brief dye position used in dye modeling publication
	core::Size get_akash_dyepos() const;
	/// @brief used for unbound mode
	bool get_unbound_mode() const;
	/// @brief used to test anchoring via constraints
	bool get_anchor_via_constraints() const;
	/// @brief VDW weight in centroid scorefunction
	core::Real get_VDW_weight() const;
	/// @brief chainbreak weight in fullatom scorefunction
	core::Real get_chainbreak_weight() const;
	/// @brief allow anchor to repack
	bool get_allow_anchor_repack() const;
	/// @brief resfile for design
	std::string const & get_resfile_1() const;
	/// @brief later-stage resfile if desired
	std::string const & get_resfile_2() const;
	// /// @brief whether to automatically initialize from the options system; defaults to true
	//bool get_autoinitialize() const;
	/// @brief copy of cmdline option loop_file
	std::string const & get_loop_file() const;
	/// @brief copy of cmdline option frag3
	std::string const & get_frag3() const;
	/// @brief do not use fragments?
	bool get_no_frags() const;
	/// @brief special anchor_noise_constraints_mode
	bool get_anchor_noise_constraints_mode() const;
	/// @brief special super_secret_fixed_interface_mode
	bool get_super_secret_fixed_interface_mode() const;

	/// @brief read options from the option system
	void read_options();

	/// @brief get string name for interface_calc_
	std::string const & interface_calc() const;
	/// @brief get string name for neighborhood_calc_
	std::string const & neighborhood_calc() const;


private:
	////////////////////////private functions generate internal data from input and defaults/////////////////
	/// @brief constructor subunit; generates movemaps from loops, including master movemap_all_
	void setup_movemaps();
	/// @brief setup_movemaps subunit
	void set_movemap( core::kinematics::MoveMapOP movemap, core::Size seqpos, bool omega = false );
	/// @brief setup_movemaps subunit; boolean controls fullatom ignoring of fixing movemap if constraints exist
	void fix_anchor( core::kinematics::MoveMapOP movemap, bool const centroid );
	/// @brief constructor subunit, determines which loop is anchor's loop
	void locate_anchor_loop();
	/// @brief constructor subunit, rearranges multiple-loops input structure into internal data structure
	void input_loops_into_tuples( protocols::loops::Loops const & loops );
	/// @brief determined which pointers are unset, and sets them if possible
	void set_unset_defaults();
	/// @brief set_unset_defaults subunit, sets default scorefunctions as necessary
	void set_unset_scorefunctions();
	/// @brief set_unset_defaults subunit, sets default packertask factory as necessary
	void set_unset_packertask_factory();
	/// @brief randomly reset just one cutpoint; used by pick_new_cutpoints
	core::Size pick_new_cutpoint( core::Size const loopstart, core::Size const loopend );

	/// @brief the anchor itself
	protocols::anchored_design::AnchorCOP anchor_;
	/// @brief location of the anchored loop within loop_mm_tuples vector
	core::Size anchor_loop_index_;
	/// @brief movemap allowing all loops to move; fullatom phase
	core::kinematics::MoveMapOP movemap_fa_all_;
	/// @brief movemap allowing all loops to move; centroid phase
	core::kinematics::MoveMapOP movemap_cen_all_;
	/// @brief pairs of loops and movemaps to move those loops for loop closure/smallmover; fullatom versions
	Loop_mm_tuples loops_and_fa_mms_;
	/// @brief pairs of loops and movemaps to move those loops for loop closure/smallmover; centroid versions
	Loop_mm_tuples loops_and_cen_mms_;
	/// @brief a copy of the input loops object, in case a Loops type is needed
	protocols::loops::Loops loops_;

	/// @brief fragments if we've got them
	core::fragment::FragSetOP fragset_;
	/// @brief PackerTask factory
	core::pack::task::TaskFactoryOP task_factory_;
	/// @brief second TaskFactory for more rotamers late in refinement
	core::pack::task::TaskFactoryOP late_factory_;
	/// @brief fullatom scorefunction
	core::scoring::ScoreFunctionOP fullatom_scorefunction_;
	/// @brief centroid scorefunction
	core::scoring::ScoreFunctionOP centroid_scorefunction_;
	/// @brief centroid scorefunction with no-derivative functions off (for minimizing)
	core::scoring::ScoreFunctionOP centroid_scorefunction_min_;

	/// @brief calculator name
	std::string const interface_calc_;
	/// @brief calculator name
	std::string const neighborhood_calc_;

	//option system replacement internals
	/// @brief dye position used in dye modeling publication
	core::Size akash_dyepos_;
	/// @brief used for unbound mode
	bool unbound_mode_;
	/// @brief used to test anchoring via constraints
	bool anchor_via_constraints_;
	/// @brief VDW weight in centroid scorefunction
	core::Real VDW_weight_;
	/// @brief chainbreak weight in fullatom scorefunction
	core::Real chainbreak_weight_;
	/// @brief allow anchor to repack
	bool allow_anchor_repack_;
	/// @brief resfile for design
	std::string resfile_1_;
	/// @brief later-stage resfile if desired
	std::string resfile_2_;
	// /// @brief whether to automatically initialize from the options system; defaults to true
	//bool autoinitialize_;
	/// @brief copy of loop_file cmdline option
	std::string loop_file_;
	/// @brief copy of frag3 cmdline option
	std::string frag3_;
	/// @brief use no fragments
	bool no_frags_;
	/// @brief special anchor_noise_constraints_mode
	bool anchor_noise_constraints_mode_;
	/// @brief special super_secret_fixed_interface_mode
	bool super_secret_fixed_interface_mode_;

};//end AnchorMoversData class

}//AnchoredDesign
}//protocols

#endif //INCLUDED_protocols_AnchoredDesign_AnchorMoversData_HH
