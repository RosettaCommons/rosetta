// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file   protocols/antibody_legacy/CDRH3Modeler.hh
/// @brief
/// @author Aroop Sircar (aroopsircar@yahoo.com)


#ifndef INCLUDED_protocols_antibody_legacy_CDRH3Modeler_hh
#define INCLUDED_protocols_antibody_legacy_CDRH3Modeler_hh


// Rosetta headers
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/FragData.fwd.hh>
#ifdef WIN32
#include <core/fragment/FragData.hh> // WIN32 INCLUDE
#include <core/fragment/FragSet.hh>
#endif
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/antibody_legacy/AntibodyClass.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/Mover.hh>


// ObjexxFCL Headers

// C++ Headers

// Utility Headers
#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace antibody_legacy {

class CDRH3Modeler;
typedef utility::pointer::owning_ptr< CDRH3Modeler > CDRH3ModelerOP;
typedef utility::pointer::owning_ptr< const CDRH3Modeler > CDRH3ModelerCOP;

//////////////////////////////////////////////////////////////////////////
/// @brief Ab initio modeling of CDR H3 loop
/// @details
class CDRH3Modeler : public moves::Mover {
public:
	/// @brief default constructor
	CDRH3Modeler( utility::vector1<core::fragment::FragSetOP> cdr_h3_frags);

	/// @brief default destructor
	~CDRH3Modeler();

	/// @brief enable CDR H3 loop building
	inline void model_h3( bool setting ) {
		do_h3_modeling_ = setting;
	}

	/// @brief enable benchmark mode
	inline void enable_benchmark_mode( bool setting ) {
		benchmark_ = setting;
	}

	/// @brief enable camelid modeling mode
	inline void set_camelid( bool setting ) {
		is_camelid_ = setting;

		if( is_camelid_ ) {
			H3_filter_ = false;
			snug_fit_ = false;
			docking_local_refine_ = false;
		}
	}

	/// @brief set centroid mode loop building
	inline void set_centroid_loop_building( bool setting ) {
		apply_centroid_mode_ = setting;
	}

	/// @brief set fullatom mode loop building
	inline void set_fullatom_loop_building( bool setting ) {
		apply_fullatom_mode_ = setting;
	};

	void set_default();
	virtual void apply( core::pose::Pose & pose_in );
	virtual std::string get_name() const;

	/// @brief Build centroid mode CDR H3 loop
	void build_centroid_loop();

	/// @brief Build fullatom mode CDR H3 loop
	void build_fullatom_loop();

	/////////////////////////////////////////////////////////////////////////
	/// @brief set scorefunction for low resolution of CDR H3 modeling
	/// @details
	void set_lowres_score_func(
	    core::scoring::ScoreFunctionOP lowres_scorefxn );

	/////////////////////////////////////////////////////////////////////////
	/// @brief set scorefunction for high resolution of CDR H3 modeling
	/// @details
	void set_highres_score_func(
	    core::scoring::ScoreFunctionOP highres_scorefxn
	);

	/// @brief insert C-terminal fragments
	void antibody_modeling_insert_ter();

	/// @brief store CDR H3 C-terminal fragments
	void store_H3_cter_fragment(
	    utility::vector1< core::fragment::FragData > & base_library_in
	);

	/// @brief return false if any cdr cutpoint is broken
	bool cutpoints_separation();

	// Compute the separation at the cutpoint. The N-C distance of the
	// peptide bond which should be formed at the cutpoint. A closed loop is
	// assumed to have a gap < 1.9 Ang
	core::Real cutpoint_separation(
	    core::pose::Pose & pose_in,
	    Size cutpoint
	);

	void scored_frag_close(
	    core::pose::Pose & pose_in,
	    loops::Loop const trimmed_cdr_h3
	);

	bool CDR_H3_filter(
	    const core::pose::Pose & pose_in,
	    core::Size const loop_begin,
	    core::Size const size,
	    char const light_chain = 'L' );

	void loop_fa_relax(
	    core::pose::Pose & pose_in,
	    core::Size const loop_begin,
	    core::Size const loop_end
	);

	void loop_centroid_relax(
	    core::pose::Pose & pose_in,
	    core::Size const loop_begin,
	    core::Size const loop_end );

	void setup_packer_task( core::pose::Pose & pose_in );

	// CDR H3 C-terminal fragments
	utility::vector1< core::fragment::FragData > H3_base_library;

private:
	core::pose::Pose template_pose_;
	core::pose::Pose start_pose_;
	/// @brief Number of ADDITIONAL residues modeled from H3_CTERM
	///        These residues range from H:n-2,n-1,n,n+1 of H3
	core::Size base_;
	core::Size c_ter_stem_;
	// constraints
	core::Real cen_cst_;
	core::Real high_cst_;
	// true if command line option "model_h3" is set
	// enables modeling of CDR H3 loop
	bool do_h3_modeling_;
	core::scoring::ScoreFunctionOP fa_scorefxn_;

	utility::vector1< core::fragment::FragSetOP > cdr_h3_frags_;

	// score functions
	core::scoring::ScoreFunctionOP lowres_scorefxn_;
	core::scoring::ScoreFunctionOP highres_scorefxn_;

	/// @brief benchmark flag
	bool benchmark_;
	/// @brief Centroid mode loop building
	bool apply_centroid_mode_;
	/// @brief Fullatom mode loop building
	bool apply_fullatom_mode_;
	/// @brief flag indicating that current loop being modeled is CDR H3
	bool current_loop_is_H3_;
	/// @brief actually enables H3 filter for H3 operations
	bool H3_filter_;
	/// @brief build H3 only
	bool antibody_build_;
	/// @brief refine H3 only
	bool antibody_refine_;
	/// @brief lower amplitude during base relaxation
	bool min_base_relax_;
	/// @brief use random cutpoints for h3 modeling
	bool h3_random_cut_;
	/// @brief cutpoint whose separation is computed in scorefile
	Size decoy_loop_cutpoint_;
	/// @brief enable docking local refine of LH chains & simultaneous H3 min
	bool snug_fit_;
	/// @brief loop_building in docking
	bool loops_flag_;
	bool docking_local_refine_;
	/// @brief insert fragment in docking
	bool dle_flag_;
	/// @brief just refine input loop
	bool refine_input_loop_;
	/// @brief number of flanking residues:default 5
	core::Size h3_flank_;
	/// @brief relax flanking regions of h3
	bool flank_relax_;
	/// @brief freeze h3 during all cdr relax and local refine
	bool freeze_h3_;
	/// @brief is camelid antibody without light chain
	bool is_camelid_;
	/// @brief size of loop above which 9mer frags are used
	core::Size cutoff_9_; // default 16
	/// @brief size of loop above which 3mer frags are used
	core::Size cutoff_3_; // default 6

	Antibody antibody_in_;

	//packer task
	core::pack::task::TaskFactoryOP tf_;
	core::pack::task::TaskFactoryOP init_task_factory_;

}; // class CDRH3Modeler

// read CDR H3 C-terminal fragments (size: 4)
void read_H3_cter_fragment(
    Antibody & antibody_in,
    utility::vector1< core::fragment::FragData > & H3_base_library,
    bool is_camelid
);

void simple_one_loop_fold_tree(
    core::pose::Pose & pose,
    loops::Loop const & loop
);

void simple_fold_tree(
    core::pose::Pose & pose_in,
    core::Size jumppoint1,
    core::Size cutpoint,
    core::Size jumppoint2
);

} // moves
} // protocols


#endif
