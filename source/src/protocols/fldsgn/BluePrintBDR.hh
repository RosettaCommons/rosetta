// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/fldsgn/BluePrintBDR.hh
/// @brief
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_fldsgn_BluePrintBDR_hh
#define INCLUDED_protocols_fldsgn_BluePrintBDR_hh

// unit headers
#include <protocols/fldsgn/BluePrintBDR.fwd.hh>

// type headers
#include <core/types.hh>

// project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/dssp/Dssp.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>
#include <protocols/forge/build/BuildInstruction.fwd.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/components/VarLengthBuild.fwd.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/InvrotTree.fwd.hh>
// AUTO-REMOVED #include <protocols/forge/constraints/SheetConstraintsRCG.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>



namespace protocols {
namespace fldsgn {


class BluePrintBDR : public protocols::moves::Mover {


private: // typedefs


	typedef protocols::moves::Mover Super;


public: // typedefs


	typedef std::string String;

	typedef core::Real Real;
	typedef core::Size Size;

	typedef core::scoring::dssp::Dssp Dssp;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;

	typedef protocols::jd2::parser::BluePrint BluePrint;
	typedef protocols::jd2::parser::BluePrintOP BluePrintOP;
	typedef protocols::forge::build::BuildInstructionOP BuildInstructionOP;
	typedef protocols::forge::build::BuildManager BuildManager;
	typedef protocols::forge::build::Interval Interval;
	typedef protocols::forge::components::VarLengthBuild VarLengthBuild;
	typedef protocols::forge::components::VarLengthBuildOP VarLengthBuildOP;
	typedef protocols::moves::MoverOP MoverOP;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public: // construct/destruct


	/// @brief default constructor
	BluePrintBDR();

	/// @brief value constructor
	BluePrintBDR( String const & filename, bool const ss_from_blueprint=true );

	/// @brief value constructor
	BluePrintBDR( BluePrintOP const & blueprintOP, bool const ss_from_blueprint=true );

	/// @brief copy constructor
	BluePrintBDR( BluePrintBDR const & rval );

	/// @brief default destructor
	virtual	~BluePrintBDR();


private: // disallow assignment


	/// @brief copy assignment
	/// @remarks Mover base class prevents this from working properly...
	BluePrintBDR & operator =( BluePrintBDR const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	MoverOP clone() const;

	/// @brief create this type of object
	virtual
	MoverOP fresh_instance() const;


public: // accessors


	/// @brief the centroid level score function, default "fldsgn_cen"
	ScoreFunction const & scorefunction() const;

	////
	Size instruction_size() const { return manager_.size(); }

	/// @brief the loop mover string to use during centroid build
	///  (default "RemodelLoopMover")
	/// @remarks set to either a string the create_loop_mover() LoopMoverFactory
	///  recognizes or the "RemodelLoopMover"
	inline String const & loop_mover_str() const {	return loop_mover_str_; }

	/// @brief use full-mer fragments when building the loop? (default false)
	inline bool use_fullmer() const {	return use_fullmer_; }

	/// @brief the number of fragments to pick at each position
	inline Size num_fragpick() const { return num_fragpick_; }

	/// @brief use sequence biased fragments when building the loop? (default false)
	inline bool use_sequence_bias() const {	return use_sequence_bias_; }

	/// @brief use abego biased fragments when building the loop? (default false)
	inline bool use_abego_bias() const { return use_abego_bias_; }

	/// @brief the maximum allowed linear chainbreak (default 0.07)
	inline Real max_linear_chainbreak() const {	return max_linear_chainbreak_; }

	/// @brief define secondary structrue by blueprint
	inline bool is_initialized() const { return initialized_;	}


public: // mutators


	/// @brief add instruction to the manager of this BluePrintBDR (no copy)
	/// @param[in] bi BuildInstruction
	/// @param[in] aa_during_design_refine The allowed amino acid sequence
	///  during design.  Only applicable to BuildInstructions like
	///  SegmentRebuild and SegmentInsert.  Make sure the length of this
	///  string matches up properly, and remember to use any special characters,
	///  e.g. the insertion character for SegmentInsert
	void add_instruction(	BuildInstructionOP bi );

	/// @brief create directed dependency between two instructions
	void create_directed_dependency( BuildInstructionOP u, BuildInstructionOP v	);

	/// @brief use full-mer fragments when building the loop?
	inline void use_fullmer( bool const flag ) {	use_fullmer_ = flag; }

	/// @brief use sequence biased fragments when building the loop? (default false)
	inline void use_sequence_bias( bool const flag ) { use_sequence_bias_ = flag; }

	/// @brief use abego biased fragments when building the loop? (default false)
	inline void use_abego_bias( bool const flag ) { use_abego_bias_ = flag; }

	/// @brief the maximum allowed linear chainbreak
	inline void max_linear_chainbreak( Real const threshold ) { max_linear_chainbreak_ = threshold; }

	/// @brief the loop mover string to use during centroid build
	/// @remarks set to either a string the create_loop_mover() LoopMoverFactory
	///  recognizes or the "RemodelLoopMover"
	inline void loop_mover_str( String const & loop_mover_str ) { loop_mover_str_ = loop_mover_str;	}

	/// @brief the number of fragments to pick at each position
	inline void num_fragpick( Size const num ){ num_fragpick_ = num;	}

	/// @brief define secondary structrue by blueprint
	inline void ss_from_blueprint( bool const flag ){ ss_from_blueprint_ = flag;	}

	/// @brief set the centroid level score function
	void scorefunction( ScoreFunction const & sfx );

	/// @brief set the centroid level score function
	void scorefunction( ScoreFunctionOP sfx );

	/// @brief set blueprint file by filename
	void set_blueprint( String const & filename );

	/// @brief set blueprint file
	void set_blueprint( BluePrintOP const & blp );

	/// @brief set constraints between N- and C- Ca atoms
	void set_constraints_NtoC( Real const & weight );

	/// @brief set constraint file
	void set_constraint_file( String const & constraint_file );

	/// @brief dump pdb when this protocol failed
	void dump_pdb_when_fail( String const & dump_pdb_when_fail );

	/// @brief set list of remodel constraint generators
	void set_rcgs( utility::vector1< protocols::forge::remodel::RemodelConstraintGeneratorOP > const & rcgs );

public: // virtual main methods


	/// @brief apply defined moves to given Pose
	virtual
	void apply( Pose & pose );

	virtual
	std::string get_name() const;


public: //parser


	/// @brief parse xml file
	void parse_my_tag( TagCOP tag,
										 basic::datacache::DataMap & data,
										 Filters_map const &,
										 Movers_map const &,
										 Pose const & );


private: // protocol methods


	/// @brief run the centroid level build stage
	/// @return true if loop closed, false otherwise
	bool centroid_build( Pose & pose );

	/// @brief set instruction by blueprint
	bool set_instruction_blueprint( Pose const & pose );

	/// @brief set up folding aroung theozyme
	/// in a vlb
	void
	setup_invrot_tree_in_vlb( VarLengthBuild & vlb, Pose & pose ) const;


private: // data

	/// @brief bluerprint file for setting build instruction
	BluePrintOP blueprint_;

	/// @brief the BuildManager
	BuildManager manager_;

	/// @brief the centroid scorefunction to use, default "remodel_cen"
	ScoreFunctionOP sfx_;

	/// @brief the loop mover string to use during centroid build
	///  (default "RemodelLoopMover")
	/// @remarks set to either a string the create_loop_mover() LoopMoverFactory
	///  recognizes or the "RemodelLoopMover"
	String loop_mover_str_;

	/// @brief use full-mer fragments when building the loop? (default false)
	bool use_fullmer_;

	/// @brief the number of fragments to pick at each position
	Size num_fragpick_;

	/// @brief use sequence biased fragments when building the loop? (default false)
	bool use_sequence_bias_;

	/// @brief use abego biased fragments when building the loop? (default false)
	bool use_abego_bias_;

	/// @brief the maximum allowed linear chainbreak (default 0.07)
	Real max_linear_chainbreak_;

	/// @brief vlb_ is intialized or not
	bool initialized_;

	/// @brief define secondary structure by blueprint
	bool ss_from_blueprint_;

	/// @brief add constraints between N- and C- terminal Ca atoms
	Real constraints_NtoC_;

	/// @brief add constraints betwee Ca atoms in sheet
	Real constraints_sheet_;

	/// @brief constraint file
	String constraint_file_;

	/// @brief dump pdb when this protocol failed
	String dump_pdb_when_fail_;

	/// @brief loop mover for loop building, default RemodelLoopMover
	String loop_mover_;

	/// @brief number of allowed_closure_attempts_ of RemodelLoopMover
	Size rmdl_attempts_;

	/// @brief Entire sequence except for rebuilding regions become poly Val
	bool use_poly_val_;

	///@brief if we a fold tree has already been set up prior
	bool tell_vlb_to_not_touch_fold_tree_;

	/// @brief in case we're folding up around a ligand
	protocols::toolbox::match_enzdes_util::InvrotTreeOP invrot_tree_;
	toolbox::match_enzdes_util::EnzConstraintIOCOP enzcst_io_;

	/// @brief User-specified constraint generators which will be called from the VarLengthBuild
	utility::vector1< protocols::forge::remodel::RemodelConstraintGeneratorOP > rcgs_;

private: // per-stage movers


	/// @brief VLB for centroid_build
	/// @remarks Store it here instead of re-instantiation each centroid_build()
	///  call so fragment caching works.
	/// @warning Remember to set this to null if the BuildManager changes.
	VarLengthBuildOP vlb_;


};


} // namespace fldsgn
} // namespace protocols


#endif /* INCLUDED_protocols_forge_components_BluePrintBDR_HH */
