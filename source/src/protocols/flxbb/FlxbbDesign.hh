// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file protocols/flxbb/FlxbbDesign.hh
/// @brief perform cycles of design and relax with filter
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

#ifndef INCLUDED_protocols_flxbb_FlxbbDesign_hh
#define INCLUDED_protocols_flxbb_FlxbbDesign_hh

// Unitt Header
#include <protocols/flxbb/FlxbbDesign.fwd.hh>

// Package Headers
#include <protocols/jd2/parser/BluePrint.fwd.hh>
#include <protocols/flxbb/DesignTask.fwd.hh>
#include <protocols/flxbb/FilterStructs.fwd.hh>
#ifdef __clang__
#include <protocols/flxbb/FilterStructs.hh>
#endif

// Project Headers
#include <core/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <utility/tag/Tag.fwd.hh>

#include <core/pack/task/operation/TaskOperation.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#ifdef WIN32
	#include <core/pack/task/operation/TaskOperation.hh>
	#include <protocols/flxbb/DesignTask.hh>
#endif

namespace protocols {
namespace flxbb{


///////////////////////////////////////////////////////////////////////////////////////////////////////
class FlxbbDesign: public protocols::moves::Mover {
public:


	typedef protocols::moves::Mover Super;
	typedef std::string String;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::kinematics::MoveMapOP MoveMapOP;
	typedef core::pose::Pose Pose;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;

	typedef protocols::moves::MoverOP MoverOP;
	typedef protocols::jd2::parser::BluePrint BluePrint;
	typedef protocols::jd2::parser::BluePrintOP BluePrintOP;
	typedef protocols::flxbb::DesignTaskOP DesignTaskOP;
	typedef protocols::flxbb::DesignTaskSet DesignTaskSet;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public: // constructor/destructor


	/// @brief default constructor
	FlxbbDesign();

	/// @brief value constructor
	FlxbbDesign(
		ScoreFunctionOP const sfxnd,
		ScoreFunctionOP const sfxnr,
		Size const ncycle = 3,
		String const layer_mode = "",
		bool const use_origseq_for_not_dsgned_layer = true,
		bool const no_relax = false );

	/// @brief copy constructor
	FlxbbDesign( FlxbbDesign const & rval );

	/// @brief destructor
	virtual ~FlxbbDesign();


public: // virtual constructors


	/// @brief clone this object
	virtual
	MoverOP clone() const;


	/// @brief create this type of object
	virtual
	MoverOP fresh_instance() const;


public:


	/// @brief read parameters
	void read_options();

	/// @brief initialize filter parameters
	void initialize_filter( Size const filter_trial, String const & filter_type );

	/// @brief register options
	static void register_options();


public: // virtual main operation


	/// @brief mover apply
	virtual void apply( Pose & pose );

	virtual std::string get_name() const;

public:// mutators


	/// @brief set the score function for fixbb design
	void set_scorefxn_design( ScoreFunctionOP const scorefxn );

	/// @brief set the score function for relax
	void set_scorefxn_relax( ScoreFunctionOP const scorefxn );

	/// @brief the number of cycles of fixbb and relax
	void set_ncycles( Size const ncycles );

	/// @brief set layer mode
	void set_layer_mode( String const & layer_mode );

	/// @brief use original sequence for not designed region in layer_mode,
	/// otherwise the residues of the region turned into Ala ( default true )
	void use_origseq_for_not_dsgned_layer( bool const use );

	/// @brief relax is not performed after design (default false)
	void no_relax( bool const no_relax );

	/// @brief design protocol will not be performed (default false)
	void no_design( bool const no_design );

	/// @brief set BluePrintOP
	void set_blueprint( BluePrintOP const blueprint );

	/// @brief set weight of constraints_sheet which constrains between Ca atoms in beta-sheets
	/// if this weight > 0.0, the constraints is applied ( default = -1.0 )
	void set_weight_constraints_sheet( Real const value );

	/// @brief set weight of constraints_NtoC which constrain between Ca atoms of C- and N-terminal
	/// if this weight > 0.0, the constraints is applied ( default = -1.0 )
	void set_weight_constraints_NtoC( Real const value );

	/// @brief set movemap for relax
	void movemap_from_blueprint( bool const value );

	/// @brief set movemap for relax
	void set_movemap( MoveMapOP const movemap );


public:// mutators relevant to the DesignTaskSet


	/// @brief set DesignTaskSet
	/// Once you set DesignTaskSet, this controls almost every setups of this class
	void set_design_taskset( DesignTaskSet const & design_taskset );

	/// @brief add the design task which cotrols the iteration of fixbb design and relax
	void add_design_task( DesignTaskOP const design_task );

	/// @brief filtering during design
	void set_filter_during_design( FilterStructsOP const filter_during_design );

	/// @brief clear DesignTaskSet
	void clear_design_taskset();


public:// parser

	virtual void parse_my_tag( TagCOP const tag,
														 basic::datacache::DataMap & data,
														 Filters_map const &,
														 Movers_map const &,
														 Pose const & );


private:// helper functions


	/// @brief
	DesignTaskSet build_design_taskset( Pose const & pose );


	/// @brief make movemap from blueprint for relax
	MoveMapOP get_movemap_from_blueprint( Pose const & pose, BluePrintOP const blueprint );


private:


	/// @brief score function for design
	ScoreFunctionOP scorefxn_design_;

	/// @brief score function for relax
	ScoreFunctionOP scorefxn_relax_;

	/// @brief the number of cycles of fixbb and relax
	Size nflxbb_cycles_;

	/// @brief perform fixbb in layer mode
	String layer_mode_;

	/// @brief use original sequence for not designed region in layer_mode, otherwise the residues of the region turned into Ala
	bool use_origseq_for_not_dsgned_layer_;

	/// @brief relax is not performed after design (default false)
	bool no_relax_;

	/// @brief design protocol will not be performed (default false)
	bool no_design_;

	/// @brief control the method of design with filter and relax
	DesignTaskSet design_taskset_;

	/// @brief filter during design
	FilterStructsOP filter_during_design_;

	/// @brief blueprint for setting up constraints of beta sheet
	BluePrintOP blueprint_;

	/// @brief resfile
	String resfile_;

	/// @brief weight of constraints_sheet which constrains between Ca atoms in beta-sheets
	/// if this weight > 0.0, the constraints is applied ( default = -1.0 )
  Real constraints_sheet_;

	/// @brief weight of constraints_NtoC which constrains between Ca atoms of C- and N-terminal
	/// if this weight > 0.0, the constraints is applied ( default = -1.0 )
	Real constraints_NtoC_;

	/// @brief set movemap from blueprint ( default false )
	bool movemap_from_blueprint_;

	/// @brief mavemap for relax ( defalut NULL )
	MoveMapOP movemap_;

	/// @brief task operations
	utility::vector1< TaskOperationOP > task_operations_;

	/// @brief use fast relax or not
	bool use_fast_relax_;

	/// @brief clear all residues before design
	bool clear_all_residues_;

	// constraint of backbone to fixbb-designed structure during relax
	bool relax_constraint_to_design_;

	// @brief Exclude aromatic chi2 rotamers, of which angles are around 0
	bool limit_aroma_chi2_;


};

///////////////////////////////////////////////////////////////////////////////////////////////////////
class FlxbbDesignPack: public protocols::simple_moves::PackRotamersMover {
public:

	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef protocols::flxbb::FilterStructsOP FilterStructsOP;

public:

	FlxbbDesignPack();

	FlxbbDesignPack(
  	ScoreFunctionCOP scorefxn,
		PackerTaskCOP task,
		FilterStructsOP filter=0 );

	virtual ~FlxbbDesignPack();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:

	FilterStructsOP filter_;

};


} // namespace flxbb
} // namespace protocols

#endif
