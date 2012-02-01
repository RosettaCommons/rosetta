// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file 
/// @brief 
/// @author Javier Castellanos ( javiercv@uw.edu )

#ifndef INCLUDED_devel_constrained_sequence_ConstrainedDesign_HH
#define INCLUDED_devel_constrained_sequence_ConstrainedDesign_HH
// Package headers
#include <devel/constrained_sequence_design/ConstraintManager.fwd.hh>
#include <devel/constrained_sequence_design/SequenceConstraint.fwd.hh>
#include <devel/constrained_sequence_design/SequenceConstraintSet.fwd.hh>
// project headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/relax/RelaxProtocolBase.fwd.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

namespace devel {
namespace constrained_sequence_design {

class ConstrainedDesign: public protocols::moves::Mover {
public:
	typedef core::pose::Pose Pose;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::pack::task::TaskFactoryOP TaskFactoryOP;
	typedef protocols::moves::MoverOP MoverOP;
	typedef protocols::jd2::parser::BluePrintOP BluePrintOP;
	typedef protocols::relax::RelaxProtocolBaseOP RelaxProtocolBaseOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef protocols::moves::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;
	typedef utility::tag::TagPtr TagPtr;

public:
  ConstrainedDesign();
	ConstrainedDesign(const ConstrainedDesign& rval);

  // --- virtual functions from mover ---
  virtual std::string get_name() const { return "ConstrainedDesign"; }
  virtual void apply(Pose& pose);

	// --- virtual copy constructors 
	virtual MoverOP clone() const;


	/// @brief create this type of object
	virtual MoverOP fresh_instance() const;

	// --- setters ---
	/// @brief set the relax score function
	void set_relax_scorefxn(ScoreFunctionOP scorefxn);
	/// @brief set the relax mover
	void set_relax_mover(RelaxProtocolBaseOP relax_mover);

	virtual void parse_my_tag( TagPtr const tag,
														 DataMap & data,
														 Filters_map const &,
														 Movers_map const &,
														 Pose const & );

private:
	/// @brief score function for design
	ScoreFunctionOP scorefxn_design_;

	/// @brief relax mover
	RelaxProtocolBaseOP relax_;

	/// @brief the number of cycles of fixbb and relax
	Size nflxbb_cycles_;

	// @brief N-Termini C-Termini constraint
	Real constraints_NtoC_;

	// @brief beta sheet constraint
	Real constraints_sheet_;

	// @brief clear all residues
	bool clear_all_residues_;

	/// @brief blueprint for setting up constraints of beta sheet
	BluePrintOP blueprint_;

	/// @brief resfile, the behaivour of the operations differs
	// due to the sequence constraints.
	// AUTO		: 
	// PICKAA	: 
	// NATRO	:
	// POLAR	:
	// APOLAR	:
	std::string resfile_;

	/// @brief number of cycles of fix backbone design.
	/// Between the cycles the structure is minimized and
	/// the packertask updated to adjust the new sequence 
	/// to the sequence constraints.
	Size ncycles_;
	
	/// @brief number of cycles to run before updating
	/// the packertask with the sequence constriants.
	Size update_packertask_step_;
	
	/// @brief number of cycles to run before 
	/// relaxing the structure.
	Size relax_step_;

	/// @brief TaskFactory with only constraint_manager_ as TaskOperation.
	TaskFactoryOP task_factory_;	

	/// @brief SequenceConstraints
	SequenceConstraintSetOP constraints_;	
};

} // devel
} // constrained_sequence_design
#endif
