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

// Unit headers
#include <devel/constrained_sequence_design/ConstrainedDesign.hh>
#include <devel/constrained_sequence_design/ConstrainedDesignMoverCreator.hh>

// Types header
#include <core/types.hh>

// Package headers
#include <devel/constrained_sequence_design/ConstraintManager.hh>
#include <devel/constrained_sequence_design/SequenceConstraintSet.hh>
#include <devel/constrained_sequence_design/SequenceConstraintLoader.hh>

// Project headers

#include <protocols/flxbb/utility.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/FastRelax.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/simple_moves/MakePolyXMover.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/moves/DataMap.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

// STL header
#include <algorithm>

static basic::Tracer TR("devel.constrained_sequence_design.ConstrainedDesign");

namespace devel {
namespace constrained_sequence_design {

using namespace core;
using namespace utility;
typedef protocols::moves::MoverOP MoverOP;

std::string
ConstrainedDesignMoverCreator::keyname() const
{
	return ConstrainedDesignMoverCreator::mover_name();
}

MoverOP
ConstrainedDesignMoverCreator::create_mover() const {
	return new ConstrainedDesign;
}

std::string
ConstrainedDesignMoverCreator::mover_name()
{
	return "ConstrainedDesign";
}

// @brief default constructor
ConstrainedDesign::ConstrainedDesign():
	protocols::moves::Mover( "ConstrainedDesign" ),
	scorefxn_design_( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::STANDARD_WTS, core::scoring::SCORE12_PATCH ) ),
	relax_( new protocols::relax::FastRelax(core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::STANDARD_WTS, core::scoring::SCORE12_PATCH ), 5) ) ,
	blueprint_( NULL ),
	resfile_( "" ),
	constraints_sheet_( -1.0 ),
	constraints_NtoC_( -1.0 ),
	clear_all_residues_( false ),
	ncycles_( 40 ),
  relax_step_( 10 ),
	update_packertask_step_( 5 ),
	task_factory_( new core::pack::task::TaskFactory() )
{ } 

/// @brief create this type of object
ConstrainedDesign::ConstrainedDesign(const ConstrainedDesign& rval):
		protocols::moves::Mover(rval),
		scorefxn_design_(rval.scorefxn_design_),
		relax_(rval.relax_),
		blueprint_(rval.blueprint_),
		resfile_(rval.resfile_),
		constraints_sheet_(rval.constraints_sheet_),
		constraints_NtoC_(rval.constraints_NtoC_),
		clear_all_residues_(rval.clear_all_residues_),
		ncycles_(rval.ncycles_),
		relax_step_(rval.relax_step_),
		update_packertask_step_(rval.update_packertask_step_),
		task_factory_( rval.task_factory_ ),
		constraints_( rval.constraints_ )
{ } 

/// @brief clone this object
MoverOP 
ConstrainedDesign::clone() const {
		return new ConstrainedDesign( *this );
}

/// @brief create this type of object
MoverOP 
ConstrainedDesign::fresh_instance() const {
		return new ConstrainedDesign();
}

/// @brief set the relax score function
void 
ConstrainedDesign::set_relax_scorefxn(ScoreFunctionOP scorefxn){
	relax_->set_scorefxn(scorefxn);
}

/// @brief set the relax mover
void
ConstrainedDesign::set_relax_mover(RelaxProtocolBaseOP relax_mover) {
		relax_ = relax_mover;
}

void
ConstrainedDesign::apply(Pose& pose) { 
	using core::util::switch_to_residue_type_set;
	using core::scoring::ScoreFunctionCOP;
	using core::scoring::constraints::ConstraintSetOP;
	using core::scoring::constraints::ConstraintSet;
	using core::pack::task::PackerTaskOP;
	using protocols::moves::MoverOP;
	using protocols::flxbb::constraints_sheet;
	using protocols::flxbb::constraints_NtoC;
	using protocols::simple_moves::PackRotamersMover;
	using protocols::relax::FastRelax;
	using protocols::simple_moves::MakePolyXMover;
	using protocols::relax::RelaxProtocolBaseCOP;

	// set pose to fullatom
	if( ! pose.is_fullatom() ){
		core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
	}
	// set sequence constraints
	ConstraintManagerOP constraint_manager = new ConstraintManager(pose, constraints_);
	for(SequenceConstraintSet::iterator it = constraints_->begin(); it != constraints_->end(); ++it)
		constraint_manager->add_constraint(*it);
	task_factory_->clear();
	task_factory_->push_back(constraint_manager);


	// set constraints
	ConstraintSetOP cstset = new ConstraintSet;
	// set weight of constraints
	if( constraints_NtoC_ > 0.0 || constraints_sheet_ > 0.0 ){
		RelaxProtocolBaseCOP rpb(relax_);
//		ScoreFunction &  scorefxn = rpb->get_scorefxn();
		Real cst_weight( rpb->get_scorefxn()->get_weight( core::scoring::atom_pair_constraint ) );
		runtime_assert( cst_weight > 0.0 );
	}

	// constraints in beta-sheet
	if( constraints_sheet_ > 0.0 ){
		if( blueprint_ ){
			cstset->add_constraints( constraints_sheet( pose, blueprint_, constraints_sheet_ ) );
		}else{
			cstset->add_constraints( constraints_sheet( pose, constraints_sheet_ ) );
		}
	}
	// constraints between N and C
	if( constraints_NtoC_ > 0.0 ){
		cstset->add_constraints( constraints_NtoC( pose, constraints_NtoC_ ) );
	}
	// attach constraints to pose
	pose.add_constraints( cstset->get_all_constraints() );

	// make pose to all ala
	// altert !! we might want to have the functionality to keep_disulfide or not
	if( clear_all_residues_ ) {
		MakePolyXMover  bap( "ALA", false/*kee_pro*/, true /*keep_gly*/, false /*keep_disulfide_cys*/ );
		bap.apply( pose );
	}

	// do _ncycles cycles of fix backbone design.
	for(Size ncycle = 1; ncycle <= ncycles_; ++ncycle) {
	  TR << "cycle = " << ncycle <<  std::endl;
		PackerTaskOP task = task_factory_->create_task_and_apply_taskoperations(pose);

		// relax every relax_step_ number of steps
		if( ! (ncycle % relax_step_)) {
				TR << "relaxing the pose" << std::endl; 
				relax_->apply(pose);
		}

		// run fixbb design
		PackRotamersMover packer( scorefxn_design_, task, 1 );
		packer.apply(pose);
		
	}

} // ConstraintedDesign::apply


/// @brief parse xml
void
ConstrainedDesign::parse_my_tag(
	TagPtr const tag,
	DataMap & data,
	Filters_map const &,
	Movers_map const &,
	Pose const & pose )
{
	using std::string;
	using core::scoring::ScoreFunctionOP;
	string const sfxn_design ( tag->getOption<string>( "sfxn_design", "score12" ));
	string const sfxn_relax  ( tag->getOption<string>( "sfxn_relax", "score12" ) );
	ScoreFunctionOP scorefxn_design = data.get< ScoreFunction * >( "scorefxns", sfxn_design );
	ScoreFunctionOP scorefxn_relax  = data.get< ScoreFunction * >( "scorefxns", sfxn_relax  );
	/// clear amino acid sequence at the very beginning of FlxbbDesign
	clear_all_residues_ = tag->getOption<bool>( "clear_all_residues", 0 );

	// resfile for fixbb
	resfile_ = tag->getOption<string>( "resfile", "" );

	// constraint N- and C- terminal
	constraints_NtoC_ = tag->getOption<Real>( "constraints_NtoC", -1.0 );

	// constraint N- and C- terminal
	constraints_sheet_ = tag->getOption<Real>( "constraints_sheet", -1.0 );

	// get the constriant and load into memory
	string const cset_name  ( tag->getOption<string>( "constraint_set") );
	TR << cset_name << std::endl;
	constraints_ = new SequenceConstraintSet(*data.get< SequenceConstraintSet * >("sequence_constraints", cset_name));

} /// parse_my_tag

} // devel
} // constrained_sequence_design

