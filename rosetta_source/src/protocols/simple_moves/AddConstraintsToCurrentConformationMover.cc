// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#include <protocols/simple_moves/AddConstraintsToCurrentConformationMover.hh>
#include <protocols/simple_moves/AddConstraintsToCurrentConformationMoverCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/SOGFunc.hh>
#include <core/scoring/constraints/ScalarWeightedFunc.hh>

// task operation
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <numeric/xyzVector.hh>

// utility
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

// option
#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

static basic::Tracer TR( "protocols.simple_moves.AddConstraintsToCurrentConformationMover" );

namespace protocols {
namespace simple_moves {

using namespace core;
using namespace basic::options;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;
using namespace constraints;
	
AddConstraintsToCurrentConformationMover::AddConstraintsToCurrentConformationMover(){
	use_distance_cst_ = false;
	max_distance_ = 12.0;
	coord_dev_ = 1.0;
	bound_width_ = 0.;
	min_seq_sep_ = 8;
	cst_weight_ = 1.0;
}

AddConstraintsToCurrentConformationMover::~AddConstraintsToCurrentConformationMover(){}
		
void AddConstraintsToCurrentConformationMover::apply( core::pose::Pose & pose )
{
	using namespace conformation;
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose::datacache;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::id;
	using namespace protocols::moves;
	using namespace core::scoring;
	
	if (!use_distance_cst_) {
		// this is not quite right without adding a virtual residue
		/*
		if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
			core::pose::addVirtualResAsRoot(pose);
		}
		 */
	
	// TR << pose.fold_tree() << std::endl;
	Size nres = pose.total_residue();

	// find anchor residue
	numeric::xyzVector<core::Real> sum_xyz(0.0);
	numeric::xyzVector<core::Real> anchor_xyz(0.0);
	core::Real natom = 0.0;
	for ( Size ires = 1; ires <= nres; ++ires ) {
		if ( pose.residue_type(ires).has("CA") ) {
			Size iatom = pose.residue_type(ires).atom_index("CA");
			sum_xyz += pose.residue(ires).xyz(iatom);
			natom += 1.;
		}
		if (natom > 1e-3) {
			anchor_xyz = sum_xyz / natom;
		}
	}
	core::Real min_dist2 = 1e9;
	Size best_anchor = 0;
	for ( Size ires = 1; ires <= nres; ++ires ) {
		if ( pose.residue_type(ires).has("CA") ) {
			Size iatom = pose.residue_type(ires).atom_index("CA");
			core::Real dist2 = pose.residue(ires).xyz(iatom).distance_squared(anchor_xyz);
			if (dist2 < min_dist2) {
				min_dist2 = dist2;
				best_anchor = ires;
			}
		}
	}
	
	if (best_anchor == 0) return;
	Size best_anchor_atom = pose.residue_type(best_anchor).atom_index("CA");

	//Real const coord_sdev( option[ OptionKeys::relax::coord_cst_stdev ] );
	for ( Size ires = 1; ires <= nres; ++ires ) {
		Size iatom;
		if ( pose.residue_type(ires).has("CA") ) {
			iatom = pose.residue_type(ires).atom_index("CA");
		}
		else if ( pose.residue_type(ires).is_DNA() ) {
			iatom = 0;
		}
		else {
			continue;
		}
		
		if ( bound_width_ < 1e-3 ) {
			if (iatom != 0) {
				pose.add_constraint( new CoordinateConstraint(
															  AtomID(iatom,ires), AtomID(best_anchor_atom,best_anchor), pose.residue(ires).xyz(iatom),
															  new HarmonicFunc( 0.0, coord_dev_ ) ) );
				TR.Debug << "Constraint added to residue " << ires << ", atom " << iatom << std::endl;
			}
			else {
				for (iatom = 1; iatom <= pose.residue_type(ires).nheavyatoms(); ++iatom) {
					pose.add_constraint( new CoordinateConstraint(
																  AtomID(iatom,ires), AtomID(best_anchor_atom,best_anchor), pose.residue(ires).xyz(iatom),
																  new HarmonicFunc( 0.0, coord_dev_ ) ) );
				}
			}
		} else {
			//Real const cst_width( option[ OptionKeys::relax::coord_cst_width ]() );
			if (iatom != 0) {
			pose.add_constraint( new CoordinateConstraint(
								  AtomID(iatom,ires), AtomID(best_anchor_atom,best_anchor), pose.residue(ires).xyz(iatom),
								  new BoundFunc( 0, bound_width_, coord_dev_, "xyz" )) );
			TR << "Constraint added to residue " << ires << ", atom " << iatom << std::endl;
			}
			else {
				for (iatom = 1; iatom <= pose.residue_type(ires).nheavyatoms(); ++iatom) {
				pose.add_constraint( new CoordinateConstraint(
															  AtomID(iatom,ires), AtomID(best_anchor_atom,best_anchor), pose.residue(ires).xyz(iatom),
															  new BoundFunc( 0, bound_width_, coord_dev_, "xyz" )) );
				}				
			}
		}
	}
	}
	else {
		// distance constraints
		for (Size ires=1; ires<=pose.total_residue(); ++ires) {
			if ( ! pose.residue_type(ires).has("CA") ) continue;
			core::Size iatom = pose.residue_type(ires).atom_index("CA");
			
			for (Size jres=ires+min_seq_sep_; jres<=pose.total_residue(); ++jres) {
				if ( ! pose.residue_type(jres).has("CA") ) continue;
				core::Size jatom = pose.residue_type(jres).atom_index("CA");
				
				core::Real dist = pose.residue(ires).xyz(iatom).distance( pose.residue(jres).xyz(jatom) );
				if ( dist <= max_distance_ ) {
					TR.Debug << "Constraint added to residue " << ires << ", atom " << iatom << " and residue " << jres << ", atom " << jatom << " with weight " << cst_weight_ << std::endl;

					pose.add_constraint(
										new core::scoring::constraints::AtomPairConstraint( core::id::AtomID(iatom,ires),
																						   core::id::AtomID(jatom,jres), 
																						   new core::scoring::constraints::ScalarWeightedFunc( cst_weight_, new core::scoring::constraints::SOGFunc( dist, coord_dev_ )  )
																						   )
										);
				}
			}
		}

	}
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
AddConstraintsToCurrentConformationMover::parse_my_tag(
								TagPtr const tag,
								moves::DataMap & datamap,
								Filters_map const & filters,
								moves::Movers_map const & movers,
								Pose const & pose
								)
{
	if ( tag->hasOption("use_distance_cst") ) {
		use_distance_cst_ = tag->getOption<bool>("use_distance_cst");
	}
	if ( tag->hasOption("max_distance") ) {
		max_distance_ = tag->getOption<core::Real>("max_distance");
	}
	if ( tag->hasOption("coord_dev") ) {
		coord_dev_ = tag->getOption<core::Real>("coord_dev");
	}
	if ( tag->hasOption("bound_width") ) {
		bound_width_ = tag->getOption<core::Real>("bound_width");
	}
	if ( tag->hasOption("min_seq_sep") ) {
		min_seq_sep_ = tag->getOption<core::Size>("min_seq_sep");
	}
	if ( tag->hasOption("cst_weight") ) {
		cst_weight_ = tag->getOption<core::Real>("cst_weight");
	}

	parse_task_operations( tag, datamap, filters, movers, pose );
}

void
AddConstraintsToCurrentConformationMover::parse_task_operations(
										 TagPtr const tag,
										 moves::DataMap const & datamap,
										 Filters_map const &,
										 moves::Movers_map const &,
										 Pose const &
										 )
{
	TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == 0) return;
	task_factory( new_task_factory );
}

void AddConstraintsToCurrentConformationMover::task_factory( TaskFactoryCOP tf )
{
	runtime_assert( tf );
	task_factory_ = tf;
}
	
moves::MoverOP AddConstraintsToCurrentConformationMover::clone() const { return new AddConstraintsToCurrentConformationMover( *this ); }
moves::MoverOP AddConstraintsToCurrentConformationMover::fresh_instance() const { return new AddConstraintsToCurrentConformationMover; }

protocols::moves::MoverOP
AddConstraintsToCurrentConformationMoverCreator::create_mover() const {
	return new AddConstraintsToCurrentConformationMover;
}
	
std::string
AddConstraintsToCurrentConformationMoverCreator::keyname() const
{
	return AddConstraintsToCurrentConformationMoverCreator::mover_name();
}
	
std::string
AddConstraintsToCurrentConformationMoverCreator::mover_name()
{
	return "AddConstraintsToCurrentConformationMover";
}
	
std::string
AddConstraintsToCurrentConformationMover::get_name() const {
	return "AddConstraintsToCurrentConformationMover";
}
	
} // moves
} // protocols
