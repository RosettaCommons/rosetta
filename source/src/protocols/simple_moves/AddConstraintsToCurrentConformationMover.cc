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

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/SOGFunc.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>

// task operation
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <numeric/xyzVector.hh>

// utility
#include <utility/tag/Tag.hh>
#include <utility/fixedsizearray1.hh>
#include <basic/Tracer.hh>

// option
#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.AddConstraintsToCurrentConformationMover" );

namespace protocols {
namespace simple_moves {

using namespace core;
using namespace basic::options;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;
using namespace constraints;

AddConstraintsToCurrentConformationMover::AddConstraintsToCurrentConformationMover() {
	has_task_factory_= false;
	use_distance_cst_ = false;
	CA_only_ = true;
	bb_only_ = false;
	inter_chain_ = true;
	max_distance_ = 12.0;
	coord_dev_ = 1.0;
	bound_width_ = 0.;
	min_seq_sep_ = 8;
	cst_weight_ = 1.0;
}

AddConstraintsToCurrentConformationMover::~AddConstraintsToCurrentConformationMover() {}

bool AddConstraintsToCurrentConformationMover::residue_to_constrain(Size const & i) const {
	if ( has_task_factory_ ) {
		return (!task_->residue_task(i).being_designed()) && (!task_->residue_task(i).being_packed());
	}
	return true;
}

void AddConstraintsToCurrentConformationMover::apply( core::pose::Pose & pose ) {
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

	if ( has_task_factory_ ) {
		task_ = task_factory_->create_task_and_apply_taskoperations( pose );
	}

	if ( !use_distance_cst_ ) {
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
			if ( natom > 1e-3 ) {
				anchor_xyz = sum_xyz / natom;
			}
		}
		core::Real min_dist2 = 1e9;
		Size best_anchor = 0;
		for ( Size ires = 1; ires <= nres; ++ires ) {
			if ( pose.residue_type(ires).has("CA") ) {
				Size iatom = pose.residue_type(ires).atom_index("CA");
				core::Real dist2 = pose.residue(ires).xyz(iatom).distance_squared(anchor_xyz);
				if ( dist2 < min_dist2 ) {
					min_dist2 = dist2;
					best_anchor = ires;
				}
			}
		}

		if ( best_anchor == 0 ) return;
		Size best_anchor_atom = pose.residue_type(best_anchor).atom_index("CA");

		for ( Size ires = 1; ires <= nres; ++ires ) {
			Size iatom_start=1, iatom_stop=pose.residue(ires).nheavyatoms();
			if ( pose.residue_type(ires).is_DNA() ) {
				iatom_stop = 0;
			} else if ( pose.residue_type(ires).is_protein() ) {
				if ( CA_only_ && pose.residue_type(ires).has("CA") ) {
					iatom_start = iatom_stop = pose.residue_type(ires).atom_index("CA");
				} else if ( bb_only_ ) {
					iatom_stop = pose.residue_type(ires).last_backbone_atom();
				}
			} else {
				continue;
			}

			if ( residue_to_constrain(ires) ) {
				for ( Size iatom = iatom_start; iatom <= iatom_stop; ++iatom ) {
					pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new CoordinateConstraint(
						AtomID(iatom,ires), AtomID(best_anchor_atom,best_anchor), pose.residue(ires).xyz(iatom), cc_func_ ) ) ) );
					TR.Debug << "coordinate constraint added to residue " << ires << ", atom " << iatom << std::endl;
				}
			}
		}//loop through residues
	} else { //cartesian
		// distance constraints
		Size nres;
		if ( core::pose::symmetry::is_symmetric(pose) ) {
			core::conformation::symmetry::SymmetryInfoCOP symm_info;
			core::conformation::symmetry::SymmetricConformation & SymmConf (
				dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
			symm_info = SymmConf.Symmetry_Info();
			nres = symm_info->num_independent_residues();
		} else {
			nres=pose.total_residue();
		}

		for ( Size ires=1; ires<=nres; ++ires ) {
			if ( pose.residue(ires).aa() == core::chemical::aa_vrt ) continue;
			if ( !residue_to_constrain(ires) ) continue;
			utility::fixedsizearray1<core::Size,2> iatoms(0); // both are 0
			if ( pose.residue_type(ires).has("CA")  ) iatoms[1] = pose.residue_type(ires).atom_index("CA");
			if ( !CA_only_ ) iatoms[2] = pose.residue_type(ires).nbr_atom();

			for ( Size jres=ires+min_seq_sep_; jres<=pose.total_residue(); ++jres ) {
				if ( pose.residue(jres).aa() == core::chemical::aa_vrt ) continue;
				if ( !inter_chain_ && pose.chain(ires)!=pose.chain(jres) ) continue;
				if ( !residue_to_constrain(jres) ) continue;
				utility::fixedsizearray1<core::Size,2> jatoms(0);
				if ( pose.residue_type(jres).has("CA") )           jatoms[1] = pose.residue_type(jres).atom_index("CA");
				if ( !CA_only_ ) jatoms[2] = pose.residue_type(jres).nbr_atom();

				for ( utility::fixedsizearray1<core::Size,2>::const_iterator iiatom = iatoms.begin(); iiatom != iatoms.end(); ++iiatom ) {
					for ( utility::fixedsizearray1<core::Size,2>::const_iterator jjatom = jatoms.begin(); jjatom != jatoms.end(); ++jjatom ) {
						Size const &iatom(*iiatom), &jatom(*jjatom);
						if ( iatom==0 || jatom==0 ) continue;

						core::Real dist = pose.residue(ires).xyz(iatom).distance( pose.residue(jres).xyz(jatom) );
						if ( dist > max_distance_ ) continue;

						pose.add_constraint(
							scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::AtomPairConstraint( core::id::AtomID(iatom,ires), core::id::AtomID(jatom,jres), core::scoring::func::FuncOP( new core::scoring::func::ScalarWeightedFunc( cst_weight_, core::scoring::func::FuncOP( new core::scoring::func::SOGFunc( dist, coord_dev_ ) ) ) ) ) ) ) );
						TR.Debug << "atom_pair_constraint added to residue " << ires << ", atom " << iatom << " and residue " << jres << ", atom " << jatom << " with weight " << cst_weight_ << std::endl;
					}
				}
			}
		}
	}

}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
AddConstraintsToCurrentConformationMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	Filters_map const & filters,
	moves::Movers_map const & movers,
	Pose const & pose
) {
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
	if ( bound_width_ < 1e-3 ) cc_func_ = core::scoring::func::FuncOP( new core::scoring::func::HarmonicFunc(0.0,coord_dev_) );
	else cc_func_ = core::scoring::func::FuncOP( new BoundFunc(0,bound_width_,coord_dev_,"xyz") );

	if ( tag->hasOption("min_seq_sep") ) {
		min_seq_sep_ = tag->getOption<core::Size>("min_seq_sep");
	}
	if ( tag->hasOption("cst_weight") ) {
		cst_weight_ = tag->getOption<core::Real>("cst_weight");
	}

	if ( tag->hasOption( "task_operations" ) ) {
		TR << "WARNING: task_operations only active for proteins" << std::endl;
		has_task_factory_=true;
		parse_task_operations( tag, datamap, filters, movers, pose );
	}

	if ( tag->hasOption("CA_only") ) {
		CA_only_ = tag->getOption<bool>("CA_only",true);
	}

	if ( tag->hasOption("bb_only") ) {
		bb_only_ = tag->getOption<bool>("bb_only",false);
	}

	if ( tag->hasOption("inter_chain_") ) {
		inter_chain_ = tag->getOption<bool>("inter_chain",true);
	}

}

void
AddConstraintsToCurrentConformationMover::parse_task_operations(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const &
) {
	TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == 0 ) return;
	task_factory( new_task_factory );
}

void AddConstraintsToCurrentConformationMover::task_factory( TaskFactoryCOP tf ) {
	runtime_assert( tf != 0 );
	task_factory_ = tf;
}

moves::MoverOP AddConstraintsToCurrentConformationMover::clone() const {
	return moves::MoverOP( new AddConstraintsToCurrentConformationMover( *this ) );
}
moves::MoverOP AddConstraintsToCurrentConformationMover::fresh_instance() const {
	return moves::MoverOP( new AddConstraintsToCurrentConformationMover );
}

protocols::moves::MoverOP
AddConstraintsToCurrentConformationMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddConstraintsToCurrentConformationMover );
}

std::string
AddConstraintsToCurrentConformationMoverCreator::keyname() const {
	return AddConstraintsToCurrentConformationMoverCreator::mover_name();
}

std::string
AddConstraintsToCurrentConformationMoverCreator::mover_name() {
	return "AddConstraintsToCurrentConformationMover";
}

std::string
AddConstraintsToCurrentConformationMover::get_name() const {
	return "AddConstraintsToCurrentConformationMover";
}

} // moves
} // protocols
