// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/GeometryFilter.cc
/// @brief
/// @author Possu Huang (possu@u.washington.edu)
/// @author Lei Shi (shilei@u.washington.edu)
// Project Headers

#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/types.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <map>
#include <cmath>
#include <numeric/random/random.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/simple_filters/DdgFilter.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_filters/GeometryFilter.hh>
#include <protocols/simple_filters/GeometryFilterCreator.hh>
#include <protocols/simple_moves/ddG.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <string>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.GeometryFilter" );

protocols::filters::FilterOP
GeometryFilterCreator::create_filter() const { return new GeometryFilter; }

std::string
GeometryFilterCreator::keyname() const { return "Geometry"; }

GeometryFilter::~GeometryFilter(){}

void
GeometryFilter::parse_my_tag( utility::tag::TagPtr const tag, moves::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	omega_cutoff_ = tag->getOption<core::Real>( "omega", 165 );
	cart_bonded_cutoff_ = tag->getOption<core::Real>( "cart_bonded", 20 );
	filename_ = tag->getOption< std::string >( "cstfile", "none" );
	cst_cutoff_ = tag->getOption< core::Real >( "cst_cutoff", 10000 );
	start_ = tag->getOption< core::Size>( "start", 1 );
	end_ = tag->getOption< core::Size >( "end", 100000 );
}

bool
GeometryFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const dist( compute( pose ) );
	bool status = false;
	if( dist == 1 ) {
		TR << "passing." << std::endl;
		status = true;
	}
	else TR << "failing." << std::endl;
	return status;
}

void
GeometryFilter::report( std::ostream & /*out*/, core::pose::Pose const & pose ) const {
	/*core::Real const dist =*/ compute( pose );
}

core::Real
GeometryFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const dist( compute( pose ) );
	return( dist );
}


// take per-residue cart_bonded score for evaluating out-liers
core::Real
GeometryFilter::compute( core::pose::Pose const & pose ) const {
	using namespace core::scoring;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;
	using namespace protocols::simple_filters;
	using namespace protocols::simple_moves;
	using core::scoring::ScoreType;

	core::pose::Pose copy_pose;
	if (is_symmetric( pose )){
			extract_asymmetric_unit( pose, copy_pose, false);
	}
	else {
		copy_pose = pose;
	}

	copy_pose.update_residue_neighbors();

	// scoring is necessary for Interface to work reliably
	core::scoring::ScoreFunctionOP scorefxn( core::scoring::getScoreFunction() );
	scorefxn->set_weight( cart_bonded, 1.0);
	(*scorefxn)(copy_pose);

	TR << "Scan residues between " << std::max(start_,Size(1)) << " and " << std::min(copy_pose.total_residue(),end_) << std::endl;
	for ( Size resnum = std::max(start_,Size(1)); resnum < std::min(copy_pose.total_residue(),end_); resnum++ ){

				if ( copy_pose.fold_tree().is_cutpoint( resnum+1 ) || copy_pose.fold_tree().is_jump_point( resnum+1 ) ) continue;

				core::Real const weight( (*scorefxn)[ core::scoring::ScoreType( cart_bonded ) ] );
				core::Real const score( copy_pose.energies().residue_total_energies( resnum )[ core::scoring::ScoreType( cart_bonded ) ]);
				core::Real weighted_score = weight * score ;

				// also gets omega values:
				core::Real omega = copy_pose.omega(resnum);

				TR << "residue " << resnum << " name " << copy_pose.residue( resnum ).name3() << " cart_bonded term: " << weighted_score ;
				TR << " omega angle: " << omega << std::endl;
				
//				copy_pose.fold_tree()
				if ( (std::abs(omega) > 180-omega_cutoff_ && std::abs(omega) < omega_cutoff_ ) && copy_pose.residue( resnum+1 ).name3() == "PRO" )  return(0);
				if (std::abs(omega) < omega_cutoff_ && copy_pose.residue( resnum+1 ).name3() != "PRO" )  return(0);
				if (weighted_score >= cart_bonded_cutoff_ ) return(0);
				//if (weighted_score >= cart_bonded_cutoff_ && ( resnum != 1 || resnum != copy_pose.total_residue() ) ) return(0);
	}

	//check the last residues
	Size resnum=std::min(copy_pose.total_residue(),end_);
  core::Real const weight( (*scorefxn)[ core::scoring::ScoreType( cart_bonded ) ] );
  core::Real const score( copy_pose.energies().residue_total_energies( resnum )[ core::scoring::ScoreType( cart_bonded ) ]);
  core::Real weighted_score = weight * score ;
	TR << "residue " << resnum << " name " << copy_pose.residue( resnum ).name3() << " cart_bonded term: " << weighted_score ;
        TR << " omega angle: NA" << std::endl;
	if (weighted_score >= cart_bonded_cutoff_) return(0);

	if (filename_ != "none"){
		ConstraintSetMoverOP cst_set_mover = new ConstraintSetMover();
		cst_set_mover->constraint_file( filename_ );
		cst_set_mover->apply( copy_pose );

		// only evaluate atom_pair constraints for now
		scorefxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
		(*scorefxn)( copy_pose );
		ScoreTypeFilter const constraint_filter( scorefxn , atom_pair_constraint, cst_cutoff_ );
    bool CScore(constraint_filter.apply( copy_pose ));
    if (!CScore){
			return(0);
		}
	}

	return(1);
}


}
}
