// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <basic/datacache/DataMap.hh>
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

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.GeometryFilter" );

protocols::filters::FilterOP
GeometryFilterCreator::create_filter() const { return protocols::filters::FilterOP( new GeometryFilter ); }

std::string
GeometryFilterCreator::keyname() const { return "Geometry"; }

GeometryFilter::GeometryFilter() :
	filters::Filter( "GeometryFilter" ),
	omega_cutoff_( 165.0 ),
	cart_bonded_cutoff_( 20.0 ),
	filename_( "none" ),
	cst_cutoff_( 10000.0 ),
	start_( 1 ),
	end_( 100000 ),
	count_bad_residues_( false ),
	selector_()
{}

GeometryFilter::~GeometryFilter() = default;

void
GeometryFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	omega_cutoff_ = tag->getOption<core::Real>( "omega", omega_cutoff_ );
	cart_bonded_cutoff_ = tag->getOption<core::Real>( "cart_bonded", cart_bonded_cutoff_ );
	filename_ = tag->getOption< std::string >( "cstfile", filename_ );
	cst_cutoff_ = tag->getOption< core::Real >( "cst_cutoff", cst_cutoff_ );
	start_ = tag->getOption< core::Size>( "start", start_ );
	end_ = tag->getOption< core::Size >( "end", end_ );
	count_bad_residues_ = tag->getOption< bool >( "count_bad_residues", count_bad_residues_ );
	selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, data );
}

bool
GeometryFilter::apply( core::pose::Pose const & pose ) const {
	core::Size const bad_residues = compute( pose );
	bool status = false;
	if ( bad_residues == 0 ) {
		TR << "passing." << std::endl;
		status = true;
	} else TR << "failing." << std::endl;
	return status;
}

void
GeometryFilter::report( std::ostream & /*out*/, core::pose::Pose const & pose ) const {
	/*core::Real const dist =*/ compute( pose );
}

core::Real
GeometryFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Size const bad_residues = compute( pose );
	if ( count_bad_residues_ ) {
		return core::Real( bad_residues );
	} else {
		if ( bad_residues == 0 ) {
			return core::Real( 1.0 );
		} else {
			return core::Real( 0.0 );
		}
	}
}


// take per-residue cart_bonded score for evaluating out-liers
core::Size
GeometryFilter::compute( core::pose::Pose const & pose ) const {
	using namespace core::scoring;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;
	using namespace protocols::simple_filters;
	using namespace protocols::simple_moves;
	using core::scoring::ScoreType;

	core::pose::Pose copy_pose;
	if ( is_symmetric( pose ) ) {
		extract_asymmetric_unit( pose, copy_pose, false);
	} else {
		copy_pose = pose;
	}

	copy_pose.update_residue_neighbors();

	core::Size bad_residues = 0;

	// scoring is necessary for Interface to work reliably
	core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
	scorefxn->set_weight( cart_bonded, 1.0);
	(*scorefxn)(copy_pose);

	// find subset of residues to scan
	core::select::residue_selector::ResidueSubset scan_me;
	if ( selector_ ) {
		scan_me = selector_->apply( copy_pose );
	} else {
		scan_me = core::select::residue_selector::ResidueSubset( copy_pose.total_residue(), true );
	}

	core::Size const start = std::max( start_, Size( 1 ) );
	core::Size const stop = std::min( copy_pose.total_residue(), end_ );
	TR << "Scan residues between " << start << " and " << stop << std::endl;
	for ( Size resnum = start; resnum < stop; resnum++ ) {
		if ( !scan_me[ resnum ] ) continue;

		if ( copy_pose.fold_tree().is_cutpoint( resnum+1 ) || copy_pose.fold_tree().is_jump_point( resnum+1 ) ) continue;
		// TL: skip checking omega values at the end of chains
		//     they are zero and meaningless, but will cause this filter to fail
		if ( copy_pose.chain( resnum ) != copy_pose.chain( resnum+1 ) ) continue;

		// TL: skip if the current residue is not a protein residue, prevents filtering on bogus numbers
		if ( !copy_pose.residue( resnum ).is_protein() ) continue;

		core::Real const weight( (*scorefxn)[ core::scoring::ScoreType( cart_bonded ) ] );
		core::Real const score( copy_pose.energies().residue_total_energies( resnum )[ core::scoring::ScoreType( cart_bonded ) ]);
		core::Real weighted_score = weight * score ;

		// also gets omega values:
		core::Real omega = copy_pose.omega(resnum);

		TR << "residue " << resnum << " name " << copy_pose.residue( resnum ).name3() << " cart_bonded term: " << weighted_score ;
		TR << " omega angle: " << omega << std::endl;

		if ( (std::abs(omega) > 180-omega_cutoff_ && std::abs(omega) < omega_cutoff_ ) && copy_pose.residue( resnum+1 ).name3() == "PRO" )  {
			TR << "omega " << resnum <<" "<< copy_pose.residue( resnum+1 ).name3() << " fail " << std::endl;
			++bad_residues;
		}

		if ( std::abs(omega) < omega_cutoff_ && copy_pose.residue( resnum+1 ).name3() != "PRO" )  {
			TR << "omega " << resnum <<" "<< copy_pose.residue( resnum+1 ).name3() << " fail " << std::endl;
			++bad_residues;
		}

		if ( weighted_score >= cart_bonded_cutoff_ ) {
			TR << "cart_bond " << resnum <<" "<< copy_pose.residue( resnum+1 ).name3() << " fail " << std::endl;
			++bad_residues;
		}
	}

	//check the last residues
	Size resnum = std::min( copy_pose.total_residue(), end_ );
	core::Real const weight( (*scorefxn)[ core::scoring::ScoreType( cart_bonded ) ] );
	core::Real const score( copy_pose.energies().residue_total_energies( resnum )[ core::scoring::ScoreType( cart_bonded ) ]);
	core::Real weighted_score = weight * score ;
	TR << "residue " << resnum << " name " << copy_pose.residue( resnum ).name3() << " cart_bonded term: " << weighted_score ;
	TR << weighted_score  << std::endl;

	if ( filename_ != "none" ) {
		TR << "Evaluate constraint energy?" << std::endl;
		ConstraintSetMoverOP cst_set_mover( new ConstraintSetMover() );
		cst_set_mover->constraint_file( filename_ );
		cst_set_mover->apply( copy_pose );

		// only evaluate atom_pair constraints for now
		scorefxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
		(*scorefxn)( copy_pose );
		ScoreTypeFilter const constraint_filter( scorefxn , atom_pair_constraint, cst_cutoff_ );
		bool CScore(constraint_filter.apply( copy_pose ));
		if ( !CScore ) {
			++bad_residues;
		}
	}

	return bad_residues;
}


}
}
