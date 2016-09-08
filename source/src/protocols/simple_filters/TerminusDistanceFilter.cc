// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/TerminusDistanceFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)
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
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/types.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <map>
#include <numeric/random/random.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/simple_filters/DdgFilter.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_filters/TerminusDistanceFilter.hh>
#include <protocols/simple_filters/TerminusDistanceFilterCreator.hh>
#include <protocols/simple_moves/ddG.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <string>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.TerminusDistanceFilter" );

protocols::filters::FilterOP
TerminusDistanceFilterCreator::create_filter() const { return protocols::filters::FilterOP( new TerminusDistanceFilter ); }

std::string
TerminusDistanceFilterCreator::keyname() const { return "TerminusDistance"; }

TerminusDistanceFilter::~TerminusDistanceFilter()= default;

void
TerminusDistanceFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	jump_num_ = tag->getOption<core::Size>( "jump_number", 1 );
	distance_ = tag->getOption<core::Size>( "distance", 5 );

	TR<<"Distance From Terminus filter over jump number " << jump_num_ << " with cutoff " << distance_ << std::endl;
}

bool
TerminusDistanceFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const dist( compute( pose ) );
	TR<<"near terminus: "<<dist<<". " ;
	bool const status = (dist <= distance_) ? (false) : (true);
	if ( status ) TR << "passing." << std::endl;
	else TR << "failing." << std::endl;
	return status;
}

void
TerminusDistanceFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const dist( compute( pose ) );
	out<<"near terminus: "<< dist<<'\n';
}

core::Real
TerminusDistanceFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const dist( compute( pose ) );
	return( dist );
}

core::Real
TerminusDistanceFilter::compute( core::pose::Pose const & pose ) const {
	core::pose::Pose copy_pose = pose;
	runtime_assert( copy_pose.num_jump() >= jump_num_ );

	// scoring is necessary for Interface to work reliably
	core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
	(*scorefxn)(copy_pose);

	protocols::scoring::Interface iface( jump_num_ );
	iface.distance( 8 );
	iface.calculate( copy_pose );
	core::Real min_dist(1000);

	for ( core::Size i=1; i <= pose.size(); ++i ) {
		core::Real dist(1000);
		if ( !iface.is_interface( i ) ) continue; // keep going if we're not at the interface

		core::Size const chain = copy_pose.residue( i ).chain();
		core::Size const N_dist = i - copy_pose.conformation().chain_begin( chain );
		core::Size const C_dist = copy_pose.conformation().chain_end( chain ) - i;
		dist = ( N_dist <= C_dist ) ? ( N_dist) : ( C_dist ) ;
		if ( ( N_dist < distance_ ) || ( C_dist < distance_ ) ) {
			return dist;
		}
		min_dist = ( dist < min_dist ) ? (dist) : (min_dist);
	}
	return( min_dist );
}


}
}
