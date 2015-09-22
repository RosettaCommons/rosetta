// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ResiduesInInterfaceFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#include <protocols/simple_filters/ResiduesInInterfaceFilter.hh>
#include <protocols/simple_filters/ResiduesInInterfaceFilterCreator.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/scoring/Interface.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreType.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <utility/tag/Tag.hh>
//#include <protocols/moves/ResidueMover.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <basic/MetricValue.hh>
#include <numeric/random/random.hh>
#include <core/chemical/AtomType.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>

#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

//Objectxxxx header
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>

// Utility Headers

// Unit Headers
#include <protocols/simple_moves/ddG.hh>
//#include <protocols/protein_interface_design/design_utils.hh>

// C++ headers
#include <map>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>


namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer residues_in_interface_tracer( "protocols.simple_filters.ResiduesInInterfaceFilter" );

protocols::filters::FilterOP
ResiduesInInterfaceFilterCreator::create_filter() const { return protocols::filters::FilterOP( new ResiduesInInterfaceFilter ); }

std::string
ResiduesInInterfaceFilterCreator::keyname() const { return "ResInInterface"; }

ResiduesInInterfaceFilter::~ResiduesInInterfaceFilter() {}

/// @details a utility for docking_filter
core::Size ResiduesInInterfaceFilter::compute( core::pose::Pose const & pose ) const
{
	Size interface_counter( 0 );

	runtime_assert( rb_jump_ <= pose.num_jump() );

	core::pose::Pose temp_pose( pose );
	core::scoring::ScoreFunctionOP scorefxn(
		core::scoring::get_score_function() );
	(*scorefxn)(temp_pose);

	protocols::scoring::Interface interface(rb_jump_);
	interface.calculate( temp_pose );

	for ( Size i = 1; i <= temp_pose.total_residue(); i++ ) {
		if ( !temp_pose.residue(i).is_protein() ) continue;
		if ( interface.is_interface( i ) ) interface_counter++;
	}
	return( interface_counter );
}

void
ResiduesInInterfaceFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	residues_in_interface_threshold_ = tag->getOption<core::Size>( "residues", 20 );
	rb_jump_ = tag->getOption<core::Size>( "jump_number", 1 );

	residues_in_interface_tracer<<"residues in interface filter over jump number " << rb_jump_ << " with threshold "<<residues_in_interface_threshold_<<std::endl;
}

bool
ResiduesInInterfaceFilter::apply( core::pose::Pose const & pose ) const {
	core::Size const interface_res( compute( pose ));
	residues_in_interface_tracer<<"There are "<<interface_res<<" residues in the interface.";
	if ( interface_res <= residues_in_interface_threshold_ ) {
		residues_in_interface_tracer<<" Breaking out."<<std::endl;
		return( false );
	} else residues_in_interface_tracer<<std::endl;
	return( true );
}

void
ResiduesInInterfaceFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Size const interface_res( compute( pose ));
	out<<"residues_in_interface "<<interface_res<<'\n';
}

core::Real
ResiduesInInterfaceFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Size const interface_res( compute( pose ));
	return( (core::Real) interface_res );
}

}
}
