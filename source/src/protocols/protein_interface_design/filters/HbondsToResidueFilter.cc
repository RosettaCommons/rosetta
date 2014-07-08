// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/HbondsToResidueFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)
#include <protocols/protein_interface_design/filters/HbondsToResidueFilter.hh>
#include <protocols/protein_interface_design/filters/HbondsToResidueFilterCreator.hh>


// Project Headers
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
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreType.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
//#include <protocols/moves/ResidueMover.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <basic/MetricValue.hh>
#include <numeric/random/random.hh>
// AUTO-REMOVED #include <boost/foreach.hpp>
#include <core/chemical/AtomType.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

//Objectxxxx header
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>

// Utility Headers

// Unit Headers
#include <protocols/simple_moves/ddG.hh>
#include <protocols/protein_interface_design/design_utils.hh>

// C++ headers
#include <map>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>

//Auto Headers
#include <protocols/simple_filters/DdgFilter.hh>


using namespace core;
using namespace core::scoring;
using namespace ObjexxFCL::format;

static numeric::random::RandomGenerator RG( 140845 ); // <- Magic number, do not change it!!!

namespace protocols {
namespace protein_interface_design {
namespace filters {

static basic::Tracer TR( "protocols.protein_interface_design.filters.HbondsToResidueFilter" );
using core::pose::Pose;

protocols::filters::FilterOP
HbondsToResidueFilterCreator::create_filter() const { return new HbondsToResidueFilter; }

std::string
HbondsToResidueFilterCreator::keyname() const { return "HbondsToResidue"; }

bool
HbondsToResidueFilter::apply( Pose const & pose ) const {
	core::Size hbonded_res( compute( pose ) );
	TR<<"found "<<hbonded_res<< " hbond to target residue " << resnum_;
	if( hbonded_res >= partners_ ) {
		TR << ". passing." << std::endl;
		return( true );
	}
	else {
		TR << ". failing." << std::endl;
		return( false );
	}
}

void
HbondsToResidueFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & pose )
{
	partners_ = tag->getOption<core::Size>( "partners" );
	energy_cutoff_ = tag->getOption<core::Real>( "energy_cutoff", -0.5 );
	bb_bb_ = tag->getOption<bool>( "bb_bb", 0 );
	backbone_ = tag->getOption<bool>( "backbone", 0 );
	sidechain_ = tag->getOption<bool>( "sidechain", 1 );
	resnum_ = core::pose::get_resnum( tag, pose );

	TR<<"Hbonds to residue filter for resnum "<<resnum_<<" with "<<partners_<<" hbonding partners"<<std::endl;
}

void
HbondsToResidueFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Size hbonded_res( compute( pose ) );

	out<<"Number of residues hbonded to "<<resnum_<< " is " << hbonded_res <<'\n';
}

core::Real
HbondsToResidueFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Size hbonded_res( compute( pose ) );
	return( hbonded_res );
}

core::Size
HbondsToResidueFilter::compute( Pose const & pose ) const {
	using core::Size;

	core::pose::Pose temp_pose( pose );
	core::scoring::ScoreFunctionOP scorefxn(get_score_function() );
	(*scorefxn)(temp_pose);
	/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( temp_pose );

	Size const chain2begin( temp_pose.conformation().chain_begin( 2 ) );
	Size partner_begin, partner_end;
	if( resnum_ >= chain2begin ) {
		partner_begin = 1; partner_end = chain2begin-1;
	}
	else {
		partner_begin = chain2begin; partner_end = temp_pose.total_residue();
	}
	std::set<Size> binders;
	for( Size i=partner_begin; i<=partner_end; ++i ) binders.insert( i );

	std::list< Size> hbonded_res( hbonded( temp_pose, resnum_, binders, backbone_, sidechain_, energy_cutoff_, bb_bb_) );

	return( hbonded_res.size() );
}

HbondsToResidueFilter::~HbondsToResidueFilter() {}

}
} // protein_interface_design
} // devel
