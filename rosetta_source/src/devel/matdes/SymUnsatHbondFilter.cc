// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/matdes/SymUnsatHbondFilter.cc
///
/// @brief Determine the number of unsatified hydrogen bonds at an interface
/// Works with both symmetric and asymmetric poses. 
///
/// @detailed Note: the BuriedUnsatHbondFilter does not work properly for symmetric poses. 
/// Options: “jump”: defaults to 1 and “threshold”: defaults to 20. 
/// Takes the current pose, gets the dofs, uses the RigidBodyTransMover to translate the pose into its unbound state 
/// (Note: does not repack in unbound state), goes through every heavy atom in the asymmetric unit and 
/// finds cases where a polar is considered buried in the bound state, but not in the unbound state. 
/// The output includes the number of unsatisfied hydrogen bonds, the specific residues and atoms that 
/// are unsatisfied, and a formatted string for easy selection in pymol.
/// 
/// example usage: <SymUnsatHbonds name=uhb jump=1 cutoff=20 />
///
/// @last_modified April 21 2012
/// 
/// @author Jacob Bale (balej@u.washington.edu)

// Unit Headers
#include <devel/matdes/SymUnsatHbondFilter.hh>
#include <devel/matdes/SymUnsatHbondFilterCreator.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <ObjexxFCL/format.hh>
#include <core/pose/symmetry/util.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <basic/MetricValue.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/id/AtomID_Map.hh>
#include <basic/MetricValue.hh>

// C++ Headers

namespace devel {
namespace matdes {

static basic::Tracer TR( "devel.matdes.SymUnsatHbondFilter" );

protocols::filters::FilterOP
SymUnsatHbondFilterCreator::create_filter() const 
{ 
return new SymUnsatHbondFilter; 
}

std::string
SymUnsatHbondFilterCreator::keyname() const { 
return "SymUnsatHbonds"; 
}

SymUnsatHbondFilter::SymUnsatHbondFilter( core::Size const upper_threshold, core::Size const jump_num ):
	Filter( "SymUnsatHbonds" ),
	upper_threshold_( upper_threshold ),
	jump_num_( jump_num )
{}

SymUnsatHbondFilter::~SymUnsatHbondFilter() {}

// Find residues with buried polar atoms.
core::Real
SymUnsatHbondFilter::compute( core::pose::Pose const & pose ) const
{
	Size nres_asymmetric_unit;	
	if(core::conformation::symmetry::is_symmetric( pose.conformation() )){
		core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
  	nres_asymmetric_unit = symm_info->num_independent_residues();
	} else {
	nres_asymmetric_unit = pose.n_residue();	
	}

	core::pose::Pose bound = pose;
  core::pose::Pose unbound = bound;

	int sym_aware_jump_id = core::pose::symmetry::get_sym_aware_jump_num( unbound, jump_num_ ); // JB 120420
	protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( unbound, sym_aware_jump_id ) );
	translate->step_size( 1000.0 );
	translate->apply( unbound );

//  Uncomment to verify that symmetric pose is being generated and unbound properly.
//	bound.dump_pdb("bound.pdb");
//	unbound.dump_pdb("unbound.pdb");

  core::pose::metrics::PoseMetricCalculatorOP unsat_calc_bound = new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("default", "default");
  basic::MetricValue< core::id::AtomID_Map<bool> > bound_Amap;
  unsat_calc_bound->get("atom_bur_unsat", bound_Amap, bound);

  core::pose::metrics::PoseMetricCalculatorOP unsat_calc_unbound = new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("default", "default");
  basic::MetricValue< core::id::AtomID_Map<bool> > unbound_Amap;
  unsat_calc_unbound->get("atom_bur_unsat", unbound_Amap, unbound);

	core::id::AtomID_Map<bool> bound_am = bound_Amap.value();
	core::id::AtomID_Map<bool> unbound_am = unbound_Amap.value();
	Size buried_unsat_polars = 0;
	std::string select_buried_unsat_polars("select buried_unsat_polars, (");
	for (Size ir=1; ir<=nres_asymmetric_unit; ir++) {
		Size flag = 0;
		for (Size ia=1; ia<=bound.residue(ir).nheavyatoms(); ia++) {
			if (bound_am[core::id::AtomID(ia,ir)] != unbound_am[core::id::AtomID(ia,ir)]) {
				buried_unsat_polars++;
				if (flag == 0) {
					TR << "buried unsat polar(s): " << bound.residue(ir).name3() << ir << "\t" << bound.residue(ir).atom_name(ia);
					select_buried_unsat_polars.append("resi " + ObjexxFCL::string_of(ir) + " and name " + bound.residue(ir).atom_name(ia) + "+ ");   
					flag = 1;
				} else {
					TR << "," << bound.residue(ir).atom_name(ia);
				}
			}
		}
		if (flag) TR << std::endl;
	}
	select_buried_unsat_polars.erase(select_buried_unsat_polars.end()-2,select_buried_unsat_polars.end());
	TR << select_buried_unsat_polars << ") and ";
	TR << std::endl;
	return buried_unsat_polars;
}


bool
SymUnsatHbondFilter::apply( core::pose::Pose const & pose ) const 
{
	core::Real const unsat_hbonds( compute( pose ) );

	TR<<"# unsatisfied hbonds: "<<unsat_hbonds<<". ";
	if( unsat_hbonds <= upper_threshold_ ){
		TR<<"passing."<<std::endl;
		return true;
	}
	else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

void
SymUnsatHbondFilter::report( std::ostream & out, core::pose::Pose const & pose ) const 
{
	core::Real const unsat_hbonds( compute( pose ));
	out<<"# unsatisfied hbonds: "<< unsat_hbonds<<'\n';
}

core::Real
SymUnsatHbondFilter::report_sm( core::pose::Pose const & pose ) const 
{
	core::Real const unsat_hbonds( compute( pose ));
	return( unsat_hbonds );
}

void
SymUnsatHbondFilter::parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap & /*datamap*/, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	jump_num_ = tag->getOption<core::Size>( "jump", 1 );
	upper_threshold_ = tag->getOption<core::Size>( "cutoff", 20 );

	TR<<"Buried Unsatisfied Hbond filter over jump number " << jump_num_ << " with cutoff " << upper_threshold_ << std::endl;
}

} //namespace matdes
} //namespace devel
