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
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <basic/MetricValue.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/id/AtomID_Map.hh>
#include <basic/MetricValue.hh>

// C++ Headers

namespace devel {
namespace matdes {

static basic::Tracer TR( "devel.matdes.SymUnsatHbondFilter" );


// @brief default constructor
SymUnsatHbondFilter::SymUnsatHbondFilter():
  upper_threshold_( 20 ),
  jump_num_( 1 ),
  sym_dof_names_( "" ),
	verbose_( 0 ),
	write2pdb_( 0 )
{}

// @brief constructor with arguments 
SymUnsatHbondFilter::SymUnsatHbondFilter( core::Size const upper_cutoff, core::Size const jump, std::string const sym_dofs, bool verb, bool write ):
	Filter( "SymUnsatHbonds" ),
	upper_threshold_( upper_cutoff ),
	jump_num_( jump ),
	sym_dof_names_( sym_dofs ),
	verbose_( verb ),
	write2pdb_( write )
{}

// @brief copy constructor
SymUnsatHbondFilter::SymUnsatHbondFilter( SymUnsatHbondFilter const & rval ):
	Super( rval ),
	upper_threshold_( rval.upper_threshold_ ),
	jump_num_( rval.jump_num_ ),
	sym_dof_names_( rval.sym_dof_names_ ),
	verbose_( rval.verbose_ ),
	write2pdb_( rval.write2pdb_ )
{}

protocols::filters::FilterOP
SymUnsatHbondFilter::fresh_instance() const{
  return new SymUnsatHbondFilter();
}

protocols::filters::FilterOP
SymUnsatHbondFilter::clone() const{
  return new SymUnsatHbondFilter( *this );
}

// @brief getters
core::Size SymUnsatHbondFilter::upper_threshold() const { return upper_threshold_; }
core::Size SymUnsatHbondFilter::jump_num() const { return jump_num_; }
std::string SymUnsatHbondFilter::sym_dof_names() const { return sym_dof_names_; }
bool SymUnsatHbondFilter::verbose() const { return verbose_; }
bool SymUnsatHbondFilter::write2pdb() const { return write2pdb_; }

// @brief setters
void SymUnsatHbondFilter::upper_threshold( core::Size const upper_cutoff ) { upper_threshold_ = upper_cutoff; }
void SymUnsatHbondFilter::jump_num( core::Size const jump ) { jump_num_ = jump; }
void SymUnsatHbondFilter::sym_dof_names( std::string const sym_dofs ) { sym_dof_names_ = sym_dofs; }
void SymUnsatHbondFilter::verbose( bool const verb ) { verbose_ = verb; }
void SymUnsatHbondFilter::write2pdb( bool const write ) { write2pdb_ = write; }

// Find residues with buried polar atoms.
core::Real
SymUnsatHbondFilter::compute( core::pose::Pose const & pose, bool const & verb, bool const & write ) const
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

	int sym_aware_jump_id = 0;
	if ( sym_dof_names() != "" ) {
		utility::vector1<std::string> sym_dof_name_list;
		sym_dof_name_list = utility::string_split( sym_dof_names() , ',' );
		for (Size i = 1; i <= sym_dof_name_list.size(); i++) {
			sym_aware_jump_id = core::pose::symmetry::sym_dof_jump_num( pose, sym_dof_name_list[i] );
			protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( unbound, sym_aware_jump_id ) );
			translate->step_size( 1000.0 );
			translate->apply( unbound );
		}
	} else {
		sym_aware_jump_id = core::pose::symmetry::get_sym_aware_jump_num( pose, jump_num() );
		protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( unbound, sym_aware_jump_id ) );
		translate->step_size( 1000.0 );
		translate->apply( unbound );
	}

	//  Uncomment to verify that symmetric pose is being generated and unbound properly.
	//	bound.dump_pdb("bound.pdb");
	//	unbound.dump_pdb("unbound.pdb");

	//fpd we need to score pose here!!
	core::scoring::ScoreFunctionOP scorehbond;
	if (core::conformation::symmetry::is_symmetric( pose.conformation() )) {
		scorehbond = new core::scoring::symmetry::SymmetricScoreFunction( );
	} else {
		scorehbond = new core::scoring::ScoreFunction( );
	}
	scorehbond->set_weight( core::scoring::hbond_lr_bb, 1.0 );
	scorehbond->set_weight( core::scoring::hbond_sr_bb, 1.0 );
	scorehbond->set_weight( core::scoring::hbond_bb_sc, 1.0 );
	scorehbond->set_weight( core::scoring::hbond_sc, 1.0 );

	(*scorehbond)(bound);
  core::pose::metrics::PoseMetricCalculatorOP unsat_calc_bound = new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("default", "default");
  basic::MetricValue< core::id::AtomID_Map<bool> > bound_Amap;
  unsat_calc_bound->get("atom_bur_unsat", bound_Amap, bound);

	(*scorehbond)(unbound);
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
					select_buried_unsat_polars.append("resi " + ObjexxFCL::string_of(ir) + " and name " + bound.residue(ir).atom_name(ia) + "+ ");   
					if ( verb ) { 
						TR << "buried unsat polar(s): " << bound.residue(ir).name3() << ir << "\t" << bound.residue(ir).atom_name(ia);
					}
					if ( write ) { 
						write_to_pdb( bound.residue(ir).name3(), ir, bound.residue(ir).atom_name(ia) ); 
					}
					flag = 1;
				} else {
					if ( verb ) {
						TR << "," << bound.residue(ir).atom_name(ia);
					}
				}
			}
		}
		if ( verb ) {
			if (flag) { 
				TR << std::endl;
			}
		}
	}
	select_buried_unsat_polars.erase(select_buried_unsat_polars.end()-2,select_buried_unsat_polars.end());
	if ( verb ) {
		TR << select_buried_unsat_polars << ") and " << protocols::jd2::JobDistributor::get_instance()->current_output_name() ;
		TR << std::endl;
	}
	if ( write ) { 
		write_pymol_string_to_pdb( select_buried_unsat_polars ); 
	}
	return buried_unsat_polars;
}

void SymUnsatHbondFilter::write_to_pdb( std::string const residue_name, core::Size const residue, std::string const atom_name ) const
{

	protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
	std::string filter_name = this->name();
	std::string user_name = this->get_user_defined_name();
	std::string unsat_pols_string = filter_name + " " + user_name + ": " + residue_name + ObjexxFCL::string_of(residue) + " " + atom_name ;
	job->add_string(unsat_pols_string);
}

void SymUnsatHbondFilter::write_pymol_string_to_pdb( std::string const pymol_selection ) const
{

	protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
	std::string filter_name = this->name();
	std::string user_name = this->get_user_defined_name();
	std::string pymol_string = filter_name + " " + user_name + ": " + pymol_selection + ") and " + protocols::jd2::JobDistributor::get_instance()->current_output_name();
	job->add_string(pymol_string);
}

bool
SymUnsatHbondFilter::apply( core::pose::Pose const & pose ) const 
{
	core::Real const unsat_hbonds( compute( pose, verbose(), write2pdb() ) );

	TR<<"# unsatisfied hbonds: "<<unsat_hbonds<<". ";
	if( unsat_hbonds <= upper_threshold() ){
		TR<<"passing."<<std::endl;
		return true;
	}
	else {
		TR<<"failing."<<std::endl;
		return false;
	}
}


// @brief parse xml
void
SymUnsatHbondFilter::parse_my_tag( utility::tag::TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	upper_threshold( tag->getOption<core::Size>( "cutoff", 20 ) );
	jump_num( tag->getOption<core::Size>( "jump", 1 ) );
	sym_dof_names( tag->getOption< std::string >( "sym_dof_names" , "" ) );
	verbose( tag->getOption< bool >( "verbose", 0 ) );
	write2pdb( tag->getOption< bool >("write2pdb", 0) );

	TR<<"Buried Unsatisfied Hbond filter over jump number " << jump_num() << " with cutoff " << upper_threshold() << std::endl;
}

void SymUnsatHbondFilter::parse_def( utility::lua::LuaObject const & def,
				utility::lua::LuaObject const & ,
				utility::lua::LuaObject const & ) {
	upper_threshold( def["cutoff"] ? def["cutoff"].to<core::Size>() : 20 );
	jump_num( def["jump"] ? def["jump"].to<core::Size>() : 1 );
	verbose( def["verbose"] ? def["verbose"].to<bool>() : false );
	write2pdb( def["write2pdb"] ? def["write2pdb"].to<bool>() : false );

	TR<<"Buried Unsatisfied Hbond filter over jump number " << jump_num() << " with cutoff " << upper_threshold() << std::endl;
}
core::Real
SymUnsatHbondFilter::report_sm( core::pose::Pose const & pose ) const 
{
	core::Real const unsat_hbonds( compute( pose, false, false ));
	return( unsat_hbonds );
}

void
SymUnsatHbondFilter::report( std::ostream & out, core::pose::Pose const & pose ) const 
{
	core::Real const unsat_hbonds( compute( pose, false, false ));
	out<<"# unsatisfied hbonds: "<< unsat_hbonds<<'\n';
}

protocols::filters::FilterOP 
SymUnsatHbondFilterCreator::create_filter() const { return new SymUnsatHbondFilter; }

std::string 
SymUnsatHbondFilterCreator::keyname() const { return "SymUnsatHbonds"; }

} //namespace matdes
} //namespace devel
