// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/matdes/SymUnsatHbondFilter.cc
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
/// @last_modified April 29 2014
/// 
/// @author Jacob Bale (balej@u.washington.edu)

// Unit Headers
#include <protocols/matdes/SymUnsatHbondFilter.hh>
#include <protocols/matdes/SymUnsatHbondFilterCreator.hh>

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
#include <core/pose/PDBInfo.hh>
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
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>

// C++ Headers

namespace protocols {
namespace matdes {

static thread_local basic::Tracer TR( "protocols.matdes.SymUnsatHbondFilter" );


// @brief default constructor
SymUnsatHbondFilter::SymUnsatHbondFilter():
  upper_threshold_( 20 ),
  jump_num_( 1 ),
  sym_dof_names_( "" ),
	verbose_( 0 ),
	write2pdb_( 0 ),
  mode_( "bound_vs_unbound" ),
	compare_to_ref_( 0 ),
  reference_pose_( /* NULL */ )
{}

// @brief constructor with arguments 
SymUnsatHbondFilter::SymUnsatHbondFilter( core::Size const upper_cutoff, core::Size const jump, std::string const sym_dofs, bool verb, bool write, std::string const mode, bool compare_to_ref, core::pose::PoseOP reference_pose ):
	Filter( "SymUnsatHbonds" ),
	upper_threshold_( upper_cutoff ),
	jump_num_( jump ),
	sym_dof_names_( sym_dofs ),
	verbose_( verb ),
	write2pdb_( write ),
	mode_( mode ),
	compare_to_ref_( compare_to_ref ),
	reference_pose_( reference_pose )
{}

// @brief copy constructor
SymUnsatHbondFilter::SymUnsatHbondFilter( SymUnsatHbondFilter const & rval ):
	Super( rval ),
	upper_threshold_( rval.upper_threshold_ ),
	jump_num_( rval.jump_num_ ),
	sym_dof_names_( rval.sym_dof_names_ ),
	verbose_( rval.verbose_ ),
	write2pdb_( rval.write2pdb_ ),
	mode_( rval.mode_ ),
	compare_to_ref_( rval.compare_to_ref_ ),
	reference_pose_( rval.reference_pose_ )
{}

protocols::filters::FilterOP
SymUnsatHbondFilter::fresh_instance() const{
  return protocols::filters::FilterOP( new SymUnsatHbondFilter() );
}

protocols::filters::FilterOP
SymUnsatHbondFilter::clone() const{
  return protocols::filters::FilterOP( new SymUnsatHbondFilter( *this ) );
}

// @brief getters
core::Size SymUnsatHbondFilter::upper_threshold() const { return upper_threshold_; }
core::Size SymUnsatHbondFilter::jump_num() const { return jump_num_; }
std::string SymUnsatHbondFilter::sym_dof_names() const { return sym_dof_names_; }
bool SymUnsatHbondFilter::verbose() const { return verbose_; }
bool SymUnsatHbondFilter::write2pdb() const { return write2pdb_; }
std::string SymUnsatHbondFilter::mode() const { return mode_; }
core::pose::PoseOP SymUnsatHbondFilter::reference_pose() const{ return reference_pose_; }
bool SymUnsatHbondFilter::compare_to_ref() const { return compare_to_ref_; }

// @brief setters
void SymUnsatHbondFilter::upper_threshold( core::Size const upper_cutoff ) { upper_threshold_ = upper_cutoff; }
void SymUnsatHbondFilter::jump_num( core::Size const jump ) { jump_num_ = jump; }
void SymUnsatHbondFilter::sym_dof_names( std::string const sym_dofs ) { sym_dof_names_ = sym_dofs; }
void SymUnsatHbondFilter::verbose( bool const verb ) { verbose_ = verb; }
void SymUnsatHbondFilter::write2pdb( bool const write ) { write2pdb_ = write; }
void SymUnsatHbondFilter::mode( std::string const mode ) { mode_ = mode; }
void SymUnsatHbondFilter::reference_pose( core::pose::PoseOP reference_pose ){ reference_pose_ = reference_pose; }
void SymUnsatHbondFilter::compare_to_ref( bool const compare_to_ref ) { compare_to_ref_ = compare_to_ref; }

// Find residues with buried polar atoms.
core::Real
SymUnsatHbondFilter::compute( core::pose::Pose const & pose, bool const & verb, bool const & write ) const
{
	Size nres_asymmetric_unit, ref_nres_asymmetric_unit;	
	if(core::conformation::symmetry::is_symmetric( pose.conformation() )){
		core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
  	nres_asymmetric_unit = symm_info->num_independent_residues();
	} else {
		nres_asymmetric_unit = pose.n_residue();	
	}

	core::pose::Pose refp;
	if( compare_to_ref_ ) {
		refp = *reference_pose_;
		if(core::conformation::symmetry::is_symmetric( refp.conformation() )){
			core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(refp);
			ref_nres_asymmetric_unit = symm_info->num_independent_residues();
		} else {
			ref_nres_asymmetric_unit = refp.n_residue();	
		}
		if( nres_asymmetric_unit != ref_nres_asymmetric_unit )
    	utility_exit_with_message( "Reference pose and current pose have a different number of residues" );
	}

	core::pose::Pose bound = pose;
  core::pose::Pose unbound = bound;
  core::pose::Pose unbound_refp = refp;

	int sym_aware_jump_id = 0;
	if ( sym_dof_names() != "" ) {
		utility::vector1<std::string> sym_dof_name_list;
		sym_dof_name_list = utility::string_split( sym_dof_names() , ',' );
		for (Size i = 1; i <= sym_dof_name_list.size(); i++) {
			sym_aware_jump_id = core::pose::symmetry::sym_dof_jump_num( pose, sym_dof_name_list[i] );
			protocols::rigid::RigidBodyTransMoverOP translate_unbound( new protocols::rigid::RigidBodyTransMover( unbound, sym_aware_jump_id ) );
			translate_unbound->step_size( 1000.0 );
			translate_unbound->apply( unbound );
		}
	} else {
		sym_aware_jump_id = core::pose::symmetry::get_sym_aware_jump_num( pose, jump_num() );
		protocols::rigid::RigidBodyTransMoverOP translate_unbound( new protocols::rigid::RigidBodyTransMover( unbound, sym_aware_jump_id ) );
		translate_unbound->step_size( 1000.0 );
		translate_unbound->apply( unbound );
	}

	//fpd we need to score pose here
	core::scoring::ScoreFunctionOP scorehbond, scorehbond_refp;
	if (core::conformation::symmetry::is_symmetric( pose.conformation() )) {
		scorehbond = core::scoring::ScoreFunctionOP( new core::scoring::symmetry::SymmetricScoreFunction( ) );
	} else {
		scorehbond = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction( ) );
	}
	if( compare_to_ref_ ) {
		if (core::conformation::symmetry::is_symmetric( refp.conformation() )) {
			scorehbond_refp = core::scoring::ScoreFunctionOP( new core::scoring::symmetry::SymmetricScoreFunction( ) );
		} else {
			scorehbond_refp = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction( ) );
		}
	}

	(*scorehbond)(bound);
  core::pose::metrics::PoseMetricCalculatorOP unsat_calc_bound( new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("default", "default") );
  basic::MetricValue< core::id::AtomID_Map<bool> > bound_Amap;
  unsat_calc_bound->get("atom_bur_unsat", bound_Amap, bound);

	(*scorehbond)(unbound);
  core::pose::metrics::PoseMetricCalculatorOP unsat_calc_unbound( new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("default", "default") );
  basic::MetricValue< core::id::AtomID_Map<bool> > unbound_Amap;
  unsat_calc_unbound->get("atom_bur_unsat", unbound_Amap, unbound);

  basic::MetricValue< core::id::AtomID_Map<bool> > unbound_refp_Amap;
	if( compare_to_ref_ ) {
		(*scorehbond_refp)(unbound_refp);
 		core::pose::metrics::PoseMetricCalculatorOP unsat_calc_unbound_refp( new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("default", "default") );
 		unsat_calc_unbound_refp->get("atom_bur_unsat", unbound_refp_Amap, unbound_refp);
	}
	// Modes: bound_vs_unbound, unbound_design_vs_reference, unbound_mutated_sidechains, all
	// For bound_vs_unbound: Report buriedunsatpolars that exist in bound state, but not in unbound state. Currently no repacking is done in either state.
	// For unbound_design_vs_reference: For bb polars and non-mutated sidechains, compare same atoms design vs native. 
	// For mutated sidechains in the unbound designed state: Simply report any buriedunsatpolars and don't compare back to native.
	// TODO JBB: Add functionality to allow a relax mover to be applied if desired by the user.
	core::id::AtomID_Map<bool> bound_am = bound_Amap.value();
	core::id::AtomID_Map<bool> unbound_am = unbound_Amap.value();
	core::id::AtomID_Map<bool> unbound_refp_am;
	if( compare_to_ref_ ) unbound_refp_am = unbound_refp_Amap.value();
	Size uhb_bound_vs_unbound = 0;
	Size uhb_unbound_design_vs_reference = 0;
	Size uhb_unbound_mutated_sidechains = 0;
	Size uhb_all = 0;
	std::string select_uhb_bound_vs_unbound("select uhb_bound_vs_unbound, (");
	std::string select_uhb_unbound_design_vs_reference("select uhb_unbound_design_vs_reference, (");
	std::string select_uhb_unbound_mutated_sidechains("select uhb_unbound_mutated_sidechains, (");
	std::string select_uhb_all("select uhb_all, (");
	for (Size ir=1; ir<=nres_asymmetric_unit; ir++) {
		core::Size output_resi;
		for (Size ia=1; ia<=bound.residue(ir).nheavyatoms(); ia++) {
			if((mode_ == "all") || (mode_ == "bound_vs_unbound")){
				if (bound_am[core::id::AtomID(ia,ir)] && !unbound_am[core::id::AtomID(ia,ir)]) {
					uhb_all++;
					uhb_bound_vs_unbound++;
					if ( write ) {
						write_to_pdb( pose, "bound_vs_unbound", bound.residue(ir).name3(), ir, bound.residue(ir).atom_name(ia) );
						write_to_pdb( pose, "all", bound.residue(ir).name3(), ir, bound.residue(ir).atom_name(ia) );
					}
					output_resi = ir;
					if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
						output_resi = pose.pdb_info()->number( ir );
					}
					select_uhb_bound_vs_unbound.append("resi " + ObjexxFCL::string_of(output_resi) + " and name " + bound.residue(ir).atom_name(ia) + "+ ");   
					select_uhb_all.append("resi " + ObjexxFCL::string_of(output_resi) + " and name " + bound.residue(ir).atom_name(ia) + "+ ");   
				}
			} 
			if ((mode_ == "all") || (mode_ == "unbound_design_vs_reference")) {
				if ((unbound.residue(ir).name3() == unbound_refp.residue(ir).name3()) || (unbound.residue(ir).atom_is_backbone(ia))) {
					if (unbound_am[core::id::AtomID(ia,ir)] && !unbound_refp_am[core::id::AtomID(ia,ir)]) {
						uhb_all++;
						uhb_unbound_design_vs_reference++;
						if ( write ) {
							write_to_pdb( pose, "unbound_design_vs_reference", bound.residue(ir).name3(), ir, bound.residue(ir).atom_name(ia) );
							write_to_pdb( pose, "all", bound.residue(ir).name3(), ir, bound.residue(ir).atom_name(ia) );
						}
						output_resi = ir;
						if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
							output_resi = pose.pdb_info()->number( ir );
						}
						select_uhb_unbound_design_vs_reference.append("resi " + ObjexxFCL::string_of(output_resi) + " and name " + bound.residue(ir).atom_name(ia) + "+ ");   
						select_uhb_all.append("resi " + ObjexxFCL::string_of(output_resi) + " and name " + bound.residue(ir).atom_name(ia) + "+ ");   
					}
				}
			} 
			if ((mode_ == "all") || (mode_ == "unbound_mutated_sidechains")) {
				if ((!unbound.residue(ir).atom_is_backbone(ia)) && (unbound.residue(ir).name3() != unbound_refp.residue(ir).name3())) {
					if (unbound_am[core::id::AtomID(ia,ir)]) {
						uhb_all++;
						uhb_unbound_mutated_sidechains++;
						if ( write ) {
							write_to_pdb( pose, "unbound_mutated_sidechains", bound.residue(ir).name3(), ir, bound.residue(ir).atom_name(ia) );
							write_to_pdb( pose, "all", bound.residue(ir).name3(), ir, bound.residue(ir).atom_name(ia) );
						}
						output_resi = ir;
						if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
							output_resi = pose.pdb_info()->number( ir );
						}
						select_uhb_unbound_mutated_sidechains.append("resi " + ObjexxFCL::string_of(output_resi) + " and name " + bound.residue(ir).atom_name(ia) + "+ ");   
						select_uhb_all.append("resi " + ObjexxFCL::string_of(output_resi) + " and name " + bound.residue(ir).atom_name(ia) + "+ ");   
					}
				}
			}
		}
	}
	// Write out the info to the tracer and/or pdb.
	if ((mode_ == "all") || (mode_ == "bound_vs_unbound")) {
		select_uhb_bound_vs_unbound.erase(select_uhb_bound_vs_unbound.end()-2,select_uhb_bound_vs_unbound.end());
		if ( verb ) {
			TR << select_uhb_bound_vs_unbound << ") and " << protocols::jd2::JobDistributor::get_instance()->current_output_name() ;
			TR << std::endl;
		}
		if ( write ) { 
			write_pymol_string_to_pdb( select_uhb_bound_vs_unbound ); 
		}
	}
	if ((mode_ == "all") || (mode_ == "unbound_design_vs_reference")) {
		select_uhb_unbound_design_vs_reference.erase(select_uhb_unbound_design_vs_reference.end()-2,select_uhb_unbound_design_vs_reference.end());
		if ( verb ) {
			TR << select_uhb_unbound_design_vs_reference << ") and " << protocols::jd2::JobDistributor::get_instance()->current_output_name() ;
			TR << std::endl;
		}
		if ( write ) { 
			write_pymol_string_to_pdb( select_uhb_unbound_design_vs_reference ); 
		}
	}
	if ((mode_ == "all") || (mode_ == "unbound_mutated_sidechains")) {
		select_uhb_unbound_mutated_sidechains.erase(select_uhb_unbound_mutated_sidechains.end()-2,select_uhb_unbound_mutated_sidechains.end());
		if ( verb ) {
			TR << select_uhb_unbound_mutated_sidechains << ") and " << protocols::jd2::JobDistributor::get_instance()->current_output_name() ;
			TR << std::endl;
		}
		if ( write ) { 
			write_pymol_string_to_pdb( select_uhb_unbound_mutated_sidechains ); 
		}
	}
	if ((mode_ == "all")) {
		select_uhb_all.erase(select_uhb_all.end()-2,select_uhb_all.end());
		if ( verb ) {
			TR << select_uhb_all << ") and " << protocols::jd2::JobDistributor::get_instance()->current_output_name() ;
			TR << std::endl;
		}
		if ( write ) { 
			write_pymol_string_to_pdb( select_uhb_all ); 
		}
	}
	if (mode_ == "bound_vs_unbound") return uhb_bound_vs_unbound;
	else if (mode_ == "unbound_design_vs_reference") return uhb_unbound_design_vs_reference;
	else if (mode_ == "unbound_mutated_sidechains") return uhb_unbound_mutated_sidechains;
	else return uhb_all;
}

void SymUnsatHbondFilter::write_to_pdb( core::pose::Pose const & pose, std::string const mode, std::string const residue_name, core::Size const residue, std::string const atom_name ) const
{

	protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
	std::string filter_name = this->name();
	std::string user_name = this->get_user_defined_name();
	core::Size output_resi = residue;
	if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
		output_resi = pose.pdb_info()->number( residue );
	}
	std::string unsat_pols_string = filter_name + " " + user_name + " " + mode + " mode: " + residue_name + ObjexxFCL::string_of(output_resi) + " " + atom_name ;
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
	core::Real const unsat_hbonds( compute( pose, false, false ) );

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
SymUnsatHbondFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	upper_threshold( tag->getOption<core::Size>( "cutoff", 20 ) );
	jump_num( tag->getOption<core::Size>( "jump", 1 ) );
	sym_dof_names( tag->getOption< std::string >( "sym_dof_names" , "" ) );
	verbose( tag->getOption< bool >( "verbose", 0 ) );
	write2pdb( tag->getOption< bool >("write2pdb", 0) );
	mode( tag->getOption< std::string >( "mode" , "bound_vs_unbound" ) );
	if( (mode_ == "unbound_design_vs_reference") || mode_ == "unbound_mutated_sidechains" || (mode_ == "all") ) {
		compare_to_ref(true);
		if( tag->getOption< bool >( "use_native", false )) {
		  if( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		    std::string const reference_pdb = basic::options::option[ basic::options::OptionKeys::in::file::native ]();
		    core::pose::PoseOP temp_pose( new core::pose::Pose );
		    core::import_pose::pose_from_pdb( *temp_pose, reference_pdb );
		    reference_pose_ = temp_pose;
		  } else {
		    utility_exit_with_message("Native PDB not specified on command line.");
		  }
		} else if ( tag->hasOption("reference_name") ){
		  std::string reference_pose_name = tag->getOption< std::string >( "reference_name", "" );
		  reference_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag,data );
		} else if( tag->hasOption("reference_pdb") ){
		  std::string reference_pdb_filename( tag->getOption< std::string >( "reference_pdb", "" ) );
		  reference_pose_ = core::import_pose::pose_from_pdb( reference_pdb_filename );
		} else {
		  utility_exit_with_message("No valid reference pdb or pose specified for SymUnsatHbondFilter.");
		}
	} else if( mode_ != "bound_vs_unbound" ){
		utility_exit_with_message("Invalid mode specified for SymUnsatHbondFilter; valid options are: bound_vs_unbound, unbound_design_vs_reference, unbound_mutated_sidechains or all");
	}
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
	core::Real const unsat_hbonds( compute( pose, verbose(), write2pdb() ));
	out<<"# unsatisfied hbonds: "<< unsat_hbonds<<'\n';
}

protocols::filters::FilterOP 
SymUnsatHbondFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SymUnsatHbondFilter ); }

std::string 
SymUnsatHbondFilterCreator::keyname() const { return "SymUnsatHbonds"; }

} //namespace matdes
} //namespace protocols
