// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/BuriedUnsatHbondFilter.cc
/// @brief filter to report buried polar atoms that are not satisfied by hydrogen bonds
/// @details significantly updated in 2017 to have more generous definition of h-bonds
///  (previously, legit h-bonds were excluded because of sfxn exceptions); users can now choose
///  between legacy SASA and VSASA for burial; poses with more than 3 chains now supported; the
///  way unsats are counted and reported is now different (before, all unsats were counted as equal,
///  which is misleading); users can choose different reporting behaviours; the Filter is now Symmetry
///  compatible, and users can pass sym_dof_names as in Jacob Bale's earlier modified filter; users can
///  still use ddG-style behavior if desired; legacy options can be restored by setting legacy=true,
///  but this is only recommended for benchmarking purposes
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)
/// @author Scott Boyken (sboyken@gmail.com), major updates and refactoring

#include <protocols/simple_filters/BuriedUnsatHbondFilter.hh>
#include <protocols/simple_filters/BuriedUnsatHbondFilterCreator.hh>

#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/MetricValue.hh>
#include <protocols/simple_pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/scoring/InterfaceInfo.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>
#include <basic/options/keys/bunsat_calc2.OptionKeys.gen.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer buried_unsat_hbond_filter_tracer( "protocols.simple_filters.BuriedUnsatHbondFilter" );

BuriedUnsatHbondFilter::BuriedUnsatHbondFilter() :
	Filter( "BuriedUnsatHbonds" ),
	name_of_sasa_calc_( "default" ),
	sym_dof_names_( "" ),
	legacy_options_( false ),
	generous_hbonds_( true ),
	legacy_counting_( false ),
	use_vsasa_( true ),
	ignore_surface_res_( false ),
	ignore_bb_heavy_unsats_( false ),
	use_sc_neighbors_( false ),
	only_interface_( false ),
	use_ddG_style_( false ),
	print_out_info_to_pdb_( false ),
	use_reporter_behavior_( true ), // default is now report all heavy unsats
	use_hbnet_behavior_( false ),
	report_all_unsats_( false ), // options from TJ
	report_all_heavy_atom_unsats_( false ), // IF ALL REORTER OPTIONS ARE FALSE, THIS BECOMES TRUE AND DEFAULT REPORTS ALL HEAVY UNSATS
	report_sc_heavy_atom_unsats_( false ),
	report_bb_heavy_atom_unsats_( false ),
	report_nonheavy_unsats_( false ),
	probe_radius_( basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius] ),
	burial_cutoff_( basic::options::option[basic::options::OptionKeys::pose_metrics::atomic_burial_cutoff] ),
	residue_surface_cutoff_( 45.0 ),
	upper_threshold_( 20 ),
	jump_num_( 1 ),
	residue_selector_( /* NULL */ ),
	sfxn_( /* NULL */ ),
	task_factory_( /* NULL */ )
{
	if ( use_vsasa_ ) {
		burial_cutoff_ = basic::options::option[basic::options::OptionKeys::bunsat_calc2::sasa_burial_cutoff];
		residue_surface_cutoff_ = 20.0;
	}
	sfxn_ = core::scoring::get_score_function();
}

// if calling from code, specifying reporter options using setter functions, if desired
BuriedUnsatHbondFilter::BuriedUnsatHbondFilter( core::Size const upper_threshold,
	core::select::residue_selector::ResidueSelectorOP residue_selector /* NULL */,
	bool const use_legacy_options /* false */,
	bool const generous_hbonds /* true */,
	bool const use_vsasa /* true */,
	bool const ignore_surface_res /* false */,
	bool const ignore_bb_heavy_unsats /* false */,
	bool const use_sc_neighbors /* false */,
	bool const only_interface /* false */,
	bool const use_ddG_style /* false */,
	bool const use_hbnet_behavior /* false */
) :
	Filter( "BuriedUnsatHbonds" ),
	name_of_sasa_calc_( "default" ),
	sym_dof_names_( "" ),
	legacy_options_( use_legacy_options ),
	generous_hbonds_( generous_hbonds ),
	legacy_counting_( false ),
	use_vsasa_( use_vsasa ),
	ignore_surface_res_( ignore_surface_res ),
	ignore_bb_heavy_unsats_( ignore_bb_heavy_unsats ),
	use_sc_neighbors_( use_sc_neighbors ),
	only_interface_( only_interface ),
	use_ddG_style_( use_ddG_style ),
	print_out_info_to_pdb_( false ),
	use_reporter_behavior_( true ), // default is now report all heavy unsats
	use_hbnet_behavior_( use_hbnet_behavior ),
	report_all_unsats_( false ), // options from TJ
	report_all_heavy_atom_unsats_( false ), // IF ALL REORTER OPTIONS ARE FALSE, THIS BECOMES TRUE AND DEFAULT REPORTS ALL HEAVY UNSATS
	report_sc_heavy_atom_unsats_( false ),
	report_bb_heavy_atom_unsats_( false ),
	report_nonheavy_unsats_( false ),
	probe_radius_( basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius] ),
	burial_cutoff_( basic::options::option[basic::options::OptionKeys::pose_metrics::atomic_burial_cutoff] ),
	residue_surface_cutoff_( 45.0 ),
	upper_threshold_( upper_threshold ),
	jump_num_( 1 ),
	residue_selector_( residue_selector ),
	sfxn_( /* NULL */ ),
	task_factory_( /* NULL */ )
{
	if ( use_vsasa ) {
		burial_cutoff_ = basic::options::option[basic::options::OptionKeys::bunsat_calc2::sasa_burial_cutoff];
		residue_surface_cutoff_ = 20.0;
	}
	sfxn_ = core::scoring::get_score_function();
	if ( legacy_options_ ) {
		set_legacy_options();
	}
}

///@brief copy constructor
BuriedUnsatHbondFilter::BuriedUnsatHbondFilter( BuriedUnsatHbondFilter const & rval ):
	Filter( rval ),
	name_of_sasa_calc_( rval.name_of_sasa_calc_ ),
	sym_dof_names_( rval.sym_dof_names_ ),
	legacy_options_( rval.legacy_options_ ),
	generous_hbonds_( rval.generous_hbonds_ ),
	legacy_counting_( rval.legacy_counting_ ),
	use_vsasa_( rval.use_vsasa_ ),
	ignore_surface_res_( rval.ignore_surface_res_ ),
	ignore_bb_heavy_unsats_( rval.ignore_bb_heavy_unsats_ ),
	use_sc_neighbors_( rval.use_sc_neighbors_ ),
	only_interface_( rval.only_interface_ ),
	use_ddG_style_( rval.use_ddG_style_ ),
	print_out_info_to_pdb_( rval.print_out_info_to_pdb_ ),
	use_reporter_behavior_( rval.use_reporter_behavior_ ), // options from TJ
	use_hbnet_behavior_( rval.use_hbnet_behavior_ ),
	report_all_unsats_( rval.report_all_unsats_ ),
	report_all_heavy_atom_unsats_( rval.report_all_heavy_atom_unsats_ ),
	report_sc_heavy_atom_unsats_( rval.report_sc_heavy_atom_unsats_ ),
	report_bb_heavy_atom_unsats_( rval.report_bb_heavy_atom_unsats_ ),
	report_nonheavy_unsats_( rval.report_nonheavy_unsats_ ),
	probe_radius_( rval.probe_radius_ ),
	burial_cutoff_( rval.burial_cutoff_ ),
	residue_surface_cutoff_( rval.residue_surface_cutoff_ ),
	upper_threshold_( rval.upper_threshold_ ),
	jump_num_( rval.jump_num_ ),
	residue_selector_( rval.residue_selector_ ),
	sfxn_( rval.sfxn_ ),
	task_factory_( rval.task_factory_ )
{}

BuriedUnsatHbondFilter::~BuriedUnsatHbondFilter()= default;

void
BuriedUnsatHbondFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & pose )
{
	legacy_options_ = tag->getOption<bool>( "use_legacy_options", false );
	generous_hbonds_ = tag->getOption<bool>( "generous_hbonds", true );
	use_vsasa_ = tag->getOption<bool>( "use_vsasa", true );
	ignore_surface_res_ = tag->getOption<bool>( "ignore_surface_res", false );
	ignore_bb_heavy_unsats_ = tag->getOption<bool>( "ignore_bb_heavy_unsats", false );
	use_sc_neighbors_ = tag->getOption<bool>( "use_sc_neighbors", false );
	use_ddG_style_ = tag->getOption<bool>( "use_ddG_style", false );
	print_out_info_to_pdb_ = tag->getOption<bool>( "print_out_info_to_pdb", false );
	only_interface_ = tag->getOption<bool>( "only_interface", false );

	if ( tag->hasOption( "probe_radius" ) ) name_of_sasa_calc_="nondefault"; // ensure that probe radius gets updated in Unsat Calc and SASA Calc

	probe_radius_ = tag->getOption<core::Real>( "probe_radius", basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius] ); // default is 1.4
	burial_cutoff_ = tag->getOption<core::Real>( "burial_cutoff", basic::options::option[basic::options::OptionKeys::pose_metrics::atomic_burial_cutoff] ); // default is 0.3
	residue_surface_cutoff_ = tag->getOption<core::Real>( "residue_surface_cutoff", 45.0 );
	jump_num_ = tag->getOption<core::Size>( "jump_number", 1 );
	upper_threshold_ = tag->getOption<core::Size>( "cutoff", 20 );

	// if user does not specificy and vsasa=true, set vsasa default for burial cutoff
	if ( use_vsasa_ ) {
		if ( !( tag->hasOption("burial_cutoff") ) ) burial_cutoff_ = basic::options::option[basic::options::OptionKeys::bunsat_calc2::sasa_burial_cutoff];
		if ( !( tag->hasOption("residue_surface_cutoff") ) ) residue_surface_cutoff_ = 20.0;
	}

	if ( use_sc_neighbors_ ) {
		// SASA uses atom cutoff; us_sc_neighbor uses residue cutoff
		// TODO need better solution but this catches cases where users forgot to update cutoff for sc_neighbor (rather than sasa) values
		if ( burial_cutoff_ < 1.0 ) burial_cutoff_ = 4.4;
		if ( residue_surface_cutoff_ > 7.0 ) burial_cutoff_ = 2.0;
	}

	std::string const scorefxn_key( protocols::rosetta_scripts::get_score_function_name(tag) );
	if ( datamap.has( "scorefxns", scorefxn_key ) ) {
		sfxn_ = datamap.get_ptr<core::scoring::ScoreFunction>( "scorefxns", scorefxn_key );
	} else {
		sfxn_ = core::scoring::get_score_function();
	}
	if ( tag->hasOption("residue_selector") ) {

		std::string const selector_name ( tag->getOption< std::string >( "residue_selector" ) );
		core::select::residue_selector::ResidueSelectorCOP selector;
		try {
			selector = datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name );
		} catch ( utility::excn::Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from AddCompositionConstraintMover::parse_tag()\n" + e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message );
		}
		runtime_assert( selector );
		residue_selector_ = selector->clone();
	}
	if ( tag->hasOption("task_operations") ) {
		buried_unsat_hbond_filter_tracer << " Will use task operations to tell the filter where to look for unsats. " << std::endl;
		buried_unsat_hbond_filter_tracer << " WARNING: it is now recommended to use residue_selector option instead of task_operations" << std::endl;
		task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	}

	// options from TJ Brunette for customizing reporter behavior
	// TODO Need better way of ensuring that is one reporter behavior is set to true, all others become false;
	use_reporter_behavior_ = tag->getOption<bool>( "use_reporter_behavior", true );
	use_hbnet_behavior_ = tag->getOption<bool>( "use_hbnet_behavior", false );
	if ( use_hbnet_behavior_ ) use_reporter_behavior_ = false;
	report_all_heavy_atom_unsats_ = tag->getOption<bool>( "report_all_heavy_atom_unsats", false ); // IF ALL ARE FALSE, DEFAULT IS THIS BECOMES TRUE
	report_all_unsats_ = tag->getOption<bool>( "report_all_unsats", false );
	report_sc_heavy_atom_unsats_ = tag->getOption<bool>( "report_sc_heavy_atom_unsats", false );
	report_bb_heavy_atom_unsats_ = tag->getOption<bool>( "report_bb_heavy_atom_unsats", false );
	report_nonheavy_unsats_ = tag->getOption<bool>( "report_nonheavy_unsats", false );

	if ( use_reporter_behavior_ && !report_all_heavy_atom_unsats_ && !report_sc_heavy_atom_unsats_ && !report_bb_heavy_atom_unsats_ && !report_nonheavy_unsats_ && !report_all_unsats_ ) {
		buried_unsat_hbond_filter_tracer << " WARNING! use_reporter_behavior=true, need to set one behavior to true; will use default report_all_heavy_atom_unsats: " << std::endl;
	}
	if ( use_reporter_behavior_ && ( report_sc_heavy_atom_unsats_ || report_bb_heavy_atom_unsats_ || report_nonheavy_unsats_ || report_all_unsats_ ) ) {
		report_all_heavy_atom_unsats_ = false;
	}

	sym_dof_names( tag->getOption< std::string >( "sym_dof_names" , "" ) );
	if ( sym_dof_names_ != "" ) {
		buried_unsat_hbond_filter_tracer << " you set sym_dof_names, which means will use ddG-style calculation for unsats; if you do not want this, pass a residue selector instead of defining symdofs" << std::endl;
		use_ddG_style_ = true;
	}
	if ( use_ddG_style_ && core::conformation::symmetry::is_symmetric( pose.conformation() ) && sym_dof_names_ == "" ) {
		buried_unsat_hbond_filter_tracer << " WARNING! using ddG_style, and your Pose is symmetric but you didn't define your sym_dof_names?  Is that want you want? setting only_interface=true to be safe" << std::endl;
		only_interface_ = true;
	}
	if ( legacy_options_ ) {
		set_legacy_options();
	}
}

bool
BuriedUnsatHbondFilter::apply( core::pose::Pose const & pose ) const {

	core::Real const unsat_hbonds( compute( pose ) );

	buried_unsat_hbond_filter_tracer<<"# unsatisfied hbonds: "<<unsat_hbonds<<". ";
	if ( unsat_hbonds <= upper_threshold_ ) {
		buried_unsat_hbond_filter_tracer<<"passing." <<std::endl;
		return true;
	} else {
		buried_unsat_hbond_filter_tracer<<"failing."<<std::endl;
		return false;
	}
}

void
BuriedUnsatHbondFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const unsat_hbonds( compute( pose ));
	out<<"# unsatisfied hbonds: "<< unsat_hbonds<<'\n';
}

core::Real
BuriedUnsatHbondFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const unsat_hbonds( compute( pose ));
	return( unsat_hbonds );
}

core::Real
BuriedUnsatHbondFilter::compute( core::pose::Pose const & pose ) const {

	buried_unsat_hbond_filter_tracer << "/////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
	if ( legacy_options_ ) {
		buried_unsat_hbond_filter_tracer << "    WARNING: USING LEGACY OPTIONS; ALL UNSATS TREATED AS EQUIVALENT: only recommended for benchmarking" << std::endl;
	} else if ( use_hbnet_behavior_ ) {
		buried_unsat_hbond_filter_tracer << "    USING HBNet BEHAVIOR: no heavy-atom donor/acc unsats allowed (will return 9999); if no heavy unsats, counts Hpol buried unsats" << std::endl;
	} else if ( use_reporter_behavior_ ) {
		if ( report_all_unsats_ || report_sc_heavy_atom_unsats_ ||  report_bb_heavy_atom_unsats_ || report_nonheavy_unsats_ ) {
			buried_unsat_hbond_filter_tracer << "    USER HAS SPECIFIED CUSTOM REPORTING BEHAVIOR FOR # UNSATS" << std::endl;
		} else {
			buried_unsat_hbond_filter_tracer << "    USING DEFAULT BEHAVIOR: filter will report total number of heavy-atom donor/acceptor buried unsats" << std::endl;
		}
	}
	buried_unsat_hbond_filter_tracer << std::endl << "/////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;

	Size total_res( pose.size() );

	// CHECK FOR SYMMETRY: if symmetric, only count within ASU, as is standard throughout most of Rosetta
	//   h-bonds across the symmetric interface are still considered and can satisfy ASU polar atoms
	bool symmetric( ( core::conformation::symmetry::is_symmetric( pose.conformation() ) ) ? true : false );
	if ( symmetric ) {
		buried_unsat_hbond_filter_tracer << " DETECTED THAT POSE IS SYMMETRIC:  " << std::endl;
		buried_unsat_hbond_filter_tracer << "    if symmetric pose and only_interface=true (default for symmetric case), then will only look at symmetric interface residues " << std::endl;
		buried_unsat_hbond_filter_tracer << "    if symmetric pose and only_interface=false (set explicitly in XML), then will add up all unsats in ASU " << std::endl;
		core::conformation::symmetry::SymmetricConformation const & SymmConf(dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation()));
		core::conformation::symmetry::SymmetryInfoCOP symm_info = SymmConf.Symmetry_Info();
		total_res = symm_info->num_independent_residues();
	}
	std::set< core::Size > region_to_calculate = std::set< core::Size >(); /* NULL */

	// if users want to specificy task operations to tell the filter where to look (legacy behavior, now recommended to use residue_selector)
	if ( task_factory() != NULL ) {
		buried_unsat_hbond_filter_tracer << " LOOKING FOR UNSATS ONLY AT RESIDUES DEFINED BY YOUR task_operations: " << std::endl;
		buried_unsat_hbond_filter_tracer << "   NOTE: it's now recommonded to use residue_selector option instead of this option" << std::endl;
		utility::vector1< core::Size > const selected_residues( protocols::rosetta_scripts::residue_packer_states( pose, task_factory(), true/*designable*/, true/*packable*/ ) );

		for ( core::Size const sr : selected_residues ) {
			region_to_calculate.insert(sr);
		}
	} else if ( residue_selector_ ) {
		buried_unsat_hbond_filter_tracer << " LOOKING FOR UNSATS ONLY AT RESIDUES DEFINED BY YOUR residue_selector: " << std::endl;
		core::select::residue_selector::ResidueSubset selection = residue_selector_->apply( pose );
		for ( Size r = 1; r <= selection.size(); ++r ) {
			if ( selection[ r ] ) {
				region_to_calculate.insert(r);
			}
		}
	} else if ( only_interface_ && pose.num_chains() > 1 ) { // if more than 1 chain, focus on interface

		buried_unsat_hbond_filter_tracer << " LOOKING FOR UNSATS ONLY AT INTERFACE RESIDUES: " << std::endl;

		utility::vector1< core::Size > jump_nums;

		// if symmetric pose and only_interface=true, then will only look at symmetric interface residues
		// if symmetric pose and only_interface=false, then will add up all unsats in ASU
		if ( symmetric ) {
			utility::vector1<std::string> sym_dof_name_list = core::pose::symmetry::sym_dof_names( pose );
			for ( Size i = 1; i <= sym_dof_name_list.size(); i++ ) {
				Size sym_aware_jump_id(core::pose::symmetry::sym_dof_jump_num(pose,sym_dof_name_list[i]));
				jump_nums.push_back(sym_aware_jump_id);
			}
		} else {
			for ( Size i=1; i <= pose.num_jump(); ++i ) {
				jump_nums.push_back(i);
			}
		}
		for ( const auto & j : jump_nums ) {
			protocols::scoring::InterfaceOP interf = protocols::scoring::InterfaceOP( new protocols::scoring::Interface( j ) );
			interf->distance( 8.0 );
			interf->calculate( pose ); //selects residues: sq dist of nbr atoms < interf_dist_sq (8^2)

			for ( Size resnum = 1; resnum <= total_res; resnum++ ) {
				if ( interf->is_interface( resnum ) ) {
					region_to_calculate.insert( resnum );
				}
			}
		}
	} else { // DEFAULT, total_res (ASU for symmetric case)
		for ( core::Size resnum=1; resnum <= total_res; resnum++ ) {
			region_to_calculate.insert( resnum );
		}
	}
	core::pose::Pose input_pose( pose );
	( *sfxn_ )( input_pose ); // safeguard to ensure Pose is scored before counting h-bonds and checking unsats

	using namespace protocols::simple_pose_metric_calculators;

	std::string name_of_hbond_calc = ( generous_hbonds_ ) ? "default" : "legacy";

	BuriedUnsatisfiedPolarsCalculator calc_bound( name_of_sasa_calc_, name_of_hbond_calc, region_to_calculate, burial_cutoff_, probe_radius_, residue_surface_cutoff_, generous_hbonds_, legacy_counting_, use_vsasa_, use_sc_neighbors_, ignore_surface_res_ );
	basic::MetricValue< core::Size > mv_all_heavy, mv_bb_heavy, mv_countable_nonheavy, mv_all_unsat;
	basic::MetricValue< core::id::AtomID_Map< bool > > mv_unsat_map;
	basic::MetricValue< core::id::AtomID_Map< bool > > mv_unbound_unsat_map;

	calc_bound.get("all_heavy_unsats", mv_all_heavy, input_pose);
	calc_bound.get("bb_heavy_unsats", mv_bb_heavy, input_pose);
	calc_bound.get("countable_nonheavy_unsats", mv_countable_nonheavy, input_pose);
	calc_bound.get("atom_bur_unsat", mv_unsat_map, input_pose);
	calc_bound.get("all_bur_unsat_polars", mv_all_unsat, input_pose);

	core::Size all_heavy_atom_unsats( mv_all_heavy.value() );
	core::Size bb_heavy_atom_unsats( mv_bb_heavy.value() );
	core::Size sc_heavy_atom_unsats( mv_all_heavy.value() - mv_bb_heavy.value() );
	core::Size countable_nonheavy_unsats( mv_countable_nonheavy.value() );
	core::Size all_unsats( mv_all_unsat.value() );

	buried_unsat_hbond_filter_tracer << "  buried unsats in input pose: " << std::endl;
	if ( legacy_counting_ ) {
		buried_unsat_hbond_filter_tracer << "  all_unsats = " << all_unsats << std::endl;
	} else {
		buried_unsat_hbond_filter_tracer << "  all_heavy_atom_unsats = " << all_heavy_atom_unsats << std::endl << "  bb_heavy_atom_unsats = " << bb_heavy_atom_unsats << std::endl << "  sc_heavy_atom_unsats = " << sc_heavy_atom_unsats << std::endl << "  countable_nonheavy_unsats = " << countable_nonheavy_unsats << std::endl;
	}
	// DDG style separate and compare bound to unbound
	bool ddG_was_computed( false );
	if ( use_ddG_style_ && jump_num_ > 0 && pose.num_chains() > 1 ) { // UPDATED TO USE sym_dofs for Symmetry like the SymUnsatFilter does

		buried_unsat_hbond_filter_tracer << " use_ddG_style=true: Using ddG style calculation ( will substract unsats in unbound state from those in bound state ): " << std::endl;
		runtime_assert_msg( symmetric || pose.num_chains() < 4, "ERROR: use_ddG_style not compatible with symmetry or poses with > 3 chains" );
		runtime_assert( jump_num_ <= pose.num_jump() );
		core::pose::Pose unbound( pose );
		core::Real const unbound_dist = 1000.0;

		Size sym_aware_jump_id = 0;
		if ( symmetric ) {
			if ( sym_dof_names() != "" ) {
				utility::vector1<std::string> sym_dof_name_list;
				sym_dof_name_list = utility::string_split( sym_dof_names() , ',' );
				for ( Size i = 1; i <= sym_dof_name_list.size(); i++ ) {
					sym_aware_jump_id = core::pose::symmetry::sym_dof_jump_num( pose, sym_dof_name_list[i] );
					protocols::rigid::RigidBodyTransMoverOP translate_unbound( new protocols::rigid::RigidBodyTransMover( unbound, sym_aware_jump_id ) );
					translate_unbound->step_size( unbound_dist );
					translate_unbound->apply( unbound );
				}
			} else {
				sym_aware_jump_id = core::pose::symmetry::get_sym_aware_jump_num( pose, jump_num_ );
				protocols::rigid::RigidBodyTransMoverOP translate_unbound( new protocols::rigid::RigidBodyTransMover( unbound, sym_aware_jump_id ) );
				translate_unbound->step_size( unbound_dist );
				translate_unbound->apply( unbound );
			}
		} else {
			protocols::rigid::RigidBodyTransMover trans_mover( unbound, jump_num_ );
			trans_mover.trans_axis( trans_mover.trans_axis() );
			trans_mover.step_size(unbound_dist);
			trans_mover.apply( unbound );
		}
		unbound.update_residue_neighbors();
		(*sfxn_)(unbound ); // score the new pose, or we get assertion error.

		BuriedUnsatisfiedPolarsCalculator calc_unbound( name_of_sasa_calc_, name_of_hbond_calc, region_to_calculate, burial_cutoff_, probe_radius_, residue_surface_cutoff_, generous_hbonds_, legacy_counting_, use_vsasa_, use_sc_neighbors_, ignore_surface_res_ );
		basic::MetricValue< core::Size > mv_all_heavy_unbound, mv_bb_heavy_unbound, mv_countable_nonheavy_unbound, mv_all_unsat_unbound;
		calc_unbound.get("all_heavy_unsats", mv_all_heavy_unbound, unbound);
		calc_unbound.get("bb_heavy_unsats", mv_bb_heavy_unbound, unbound);
		calc_unbound.get("countable_nonheavy_unsats", mv_countable_nonheavy_unbound, unbound);
		calc_unbound.get("atom_bur_unsat", mv_unbound_unsat_map, unbound);
		calc_unbound.get("all_bur_unsat_polars", mv_all_unsat_unbound, unbound);

		all_heavy_atom_unsats = all_heavy_atom_unsats - mv_all_heavy_unbound.value();
		bb_heavy_atom_unsats = bb_heavy_atom_unsats - mv_bb_heavy_unbound.value();
		sc_heavy_atom_unsats = sc_heavy_atom_unsats - ( mv_all_heavy_unbound.value() - mv_bb_heavy_unbound.value() );
		countable_nonheavy_unsats = countable_nonheavy_unsats - mv_countable_nonheavy_unbound.value();
		all_unsats = all_unsats - mv_all_unsat_unbound.value();

		buried_unsat_hbond_filter_tracer << "  AFTER SUBTRACTING UNBOUND STATE: " << std::endl;
		if ( legacy_counting_ ) {
			buried_unsat_hbond_filter_tracer << "  all_unsats = " << all_unsats << std::endl;
		} else {
			buried_unsat_hbond_filter_tracer << "  all_heavy_atom_unsats = " << all_heavy_atom_unsats << std::endl << "  bb_heavy_atom_unsats = " << bb_heavy_atom_unsats << std::endl << "  sc_heavy_atom_unsats = " << sc_heavy_atom_unsats << std::endl << "  countable_nonheavy_unsats = " << countable_nonheavy_unsats << std::endl;
		}
		ddG_was_computed = true;
	}

	std::ostringstream oss;
	if ( print_out_info_to_pdb_ ) {
		std::string filter_name = this->name();
		std::string user_name = this->get_user_defined_name();
		oss << std::endl << filter_name << " " << user_name + ": " << std::endl;
	}
	for ( core::Size r = 1; r <= mv_unsat_map.value().size(); ++r ) {
		for ( core::Size a = 1; a <= mv_unsat_map.value().n_atom(r); ++a ) {
			if ( ignore_bb_heavy_unsats_ && pose.residue(r).atom_is_backbone(a) && a <= pose.residue(r).nheavyatoms() /* don't want backbone H's */ ) continue;
			if ( mv_unsat_map.value()(r,a) ) {
				if ( ddG_was_computed && mv_unbound_unsat_map.value()(r,a) ) continue; // for ddG, if it's also Unsat in the Unbound case, we don't care
				std::string unsat_type = ( pose.residue(r).atom_is_polar_hydrogen(a) ) ? "Hpol" : "HEAVY";
				std::string temp_str = "      Unsatisfied " + unsat_type + " polar atom at residue " + utility::to_string( r ) + ": " + pose.residue( r ).name3() + " " + pose.residue(r).atom_name(a);
				buried_unsat_hbond_filter_tracer << temp_str << std::endl;
				if ( print_out_info_to_pdb_ ) oss << temp_str << std::endl;
			}
		}
	}
	if ( print_out_info_to_pdb_ ) {
		protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
		oss << "  all_heavy_atom_unsats = " << all_heavy_atom_unsats << std::endl << "  bb_heavy_atom_unsats = " << bb_heavy_atom_unsats << std::endl << "  sc_heavy_atom_unsats = " << sc_heavy_atom_unsats << std::endl << "  countable_nonheavy_unsats = " << countable_nonheavy_unsats << std::endl;
		job->add_string( oss.str() );
	}
	// TJ's reporter options
	if ( use_reporter_behavior_ ) { // now defaults to true, with report_all_heavy_atom_unsats_=true as default
		if ( report_all_heavy_atom_unsats_ ) {
			return core::Real(all_heavy_atom_unsats);
		}
		if ( report_sc_heavy_atom_unsats_ ) {
			return core::Real(sc_heavy_atom_unsats);
		}
		if ( report_bb_heavy_atom_unsats_ ) {
			return core::Real(bb_heavy_atom_unsats);
		}
		if ( report_nonheavy_unsats_ ) {
			return core::Real(countable_nonheavy_unsats);
		}
		if ( report_all_unsats_ ) {
			return core::Real(all_unsats);
		}
		std::cout << "need to change one of the reporters to true"; // shouldn't ever happen now
		return(911);
	} else if ( legacy_options_ ) {
		buried_unsat_hbond_filter_tracer << " USING LEGACY OPTIONS: ALL UNSATS TREATED AS EQUIVALENT (NOT RECOMMENDED!);  all unsats = " << ( all_unsats ) << std::endl;
		return core::Real( all_unsats );
	} else if ( use_hbnet_behavior_ ) {
		core::Size counted_heavy_atom_unsats( ( ignore_bb_heavy_unsats_ ) ? sc_heavy_atom_unsats : all_heavy_atom_unsats );
		if ( counted_heavy_atom_unsats > 0 ) {
			buried_unsat_hbond_filter_tracer << "  HEAVY ATOM UNSATS DETECTED! AUTOMATIC FAIL!" << std::endl << "  HEAVY UNSATS ";
			if ( ignore_bb_heavy_unsats_ ) buried_unsat_hbond_filter_tracer << "minus backbone heavy-atom unsats (ignore_bb_heavy_unsats=true)";
			buried_unsat_hbond_filter_tracer << " = " << counted_heavy_atom_unsats << std::endl;
			return 9999.9;
		}
		return core::Real( countable_nonheavy_unsats );
	}
	// shoudn't ever get here, but if we do, return default: reporter behavior that reports all heavy atom unsats
	return core::Real(all_heavy_atom_unsats);
}

void
BuriedUnsatHbondFilter::task_factory( core::pack::task::TaskFactoryOP tf ){
	task_factory_ = tf;
}

core::pack::task::TaskFactoryOP
BuriedUnsatHbondFilter::task_factory() const {
	return task_factory_;
}

std::string BuriedUnsatHbondFilter::name() const {
	return class_name();
}

std::string BuriedUnsatHbondFilter::class_name() {
	return "BuriedUnsatHbonds";
}

void BuriedUnsatHbondFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "use_legacy_options", xsct_rosetta_bool, "revert to legacy options (equivalent to old, original BuriedUnsat Filter; WARNING! If this is true, will overwrite all other options", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "generous_hbonds", xsct_rosetta_bool, "count all h-bonds (not just those scored by the default scorefxn in rosetta", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "use_vsasa", xsct_rosetta_bool, "use vsasa insteady of legacy sasa for burial calculation", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "ignore_surface_res", xsct_rosetta_bool, "many polar atoms on surface atoms get flagged as buried unsat becuause they are occluded by a long sidechina (e.g. Lys or Arg) that could easily move out of the way; this option ignores surface residues, as deinfed by SASA (default) or sc_neighbors if use_sc_neighbors=true", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "ignore_bb_heavy_unsats", xsct_rosetta_bool, "ignore bb heayy atom unsats when using hbnet-style behavior", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "use_sc_neighbors", xsct_rosetta_bool, "use sc_neighbors insteady of SASA for burial calculations", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "use_ddG_style", xsct_rosetta_bool, "perform ddG style calcation: the Unsats are calculated by subtracting all unsats in bound state from unbound state; this is how the original BuriedUnsatHBondsFilter works", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "only_interface", xsct_rosetta_bool, "restrict unsat search only to interface residues; if true and more than one chain it's ignored", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "print_out_info_to_pdb", xsct_rosetta_bool, "print all info to pdb file into addition to tracer", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "jump_number", xsct_non_negative_integer, "The jump over which to evaluate the filter; only applies to use_ddG_style", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "cutoff", xsct_non_negative_integer, "The upper threshold for counted buried unsat H-bonds, above which the filter fails", "20" )
		+ XMLSchemaAttribute::attribute_w_default( "probe_radius", xsct_real, "probe radius to use for SASA buriedness calculations; default is grabbed from sasa_calculator_probe_radius in options code, which defaults to 1.4", "1.4" )
		+ XMLSchemaAttribute::attribute_w_default( "burial_cutoff", xsct_real, "used to determine burial; deafault legacy SASA atomic_burial_cutoff is 0.3; default VSASA cutoff is 0.1; if use_sc_neighbors=true, default becomes 4.4 or can be user-specified to sc_neighbor cutoff that is desired", "0.3" )
		+ XMLSchemaAttribute::attribute_w_default( "residue_surface_cutoff", xsct_real, "cutoff to determine which residues are surface if ignore_surface_res=true; default is 45.0 for SASA, 20.0 for VSASA and 2.0 if use_sc_neighbors=true", "45.0" )
		+ XMLSchemaAttribute::attribute_w_default( "use_reporter_behavior",xsct_rosetta_bool,"report as filter score the type of unsat turned on; this is now TRUE by default","true")
		+ XMLSchemaAttribute::attribute_w_default( "use_hbnet_behavior",xsct_rosetta_bool,"no heavy unstas allowed (will return 9999); if no heavy unstas, will count Hpol unsats; FALSE by default; if set to true, will NOT use reporter behavior","false")
		+ XMLSchemaAttribute::attribute_w_default( "report_all_unsats",xsct_rosetta_bool,"report all unsats","false")
		+ XMLSchemaAttribute::attribute_w_default( "report_all_heavy_atom_unsats",xsct_rosetta_bool,"report all heavy atom unsats; IF ALL REORTER OPTIONS ARE FALSE, THIS BECOMES TRUE AND DEFAULT REPORTS ALL HEAVY UNSATS","false")
		+ XMLSchemaAttribute::attribute_w_default( "report_sc_heavy_atom_unsats",xsct_rosetta_bool,"report side chain heavy atom unsats","false")
		+ XMLSchemaAttribute::attribute_w_default( "report_bb_heavy_atom_unsats",xsct_rosetta_bool,"report back bone heavy atom unsats","false")
		+ XMLSchemaAttribute::attribute_w_default( "report_nonheavy_unsats",xsct_rosetta_bool,"report non heavy atom unsats","false")
		+ XMLSchemaAttribute( "sym_dof_names" , xs_string , "For multicomponent symmetry: what jump(s) used for ddG-like separation. (From Dr. Bale: For multicomponent systems, one can simply pass the names of the sym_dofs that control the master jumps. For one component systems, jump can still be used.)  IF YOU DEFIN THIS OPTION, Will use ddG-style separation for the calulation; if you do not want this, pass a residue selector instead of defining symdofs." )
		+ XMLSchemaAttribute( "residue_selector", xs_string, "residue selector that tells the filter to restrict the Unsat search to only those residues" );
	rosetta_scripts::attributes_for_get_score_function_name( attlist );
	rosetta_scripts::attributes_for_parse_task_operations( attlist );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string BuriedUnsatHbondFilterCreator::keyname() const {
	return BuriedUnsatHbondFilter::class_name();
}

filters::FilterOP
BuriedUnsatHbondFilterCreator::create_filter() const {
	return filters::FilterOP( new BuriedUnsatHbondFilter );
}

void BuriedUnsatHbondFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BuriedUnsatHbondFilter::provide_xml_schema( xsd );
}


}
}
