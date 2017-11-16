// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/TaskAwareScoreTypeFilter.cc
/// @brief Formerly devel/matdes/AverageInterfaceEnergyFilter. Calculates/filters based on a particular score type over the set of packable residues specified by the user-provided task operations. Can be calculate in either the bound or unbound state (a jump or sym_dof_names must be provided to create the unbound pose). Can operate in 1 of 3 modes: 1) if mode=total then compute function returns the total score over all packable residues and the apply function returns true if the total score is less than or equal to the user-specified threshold (false otherwise), 2) if mode=average then the compute function returns the average score over all packable residues and the apply function returns true if the average score is less than or equal to the user-specified threshold (false otherwise), and 3) if mode=individual then the compute function returns the number of individual residues that do not pass the user-specified score type threshold and the apply function returns true if all residues passed the user-specified threshold (false otherwise).
/// @author Jacob Bale (balej@uw.edu), Neil King (neilking@uw.edu)

#include <protocols/simple_filters/TaskAwareScoreTypeFilter.hh>
#include <protocols/simple_filters/TaskAwareScoreTypeFilterCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <utility/string_util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <protocols/jd2/util.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/pose/symmetry/util.hh>
#include <ObjexxFCL/format.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.TaskAwareScoreTypeFilter" );

/// @brief default ctor
TaskAwareScoreTypeFilter::TaskAwareScoreTypeFilter() :
	parent( "TaskAwareScoreType" ),
	task_factory_( /* NULL */ ),
	scorefxn_( /* NULL */ ),
	threshold_( 100000 ),
	bb_bb_( false ),
	write2pdb_( false ),
	sym_dof_names_( "" ),
	jump_( 0 ),
	mode_( "average" )
{}

core::pack::task::TaskFactoryOP
TaskAwareScoreTypeFilter::task_factory() const{
	return task_factory_;
}

void
TaskAwareScoreTypeFilter::task_factory( core::pack::task::TaskFactoryOP task_factory ){
	task_factory_ = task_factory;
}

core::scoring::ScoreFunctionOP
TaskAwareScoreTypeFilter::scorefxn() const{
	return scorefxn_;
}

void
TaskAwareScoreTypeFilter::scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
	scorefxn_ = scorefxn;
}

core::scoring::ScoreType
TaskAwareScoreTypeFilter::score_type() const{
	return score_type_;
}

void
TaskAwareScoreTypeFilter::score_type( core::scoring::ScoreType const st ){
	score_type_ = st;
}

core::Real
TaskAwareScoreTypeFilter::threshold() const{
	return threshold_;
}

void
TaskAwareScoreTypeFilter::threshold( core::Real const thresh ){
	threshold_ = thresh;
}

bool
TaskAwareScoreTypeFilter::bb_bb() const{
	return bb_bb_;
}

void
TaskAwareScoreTypeFilter::bb_bb( bool const bb ){
	bb_bb_ = bb;
}

std::string
TaskAwareScoreTypeFilter::mode() const {
	return mode_;
}

bool
TaskAwareScoreTypeFilter::write2pdb() const {
	return write2pdb_;
}

void
TaskAwareScoreTypeFilter::write2pdb( bool const write ) {
	write2pdb_ = write;
}

bool
TaskAwareScoreTypeFilter::individual_hbonds() const{
	return individual_hbonds_;
}

void
TaskAwareScoreTypeFilter::individual_hbonds( bool individual_hbonds ){
	individual_hbonds_ = individual_hbonds;
}

void
TaskAwareScoreTypeFilter::mode( std::string const & mode ) {
	if ( !((mode == "average") || (mode == "total") || (mode == "individual")) ) {
		utility_exit_with_message( "InputError: The specified mode for the TaskAwareScoreTypeFilter does not match either average, total, or individual." );
	}
	mode_ = mode;
}

bool TaskAwareScoreTypeFilter::unbound() const { return unbound_; }
void TaskAwareScoreTypeFilter::unbound( bool const unbound ) { unbound_ = unbound; }
std::string TaskAwareScoreTypeFilter::sym_dof_names() const { return sym_dof_names_; }
void TaskAwareScoreTypeFilter::sym_dof_names( std::string const & sym_dofs ) { sym_dof_names_ = sym_dofs; }
core::Size TaskAwareScoreTypeFilter::jump() const { return jump_; }
void TaskAwareScoreTypeFilter::jump( core::Size const jump ) { jump_ = jump; }
std::string TaskAwareScoreTypeFilter::score_type_name() const { return score_type_name_; }
void TaskAwareScoreTypeFilter::score_type_name( std::string const & name ) { score_type_name_ = name; }

bool
TaskAwareScoreTypeFilter::apply(core::pose::Pose const & pose ) const
{
	core::Real score(compute( pose, false ));
	if ( mode() == "individual" ) { //TODO: JBB - Modify to allow the user to specify both a per residue threshold and a threshold on the number of failing residues allowed to still pass the filter.
		if ( score == 0 ) {
			TR<<"passing."<<std::endl;
			return true;
		} else {
			TR<<"failing."<<std::endl;
			return false;
		}
	} else {
		if ( score <= threshold_ ) {
			TR<<"passing."<<std::endl;
			return true;
		} else {
			TR<<"failing."<<std::endl;
			return false;
		}
	}
}

core::Real
TaskAwareScoreTypeFilter::compute( core::pose::Pose const & pose, bool const & write ) const{
	runtime_assert( task_factory() != nullptr );
	core::pack::task::PackerTaskCOP packer_task( task_factory()->create_task_and_apply_taskoperations( pose ) );
	core::Size total_residue;
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
		total_residue = symm_info->num_independent_residues();
	} else {
		total_residue = pose.size();
	}
	std::string interface_pos("interface positions considered: ");
	core::Real tie=0;
	core::Real interface_energy = 0;
	core::Size interface_size=0;
	core::Real num_failing_resis = 0;
	std::string pymol_selection = "select " + this->get_user_defined_name() + "_pos, ";
	core::pose::Pose p = pose;

	// If unbound energies are being evaluated, create the unbound state
	if ( unbound() ) {
		int sym_aware_jump_id = 0;
		if ( sym_dof_names() != "" ) {
			utility::vector1<std::string> sym_dof_name_list;
			sym_dof_name_list = utility::string_split( sym_dof_names() , ',' );
			for ( Size i = 1; i <= sym_dof_name_list.size(); i++ ) {
				sym_aware_jump_id = core::pose::symmetry::sym_dof_jump_num( pose, sym_dof_name_list[i] );
				protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( p, sym_aware_jump_id ) );
				translate->step_size( 1000.0 );
				translate->apply( p );
			}
		} else if ( jump() != 0 ) {
			sym_aware_jump_id = core::pose::symmetry::get_sym_aware_jump_num( pose, jump() );
			protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( p, sym_aware_jump_id ) );
			translate->step_size( 1000.0 );
			translate->apply( p );
		} else {
			utility_exit_with_message( "If unbound is set to true, you must provide either sym_dof_names or a jump in order to create the unbound pose!" );
		}
	}

	core::scoring::hbonds::HBondSet hbset, hbset_ref;
	core::scoring::EnergyMap em;
	if ( score_type_name() == "hbond_bb_sc" ) {
		core::scoring::methods::EnergyMethodOptions myopt = scorefxn()->energy_method_options();
		myopt.hbond_options().decompose_bb_hb_into_pair_energies(true);
		scorefxn()->set_energy_method_options(myopt);
		scorefxn()->score( p );
		// get the HBondSet
		core::scoring::hbonds::fill_hbond_set(p,false,hbset);
	} else {
		scorefxn()->score( p );
	}
	//TODO JBB: Decide if this functionality is really needed/useful.  If so then, also consider adding the ability to do this with mode=individual.
	if ( individual_hbonds_ ) {
		for ( core::Size i=1; i<=total_residue; ++i ) {
			if ( packer_task->being_packed( i ) ) {
				bool flag = 0;
				core::Real resi_hbond_bb_sc_energy = 0;
				for ( core::Size j=1; j<=total_residue; j++ ) {
					for ( Size ihb = 1; ihb <= hbset.nhbonds(); ++ihb ) {
						core::scoring::hbonds::HBond const & hb(hbset.hbond(ihb));
						if ( hb.don_res()==i && hb.acc_res()==j ) {
							if ( !hb.don_hatm_is_protein_backbone() && hb.acc_atm_is_protein_backbone() ) {
								if ( flag == 0 ) interface_size++;
								flag=1;
								resi_hbond_bb_sc_energy += hb.energy();
								//TR << "Donor sc resi: " << i << p.residue(i).name3() << ". Acceptor bb resi: " << j << p.residue(j).name3() << " Energy = " << hb.energy() << std::endl;
							}
						}
						if ( hb.don_res()==j && hb.acc_res()==i ) {
							if ( hb.don_hatm_is_protein_backbone() && !hb.acc_atm_is_protein_backbone() ) {
								if ( flag == 0 ) interface_size++;
								flag=1;
								resi_hbond_bb_sc_energy += hb.energy();
								//TR << "Acceptor sc resi: " << i << p.residue(i).name3() << ". Donor bb resi: " << j << p.residue(j).name3() << " Energy = " << hb.energy() << std::endl;
							}
						}
					}
				}
				tie += resi_hbond_bb_sc_energy;
			}
		}
	} else {
		for ( core::Size resi=1; resi<=total_residue; ++resi ) {
			if ( packer_task->being_packed( resi ) ) {
				interface_size += 1;
				interface_pos.append(ObjexxFCL::string_of(resi) + "+");
				if ( mode() == "individual" ) {
					em = p.energies().residue_total_energies( resi );
					em *= scorefxn()->weights();
					if ( em[score_type()] <= threshold() ) {
						TR << score_type_name() << " " << pose.residue( resi ).name3() << "_" << resi << ": " << em[score_type()] << " pass" << std::endl;
					} else {
						TR << score_type_name() << " " << pose.residue( resi ).name3() << "_" << resi << ": " << em[score_type()] << " fail" << std::endl;
						if ( write ) { write_to_pdb( pose, resi, pose.residue( resi ).name3(), em[score_type()] ); }
						core::Size output_resi = resi;
						if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
							output_resi = pose.pdb_info()->number( resi );
						}
						pymol_selection.append(ObjexxFCL::string_of(output_resi) + "+");
						num_failing_resis += 1;
					}
				} else {
					interface_energy += p.energies().residue_total_energy( resi );
					em += p.energies().residue_total_energies( resi );
				}
			}
		}
		TR << pymol_selection << std::endl;
	}
	if ( mode() == "individual" ) {
		return num_failing_resis;
	} else {
		if ( score_type_name() == "total_score" ) {
			tie = interface_energy;
		} else {
			em *= scorefxn()->weights();
			tie = em[score_type()];
		}
		TR.Debug <<interface_pos<<std::endl;
		if ( mode() == "total" ) {
			return ( tie );
		} else if ( mode() == "average" ) {
			return( tie / interface_size );
		} else {
			utility_exit_with_message( "InputError: The specified mode for the TaskAwareScoreTypeFilter does not match either average, total, or individual." );
		}
	}
}

void TaskAwareScoreTypeFilter::write_to_pdb( core::pose::Pose const & pose, core::Size const residue, std::string const & residue_name, core::Real const score ) const
{

	std::string filter_name = this->name();
	std::string user_name = this->get_user_defined_name();
	core::Size output_resi = residue;
	if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
		output_resi = pose.pdb_info()->number( residue );
	}
	std::string output_string = filter_name + " " + user_name + ": " + residue_name + ObjexxFCL::string_of(output_resi) + " = " + ObjexxFCL::string_of(score);
	protocols::jd2::add_string_to_current_job(output_string);

}

core::Real
TaskAwareScoreTypeFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::Real score(compute( pose, false ));
	return( score );
}

void
TaskAwareScoreTypeFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out<<"TaskAwareScoreTypeFilter returns "<<compute( pose, write2pdb_ )<<std::endl;
}

void
TaskAwareScoreTypeFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	TR << "TaskAwareScoreTypeFilter"<<std::endl;
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	scorefxn(protocols::rosetta_scripts::parse_score_function(tag, data));
	score_type(core::scoring::score_type_from_name( tag->getOption<std::string>( "score_type", "total_score" ) ) );
	score_type_name( tag->getOption<std::string>( "score_type", "total_score" ) );
	threshold( tag->getOption< core::Real >( "threshold", 100000 ) );
	bb_bb( tag->getOption< bool >( "bb_bb", false ) );
	unbound( tag->getOption< bool >( "unbound", false ) );
	sym_dof_names( tag->getOption< std::string >( "sym_dof_names", "" ) );
	jump( tag->getOption< core::Size >( "jump", 0 ) );
	mode( tag->getOption< std::string >( "mode", "average" ) );
	write2pdb( tag->getOption< bool >("write2pdb", 0) );
	individual_hbonds(tag->getOption< bool >("individual_hbonds", false ));
	if ( !((mode() == "average") || (mode() == "total") || (mode() == "individual")) ) {
		utility_exit_with_message( "InputError: The specified mode for the TaskAwareScoreTypeFilter does not match either average, total, or individual." );
	}
	TR<<"with options scoretype: "<<score_type()<<" and threshold: "<<threshold()<<std::endl;
}

protocols::filters::FilterOP
TaskAwareScoreTypeFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new TaskAwareScoreTypeFilter() );
}

TaskAwareScoreTypeFilter::~TaskAwareScoreTypeFilter()= default;

protocols::filters::FilterOP
TaskAwareScoreTypeFilter::clone() const{
	return protocols::filters::FilterOP( new TaskAwareScoreTypeFilter( *this ) );
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP TaskAwareScoreTypeFilterCreator::create_filter() const { return protocols::filters::FilterOP( new TaskAwareScoreTypeFilter ); }

// XRW TEMP std::string
// XRW TEMP TaskAwareScoreTypeFilterCreator::keyname() const { return "TaskAwareScoreType"; }

std::string TaskAwareScoreTypeFilter::name() const {
	return class_name();
}

std::string TaskAwareScoreTypeFilter::class_name() {
	return "TaskAwareScoreType";
}

void TaskAwareScoreTypeFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	// AMW: should put in utility (also used elsewhere) but instead
	// will namespace with TaskAwareScoreTypeFilter_
	XMLSchemaRestriction mode_enumeration;
	mode_enumeration.name( "TaskAwareScoreTypeFilter_mode_choices" );
	mode_enumeration.base_type( xs_string );
	mode_enumeration.add_restriction( xsr_enumeration, "average" );
	mode_enumeration.add_restriction( xsr_enumeration, "total" );
	mode_enumeration.add_restriction( xsr_enumeration, "individual" );
	xsd.add_top_level_element( mode_enumeration );

	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "score_type" , xs_string , "Score term to use for the calculations. Defaults to total_score." , "total_score" )
		+ XMLSchemaAttribute::attribute_w_default( "threshold" , xsct_real , "Threshold for the highest permissible value for a passing design. Note: When using mode='individual' this is the score threshold that each individual residue must pass in order for the filter to pass. For example, the user could set threshold=3.5, score_type=fa_rep, and mode=individual to check if all the packable residues have a fa_rep score less than 3.5. The value reported when mode=individual is the number of failing residues." , "100000" )
		+ XMLSchemaAttribute::attribute_w_default( "bb_bb" , xsct_rosetta_bool , "See EnergyPerResidueFilter." , "false" )
		+ XMLSchemaAttribute::attribute_w_default( "unbound" , xsct_rosetta_bool , "Default is false. Set to true if you want the scores to be assessed in the unbound state. If unbound is set true, then the user must provide either the sym_dof_names or jump along which to separate the pose." , "false" )
		+ XMLSchemaAttribute( "sym_dof_names" , xs_string, "For use with the unbound option. Which sym_dofs should be used to separate the pose?" )
		+ XMLSchemaAttribute::attribute_w_default( "jump" , xsct_non_negative_integer , "For use with the unbound option. Which jump should be used to separate the pose?" , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "mode" , "TaskAwareScoreTypeFilter_mode_choices" , "This filter can operate in three modes. The options are 'total', 'average', or 'individual'. If mode=total, then the total score of the packable residues is returned and the threshold refers to that total value. If mode=average, then the average score of the packable residues is returned and the threshold refers to that average value. If mode=individual, then each packable residue is evaluated on an individual basis and the threshold refers to those individual values; in this case, individual residue must pass in order for the filter to pass. For example, the user could set threshold=3.5, score_type=fa_rep, and mode=individual to check if all the packable residues have a fa_rep score less than 3.5. The value reported from the filter when mode=individual is the number of failing residues." , "average" )
		+ XMLSchemaAttribute::attribute_w_default( "write2pdb" , xsct_rosetta_bool , "Only for use with mode=individual. Default is false. Setting to true will cause the individual scores to be output to the bottom of the pdb." , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "individual_hbonds" , xsct_rosetta_bool , "This seems to still be in development." , "false" ) ;

	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist ) ;
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist ) ;

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Takes task operations to determine the packable residues and then calculates/filters based on the scores of those residues in one of three possible modes: 1) total: the total score of all packable residues 2) average: the average score of all packable residues 3) individual: the scores of individual residues (in individual mode each residue must pass the user-defined threshold in order for the filter to pass). In each mode, the score type used to assess the score is set by the user (default is total_score).", attlist );
}

std::string TaskAwareScoreTypeFilterCreator::keyname() const {
	return TaskAwareScoreTypeFilter::class_name();
}

protocols::filters::FilterOP
TaskAwareScoreTypeFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new TaskAwareScoreTypeFilter );
}

void TaskAwareScoreTypeFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TaskAwareScoreTypeFilter::provide_xml_schema( xsd );
}


} // simple_filters
} // protocols
