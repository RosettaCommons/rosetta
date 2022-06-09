// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/SixDoFGridDockMover.cc
/// @brief Class for enumerating docked orientations (3 translations and 3 rotations) between two chains
/// @author Odessa Goudy (oda@email.unc.edu)

// Unit headers
#include <protocols/simple_moves/SixDoFGridDockMover.hh>
#include <protocols/simple_moves/SixDoFGridDockMoverCreator.hh>

#include <protocols/jd2/util.hh> //uses jd2
#include <iostream> //for array
#include <array> //for array
#include <vector> //for vector
#include <protocols/rosetta_scripts/util.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/select/util.hh> //for residue selectors, get_residues_from_subset
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/Energies.hh> // for scoring?
#include <core/scoring/ScoreFunction.hh> // for scoring?
#include <core/pose/extra_pose_info_util.hh> //for score.sc file, setPoseExtraScore function

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/datacache/DataMap.hh>

#include <utility/stream_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>
#include <utility/LexicographicalIterator.hh>
#include <utility/string_util.hh>
#include <utility/VBWrapper.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <protocols/rigid/RigidBodyMover.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.simple_moves.SixDoFGridDockMover" );

namespace protocols {
namespace simple_moves {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
SixDoFGridDockMover::SixDoFGridDockMover():
	protocols::moves::Mover( SixDoFGridDockMover::mover_name() ),
	dof_sample_index_( 0 )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
}

////////////////////////////////////////////////////////////////////////////////
/// @Brief Copy constructor
SixDoFGridDockMover::SixDoFGridDockMover( SixDoFGridDockMover const & ) = default;


////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
SixDoFGridDockMover::~SixDoFGridDockMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
SixDoFGridDockMover::apply( core::pose::Pose & pose ){
	using core::Size;

	numeric::xyzVector< core::Real > selected_pos(0., 0., 0.);
	numeric::xyzVector< core::Real > selected_pos2a(0., 0., 0.);
	numeric::xyzVector< core::Real > selected_pos2b(0., 0., 0.);

	numeric::xyzVector< core::Real > axis1(0., 0., 0.);
	numeric::xyzVector< core::Real > axis2(0., 0., 0.);
	numeric::xyzVector< core::Real > axis3(0., 0., 0.);

	if ( selector1_ == nullptr ) {
		utility_exit_with_message("selector1 does not exist in apply");
	}
	if ( selector2a_ == nullptr ) {
		utility_exit_with_message("selector2a does not exist in apply");
	}
	if ( selector2b_ == nullptr ) {
		utility_exit_with_message("selector2b does not exist in apply");
	}

	if ( selector1_ ) {

		TR << "selector1 exists in apply: " << selector1_ << std::endl;
		utility::vector1< core::Size > residue1 = core::select::get_residues_from_subset( selector1_->apply( pose ) );
		//get_residues_from_subset necessary for conversion back and forth from pdb numbering to rosetta numbering

		//below takes the vector of sizes residue1 and iterates through all of the values within residue1 vector
		//in our class, residue1 is only size 1 and contains the rosetta numbering of the residue (in example, 33)
		if ( residue1.size() != 1 ) {
			utility_exit_with_message( "ERROR: Provide one residue for residue 1" );
		}
		selected_pos = pose.residue( residue1.front() ).xyz("CA");
		TR << "selector1_ xyz " << selected_pos[0] << " " << selected_pos[1] << " " <<  selected_pos[2] << std::endl;

		utility::vector1< core::Size > residue2a;
		residue2a = core::select::get_residues_from_subset( selector2a_->apply( pose ) );
		TR << "selector2a residue number: " << residue2a <<std::endl;
		if ( residue2a.size() != 1 ) {
			utility_exit_with_message( "ERROR: Provide one residue for residue 2a" );
		}
		selected_pos2a = pose.residue( residue2a.front() ).xyz("CA");
		TR << "selector2a_ xyz " << selected_pos2a[0] << " " << selected_pos2a[1] << " " << selected_pos2a[2] << std::endl;

		utility::vector1< core::Size > residue2b;
		residue2b = core::select::get_residues_from_subset( selector2b_->apply( pose ) );
		TR << "selector2b residue number: " << residue2b <<std::endl;
		if ( residue2b.size() != 1 ) {
			utility_exit_with_message( "ERORR: Provide one residue for residue 2b" );
		}
		selected_pos2b = pose.residue( residue2b.front() ).xyz("CA");
		TR << "selector2b_ xyz " << selected_pos2b[0] << " " << selected_pos2b[1] << " " << selected_pos2b[2] << std::endl;


		//Define three othogonal axes and check angle between given residue selectors with res2a as the vertex (degree check)
		axis3 = ( selected_pos2a - selected_pos );
		core::Real dot_product = dot(axis3, selected_pos2b - selected_pos2a); //used to find angle between given residue selectors
		TR << "dot-product of axis3 and chosen axes: " << dot_product << std::endl;
		core::Real axis3_length = axis3.length();
		core::Real b_a_length = (selected_pos2b - selected_pos2a).length();
		core::Real angle_rad = acos( dot_product / ( axis3_length * b_a_length ));
		core::Real angle_deg =  angle_rad * 180.0 / 3.14159;
		TR << "angle in deg: " << angle_deg << std::endl;
		if ( degree_check_ ) { //degree check
			if ( angle_deg < 60 || angle_deg > 120 ) {
				utility_exit_with_message("ERROR: Poorly chosen residue selectors. Angle is outside the range of 60-120 degrees.");
			}
		}
		TR << "axis 3 " << axis3[0] << std::endl;
		axis1 = axis3.cross( selected_pos2b - selected_pos2a ); //creating orthogonal axes to already defined axis3
		TR << "axis 1 " << axis1[0] << std::endl;
		axis2 = axis1.cross( axis3 );
		TR << "axis 2 " << axis2[0] << std::endl;

	} else {
		utility_exit_with_message("selector1 does not exist in apply");
	}

	rigid::RigidBodyDeterministicSpinMover axis1_spin_mover ( 1, axis1, selected_pos2a, lex_position( dof_values_, dof_sample_index_, 1) ); //function arguments are jump, axis to apply the rotation, center of rotation, rotation degree
	rigid::RigidBodyDeterministicSpinMover axis2_spin_mover ( 1, axis2, selected_pos2a, lex_position( dof_values_, dof_sample_index_, 2) );
	rigid::RigidBodyDeterministicSpinMover axis3_spin_mover ( 1, axis3, selected_pos2a, lex_position( dof_values_, dof_sample_index_, 3) );

	rigid::RigidBodyTransMover axis1_trans_mover(pose,1); //function arguments are pose to apply values to, jump
	axis1_trans_mover.trans_axis( axis1 ); //set axis for trans mover
	axis1_trans_mover.step_size( lex_position( dof_values_, dof_sample_index_, 4) ); //set step size for trans mover

	rigid::RigidBodyTransMover axis2_trans_mover(pose,1);
	axis2_trans_mover.trans_axis( axis2 );
	axis2_trans_mover.step_size( lex_position( dof_values_, dof_sample_index_, 5) );

	rigid::RigidBodyTransMover axis3_trans_mover(pose,1);
	axis3_trans_mover.trans_axis( axis3 );
	axis3_trans_mover.step_size( lex_position( dof_values_, dof_sample_index_, 6) );

	axis1_spin_mover.apply(pose); //apply the rotations first, then translations - to reduce lever arm effect
	axis2_spin_mover.apply(pose);
	axis3_spin_mover.apply(pose);
	axis1_trans_mover.apply(pose);
	axis2_trans_mover.apply(pose);
	axis3_trans_mover.apply(pose);

	TR << "dof_sample_index" << dof_sample_index_ << std::endl; //prints out job name/number

	//following section adds dof values to the score file
	core::pose::setPoseExtraScore(pose, "axis1_rot", lex_position( dof_values_, dof_sample_index_, 1)); //function arguments are pose to read from, scorefile column name, specific DoF value to report
	core::pose::setPoseExtraScore(pose, "axis2_rot", lex_position( dof_values_, dof_sample_index_, 2));
	core::pose::setPoseExtraScore(pose, "axis3_rot", lex_position( dof_values_, dof_sample_index_, 3));
	core::pose::setPoseExtraScore(pose, "axis1_trans", lex_position( dof_values_, dof_sample_index_, 4));
	core::pose::setPoseExtraScore(pose, "axis2_trans", lex_position( dof_values_, dof_sample_index_, 5));
	core::pose::setPoseExtraScore(pose, "axis3_trans", lex_position( dof_values_, dof_sample_index_, 6));
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
SixDoFGridDockMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
SixDoFGridDockMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
)
{
	using core::Size;

	dof_values_.resize(6);

	std::string range_option_name;
	std::string vals_option_name;

	debug_assert(data.has("JobLabels", "nstruct_index")); //included to ensure that users can reinitialize the mover

	for ( Size ii = 1; ii <= 6; ++ii ) {
		if ( ii <= 3 ) { //this for loop gets is split so that we can store the name into a simpler variable. the first half (1-3) is for rotation DoFs, then (4-6) for translation DoFs.
			range_option_name = "range_rot_axis_" + utility::to_string(ii);
			vals_option_name = "values_rot_axis_" + utility::to_string(ii);
		} else if ( ii > 3 ) {
			range_option_name = "range_trans_axis_" + utility::to_string(ii - 3);
			vals_option_name = "values_trans_axis_" + utility::to_string(ii - 3);
		}

		if ( tag->hasOption(range_option_name) && tag->hasOption(vals_option_name) ) {
			utility_exit_with_message("ERROR: User cannot provide *both* a range option and a value option for a single DOF.");
		}
		if ( tag->hasOption( range_option_name ) ) {
			dof_values_[ii] = parse_range( tag, range_option_name );
			TR << "Search range for DOF " << ii << " : " << dof_values_[ii] << std::endl;

		} else if ( tag->hasOption( vals_option_name ) ) {
			dof_values_[ii] = utility::string_split( tag->getOption< std::string >( vals_option_name ), ',', core::Real() );
			TR << "User provided value(s) for DOF " << ii << ": " << dof_values_[ii] << std::endl;
		} else {
			utility_exit_with_message("ERROR: User must provide at least one option (either range option or value option) for each DOF to be sampled.");
		}
	}


	//Residue selectors for x, y, z coordinate system set up


	if ( tag->hasOption("dof_residue_selector_1") && tag->hasOption("dof_residue_selector_2a") && tag->hasOption("dof_residue_selector_2b") ) {
		selector1_ = core::select::residue_selector::parse_residue_selector( tag, data, "dof_residue_selector_1" );
		selector2a_ = ( core::select::residue_selector::parse_residue_selector( tag, data, "dof_residue_selector_2a" ) );
		selector2b_ = ( core::select::residue_selector::parse_residue_selector( tag, data, "dof_residue_selector_2b" ) );

		if ( selector1_ ) {
			TR << "selector1 exists in parse_my_tag! " << selector1_ << std::endl;
		} else {
			utility_exit_with_message( "ERROR: Please declare three residue selectors in the mover option." );
		}
	}

	// parse max_samples here and exit iff total search exceeds user provided max_samples
	if ( tag->hasOption("max_samples") ) {
		max_samples_ = tag->getOption< core::Size >( "max_samples" );
		TR << "User provided value for maximum samples/scope is: " << max_samples_ << std::endl;
	} else {
		max_samples_ = 10000;
		TR << "Default value for maximum samples/scope is set at: " << max_samples_ << std::endl;
	}

	if ( max_samples_ ) {
		core:: Real log_sum = log(dof_values_[1].size()); //note, caluclation compares the LOG of the samplespace to the LOG of the maximum expected output
		TR << "log(dof1)" << log_sum << std::endl;
		for ( Size ii = 2; ii <= 6; ++ii ) {
			log_sum += log(dof_values_[ii].size());
		}
		TR << "Scope of log(search space) is: " << log_sum << std::endl;
		if ( log_sum > log(max_samples_) ) {
			utility_exit_with_message( "ERROR: Based on user input DoF values, the scope of the search exceeds either the user input or default value (of 10,000) for the maximum number of samples. Exiting. Change max samples with the max_samples option." );
		}
	}

	if ( tag->hasOption("degree_check") ) {
		degree_check_ = tag->getOption< bool >( "degree_check" );
	}

	typedef utility::pointer::shared_ptr< utility::VBWrapper<Size> > WrappedSizeOP; //VBWrapper added to be compatible with jd3 and jd2, namely to store job ID data in a wrapper to work with jd3
	if ( jd2::jd2_used() &&
			! basic::options::option[basic::options::OptionKeys::run::reinitialize_mover_for_each_job] ) {

		// Consider it ok if all the dof_values_ lists have exactly 1 value.
		bool ok = true;
		for ( Size ii = 1; ii <= 6; ++ii ) {
			if ( dof_values_[ii].size() != 1 ) {
				ok = false;
				break;
			}
		}

		if ( !ok ) {
			std::string error = "ERROR: When using the SixDoFGridDockMover with the jd2 job distributor,"
				" you must use the -run:reinitialize_mover_for_each_job flag or the mover will not"
				" enumerate the full set of combinations you're looking for.\n";
			utility_exit_with_message(error);
		}
	}
	WrappedSizeOP wrapped_ind = data.get_ptr<utility::VBWrapper<Size>>("JobLabels", "nstruct_index");
	dof_sample_index_ = wrapped_ind->data();

}


void SixDoFGridDockMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	attlist + XMLSchemaAttribute( "range_trans_axis_1", xsct_real_cslist, "Comma-separated list of DoF range: minimum, maximum, step size" )
		+ XMLSchemaAttribute( "range_trans_axis_2", xsct_real_cslist, "Comma-separated list of DoF range: minimum, maximum, step size" )
		+ XMLSchemaAttribute( "range_trans_axis_3", xsct_real_cslist, "Comma-separated list of DoF range: minimum, maximum, step size" )
		+ XMLSchemaAttribute( "range_rot_axis_1", xsct_real_cslist, "Comma-separated list of DoF range: minimum, maximum, step size" )
		+ XMLSchemaAttribute( "range_rot_axis_2", xsct_real_cslist, "Comma-separated list of DoF range: minimum, maximum, step size" )
		+ XMLSchemaAttribute( "range_rot_axis_3", xsct_real_cslist, "Comma-separated list of DoF range: minimum, maximum, step size" )
		+ XMLSchemaAttribute( "values_trans_axis_1", xsct_real_cslist, "Comma-separated list of DoF values to sample" )
		+ XMLSchemaAttribute( "values_trans_axis_2", xsct_real_cslist, "Comma-separated list of DoF values to sample" )
		+ XMLSchemaAttribute( "values_trans_axis_3", xsct_real_cslist, "Comma-separated list of DoF values to sample" )
		+ XMLSchemaAttribute( "values_rot_axis_1", xsct_real_cslist, "Comma-separated list of DoF values to sample" )
		+ XMLSchemaAttribute( "values_rot_axis_2", xsct_real_cslist, "Comma-separated list of DoF values to sample" )
		+ XMLSchemaAttribute( "values_rot_axis_3", xsct_real_cslist, "Comma-separated list of DoF values to sample" )
		+ XMLSchemaAttribute::attribute_w_default( "degree_check" , xsct_rosetta_bool, "True or false to turn on or turn off degree_check", "true")
		+ XMLSchemaAttribute::attribute_w_default( "max_samples", xsct_non_negative_integer, "Mover will throw an error if the actual number of samples is greater than the user-proved value", "10000" );
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "dof_residue_selector_1",
		"Selector used to indicate which residue should be used during x, y, z coordinate system set up. This residue is generally selected on the static protein." );
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "dof_residue_selector_2a",
		"Selector used to indicate which residue should be used during x, y, z coordinate system set up. This residue is generally selected on the moving protein." );
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "dof_residue_selector_2b",
		"Selector used to indicate which residue should be used during x, y, z coordinate system set up. This residue is generally selected on the moving protein." );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "The Six Degree of Freedom Grid Dock mover utilizes the basic six degrees of freedom (DoF) to dock the protein in 3D space: translating and rotating one domain around three user-definied axes, commonly known as the x-axis, y-axis, and z-axis. The user defines both the coordinate frame and the sample sapce allowing for a user-friendly and highly customizable docking system. When complete, this mover generates output pdb files of each newly docked pose.", attlist );
}


//Following defines a function named parse_range that takes the user given minimum, maximum, and step-size to create the full sampling space. It then takes this full sample space and stores it within dof_vector (vector of vectors that stores the sample space for all six dofs). parse_range takes two inputs: a tag and a string named range_option_name. The tag parameter stores the minimum, maximum, and step-size information for all 6 dofs, and the string range_option_name specifies which dof to apply parse_range() to.
utility::vector1< core::Real >
SixDoFGridDockMover::parse_range( utility::tag::TagCOP tag, std::string range_option_name ) { //function arguments are tag containg the minimum, maximum, and step size AND a string specifying which DoF to apply to.
	utility::vector1< core::Real > dof_range_input = utility::string_split( tag->getOption< std::string >( range_option_name ), ',', core::Real() );
	if ( dof_range_input.size() != 3 ) {
		utility_exit_with_message("ERROR: User has not provided three values for the range option as: minimum, maximum, and step size.");
	}
	TR << "tag " << tag << std::endl;
	TR << "range_option_name " << range_option_name << std::endl;
	core::Real dof_range_min = dof_range_input[1];
	core::Real dof_range_max = dof_range_input[2];
	core::Real dof_range_ss = dof_range_input[3];
	if ( dof_range_ss == 0 ) {
		utility_exit_with_message("ERROR: User has entered 0 for the step size. Change step size to a non-zero value.");
	}
	if ( dof_range_ss < 0 ) {
		utility_exit_with_message("ERROR: User has entered a value less than 0 for the step size. Change step size to a positive value greater than zero.");
	}
	TR << "User provided range for DOF: " << tag->getOption< std::string >( range_option_name ) << std::endl;
	TR << "Minimum: " << dof_range_input[1] << std::endl;
	TR << "Maximum: "  << dof_range_input[2] << std::endl;
	TR << "Step Size: " << dof_range_input[3] << std::endl;
	utility::vector1< core::Real > dof_vector;
	for ( core::Real jj = dof_range_min; jj <= dof_range_max; jj += dof_range_ss ) {
		dof_vector.push_back(jj);
	}
	return dof_vector;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
SixDoFGridDockMover::fresh_instance() const
{
	return utility::pointer::make_shared< SixDoFGridDockMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
SixDoFGridDockMover::clone() const
{
	return utility::pointer::make_shared< SixDoFGridDockMover >( *this );
}

std::string SixDoFGridDockMover::get_name() const {
	return mover_name();
}

std::string SixDoFGridDockMover::mover_name() {
	return "SixDoFGridDockMover";
}

core::Size
SixDoFGridDockMover::dof_sample_index() const
{
	return dof_sample_index_;
}

void
SixDoFGridDockMover::dof_sample_index(core::Size setting)
{
	dof_sample_index_ = setting;
}

//Following defines a function named lex_position, which returns the specific dof value to apply to a pose when given a dof_values (vector of vectors), the job number, and the desired dof.
core::Real
SixDoFGridDockMover::lex_position( utility::vector1< utility::vector1 < core::Real > > dof_values, core::Size index, core::Size dof )
//function arguments are vector of vectors containing the sample space, job ID, the specific degree of freedom.
// specific degree of freedom order are:
// 1: rot axis 1, 2: rot axis 2, 3: rot axis 3, 4: trans axis 1, 5: trans axis 2, 6: trans axis 3
{
	utility::vector1< core::Size > dimsizes(6, 0);
	for ( core::Size ii = 1; ii <= 6; ++ii ) {
		dimsizes[ ii ] = dof_values[ ii ].size();
	}

	utility::LexicographicalIterator lex(dimsizes);
	std::cout << "Dim sizes:";
	for ( core::Size ii = 1; ii <= 6; ++ii ) {
		std::cout << " " << dimsizes[ii];
	}
	std::cout << std::endl;
	lex.set_position_from_index( index );
	TR << "dof_values_: " << dof_values << std::endl;
	for ( core::Size ii =1; ii <= 6; ++ii ) {
		TR << "index at " << ii << " index: " << lex [ii] << " dof value: " << dof_values[ii][ lex[ ii ] ] << std::endl;
	}
	return dof_values[ dof ][ lex[ dof ]];
}


/////////////// Creator ///////////////

protocols::moves::MoverOP
SixDoFGridDockMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< SixDoFGridDockMover >();
}

std::string
SixDoFGridDockMoverCreator::keyname() const
{
	return SixDoFGridDockMover::mover_name();
}

void SixDoFGridDockMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SixDoFGridDockMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, SixDoFGridDockMover const & mover )
{
	mover.show(os);
	return os;
}



} //simple_moves
} //protocols
