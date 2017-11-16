// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/movers/FastDesign.cc
/// @brief Similar to FastRelax, but with design and a few additional options.
/// @details Carries out alternating rounds of design (packing) and minimization, ramping up the repulsive term from a very low value from one iteration to the next.
/// @author Tom Linsky
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- Added support for D-amino acids.

//Unit Headers
#include <protocols/denovo_design/movers/FastDesign.hh>
#include <protocols/denovo_design/movers/FastDesignCreator.hh>

//Protocol Headers
#include <protocols/constraint_generator/AddConstraints.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>
#include <protocols/constraint_generator/RemoveConstraints.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>

//Core Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>

//Basic Headers
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/tag/Tag.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

//C++ Headers

#ifdef GL_GRAPHICS
#include <protocols/viewer/viewers.hh>
#endif

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

static basic::Tracer TR( "protocols.denovo_design.movers.FastDesign" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace movers {
////////////////////////////////////////////////////////////////////////////////////////////////////

// XRW TEMP std::string
// XRW TEMP FastDesignCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return FastDesign::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP FastDesignCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new FastDesign() );
// XRW TEMP }

///  ---------------------------------------------------------------------------------
///  FastDesign main code:
///  ---------------------------------------------------------------------------------

/// @brief default constructor
FastDesign::FastDesign():
	FastRelax(),
	clear_designable_residues_( false ),
	run_count_( 0 ),
	cgs_()
{
	set_enable_design( true );
	//read_script_file( "", default_repeats_ );
}

/// @brief Constructor with some options
///
FastDesign::FastDesign(
	core::scoring::ScoreFunctionOP scorefxn_in,
	core::Size standard_repeats
) :
	FastRelax( scorefxn_in, standard_repeats ),
	clear_designable_residues_(false),
	run_count_(0),
	cgs_()
{
	set_enable_design( true );
}


/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
FastDesign::~FastDesign()
{}

// XRW TEMP std::string
// XRW TEMP FastDesign::mover_name()
// XRW TEMP {
// XRW TEMP  return "FastDesign";
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP FastDesign::get_name() const
// XRW TEMP {
// XRW TEMP  return FastDesign::mover_name();
// XRW TEMP }

/// @brief return an blank version of ourselves
protocols::moves::MoverOP
FastDesign::fresh_instance() const
{
	return protocols::moves::MoverOP( new FastDesign );
}

/// Return a copy of ourselves
protocols::moves::MoverOP
FastDesign::clone() const
{
	return protocols::moves::MoverOP( new FastDesign(*this) );
}

/// @brief Create the default task factory.  Must be called before design can occur.
void
FastDesign::set_up_default_task_factory() {
	set_task_factory( create_default_task_factory() );
}


void
FastDesign::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	// make sure we create a task factory before parsing FastRelax::parse_my_tag
	// otherwise no design will occur
	set_up_default_task_factory();
	FastRelax::parse_my_tag( tag, data, filters, movers, pose );

	clear_designable_residues_ = tag->getOption< bool >( "clear_designable_residues", clear_designable_residues_ );

	// parse constraint generators
	utility::vector1< std::string > const cgs = utility::string_split( tag->getOption< std::string >( "cgs", "" ), ',' );
	for ( utility::vector1< std::string >::const_iterator cg=cgs.begin(); cg!=cgs.end(); ++cg ) {
		if ( cg->empty() ) continue;
		protocols::constraint_generator::ConstraintGeneratorCOP new_cg =
			data.get_ptr< protocols::constraint_generator::ConstraintGenerator const >( "CONSTRAINT_GENERATORS", *cg );
		if ( !new_cg ) {
			std::stringstream msg;
			msg << "FastDesign: Could not find a constraint generator named " << *cg << " in the data map.  Ensure it has been defined in an AddConstraints mover before being referenced by FastDesign."
				<< std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
		}
		cgs_.push_back( new_cg->clone() );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
FastDesign::apply( core::pose::Pose & pose )
{
	// increment run counter
	++run_count_;

	TR.Debug   << "========================== FastDesign Run: " << run_count_ << " ===================" << std::endl;
	reset_status();

	// create and print task
	core::pack::task::PackerTaskOP clear_task = get_task_factory()->create_task_and_apply_taskoperations( pose );
	clear_task->show( TR );

	// print movemap
	if ( get_movemap() && TR.Debug.visible() ) {
		TR.Debug << "Movemap: ";
		get_movemap()->show( TR.Debug );
		TR.Debug << std::endl;
	}

	// if requested, clear all residues marked as designable
	if ( clear_designable_residues_ ) {
		TR << "Clearing designable residues...";
		for ( core::Size resid=1; resid<=pose.size(); ++resid ) {
			if ( pose.residue( resid ).is_protein() && clear_task->being_designed( resid ) ) {
				if ( pose.residue( resid ).type().is_alpha_aa() ) {
					if ( pose.residue( resid ).type().is_l_aa() ) { //Note that glycine is skipped, since it's not an L-amino acid.
						// skip proline
						if ( pose.residue( resid ).type().aa() != core::chemical::aa_pro ) {
							TR << resid << "...";
							// mutate to alanine
							protocols::simple_moves::MutateResidue mut_res( resid, core::chemical::aa_ala );
							mut_res.apply( pose );
						} else {
							TR << "<" << resid << ">...";
						}
					} else if ( pose.residue( resid ).type().is_d_aa() ) {
						// skip D-proline
						if ( pose.residue( resid ).type().aa() != core::chemical::aa_dpr ) {
							TR << resid << "...";
							// mutate to D-alanine
							protocols::simple_moves::MutateResidue mut_res( resid, core::chemical::aa_dal );
							mut_res.apply( pose );
						} else {
							TR << "<" << resid << ">...";
						}
					} else { //Achiral case
						TR << "<" << resid << ">...";
					}
				}
			}
		}
		TR << std::endl;
	}

	// support for ramping reference weights
	modify_scripts_for_alternative_scorefunctions();

	FastRelax::apply( pose );

	// show final scores
	core::scoring::ScoreFunctionOP local_sfxn( get_scorefxn()->clone() );  //May need to be modified for symmetry
	core::pose::symmetry::make_score_function_consistent_with_symmetric_state_of_pose( pose, local_sfxn ); //If the pose is symmetric but the scorefunction is not, remedy this.
	local_sfxn->show( TR, pose );
	TR.flush();
}

/// @brief sets constraint weights -- used with constraint ramping
void
FastDesign::set_constraint_weight(
	core::scoring::ScoreFunctionOP local_scorefxn,
	core::scoring::EnergyMap const & full_weights,
	core::Real const weight,
	core::pose::Pose & pose ) const
{
	runtime_assert( local_scorefxn != 0 );
	if ( cgs_.size() ) {
		protocols::constraint_generator::RemoveConstraints rm_csts( cgs_, true );

		try {
			rm_csts.apply( pose );
		} catch ( protocols::constraint_generator::EXCN_RemoveCstsFailed const & e ) {
			// if removing constraints fails, we don't really care
			// they are only being removed to clean up before re-adding them below
		}

		if ( weight > 1e-6 ) {
			protocols::constraint_generator::AddConstraints( cgs_ ).apply( pose );
		}

	} else {
		local_scorefxn->set_weight( core::scoring::coordinate_constraint, full_weights[ core::scoring::coordinate_constraint ] * weight );
		local_scorefxn->set_weight( core::scoring::atom_pair_constraint, full_weights[ core::scoring::atom_pair_constraint ] * weight );
		local_scorefxn->set_weight( core::scoring::angle_constraint, full_weights[ core::scoring::angle_constraint ] * weight );
		local_scorefxn->set_weight( core::scoring::dihedral_constraint, full_weights[ core::scoring::dihedral_constraint ] * weight );
		TR << "[coordinate:atom_pair:angle:dihedral] = "
			<< local_scorefxn->get_weight( core::scoring::coordinate_constraint )
			<< " : " << local_scorefxn->get_weight( core::scoring::atom_pair_constraint )
			<< " : " << local_scorefxn->get_weight( core::scoring::angle_constraint )
			<< " : " << local_scorefxn->get_weight( core::scoring::dihedral_constraint ) << std::endl;
	}
}

core::pack::task::TaskFactoryOP
FastDesign::create_default_task_factory() const
{
	using core::pack::task::operation::TaskOperationCOP;

	core::pack::task::TaskFactoryOP local_tf( new core::pack::task::TaskFactory() );
	if ( get_task_factory() ) {
		local_tf = get_task_factory()->clone();
	} else {
		// add command line things
		local_tf->push_back(TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline() ));

		// add resfile if requested by flags
		if ( basic::options::option[ basic::options::OptionKeys::relax::respect_resfile]() &&
				basic::options::option[ basic::options::OptionKeys::packing::resfile].user() ) {
			local_tf->push_back(TaskOperationCOP( new core::pack::task::operation::ReadResfile() ));
			TR << "Using Resfile for packing step." << std::endl;
		}

		// turn off packing on residues that can't move
		// (Note that this function *doesn't* see any MoveMap info that's set via parse_my_tag().)
		if ( get_movemap() ) {
			core::pack::task::operation::PreventRepackingOP turn_off_packing( new core::pack::task::operation::PreventRepacking() );
			core::kinematics::MoveMap const & mm = *get_movemap();
			for ( core::kinematics::MoveMap::MoveMapTorsionID_Map::const_iterator mm_torsion=mm.movemap_torsion_id_begin();
					mm_torsion!=mm.movemap_torsion_id_end(); ++mm_torsion ) {
				// skip if this residue and torsion are allowed to move
				if ( mm_torsion->second ) continue;
				core::kinematics::MoveMap::MoveMapTorsionID const & torsion_id = mm_torsion->first;

				// skip if this is not chi
				if ( torsion_id.second != core::id::CHI ) continue;

				turn_off_packing->include_residue( torsion_id.first );
			}
			local_tf->push_back( turn_off_packing );
		}
	}

	//Include current rotamer by default - as before.
	local_tf->push_back( TaskOperationCOP( new core::pack::task::operation::IncludeCurrent() ) );

	if ( limit_aroma_chi2() ) {
		local_tf->push_back(TaskOperationCOP( new protocols::toolbox::task_operations::LimitAromaChi2Operation() ));
	}

	return local_tf;
}

void
FastDesign::modify_scripts_for_alternative_scorefunctions()
{
	using namespace basic::options;


	// will attempt to modify the relax script for -beta only
	// options for any other score functions could be also added below
	std::vector< std::string > filelines;

	// Now beta_nov15 as default... this is pretty dangerous for talaris
	if ( !(FastRelax::script_file_specified_) ) {
		/* // no more support for beta_nov15_patch
		if ( option[ OptionKeys::corrections::beta_nov15_patch ]() ) {
		filelines.push_back( "reference  2.01692 3.95229 -1.84895 -2.44909 1.54388 1.43603 0.25816 2.70992 -0.38208 2.00235 2.31398 -0.91852 -0.67964 -0.97481 -0.11701 -1.53805 -1.70469 2.85306 2.72731 -0.99943" ); // using "minimized context"
		filelines.push_back( "ramp_repack_min 1     0.00001  0.0"      );
		filelines.push_back( "accept_to_best"                  );
		filelines.push_back( "endrepeat "                      );
		} else
		*/
		if ( option[ OptionKeys::corrections::beta_nov16 ]() || option[ OptionKeys::corrections::beta_nov16_cart ] ) {
			TR << "Calling correction for beta_nov16, " << FastRelax::default_repeats() << " repeats..." << std::endl;
			filelines.push_back( "repeat "+ObjexxFCL::string_of(FastRelax::default_repeats()));

			filelines.push_back( "reference 0.3     3.1     -2.6     -2.55    4.8     -0.5    0.7      4.5     -1.6     4.0     3.9     -1.7     -2.0     -1.5     -1.0     -2.0    -2.0     4.0     9.0     3.7" );
			filelines.push_back( "ramp_repack_min 0.02  0.01     1.0"      );
			filelines.push_back( "reference 2.2619  4.8148  -1.6204  -1.6058  2.7602  1.0350  1.3406   2.5006  -0.6895  1.9223  2.3633  -0.3009  -4.2787   0.1077   0.0423  -0.4390 -0.7333  3.2371  4.7077  2.3379" );
			filelines.push_back( "ramp_repack_min 0.250 0.01     0.5"      );
			filelines.push_back( "reference 2.2619  4.5648  -1.6204  -1.6158  2.5602  1.1350  1.2406   2.3006  -0.7895  1.7223  2.1633  -0.3009  -4.3787   0.1077   0.0423  -0.4390 -0.7333  3.1371  4.4077  2.1379" );
			filelines.push_back( "ramp_repack_min 0.550 0.01     0.0"      );
			filelines.push_back( "reference 2.2619  4.3148  -1.6204  -1.6358  1.9602  1.4350  0.8406   1.8006  -0.8895  1.3223  1.4633  -0.3009  -4.6787  -0.1077  -0.1423  -0.5390 -0.9333  2.7371  3.7077  1.7379" );
			filelines.push_back( "ramp_repack_min 1     0.00001  0.0"      );
			filelines.push_back( "accept_to_best"                  );
			filelines.push_back( "endrepeat "                      );
		}
	}

	if ( filelines.size() > 0 ) {
		FastRelax::set_script_from_lines( filelines );
	}
}

std::string FastDesign::get_name() const {
	return mover_name();
}

std::string FastDesign::mover_name() {
	return "FastDesign";
}

void FastDesign::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaComplexTypeGeneratorOP ct_gen = FastRelax::complex_type_generator_for_fast_relax( xsd );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "clear_designable_residues", xsct_rosetta_bool, "Clear the set of designable residues?")
		+ XMLSchemaAttribute( "cgs", xs_string, "Names of previously specified constraint generators to apply during this run");

	ct_gen->add_attributes( attlist )
		.element_name( mover_name() )
		.description( "FastRelax mover used for design that can take constraint generators" )
		.write_complex_type_to_schema( xsd );


}

std::string FastDesignCreator::keyname() const {
	return FastDesign::mover_name();
}

protocols::moves::MoverOP
FastDesignCreator::create_mover() const {
	return protocols::moves::MoverOP( new FastDesign );
}

void FastDesignCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FastDesign::provide_xml_schema( xsd );
}


} // namespace movers
} // namespace denovo_design
} // namespace protocols
