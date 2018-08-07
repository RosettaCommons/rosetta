// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/BundleGridSampler.cc
/// @brief  This mover samples conformations of a helical bundle by performing a grid search in Crick parameter space.
/// @details This mover calls the MakeBundle mover (which in turn calls the MakeBundleHelix mover).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/helical_bundle/BundleGridSampler.hh>
#include <protocols/helical_bundle/BundleGridSamplerCreator.hh>
#include <protocols/helical_bundle/BundleGridSamplerHelper.hh>
#include <protocols/helical_bundle/MakeBundleHelix.hh>
#include <protocols/helical_bundle/MakeBundle.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <utility/tag/Tag.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/scoring/Energies.hh>

#include <numeric/constants.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <numeric/random/random.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <core/conformation/parametric/Parameter.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>
#include <core/scoring/rms_util.hh>
#include <core/pose/util.tmpl.hh>
#include <protocols/helical_bundle/util.hh>
#include <utility/pointer/memory.hh>

//JD2:
#include <protocols/jd2/util.hh>

// Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>

// C math
#include <cmath>

// C output headers:
#include <cstdio>
#include <sstream>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using basic::Error;
using basic::Warning;

namespace protocols {
namespace helical_bundle {

static basic::Tracer TR("protocols.helical_bundle.BundleGridSampler");

/// @brief Creator for BundleGridSampler mover.
BundleGridSampler::BundleGridSampler():
	Mover("BundleGridSampler"),
	reset_mode_(true),
	nstruct_mode_(false),
	nstruct_mode_repeats_(1),
	select_low_(true),
	n_helices_(0),
	max_samples_( 10000 ),
	make_bundle_( utility::pointer::make_shared< MakeBundle >() ),
	pre_selection_mover_(),
	pre_selection_mover_exists_(false),
	pre_selection_filter_(),
	pre_selection_filter_exists_(false),
	dump_pdbs_(false),
	pdb_prefix_("bgs_out"),
	sfxn_set_(false),
	sfxn_(),
	default_calculator_( make_bundle_->default_calculator_nonconst() )
{
	debug_assert( make_bundle_ != nullptr );
	debug_assert( default_calculator_ != nullptr );
}


/// @brief Copy constructor for BundleGridSampler mover.
BundleGridSampler::BundleGridSampler( BundleGridSampler const & src ):
	protocols::moves::Mover( src ),
	reset_mode_(src.reset_mode_),
	nstruct_mode_(src.nstruct_mode_),
	nstruct_mode_repeats_(src.nstruct_mode_repeats_),
	select_low_(src.select_low_),
	n_helices_(src.n_helices_),
	max_samples_(src.max_samples_),
	make_bundle_(  utility::pointer::dynamic_pointer_cast< MakeBundle >(src.make_bundle_->clone()) ),
	pre_selection_mover_( src.pre_selection_mover_ ), //NOTE that we're not cloning this mover, but using it straight
	pre_selection_mover_exists_(src.pre_selection_mover_exists_),
	pre_selection_filter_( src.pre_selection_filter_ ), //NOTE that we're not cloning this filter, but using it straight
	pre_selection_filter_exists_(src.pre_selection_filter_exists_),
	dump_pdbs_(src.dump_pdbs_),
	pdb_prefix_(src.pdb_prefix_),
	sfxn_set_(src.sfxn_set_),
	sfxn_(src.sfxn_), //NOTE that this is also copied without cloning
	default_calculator_( make_bundle_->default_calculator_nonconst() )
{
	debug_assert( make_bundle_ != nullptr );
	debug_assert( default_calculator_ != nullptr );
}


/// @brief Destructor for BundleGridSampler mover.
BundleGridSampler::~BundleGridSampler() = default;


/// @brief Clone operator to create a pointer to a fresh BundleGridSampler object that copies this one.
protocols::moves::MoverOP BundleGridSampler::clone() const {
	return protocols::moves::MoverOP( new BundleGridSampler( *this ) );
}


/// @brief Fresh_instance operator to create a pointer to a fresh BundleGridSampler object that does NOT copy this one.
protocols::moves::MoverOP BundleGridSampler::fresh_instance() const {
	return protocols::moves::MoverOP( new BundleGridSampler );
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
void BundleGridSampler::apply (core::pose::Pose & pose)
{
	if ( TR.visible() ) TR << "Performing Crick parameter space grid sampling." << std::endl;

	runtime_assert_string_msg( sfxn_set() && sfxn_, "In BundleGridSampler::apply(): no scorefunction has been set for this mover!" );

	core::pose::Pose lowestEpose(pose);
	core::Size const nhelices( n_helices() ); //The number of helices that have been defined.
	bool at_least_one_success(false); //Has at least one bundle been generated successfully?  False at the end if all parameter values tried result in nonsensical bundles.
	bool first_bundle(true); //Is this the first bundle generated?
	core::Real lowest_energy(0); //What's the lowest (or highest) energy encountered so far?

	std::stringstream final_report(""); //The summary of the parameters for the lowest-energy state.  This will be written out at the end.

	core::Size const total_samples( calculate_total_samples() );
	if ( TR.visible() ) TR << "A total of " << total_samples << " grid points will be sampled." << std::endl;
	runtime_assert_string_msg(total_samples <= max_samples(), "The total number of grid samples exceeds the maximum allowed!  (Note that you can increase the maximum allowed with the max_samples flag).");

	//Set up the BundleGridSamplerHelper object:
	if ( TR.visible() ) TR << "Creating BundleGridSamplerHelper object." << std::endl;
	BundleGridSamplerHelperOP sampler_helper( new BundleGridSamplerHelper );
	sampler_helper->reset();

	//Set up the DoFs over which we'll be sampling:
	for ( core::Size iparam(1); iparam <= static_cast<core::Size>( BPC_last_parameter_to_be_sampled ); ++iparam ) {
		for ( core::Size ihelix(1); ihelix<=nhelices; ++ihelix ) { //Loop through all helices.
			core::conformation::parametric::RealValuedParameterCOP curparam( utility::pointer::dynamic_pointer_cast< core::conformation::parametric::RealValuedParameter const >( make_bundle_->helix_cop(ihelix)->calculator_cop()->parameter_cop(iparam) ) );
			if ( curparam == nullptr ) continue;
			if ( curparam->can_be_sampled() && curparam->sampling_set() ) {
				sampler_helper->add_DoF( static_cast<BPC_Parameters>(iparam), ihelix, curparam->value_samples(), curparam->value_min(), curparam->value_max() );
			}
		}
	}

	//Initialize the BundleGridSamplerHelper object (i.e. perform the internal calculation that pre-generates all of the values to be sampled).
	sampler_helper->initialize_samples();

	//Loop through all grid samples
	core::Size loopstart(1);
	core::Size loopend(total_samples);
	if ( nstruct_mode() ) { //Special case: if we're just doing one set of Crick parameters per job, we need to figure out which set to do.
		if ( !protocols::jd2::jd2_used() ) {
			utility_exit_with_message(
				"In protocols::helical_bundle::BundleGridSampler::apply() function: The nstruct_mode option was used, but the current application is not using JD2.");
		}
		core::Size curjob( protocols::jd2::current_nstruct_index() );
		core::Size totaljobs( protocols::jd2::max_nstruct_index() );
		if ( curjob==0 || totaljobs==0 ) {
			utility_exit_with_message(
				"In protocols::helical_bundle::BundleGridSampler::apply() function: The nstruct_mode option was used, but invalid values were obtained for the current job index or the total number of jobs.");
		}
		if ( totaljobs < total_samples*nstruct_repeats() && TR.Warning.visible() ) {
			TR.Warning << "The BundleGridSampler mover is in nstruct mode, meaning that one set of Crick parameters will be sampled per job.  However, the total number of jobs is less than the total number of samples!  Certain sets of Crick parameters will be missed!" << std::endl ;
		}
		//The current job might be greater than the total number of samples, in which case we should wrap around:
		loopstart = ( ( (curjob-1) % total_samples) + 1 ) / nstruct_repeats();
		loopend=loopstart;
	}

	for ( core::Size i=loopstart; i<=loopend; ++i ) {
		core::pose::Pose temppose;
		if ( !reset_mode() ) temppose=pose; //If we're not resetting the pose, then we need to copy the input pose.

		//Making a copy of the MakeBundle mover:
		MakeBundleOP makebundle_copy( utility::pointer::dynamic_pointer_cast< MakeBundle>( make_bundle_->clone() ) );

		if ( nstruct_mode() ) { //If this is one-set-of-Crick-params-per-job mode
			if ( i > 1 ) {
				//If i is greater than 1, increment repeatedly until we reach the current permutation.
				//(Note that we start out on the first permutation, so there's no need to increment the first time).
				for ( core::Size j=2; j<=i; ++j ) {
					sampler_helper->increment_cur_indices(); //This is a recursive function that increments indices in a non-trivial way, and so must be called repeatedly to get the right indices.
				}
			}
		} else { //If this is NOT one-set-of-Crick-params-per-job mode
			//If i is greater than 1, increment the current permutation.
			//(Note that we start out on the first permutation, so there's no need to increment the first time).
			if ( i > 1 ) sampler_helper->increment_cur_indices();
		}

		//Set the parameters that are being varied:
		for ( core::Size j=1, jmax=sampler_helper->nDoFs(); j<=jmax; ++j ) {
			//Get an owning pointer to the helix that this DoF is part of:
			MakeBundleHelixOP curhelix( makebundle_copy->helix( sampler_helper->DoF_helix_index(j) ) );
			runtime_assert_string_msg(curhelix != nullptr, "Error in getting owning pointer to current helix in BundleGridSampler::apply() function.");

#ifdef NDEBUG
			//Release mode: static cast for speed.
			core::conformation::parametric::RealValuedParameterOP curparam( utility::pointer::static_pointer_cast< core::conformation::parametric::RealValuedParameter >( curhelix->calculator_op()->parameter( sampler_helper->DoF_type(j) ) ) );
#else
			//Debug mode: dynamic cast and check.
			core::conformation::parametric::RealValuedParameterOP curparam( utility::pointer::dynamic_pointer_cast< core::conformation::parametric::RealValuedParameter >( curhelix->calculator_op()->parameter( sampler_helper->DoF_type(j) ) ) );
			debug_assert( curparam != nullptr );
#endif
			curparam->set_value( sampler_helper->DoF_sample_value(j), true );
		}

		//Set the parameters that are copying other parameters:
		bool loopthrough_failed(false); //If we fail to copy another parameter (e.g. pitch angle), we need to know it.
		for ( core::Size j=1, jmax=n_helices(); j<=jmax; ++j ) {
			MakeBundleHelixOP curhelix( makebundle_copy->helix( j ) );
			debug_assert(curhelix != nullptr);

			BundleParametrizationCalculatorOP cur_calculator( curhelix->calculator_op() );
			debug_assert(cur_calculator != nullptr);

			for ( core::Size iparam(1); iparam <= static_cast<core::Size>(BPC_last_parameter_to_be_sampled); ++iparam ) {
				core::conformation::parametric::RealValuedParameterOP curparam( cur_calculator->real_parameter(iparam) );
				if ( curparam == nullptr ) continue;
				if ( curparam->can_be_copied() && curparam->copying_information_was_set() ) {
					debug_assert( !curparam->sampling_set() );
					debug_assert( !curparam->value_was_set() );
					core::conformation::parametric::ParametersCOP other_parameters( makebundle_copy->helix_cop( curparam->copy_from_parameters_index() )->calculator_cop()->parameters_cop() );
					debug_assert(other_parameters != nullptr);
					loopthrough_failed = curparam->copy_value_from_parameter( other_parameters->parameter_cop( iparam ), other_parameters, cur_calculator->parameters_cop() );
					if ( loopthrough_failed && TR.visible() ) {
						TR << "Failed to copy " << curparam->parameter_name() << " from helix " << curparam->copy_from_parameters_index() << " to helix " << j << "." << std::endl;
					}
				}
				if ( loopthrough_failed ) break;
			}
			if ( loopthrough_failed ) break;
		}
		if ( loopthrough_failed ) {
			if ( TR.visible() ) {
				TR << "Failed to copy at least one parameter.  Continuing to next grid point to sample." << std::endl;
				TR.flush();
			}
			continue; //If we failed to copy another parameter (e.g. pitch angle), continue to the next grid point.
		}

		//Output the parameters being attempted:
		if ( TR.visible() ) {
			TR << "Grid point " << i << ": attempting to build a helical bundle with the following parameters:" << "\n";
			TR << "Helix";
			for ( core::Size iparam(1); iparam<=static_cast<core::Size>(BPC_last_parameter_to_be_sampled); ++iparam ) {
				TR << "\t" << BundleParametrizationCalculator::parameter_name_from_enum(static_cast<BPC_Parameters>( iparam ) );
			}
			TR << "\n";
			core::Size const old_precision( TR.precision() );
			TR.precision( 4 );
			for ( core::Size j=1, jmax=n_helices(); j<=jmax; ++j ) {
				TR << j;
				MakeBundleHelixCOP curhelix( makebundle_copy->helix_cop( j ) );
				debug_assert(curhelix != nullptr);
				BundleParametrizationCalculatorCOP curcalculator( curhelix->calculator_cop() );
				debug_assert( curcalculator != nullptr );
				for ( core::Size iparam(1); iparam<=static_cast<core::Size>(BPC_last_parameter_to_be_sampled); ++iparam ) {
					core::conformation::parametric::RealValuedParameterCOP curparam( curcalculator->real_parameter_cop(iparam) );
					if ( curparam == nullptr ) continue;
					TR << "\t" << curparam->value();
				}
				TR << std::endl;
			}
			TR.precision(old_precision);
			TR.flush();
		}

		//Actually make the bundle:
		makebundle_copy->apply( temppose );

		//Check for success or failure:
		if ( makebundle_copy->last_apply_failed() ) {
			if ( TR.visible() ) {
				TR << "The current set of parameter values failed to produce a sensible helical bundle.  Continuing to next grid point to sample." << std::endl;
				TR.flush();
			}
			continue; //Go on to the next grid sample.
		} else {
			if ( TR.visible() ) {
				TR << "Bundle generated successfully." << std::endl;
				TR.flush();
			}
		}

		//Apply the preselection mover, if it exists:
		if ( preselection_mover_exists() ) {
			if ( TR.visible() ) {
				TR << "Applying preselection mover to generated helical bundle." << std::endl;
				TR.flush();
			}
			pre_selection_mover_->apply(temppose);
			if ( pre_selection_mover_->get_last_move_status() != protocols::moves::MS_SUCCESS ) {
				if ( TR.visible() ) {
					TR << "Pre-selection mover failed.  Discarding current grid sample and moving on to next." << std::endl;
					TR.flush();
				}
				continue;
			}
		}

		//Apply the preselection filter, if it exists:
		if ( preselection_filter_exists() ) {
			if ( TR.visible() ) {
				TR << "Applying filter to generated helical bundle." << std::endl;
			}
			bool const filterpassed( pre_selection_filter_->apply( temppose ) ); //Apply the filter and store the result in filterpassed.
			if ( TR.visible() ) {
				if ( filterpassed ) TR << "Filter passed!" << std::endl;
				else TR << "Filter failed!  Discarding current grid sample and moving on to next." << std::endl;
			}
			if ( !filterpassed ) continue; //Go on to the next grid sample.
		}

		at_least_one_success = true; //At this point, we've successfully generated at least one helical bundle.

		// Score the pose:
		(*sfxn_)( temppose);
		if ( first_bundle || ( selection_low() && temppose.energies().total_energy() < lowest_energy ) || ( !selection_low() && temppose.energies().total_energy() > lowest_energy) ) {
			first_bundle=false;
			lowest_energy = temppose.energies().total_energy();
			lowestEpose = temppose;
			if ( TR.visible() ) {
				char outstring[1024];
				std::string lowest_highest( "lowest" );
				if ( !selection_low() ) lowest_highest = "highest";
				sprintf( outstring, "Energy is %.4f.  This is the %s encountered so far.", temppose.energies().total_energy(), lowest_highest.c_str());
				TR << outstring << std::endl;

				{
					// Generate the summary of parameters, in case this is ultimately the lowest-energy pose:
					final_report.str("");
					final_report << "Parameters yielding the ";
					if ( selection_low() ) final_report << "lowest"; else final_report << "highest";
					final_report << "-energy bundle:" << std::endl;
					\
						final_report << "Helix";
					for ( core::Size iparam(1); iparam<=static_cast<core::Size>(BPC_last_parameter_to_be_sampled); ++iparam ) {
						final_report << "\t" << BundleParametrizationCalculator::parameter_name_from_enum( static_cast<BPC_Parameters>( iparam ) );
					}
					final_report << "\n";
					core::Size const old_precision( final_report.precision() );
					final_report.precision( 4 );
					for ( core::Size j=1, jmax=n_helices(); j<=jmax; ++j ) {
						final_report << j;
						MakeBundleHelixCOP curhelix( makebundle_copy->helix_cop( j ) );
						debug_assert(curhelix != nullptr);
						BundleParametrizationCalculatorCOP curcalculator( curhelix->calculator_cop() );
						debug_assert( curcalculator != nullptr );
						for ( core::Size iparam(1); iparam<=static_cast<core::Size>(BPC_last_parameter_to_be_sampled); ++iparam ) {
							core::conformation::parametric::RealValuedParameterCOP curparam( curcalculator->real_parameter_cop(iparam) );
							if ( curparam == nullptr ) continue;
							final_report << "\t" << curparam->value();
						}
						final_report << std::endl;
					}
					final_report.precision(old_precision);
					final_report << std::endl;
				}
			}
		} else {
			if ( TR.visible() ) {
				char outstring[1024];
				sprintf( outstring, "Energy is %.4f.", temppose.energies().total_energy());
				TR << outstring << std::endl;
			}
		}

		//Dump a PDB file, if the user has specified that this mover should do so:
		if ( pdb_output() ) {
			char outfile[1024];
			sprintf( outfile, "%s_%05lu.pdb", pdb_prefix().c_str(), static_cast<unsigned long>(i) );
			if ( TR.visible() ) TR << "Writing " << outfile << std::endl;
			temppose.dump_scored_pdb( outfile, *sfxn_ );
		}
	} //Looping through all grid samples

	if ( at_least_one_success ) {
		pose=lowestEpose;
		if ( TR.visible() ) TR << "Success!  Returning lowest-energy pose." << std::endl << final_report.str();
		set_last_move_status( protocols::moves::MS_SUCCESS );
	} else {
		if ( TR.visible() ) TR << "No parameter values returning sensible helical bundles, or helical bundles passing filters, were sampled.  Mover failed." << std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );
	}

	if ( TR.visible() ) TR.flush();
	if ( TR.Warning.visible() ) TR.Warning.flush();
	if ( TR.Debug.visible() ) TR.Debug.flush();

	return;
}

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////

/// @brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
///
void
BundleGridSampler::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const &filters,
	protocols::moves::Movers_map const &movers,
	core::pose::Pose const & /*pose*/
) {

	if ( tag->getName() != "BundleGridSampler" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This should be impossible -- the tag name does not match the mover name.");
	}

	if ( TR.visible() ) TR << "Parsing options for BundleGridSampler (\"" << tag->getOption<std::string>("name" ,"") << "\") mover." << std::endl;

	//Global options for this mover:
	runtime_assert_string_msg( tag->hasOption("scorefxn" ), "In BundleGridSampler::parse_my_tag(): A \"scorefxn\" option must be specified!");
	set_sfxn(protocols::rosetta_scripts::parse_score_function( tag, "scorefxn", data_map )->clone()); // The scorefunction.
	if ( tag->hasOption("max_samples") ) {
		auto const val( tag->getOption<core::Size>("max_samples", 10000) );
		if ( TR.visible() ) TR << "Setting maximum number of samples to " << val << "." << std::endl;
		set_max_samples(val);
	}
	if ( tag->hasOption("selection_type") ) {
		std::string const val = tag->getOption<std::string>("selection_type", "");
		runtime_assert_string_msg( val=="high" || val=="low",
			"When parsing options for the BundleGridSampler mover, could not interpret the selection_type.  This must be set to \"high\" or \"low\"." );
		if ( TR.visible() ) TR << "Setting selection type to " << val << "." << std::endl;
		if ( val=="high" ) set_selection_low(false);
		else set_selection_low(true);
	}
	if ( tag->hasOption("pre_selection_mover") ) {
		protocols::moves::MoverOP curmover = protocols::rosetta_scripts::parse_mover( tag->getOption< std::string >( "pre_selection_mover" ), movers );
		set_preselection_mover(curmover);
		if ( TR.visible() ) TR << "BundleGridSampler mover \"" << tag->getOption< std::string >("name", "") << "\" has been assigned mover \"" << tag->getOption< std::string >("pre_selection_mover") << "\" as a pre-selection mover that will be applied to all generated bundles prior to energy evaluation." << std::endl;
	}
	if ( tag->hasOption("pre_selection_filter") ) {
		protocols::filters::FilterOP curfilter = protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "pre_selection_filter" ), filters );
		set_preselection_filter(curfilter);
		if ( TR.visible() ) TR << "BundleGridSampler mover \"" << tag->getOption< std::string >("name", "") << "\" has been assigned filter \"" << tag->getOption< std::string >("pre_selection_filter") << "\" as a pre-selection filter that will be applied to all generated bundles prior to energy evaluation." << std::endl;
	}
	if ( tag->hasOption("dump_pdbs") ) {
		bool const val = tag->getOption<bool>("dump_pdbs", false);
		set_pdb_output(val);
		if ( TR.visible() ) {
			if ( val ) TR << "Setting BundleGridSampler mover \"" << tag->getOption< std::string >("name", "") << "\" to write out PDB files." << std::endl;
			else TR << "No PDB output will occur from this mover." << std::endl;
		}
	}
	if ( tag->hasOption("pdb_prefix") ) {
		std::string const val = tag->getOption<std::string>("pdb_prefix", "bgs_out");
		set_pdb_prefix(val);
		if ( TR.visible() ) TR << "Setting prefix for PDB output to " << val << "." << std::endl;
	}

	//Determine whether we'll be parsing values in degrees or in radians:
	bool const use_degrees(tag->getOption<bool>("use_degrees", false));
	make_bundle_->set_use_degrees( use_degrees );
	if ( TR.visible() ) TR << "Setting mover to interpret angles in RosettaScripts as though they were specified in " << (make_bundle_->use_degrees() ? "degrees." : "radians.") << "  Note that internally and in output, the mover uses radians." << std::endl;

	//Set symmetry options for the MakeBundle mover (symmetry, symmetry_copies):
	make_bundle_->set_symmetry_options_from_tag( tag );

	//Set reset mode for the MakeBundle mover:
	bool const resetmode( tag->getOption<bool>("reset", true) );
	if ( TR.visible() ) TR << "Setting reset mode to " << (resetmode ? "true." : "false.") << std::endl;
	set_reset_mode(resetmode);
	make_bundle_->set_reset_pose( reset_mode() );

	//Set nstruct mode and options:
	bool nstructmode = tag->getOption<bool>("nstruct_mode", false);
	if ( TR.visible() ) TR << "Setting nstruct mode to " << (nstructmode ? "true." : "false.") << "  This means that " << (nstructmode ? "each job will sample a different set of Crick parameters." : "every job will sample all sets of Crick parameters.") << std::endl;
	set_nstruct_mode(nstructmode);
	auto nstructrepeats( tag->getOption<core::Size>( "nstruct_repeats", 1 ) );
	if ( nstructrepeats<1 ) nstructrepeats=1;
	if ( TR.visible() ) TR << "Setting nstruct repeats to " << nstructrepeats << "." << std::endl;
	set_nstruct_repeats(nstructrepeats);

	if ( tag->hasOption("crick_params_file") ) {
		set_default_crick_params_file( tag->getOption<std::string>("crick_params_file") );
	}

	std::string resname_string( tag->getOption<std::string>("residue_name", "ALA") );
	utility::vector1 < std::string > resname_vect;
	parse_resnames( resname_string, resname_vect );
	make_bundle_->set_default_residue_name( resname_vect );

	//Set default helix length:
	make_bundle_->set_default_helix_length( tag->getOption< core::Size >( "helix_length", 10 ) );

	//Read default options:
	for ( core::Size i(1); i < static_cast< core::Size >( BPC_end_of_list ); ++i ) {
		core::conformation::parametric::ParameterOP curparam( default_calculator_->parameter(i) );
		if ( curparam->can_be_set() || curparam->can_be_sampled() ) {
			curparam->parse_setting( tag, curparam->can_be_set(), curparam->can_be_sampled(), false, false );
		}
	}

	//Check that at least one helix is defined:
	bool at_least_one_helix = false;

	//Parse options for each helix (sub-tags):
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	for ( auto const & branch_tag : branch_tags ) { //Loop through all sub-tags
		if ( branch_tag->getName() == "Helix" ) { //A helix has been added.  Add it, and parse its options.
			at_least_one_helix = true;

			add_helix(); //This updates this mover and the MakeBundle mover
			core::Size const helix_index(n_helices());
			if ( TR.visible() ) TR << "Added a helix." << std::endl;

			//Set crick_params_file:
			if ( branch_tag->hasOption("crick_params_file") ) {
				make_bundle_->set_crick_params_file_for_helix( branch_tag->getOption<std::string>("crick_params_file"), helix_index );
			}

			if ( branch_tag->hasOption( "residue_name" ) ) {
				std::string resname_string( branch_tag->getOption<std::string>("residue_name") );
				utility::vector1 < std::string > resname_vect;
				parse_resnames( resname_string, resname_vect );
				make_bundle_->set_residue_name_for_helix( resname_vect, helix_index );
			}

			if ( branch_tag->hasOption( "helix_length" ) ) {
				make_bundle_->set_helix_length_for_helix( branch_tag->getOption<core::Size>("helix_length"), helix_index );
			}

			for ( core::Size i(1); i < static_cast< core::Size >( BPC_end_of_list ); ++i ) {
				core::conformation::parametric::ParameterOP curparam( make_bundle_->helix(helix_index)->calculator_op()->parameter(i) );
				if ( curparam->can_be_set() || curparam->can_be_sampled() || curparam->can_be_copied() ) {
					curparam->parse_setting( branch_tag, curparam->can_be_set(), curparam->can_be_sampled(), false, curparam->can_be_copied() );
				}
			}
		} //if ( (*tag_it)->getName() == "Helix" )
	}
	runtime_assert_string_msg(at_least_one_helix, "When parsing options for the BundleGridSampler mover, found no helices!  Define at least one helix with a <Helix...> sub-tag.");

	if ( TR.visible() ) TR.flush();
} //parse_my_tag

////////////////////////////////////////////////////////////////////////////////
//          PUBLIC FUNCTIONS                                                  //
////////////////////////////////////////////////////////////////////////////////


/// @brief Set whether we're using degrees (true) or radians (false)
void BundleGridSampler::set_use_degrees( bool const use_degrees ) {
	debug_assert( make_bundle_ != nullptr ); //Should be true
	make_bundle_->set_use_degrees(use_degrees);
}

/// @brief Get whether we're using degrees (true) or radians (false)
bool BundleGridSampler::use_degrees() const {
	debug_assert( make_bundle_ != nullptr ); //Should be true
	return make_bundle_->use_degrees();
}

/// @brief Add options for a new helix
/// @details Return value is the current total number of helices after the addition.
core::Size
BundleGridSampler::add_helix() {
	increment_helix_count();
	core::Size const nhelices(n_helices());

	//Update the MakeBundle mover, too:
	make_bundle_->add_helix();
	runtime_assert_string_msg( make_bundle_->n_helices()==nhelices,
		"In protocols::helical_bundle::BundleGridSampler::add_helix() function: somehow, the MakeBundle mover is out of sync with the BundleGridSampler mover.  Consult a developer or an exorcist -- this shouldn't happen.");

	return nhelices;
}

/// @brief Set the default Crick params file.
/// @details This is used unless overridden on a helix-by-helix basis.
/// @note Triggers a read from disk!
void
BundleGridSampler::set_default_crick_params_file(
	std::string const & default_crick_file
) {
	debug_assert( make_bundle_ != nullptr );
	debug_assert( default_calculator_ != nullptr );
	make_bundle_->set_default_crick_params_file( default_crick_file );
}

/// @brief Sets the mover that will be applied to all helical bundles generated prior to energy evaluation.
/// @details Note: if this is used, there is no guarantee that the resulting geometry will still lie within the
/// parameter space.  (That is, this mover could move the backbone.)
void BundleGridSampler::set_preselection_mover ( protocols::moves::MoverOP mover )
{
	pre_selection_mover_ = mover;
	pre_selection_mover_exists_ = true;
	return;
}

/// @brief Sets the filter that will be applied to all helical bundles generated prior to energy evaluation.
/// @details See the pre_selection_filter_ private member variable for details.
void BundleGridSampler::set_preselection_filter ( protocols::filters::FilterOP filter )
{
	pre_selection_filter_ = filter;
	pre_selection_filter_exists_ = true;
	return;
}


////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////

/// @brief Is a value in a list?
///
bool BundleGridSampler::is_in_list( core::Size const val, utility::vector1 < core::Size> const &list ) const {
	core::Size const listsize(list.size());
	if ( listsize==0 ) return false;
	for ( core::Size i=1; i<=listsize; ++i ) {
		if ( list[i]==val ) return true;
	}
	return false;
}

/// @brief Calculate the number of grid points that will be sampled, based on the options set by the user.
///
core::Size
BundleGridSampler::calculate_total_samples() const {
	core::Size total_samples(1);
	core::Size const helixcount( n_helices() );

	//Loop through all helices that are defined:
	for ( core::Size i=1, imax=helixcount; i<=imax; ++i ) {
		for ( core::Size iparam(1); iparam<=static_cast<core::Size>(BPC_last_parameter_to_be_sampled); ++iparam ) {
			core::conformation::parametric::RealValuedParameterCOP curparam( utility::pointer::dynamic_pointer_cast< core::conformation::parametric::RealValuedParameter const >( make_bundle_->helix_cop(i)->calculator_cop()->parameter_cop(iparam) ) );
			if ( curparam == nullptr ) continue;
			if ( curparam->can_be_sampled() && curparam->sampling_set() ) {
				total_samples *= curparam->value_samples();
			}
		}
	}
	return total_samples;
}

std::string BundleGridSampler::get_name() const {
	return mover_name();
}

std::string BundleGridSampler::mover_name() {
	return "BundleGridSampler";
}

std::string subtag_for_bundgrid( std::string const & subtag ) {
	return "stfbundg_" + subtag;
}

void BundleGridSampler::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	BundleParametrizationCalculatorOP default_calculator( utility::pointer::make_shared< BundleParametrizationCalculator >() );

	XMLSchemaRestriction lohigh;
	lohigh.name( "BundleGridSampler_lohigh" );
	lohigh.base_type( xs_string );
	lohigh.add_restriction( xsr_enumeration, "low" );
	lohigh.add_restriction( xsr_enumeration, "high" );
	xsd.add_top_level_element( lohigh );

	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_score_function( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default( "max_samples", xsct_non_negative_integer, "Maximum number of gridpoints to be sampled.  This is provided as a user sanity check.  Set max_samples to the number of samples you think you have requested.  The BundleGridSampler will throw an error if the actual nuber of samples is greater than this.  (For example, if I accidentally tell the BundleGridSampler to sample a 10x10x10x10 parameter grid, thinking that I will get 100 samples when I will actually get 10,000, I will quickly discover my error if I have set max_samples to 100.)", "10000" )
		+ XMLSchemaAttribute("selection_type", "BundleGridSampler_lohigh", "Score criterion for selection: \"high\" or \"low\"." )
		+ XMLSchemaAttribute("pre_scoring_mover", xs_string, "A mover to apply after backbone torsions are set but before final scoring and evaluation (like a min mover or something similar)." )
		+ XMLSchemaAttribute("pre_scoring_filter", xs_string, "A filter to apply before scoring, which could help avoid wasteful scoring of bad conformations (like a bump check filter)." )
		+ XMLSchemaAttribute("pre_selection_mover", xs_string, "A mover to apply before final solution selection (like a min mover or something similar)." )
		+ XMLSchemaAttribute("pre_selection_filter", xs_string, "A filter to apply before final solution selection, which could help avoid wasteful scoring of bad conformations (like a bump check filter)." )
		+ XMLSchemaAttribute::attribute_w_default("dump_pdbs", xsct_rosetta_bool, "Write PDBs for all conformations sampled.  If false (the default), then this mover carries out no direct PDB writing..", "false" )
		+ XMLSchemaAttribute("pdb_prefix", xs_string, "A prefix to apply to all output PDBs, if the dump_pdbs option is set to \"true\"." )
		+ XMLSchemaAttribute::attribute_w_default("use_degrees", xsct_rosetta_bool, "Interpret input values as degrees, not radians.", "false" );

	add_attributes_for_make_bundle_symmetry( attlist );
	add_attributes_for_make_bundle_other_defaults( attlist );
	add_attributes_for_make_bundle_minorhelix_defaults( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default("reset", xsct_rosetta_bool, "If reset is set to \"true\" (the default), then input geometry is discarded and the BundleGridSampler builds a pose from scratch.  If \"false\", then parametric geometry is appended to the input geometry.", "true" )
		+ XMLSchemaAttribute::attribute_w_default("nstruct_mode", xsct_rosetta_bool, "If \"true\", sample a different set of mainchain torsions for each RosettaScripts job (with each successive job sampling the next gridpoint in the grid of parameter values to be sampled).  If \"false\" (the default), then each job consists of the whole mainchain sampling effort.", "false" )
		+ XMLSchemaAttribute::attribute_w_default("nstruct_repeats", xsct_non_negative_integer, "In nstruct_mode, this is the number of times each parameter gridpoint will be sampled.  This defaults to 1 (i.e. each successive RosettaScripts job goes on to the next gridpoint), but can be set higher (i.e. successive RosettaScripts jobs repeat gridpoints).", "1" );

	for ( core::Size i(1); static_cast< BPC_Parameters >(i) < BPC_end_of_list; ++i ) { //Process all defaults that can be set:
		core::conformation::parametric::ParameterCOP curparam( default_calculator->parameter_cop(i) );
		if ( curparam->can_be_set() || curparam->can_be_sampled() ) {
			curparam->provide_xsd_information( attlist, curparam->can_be_set(), false, curparam->can_be_sampled(), false );
		}
	}

	AttributeList subtag_attributes;

	add_attributes_for_helix_params( subtag_attributes );
	add_attributes_for_minor_helix_params( subtag_attributes );
	add_attributes_for_other_helix_params( subtag_attributes );

	for ( core::Size i(1); static_cast< BPC_Parameters >(i) < BPC_end_of_list; ++i ) { //Process all defaults that can be set:
		core::conformation::parametric::ParameterCOP curparam( default_calculator->parameter_cop(i) );
		if ( !curparam->global_for_parameters_set() && ( curparam->can_be_set() || curparam->can_be_copied() || curparam->can_be_sampled() ) ) {
			curparam->provide_xsd_information( subtag_attributes, curparam->can_be_set(), curparam->can_be_copied(), curparam->can_be_sampled(), false );
		}
	}

	utility::tag::XMLSchemaSimpleSubelementList ssl;
	ssl.add_simple_subelement( "Helix", subtag_attributes, "Tags describing individual helices in the bundle"/*, 0 minoccurs*/ )
		.complex_type_naming_func( & subtag_for_bundgrid );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "The BundleGridSampler is a mover that generates helical bundles using the Crick parameterization.  It can sample regular N-dimensional grids of parameter values, with efficient parallelization.", attlist, ssl );
}

std::string BundleGridSamplerCreator::keyname() const {
	return BundleGridSampler::mover_name();
}

protocols::moves::MoverOP
BundleGridSamplerCreator::create_mover() const {
	return protocols::moves::MoverOP( new BundleGridSampler );
}

void BundleGridSamplerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BundleGridSampler::provide_xml_schema( xsd );
}


} //namespace helical_bundle
} //namespace protocols
