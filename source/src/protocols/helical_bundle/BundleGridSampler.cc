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
//#include <protocols/cyclic_peptide/PeptideStubMover.hh>
//#include <numeric/crick_equations/BundleParams.hh>
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
#include <core/scoring/rms_util.hh>
#include <core/pose/util.tmpl.hh>
#include <protocols/helical_bundle/PerturbBundleOptions.fwd.hh>
#include <protocols/helical_bundle/PerturbBundleOptions.hh>
#include <protocols/helical_bundle/util.hh>

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

// XRW TEMP std::string
// XRW TEMP BundleGridSamplerCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return BundleGridSampler::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP BundleGridSamplerCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new BundleGridSampler );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP BundleGridSampler::mover_name()
// XRW TEMP {
// XRW TEMP  return "BundleGridSampler";
// XRW TEMP }


/// @brief Creator for BundleGridSampler mover.
BundleGridSampler::BundleGridSampler():
	Mover("BundleGridSampler"),
	reset_mode_(true),
	nstruct_mode_(false),
	nstruct_mode_repeats_(1),
	select_low_(true),
	n_helices_(0),
	max_samples_( 10000 ),
	default_r0_( new PerturbBundleOptions ),
	r0_(),
	default_omega0_( new PerturbBundleOptions ),
	omega0_(),
	default_delta_omega0_( new PerturbBundleOptions ),
	delta_omega0_(),
	default_delta_omega1_( new PerturbBundleOptions ),
	delta_omega1_(),
	default_delta_t_( new PerturbBundleOptions ),
	delta_t_(),
	default_z1_offset_( new PerturbBundleOptions ),
	z1_offset_(),
	default_z0_offset_( new PerturbBundleOptions ),
	z0_offset_(),
	default_epsilon_( new PerturbBundleOptions ),
	epsilon_(),
	make_bundle_( new MakeBundle ),
	pre_selection_mover_(),
	pre_selection_mover_exists_(false),
	pre_selection_filter_(),
	pre_selection_filter_exists_(false),
	dump_pdbs_(false),
	pdb_prefix_("bgs_out"),
	sfxn_set_(false),
	sfxn_()
{
	default_epsilon_->set_default_value(1.0);
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
	default_r0_(src.default_r0_->clone()),
	r0_(),
	default_omega0_(src.default_omega0_->clone()),
	omega0_(),
	default_delta_omega0_(src.default_delta_omega0_->clone()),
	delta_omega0_(),
	default_delta_omega1_(src.default_delta_omega1_->clone()),
	delta_omega1_(),
	default_delta_t_(src.default_delta_t_->clone()),
	delta_t_(),
	default_z1_offset_(src.default_z1_offset_->clone()),
	z1_offset_(),
	default_z0_offset_(src.default_z0_offset_->clone()),
	z0_offset_(),
	default_epsilon_(src.default_epsilon_->clone()),
	epsilon_(),
	make_bundle_(  utility::pointer::dynamic_pointer_cast< MakeBundle >(src.make_bundle_->clone()) ),
	pre_selection_mover_( src.pre_selection_mover_ ), //NOTE that we're not cloning this mover, but using it straight
	pre_selection_mover_exists_(src.pre_selection_mover_exists_),
	pre_selection_filter_( src.pre_selection_filter_ ), //NOTE that we're not cloning this filter, but using it straight
	pre_selection_filter_exists_(src.pre_selection_filter_exists_),
	dump_pdbs_(src.dump_pdbs_),
	pdb_prefix_(src.pdb_prefix_),
	sfxn_set_(src.sfxn_set_),
	sfxn_(src.sfxn_) //NOTE that this is also copied without cloning
{
	r0_.clear();
	omega0_.clear();
	delta_omega0_.clear();
	delta_omega1_.clear();
	delta_t_.clear();
	z1_offset_.clear();
	z0_offset_.clear();
	epsilon_.clear();
	for ( core::Size i=1,imax=src.r0_.size(); i<=imax; ++i ) r0_.push_back( src.r0_[i]->clone() );
	for ( core::Size i=1,imax=src.omega0_.size(); i<=imax; ++i ) omega0_.push_back( src.omega0_[i]->clone() );
	for ( core::Size i=1,imax=src.delta_omega0_.size(); i<=imax; ++i ) delta_omega0_.push_back( src.delta_omega0_[i]->clone() );
	for ( core::Size i=1,imax=src.delta_omega1_.size(); i<=imax; ++i ) delta_omega1_.push_back( src.delta_omega1_[i]->clone() );
	for ( core::Size i=1,imax=src.delta_t_.size(); i<=imax; ++i ) delta_t_.push_back( src.delta_t_[i]->clone() );
	for ( core::Size i=1,imax=src.z1_offset_.size(); i<=imax; ++i ) z1_offset_.push_back( src.z1_offset_[i]->clone() );
	for ( core::Size i=1,imax=src.z0_offset_.size(); i<=imax; ++i ) z0_offset_.push_back( src.z0_offset_[i]->clone() );
	for ( core::Size i=1,imax=src.epsilon_.size(); i<=imax; ++i ) epsilon_.push_back( src.epsilon_[i]->clone() );
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

	core::Size const total_samples = calculate_total_samples();
	if ( TR.visible() ) TR << "A total of " << total_samples << " grid points will be sampled." << std::endl;
	runtime_assert_string_msg(total_samples <= max_samples(), "The total number of grid samples exceeds the maximum allowed!  (Note that you can increase the maximum allowed with the max_samples flag).");

	//Set up the BundleGridSamplerHelper object:
	if ( TR.visible() ) TR << "Creating BundleGridSamplerHelper object." << std::endl;
	BundleGridSamplerHelperOP sampler_helper( new BundleGridSamplerHelper );
	sampler_helper->reset();
	for ( core::Size ihelix=1; ihelix<=nhelices; ++ihelix ) { //Loop through all defined helices.
		if ( r0(ihelix)->use_defaults() ) {
			if ( default_r0()->is_perturbable() ) {
				sampler_helper->add_DoF( bgsh_r0, ihelix, default_r0()->samples(), default_r0()->lower_value(), default_r0()->upper_value() );
			}
		} else if ( r0(ihelix)->is_perturbable() ) {
			sampler_helper->add_DoF( bgsh_r0, ihelix, r0(ihelix)->samples(), r0(ihelix)->lower_value(), r0(ihelix)->upper_value() );
		}
	}
	for ( core::Size ihelix=1; ihelix<=nhelices; ++ihelix ) { //Loop through all defined helices.
		if ( omega0(ihelix)->use_defaults() ) {
			if ( default_omega0()->is_perturbable() ) {
				sampler_helper->add_DoF( bgsh_omega0, ihelix, default_omega0()->samples(), default_omega0()->lower_value(), default_omega0()->upper_value() );
			}
		} else if ( omega0(ihelix)->is_perturbable() ) {
			sampler_helper->add_DoF( bgsh_omega0, ihelix, omega0(ihelix)->samples(), omega0(ihelix)->lower_value(), omega0(ihelix)->upper_value() );
		}
	}
	for ( core::Size ihelix=1; ihelix<=nhelices; ++ihelix ) { //Loop through all defined helices.
		if ( delta_omega0(ihelix)->use_defaults() ) {
			if ( default_delta_omega0()->is_perturbable() ) {
				sampler_helper->add_DoF( bgsh_delta_omega0, ihelix, default_delta_omega0()->samples(), default_delta_omega0()->lower_value(), default_delta_omega0()->upper_value() );
			}
		} else if ( delta_omega0(ihelix)->is_perturbable() ) {
			sampler_helper->add_DoF( bgsh_delta_omega0, ihelix, delta_omega0(ihelix)->samples(), delta_omega0(ihelix)->lower_value(), delta_omega0(ihelix)->upper_value() );
		}
	}
	for ( core::Size ihelix=1; ihelix<=nhelices; ++ihelix ) { //Loop through all defined helices.
		if ( delta_omega1(ihelix)->use_defaults() ) {
			if ( default_delta_omega1()->is_perturbable() ) {
				sampler_helper->add_DoF( bgsh_delta_omega1, ihelix, default_delta_omega1()->samples(), default_delta_omega1()->lower_value(), default_delta_omega1()->upper_value() );
			}
		} else if ( delta_omega1(ihelix)->is_perturbable() ) {
			sampler_helper->add_DoF( bgsh_delta_omega1, ihelix, delta_omega1(ihelix)->samples(), delta_omega1(ihelix)->lower_value(), delta_omega1(ihelix)->upper_value() );
		}
	}
	for ( core::Size ihelix=1; ihelix<=nhelices; ++ihelix ) { //Loop through all defined helices.
		if ( delta_t(ihelix)->use_defaults() ) {
			if ( default_delta_t()->is_perturbable() ) {
				sampler_helper->add_DoF( bgsh_delta_t, ihelix, default_delta_t()->samples(), default_delta_t()->lower_value(), default_delta_t()->upper_value() );
			}
		} else if ( delta_t(ihelix)->is_perturbable() ) {
			sampler_helper->add_DoF( bgsh_delta_t, ihelix, delta_t(ihelix)->samples(), delta_t(ihelix)->lower_value(), delta_t(ihelix)->upper_value() );
		}
	}
	for ( core::Size ihelix=1; ihelix<=nhelices; ++ihelix ) { //Loop through all defined helices.
		if ( z1_offset(ihelix)->use_defaults() ) {
			if ( default_z1_offset()->is_perturbable() ) {
				sampler_helper->add_DoF( bgsh_z1_offset, ihelix, default_z1_offset()->samples(), default_z1_offset()->lower_value(), default_z1_offset()->upper_value() );
			}
		} else if ( z1_offset(ihelix)->is_perturbable() ) {
			sampler_helper->add_DoF( bgsh_z1_offset, ihelix, z1_offset(ihelix)->samples(), z1_offset(ihelix)->lower_value(), z1_offset(ihelix)->upper_value() );
		}
	}
	for ( core::Size ihelix=1; ihelix<=nhelices; ++ihelix ) { //Loop through all defined helices.
		if ( z0_offset(ihelix)->use_defaults() ) {
			if ( default_z0_offset()->is_perturbable() ) {
				sampler_helper->add_DoF( bgsh_z0_offset, ihelix, default_z0_offset()->samples(), default_z0_offset()->lower_value(), default_z0_offset()->upper_value() );
			}
		} else if ( z0_offset(ihelix)->is_perturbable() ) {
			sampler_helper->add_DoF( bgsh_z0_offset, ihelix, z0_offset(ihelix)->samples(), z0_offset(ihelix)->lower_value(), z0_offset(ihelix)->upper_value() );
		}
	}
	for ( core::Size ihelix=1; ihelix<=nhelices; ++ihelix ) { //Loop through all defined helices.
		if ( epsilon(ihelix)->use_defaults() ) {
			if ( default_epsilon()->is_perturbable() ) {
				sampler_helper->add_DoF( bgsh_epsilon, ihelix, default_epsilon()->samples(), default_epsilon()->lower_value(), default_epsilon()->upper_value() );
			}
		} else if ( epsilon(ihelix)->is_perturbable() ) {
			sampler_helper->add_DoF( bgsh_epsilon, ihelix, epsilon(ihelix)->samples(), epsilon(ihelix)->lower_value(), epsilon(ihelix)->upper_value() );
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
			runtime_assert_string_msg(curhelix, "Error in getting owning pointer to current helix in BundleGridSampler::apply() function.");

			if ( sampler_helper->DoF_type(j) == bgsh_r0 ) curhelix->set_r0( sampler_helper->DoF_sample_value(j) );
			else if ( sampler_helper->DoF_type(j) == bgsh_omega0 ) curhelix->set_omega0( sampler_helper->DoF_sample_value(j) );
			else if ( sampler_helper->DoF_type(j) == bgsh_delta_omega0 ) curhelix->set_delta_omega0( sampler_helper->DoF_sample_value(j) );
			else if ( sampler_helper->DoF_type(j) == bgsh_delta_omega1 ) curhelix->set_delta_omega1_all( sampler_helper->DoF_sample_value(j) );
			else if ( sampler_helper->DoF_type(j) == bgsh_delta_t ) curhelix->set_delta_t( sampler_helper->DoF_sample_value(j) );
			else if ( sampler_helper->DoF_type(j) == bgsh_z1_offset ) curhelix->set_z1_offset( sampler_helper->DoF_sample_value(j) );
			else if ( sampler_helper->DoF_type(j) == bgsh_z0_offset ) curhelix->set_z0_offset( sampler_helper->DoF_sample_value(j) );
			else if ( sampler_helper->DoF_type(j) == bgsh_epsilon ) curhelix->set_epsilon( sampler_helper->DoF_sample_value(j) );
		}

		//Set the parameters that are copying other parameters:
		bool loopthrough_failed(false); //If we fail to copy another parameter (e.g. pitch angle), we need to know it.
		for ( core::Size j=1, jmax=n_helices(); j<=jmax; ++j ) {
			MakeBundleHelixOP curhelix( makebundle_copy->helix( j ) );
			runtime_assert_string_msg(curhelix, "Error in getting owning pointer to current helix in BundleGridSampler::apply() function.");

			if ( r0(j)->is_copy() ) {
				core::Size const otherhelix( r0(j)->other_helix() );
				runtime_assert_string_msg(otherhelix > 0 && otherhelix < j, "Error in getting index of helix to copy in BundleGridSampler::apply() function.");
				curhelix->set_r0( makebundle_copy->helix_cop(otherhelix)->r0() );
			}
			if ( omega0(j)->is_copy() ) {
				core::Size const otherhelix( omega0(j)->other_helix() );
				runtime_assert_string_msg(otherhelix > 0 && otherhelix < j, "Error in getting index of helix to copy in BundleGridSampler::apply() function.");
				//Special case: if we're copying the pitch angle instead of the omega0 value, we need to do some math.
				if ( omega0(j)->omega0_copies_pitch_instead() ) {
					core::Real const other_r0( makebundle_copy->helix_cop(otherhelix)->r0() );
					core::Real const other_omega0( makebundle_copy->helix_cop(otherhelix)->omega0() );
					core::Real const other_z1( makebundle_copy->helix_cop(otherhelix)->z1() );
					core::Real const other_sinalpha( other_r0*other_omega0/other_z1 );
					if ( other_sinalpha > 1 || other_sinalpha < -1 ) {
						if ( TR.visible() ) TR << "Failed to copy pitch angle.  Current parameters do not generate a sensible pitch angle for helix " << otherhelix << "." << std::endl;
						loopthrough_failed=true;
						break; //Stop looping through the helices.
					}
					//If we've got a good pitch angle, then continue:
					core::Real const other_alpha( asin(other_sinalpha) );

					core::Real const this_r0( curhelix->r0() ); //Already set above, if sampled or if copied.
					core::Real const this_z1( curhelix->z1() ); //Cannot be sampled or copied.
					/********************
					We know: tan(alpha)=2*PI*R0/P, where alpha is the pitch angle, P is the pitch (rise per turn about major helix), and R0 is the major radius.
					sin(alpha)=R0*omega0/z1
					We want: P' = P
					2*PI*RO'/tan(alpha') = 2*PI*R0/tan(alpha)
					tan(alpha) = R0/R0'*tan(alpha')
					alpha = atan(R0/R0'*tan(alpha')
					R0*omega0/z1 = sin(atan(R0/R0'*tan(asin(R0'*omega0'/z1'))))
					omega0 = z1/R0*sin(atan(R0/R0'*tan(asin(R0'*omega0'/z1'))))
					********************/
					curhelix->set_omega0( this_z1/this_r0 * sin(atan(this_r0/other_r0*tan(other_alpha))) );
				} else {
					curhelix->set_omega0( makebundle_copy->helix_cop(otherhelix)->omega0() );
				}
			}
			if ( delta_omega0(j)->is_copy() ) {
				core::Size const otherhelix( delta_omega0(j)->other_helix() );
				runtime_assert_string_msg(otherhelix > 0 && otherhelix < j, "Error in getting index of helix to copy in BundleGridSampler::apply() function.");
				curhelix->set_delta_omega0( makebundle_copy->helix_cop(otherhelix)->delta_omega0() );
			}
			if ( delta_omega1(j)->is_copy() ) {
				core::Size const otherhelix( delta_omega1(j)->other_helix() );
				runtime_assert_string_msg(otherhelix > 0 && otherhelix < j, "Error in getting index of helix to copy in BundleGridSampler::apply() function.");
				curhelix->set_delta_omega1_all( makebundle_copy->helix_cop(otherhelix)->delta_omega1_all() );
			}
			if ( delta_t(j)->is_copy() ) {
				core::Size const otherhelix( delta_t(j)->other_helix() );
				runtime_assert_string_msg(otherhelix > 0 && otherhelix < j, "Error in getting index of helix to copy in BundleGridSampler::apply() function.");
				curhelix->set_delta_t( makebundle_copy->helix_cop(otherhelix)->delta_t() );
			}
			if ( z1_offset(j)->is_copy() ) {
				core::Size const otherhelix( z1_offset(j)->other_helix() );
				runtime_assert_string_msg(otherhelix > 0 && otherhelix < j, "Error in getting index of helix to copy in BundleGridSampler::apply() function.");
				curhelix->set_z1_offset( makebundle_copy->helix_cop(otherhelix)->z1_offset() );
			}
			if ( z0_offset(j)->is_copy() ) {
				core::Size const otherhelix( z0_offset(j)->other_helix() );
				runtime_assert_string_msg(otherhelix > 0 && otherhelix < j, "Error in getting index of helix to copy in BundleGridSampler::apply() function.");
				curhelix->set_z0_offset( makebundle_copy->helix_cop(otherhelix)->z0_offset() );
			}
			if ( epsilon(j)->is_copy() ) {
				core::Size const otherhelix( epsilon(j)->other_helix() );
				runtime_assert_string_msg(otherhelix > 0 && otherhelix < j, "Error in getting index of helix to copy in BundleGridSampler::apply() function.");
				curhelix->set_epsilon( makebundle_copy->helix_cop(otherhelix)->epsilon() );
			}

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
			char outchar [1024];
			TR << "Grid point " << i << ": attempting to build a helical bundle with the following parameters:" << std::endl;
			sprintf(outchar, "Helix\tr0\tomega0\tdelta_omega0\tdelta_omega1\tdelta_t\tz1_offset\tz0_offset\tepsilon");
			TR << outchar << std::endl;
			for ( core::Size j=1, jmax=n_helices(); j<=jmax; ++j ) {
				MakeBundleHelixCOP curhelix( makebundle_copy->helix_cop( j ) );
				runtime_assert_string_msg(curhelix, "Error in getting owning pointer to current helix in BundleGridSampler::apply() function.");
				sprintf(outchar, "%lu\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f", static_cast<unsigned long>(j), curhelix->r0(), curhelix->omega0(), curhelix->delta_omega0(), curhelix->delta_omega1_all(), curhelix->delta_t(), curhelix->z1_offset(), curhelix->z0_offset(), curhelix->epsilon());
				TR << outchar << std::endl;
			}
			TR << std::endl;
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

				// Generate the summary of parameters, in case this is ultimately the lowest-energy pose:
				final_report.str("");
				final_report << "Parameters yielding the ";
				if ( selection_low() ) final_report << "lowest"; else final_report << "highest";
				final_report << "-energy bundle:" << std::endl;
				sprintf(outstring, "Helix\tr0\tomega0\tdelta_omega0\tdelta_omega1\tdelta_t\tz1_offset\tz0_offset\tepsilon");
				final_report << outstring << std::endl;
				for ( core::Size j=1, jmax=n_helices(); j<=jmax; ++j ) {
					MakeBundleHelixCOP curhelix( makebundle_copy->helix_cop( j ) );
					runtime_assert_string_msg(curhelix, "Error in getting owning pointer to current helix in BundleGridSampler::apply() function.");
					sprintf(outstring, "%lu\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f", static_cast<unsigned long>(j), curhelix->r0(), curhelix->omega0(), curhelix->delta_omega0(), curhelix->delta_omega1_all(), curhelix->delta_t(), curhelix->z1_offset(), curhelix->z0_offset(), curhelix->epsilon() );
					final_report << outstring << std::endl;
				}
				final_report << std::endl;
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


/// @brief Returns the name of this mover ("BundleGridSampler").
// XRW TEMP std::string BundleGridSampler::get_name() const{
// XRW TEMP  return "BundleGridSampler";
// XRW TEMP }

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

	//TODO: ADD SUPPORT FOR Z-OFFSET.

	if ( tag->getName() != "BundleGridSampler" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This should be impossible -- the tag name does not match the mover name.");
	}

	if ( TR.visible() ) TR << "Parsing options for BundleGridSampler (\"" << tag->getOption<std::string>("name" ,"") << "\") mover." << std::endl;

	//Global options for this mover:
	runtime_assert_string_msg( tag->hasOption("scorefxn" ), "In BundleGridSampler::parse_my_tag(): A \"scorefxn\" option must be specified!");
	set_sfxn(protocols::rosetta_scripts::parse_score_function( tag, "scorefxn", data_map )->clone()); // The scorefunction.
	if ( tag->hasOption("max_samples") ) {
		core::Size const val( tag->getOption<core::Size>("max_samples", 10000) );
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
	make_bundle_->set_use_degrees( tag->getOption<bool>("use_degrees", false) );
	if ( TR.visible() ) TR << "Setting mover to interpret angles in RosettaScripts as though they were specified in " << (make_bundle_->use_degrees() ? "degrees." : "radians.") << "  Note that internally and in output, the mover uses radians." << std::endl;

	//Set symmetry options for the MakeBundle mover (symmetry, symmetry_copies):
	make_bundle_->set_symmetry_options_from_tag( tag );

	//Set defaults for whether the mover can set bond lengths, bond angles, and dihedrals:
	make_bundle_->set_dofs_from_tag( tag );

	//Set additional defaults for the MakeBundle mover (residue_name, invert, and helix_length):
	make_bundle_->set_other_defaults_from_tag( tag );

	//Set minor helix default parameters for the MakeBundle mover (crick_params_file, omega1, and z1):
	make_bundle_->set_minorhelix_defaults_from_tag( tag );

	//Set reset mode for the MakeBundle mover:
	bool const resetmode( tag->getOption<bool>("reset", true) );
	if ( TR.visible() ) TR << "Setting reset mode to " << (resetmode ? "true." : "false.") << std::endl;
	set_reset_mode(resetmode);
	make_bundle_->set_reset_pose( reset_mode() );

	//Set nstruct mode and options:
	bool nstructmode = tag->getOption<bool>("nstruct_mode", false);
	if ( TR.visible() ) TR << "Setting nstruct mode to " << (nstructmode ? "true." : "false.") << "  This means that " << (nstructmode ? "each job will sample a different set of Crick parameters." : "every job will sample all sets of Crick parameters.") << std::endl;
	set_nstruct_mode(nstructmode);
	core::Size nstructrepeats( tag->getOption<core::Size>( "nstruct_repeats", 1 ) );
	if ( nstructrepeats<1 ) nstructrepeats=1;
	if ( TR.visible() ) TR << "Setting nstruct repeats to " << nstructrepeats << "." << std::endl;
	set_nstruct_repeats(nstructrepeats);

	// Default options applied to all helices, unless overrides are provided:
	if ( tag->hasOption("r0") ) {
		runtime_assert_string_msg(!tag->hasOption("r0_min") && !tag->hasOption("r0_max"),
			"When parsing options for the BundleGridSampler mover, found r0 defined alongside r0_min or r0_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		core::Real const val( tag->getOption<core::Real>("r0", 0.0) );
		runtime_assert_string_msg( val>=(-0.000000001), "When parsing options for the BundleGridSampler mover, found an r0 value less than 0." );
		if ( TR.visible() ) TR << "Setting default r0 value to " << val << "." << std::endl;
		default_r0()->set_default_value(val);
		default_r0()->set_perturbable(false);
		make_bundle_->set_default_r0(val);
	} else if ( tag->hasOption("r0_min") ) {
		runtime_assert_string_msg(!tag->hasOption("r0"),
			"When parsing options for the BundleGridSampler mover, found r0 defined alongside r0_min or r0_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		runtime_assert_string_msg(tag->hasOption("r0_max"),
			"When parsing options for the BundleGridSampler mover, found r0_min but no r0_max.  This does not make sense.");
		core::Real const val1( tag->getOption<core::Real>("r0_min", 0.0) );
		core::Real const val2( tag->getOption<core::Real>("r0_max", 0.0) );
		runtime_assert_string_msg( val1>=(-0.000000001), "When parsing options for the BundleGridSampler mover, found an r0_min value less than 0." );
		runtime_assert_string_msg( val2>=(-0.000000001), "When parsing options for the BundleGridSampler mover, found an r0_max value less than 0." );
		runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an r0_max value less than or equal to the r0_min value.");
		if ( TR.visible() ) TR << "Setting default r0_min and r0_max values to " << val1 << " and " << val2 << ", respectively." << std::endl;
		default_r0()->set_lower_value(val1);
		default_r0()->set_upper_value(val2);
		default_r0()->set_perturbable(true);
		make_bundle_->set_default_r0(val1);
		runtime_assert_string_msg( tag->hasOption("r0_samples"),
			"When parsing options for the BundleGridSampler mover, found r0_min and r0_max options, but no r0_samples option.  The number of r0 samples must be specified." );
		core::Size const val3( tag->getOption<core::Size>( "r0_samples", 0 ) );
		runtime_assert_string_msg( val3 > 1,
			"The number of r0 samples must be greater than 1." );
		if ( TR.visible() ) TR << "Setting default r0 samples to " << val3 << "." << std::endl;
		default_r0()->set_samples(val3);
	} else if ( tag->hasOption("r0_max") ) {
		runtime_assert_string_msg(tag->hasOption("r0_min"),
			"When parsing options for the BundleGridSampler mover, found r0_max but no r0_min.  This does not make sense.");
	}

	if ( tag->hasOption("omega0") ) {
		runtime_assert_string_msg(!tag->hasOption("omega0_min") && !tag->hasOption("omega0_max"),
			"When parsing options for the BundleGridSampler mover, found omega0 defined alongside omega0_min or omega0_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		core::Real const val( tag->getOption<core::Real>("omega0", 0.0) );
		if ( TR.visible() ) TR << "Setting default omega0 value to " << val << "." << std::endl;
		default_omega0()->set_default_value( make_bundle_->convert_angle(val) );
		default_omega0()->set_perturbable(false);
		make_bundle_->set_default_omega0( val );
	} else if ( tag->hasOption("omega0_min") ) {
		runtime_assert_string_msg(!tag->hasOption("omega0"),
			"When parsing options for the BundleGridSampler mover, found omega0 defined alongside omega0_min or omega0_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		runtime_assert_string_msg(tag->hasOption("omega0_max"),
			"When parsing options for the BundleGridSampler mover, found omega0_min but no omega0_max.  This does not make sense.");
		core::Real const val1( tag->getOption<core::Real>("omega0_min", 0.0) );
		core::Real const val2( tag->getOption<core::Real>("omega0_max", 0.0) );
		runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an omega0_max value less than or equal to the omega0_min value.");
		if ( TR.visible() ) TR << "Setting default omega0_min and omega0_max values to " << val1 << " and " << val2 << ", respectively." << std::endl;
		default_omega0()->set_lower_value( make_bundle_->convert_angle(val1) );
		default_omega0()->set_upper_value( make_bundle_->convert_angle(val2) );
		default_omega0()->set_perturbable(true);
		make_bundle_->set_default_omega0( val1 );
		runtime_assert_string_msg( tag->hasOption("omega0_samples"),
			"When parsing options for the BundleGridSampler mover, found omega0_min and omega0_max options, but no omega0_samples option.  The number of omega0 samples must be specified." );
		core::Size const val3( tag->getOption<core::Size>( "omega0_samples", 0 ) );
		runtime_assert_string_msg( val3 > 1,
			"The number of omega0 samples must be greater than 1." );
		if ( TR.visible() ) TR << "Setting default omega0 samples to " << val3 << "." << std::endl;
		default_omega0()->set_samples(val3);
	} else if ( tag->hasOption("omega0_max") ) {
		runtime_assert_string_msg(tag->hasOption("omega0_min"),
			"When parsing options for the BundleGridSampler mover, found omega0_max but no omega0_min.  This does not make sense.");
	}

	if ( tag->hasOption("delta_omega0") ) {
		runtime_assert_string_msg(!tag->hasOption("delta_omega0_min") && !tag->hasOption("delta_omega0_max"),
			"When parsing options for the BundleGridSampler mover, found delta_omega0 defined alongside delta_omega0_min or delta_omega0_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		core::Real const val( tag->getOption<core::Real>("delta_omega0", 0.0) );
		if ( TR.visible() ) TR << "Setting default delta_omega0 value to " << val << "." << std::endl;
		default_delta_omega0()->set_default_value( make_bundle_->convert_angle(val) );
		default_delta_omega0()->set_perturbable(false);
		make_bundle_->set_default_delta_omega0( val );
	} else if ( tag->hasOption("delta_omega0_min") ) {
		runtime_assert_string_msg(!tag->hasOption("delta_omega0"),
			"When parsing options for the BundleGridSampler mover, found delta_omega0 defined alongside delta_omega0_min or delta_omega0_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		runtime_assert_string_msg(tag->hasOption("delta_omega0_max"),
			"When parsing options for the BundleGridSampler mover, found delta_omega0_min but no delta_omega0_max.  This does not make sense.");
		core::Real const val1( tag->getOption<core::Real>("delta_omega0_min", 0.0) );
		core::Real const val2( tag->getOption<core::Real>("delta_omega0_max", 0.0) );
		runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an delta_omega0_max value less than or equal to the delta_omega0_min value.");
		if ( TR.visible() ) TR << "Setting default delta_omega0_min and delta_omega0_max values to " << val1 << " and " << val2 << ", respectively." << std::endl;
		default_delta_omega0()->set_lower_value( make_bundle_->convert_angle(val1) );
		default_delta_omega0()->set_upper_value( make_bundle_->convert_angle(val2) );
		default_delta_omega0()->set_perturbable(true);
		make_bundle_->set_default_delta_omega0( val1 );
		runtime_assert_string_msg( tag->hasOption("delta_omega0_samples"),
			"When parsing options for the BundleGridSampler mover, found delta_omega0_min and delta_omega0_max options, but no delta_omega0_samples option.  The number of delta_omega0 samples must be specified." );
		core::Size const val3( tag->getOption<core::Size>( "delta_omega0_samples", 0 ) );
		runtime_assert_string_msg( val3 > 1,
			"The number of delta_omega0 samples must be greater than 1." );
		if ( TR.visible() ) TR << "Setting default delta_omega0 samples to " << val3 << "." << std::endl;
		default_delta_omega0()->set_samples(val3);
	} else if ( tag->hasOption("delta_omega0_max") ) {
		runtime_assert_string_msg(tag->hasOption("delta_omega0_min"),
			"When parsing options for the BundleGridSampler mover, found delta_omega0_max but no delta_omega0_min.  This does not make sense.");
	}

	if ( tag->hasOption("delta_omega1") ) {
		runtime_assert_string_msg(!tag->hasOption("delta_omega1_min") && !tag->hasOption("delta_omega1_max"),
			"When parsing options for the BundleGridSampler mover, found delta_omega1 defined alongside delta_omega1_min or delta_omega1_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		core::Real const val( tag->getOption<core::Real>("delta_omega1", 0.0) );
		if ( TR.visible() ) TR << "Setting default delta_omega1 value to " << val << "." << std::endl;
		default_delta_omega1()->set_default_value( make_bundle_->convert_angle(val) );
		default_delta_omega1()->set_perturbable(false);
		make_bundle_->set_default_delta_omega1_all( val );
	} else if ( tag->hasOption("delta_omega1_min") ) {
		runtime_assert_string_msg(!tag->hasOption("delta_omega1"),
			"When parsing options for the BundleGridSampler mover, found delta_omega1 defined alongside delta_omega1_min or delta_omega1_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		runtime_assert_string_msg(tag->hasOption("delta_omega1_max"),
			"When parsing options for the BundleGridSampler mover, found delta_omega1_min but no delta_omega1_max.  This does not make sense.");
		core::Real const val1( tag->getOption<core::Real>("delta_omega1_min", 0.0) );
		core::Real const val2( tag->getOption<core::Real>("delta_omega1_max", 0.0) );
		runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an delta_omega1_max value less than or equal to the delta_omega1_min value.");
		if ( TR.visible() ) TR << "Setting default delta_omega1_min and delta_omega1_max values to " << val1 << " and " << val2 << ", respectively." << std::endl;
		default_delta_omega1()->set_lower_value( make_bundle_->convert_angle(val1) );
		default_delta_omega1()->set_upper_value( make_bundle_->convert_angle(val2) );
		default_delta_omega1()->set_perturbable(true);
		make_bundle_->set_default_delta_omega1_all( val1 );
		runtime_assert_string_msg( tag->hasOption("delta_omega1_samples"),
			"When parsing options for the BundleGridSampler mover, found delta_omega1_min and delta_omega1_max options, but no delta_omega1_samples option.  The number of delta_omega1 samples must be specified." );
		core::Size const val3( tag->getOption<core::Size>( "delta_omega1_samples", 0 ) );
		runtime_assert_string_msg( val3 > 1,
			"The number of delta_omega1 samples must be greater than 1." );
		if ( TR.visible() ) TR << "Setting default delta_omega1 samples to " << val3 << "." << std::endl;
		default_delta_omega1()->set_samples(val3);
	} else if ( tag->hasOption("delta_omega1_max") ) {
		runtime_assert_string_msg(tag->hasOption("delta_omega1_min"),
			"When parsing options for the BundleGridSampler mover, found delta_omega1_max but no delta_omega1_min.  This does not make sense.");
	}

	if ( tag->hasOption("delta_t") ) {
		runtime_assert_string_msg(!tag->hasOption("delta_t_min") && !tag->hasOption("delta_t_max"),
			"When parsing options for the BundleGridSampler mover, found delta_t defined alongside delta_t_min or delta_t_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		core::Real const val( tag->getOption<core::Real>("delta_t", 0.0) );
		if ( TR.visible() ) TR << "Setting default delta_t value to " << val << "." << std::endl;
		default_delta_t()->set_default_value(val);
		default_delta_t()->set_perturbable(false);
		make_bundle_->set_default_delta_t(val);
	} else if ( tag->hasOption("delta_t_min") ) {
		runtime_assert_string_msg(!tag->hasOption("delta_t"),
			"When parsing options for the BundleGridSampler mover, found delta_t defined alongside delta_t_min or delta_t_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		runtime_assert_string_msg(tag->hasOption("delta_t_max"),
			"When parsing options for the BundleGridSampler mover, found delta_t_min but no delta_t_max.  This does not make sense.");
		core::Real const val1( tag->getOption<core::Real>("delta_t_min", 0.0) );
		core::Real const val2( tag->getOption<core::Real>("delta_t_max", 0.0) );
		runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an delta_t_max value less than or equal to the delta_t_min value.");
		if ( TR.visible() ) TR << "Setting default delta_t_min and delta_t_max values to " << val1 << " and " << val2 << ", respectively." << std::endl;
		default_delta_t()->set_lower_value(val1);
		default_delta_t()->set_upper_value(val2);
		default_delta_t()->set_perturbable(true);
		make_bundle_->set_default_delta_t(val1);
		runtime_assert_string_msg( tag->hasOption("delta_t_samples"),
			"When parsing options for the BundleGridSampler mover, found delta_t_min and delta_t_max options, but no delta_t_samples option.  The number of delta_t samples must be specified." );
		core::Size const val3( tag->getOption<core::Size>( "delta_t_samples", 0 ) );
		runtime_assert_string_msg( val3 > 1,
			"The number of delta_t samples must be greater than 1." );
		if ( TR.visible() ) TR << "Setting default delta_t samples to " << val3 << "." << std::endl;
		default_delta_t()->set_samples(val3);
	} else if ( tag->hasOption("delta_t_max") ) {
		runtime_assert_string_msg(tag->hasOption("delta_t_min"),
			"When parsing options for the BundleGridSampler mover, found delta_t_max but no delta_t_min.  This does not make sense.");
	}

	if ( tag->hasOption("z1_offset") ) {
		runtime_assert_string_msg(!tag->hasOption("z1_offset_min") && !tag->hasOption("z1_offset_max"),
			"When parsing options for the BundleGridSampler mover, found z1_offset defined alongside z1_offset_min or z1_offset_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		core::Real const val( tag->getOption<core::Real>("z1_offset", 0.0) );
		if ( TR.visible() ) TR << "Setting default z1_offset value to " << val << "." << std::endl;
		default_z1_offset()->set_default_value(val);
		default_z1_offset()->set_perturbable(false);
		make_bundle_->set_default_z1_offset(val);
	} else if ( tag->hasOption("z1_offset_min") ) {
		runtime_assert_string_msg(!tag->hasOption("z1_offset"),
			"When parsing options for the BundleGridSampler mover, found z1_offset defined alongside z1_offset_min or z1_offset_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		runtime_assert_string_msg(tag->hasOption("z1_offset_max"),
			"When parsing options for the BundleGridSampler mover, found z1_offset_min but no z1_offset_max.  This does not make sense.");
		core::Real const val1( tag->getOption<core::Real>("z1_offset_min", 0.0) );
		core::Real const val2( tag->getOption<core::Real>("z1_offset_max", 0.0) );
		runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an z1_offset_max value less than or equal to the z1_offset_min value.");
		if ( TR.visible() ) TR << "Setting default z1_offset_min and z1_offset_max values to " << val1 << " and " << val2 << ", respectively." << std::endl;
		default_z1_offset()->set_lower_value(val1);
		default_z1_offset()->set_upper_value(val2);
		default_z1_offset()->set_perturbable(true);
		make_bundle_->set_default_z1_offset(val1);
		runtime_assert_string_msg( tag->hasOption("z1_offset_samples"),
			"When parsing options for the BundleGridSampler mover, found z1_offset_min and z1_offset_max options, but no z1_offset_samples option.  The number of z1_offset samples must be specified." );
		core::Size const val3( tag->getOption<core::Size>( "z1_offset_samples", 0 ) );
		runtime_assert_string_msg( val3 > 1,
			"The number of z1_offset samples must be greater than 1." );
		if ( TR.visible() ) TR << "Setting default z1_offset samples to " << val3 << "." << std::endl;
		default_z1_offset()->set_samples(val3);
	} else if ( tag->hasOption("z1_offset_max") ) {
		runtime_assert_string_msg(tag->hasOption("z1_offset_min"),
			"When parsing options for the BundleGridSampler mover, found z1_offset_max but no z1_offset_min.  This does not make sense.");
	}

	if ( tag->hasOption("z0_offset") ) {
		runtime_assert_string_msg(!tag->hasOption("z0_offset_min") && !tag->hasOption("z0_offset_max"),
			"When parsing options for the BundleGridSampler mover, found z0_offset defined alongside z0_offset_min or z0_offset_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		core::Real const val( tag->getOption<core::Real>("z0_offset", 0.0) );
		if ( TR.visible() ) TR << "Setting default z0_offset value to " << val << "." << std::endl;
		default_z0_offset()->set_default_value(val);
		default_z0_offset()->set_perturbable(false);
		make_bundle_->set_default_z0_offset(val);
	} else if ( tag->hasOption("z0_offset_min") ) {
		runtime_assert_string_msg(!tag->hasOption("z0_offset"),
			"When parsing options for the BundleGridSampler mover, found z0_offset defined alongside z0_offset_min or z0_offset_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		runtime_assert_string_msg(tag->hasOption("z0_offset_max"),
			"When parsing options for the BundleGridSampler mover, found z0_offset_min but no z0_offset_max.  This does not make sense.");
		core::Real const val1( tag->getOption<core::Real>("z0_offset_min", 0.0) );
		core::Real const val2( tag->getOption<core::Real>("z0_offset_max", 0.0) );
		runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an z0_offset_max value less than or equal to the z0_offset_min value.");
		if ( TR.visible() ) TR << "Setting default z0_offset_min and z0_offset_max values to " << val1 << " and " << val2 << ", respectively." << std::endl;
		default_z0_offset()->set_lower_value(val1);
		default_z0_offset()->set_upper_value(val2);
		default_z0_offset()->set_perturbable(true);
		make_bundle_->set_default_z0_offset(val1);
		runtime_assert_string_msg( tag->hasOption("z0_offset_samples"),
			"When parsing options for the BundleGridSampler mover, found z0_offset_min and z0_offset_max options, but no z0_offset_samples option.  The number of z0_offset samples must be specified." );
		core::Size const val3( tag->getOption<core::Size>( "z0_offset_samples", 0 ) );
		runtime_assert_string_msg( val3 > 1,
			"The number of z0_offset samples must be greater than 1." );
		if ( TR.visible() ) TR << "Setting default z0_offset samples to " << val3 << "." << std::endl;
		default_z0_offset()->set_samples(val3);
	} else if ( tag->hasOption("z0_offset_max") ) {
		runtime_assert_string_msg(tag->hasOption("z0_offset_min"),
			"When parsing options for the BundleGridSampler mover, found z0_offset_max but no z0_offset_min.  This does not make sense.");
	}

	if ( tag->hasOption("epsilon") ) {
		runtime_assert_string_msg(!tag->hasOption("epsilon_min") && !tag->hasOption("epsilon_max"),
			"When parsing options for the BundleGridSampler mover, found epsilon defined alongside epsilon_min or epsilon_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		core::Real const val( tag->getOption<core::Real>("epsilon", 1.0) );
		if ( TR.visible() ) TR << "Setting default epsilon value to " << val << "." << std::endl;
		default_epsilon()->set_default_value(val);
		default_epsilon()->set_perturbable(false);
		make_bundle_->set_default_epsilon(val);
	} else if ( tag->hasOption("epsilon_min") ) {
		runtime_assert_string_msg(!tag->hasOption("epsilon"),
			"When parsing options for the BundleGridSampler mover, found epsilon defined alongside epsilon_min or epsilon_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
		runtime_assert_string_msg(tag->hasOption("epsilon_max"),
			"When parsing options for the BundleGridSampler mover, found epsilon_min but no epsilon_max.  This does not make sense.");
		core::Real const val1( tag->getOption<core::Real>("epsilon_min", 1.0) );
		core::Real const val2( tag->getOption<core::Real>("epsilon_max", 1.0) );
		runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an epsilon_max value less than or equal to the epsilon_min value.");
		if ( TR.visible() ) TR << "Setting default epsilon_min and epsilon_max values to " << val1 << " and " << val2 << ", respectively." << std::endl;
		default_epsilon()->set_lower_value(val1);
		default_epsilon()->set_upper_value(val2);
		default_epsilon()->set_perturbable(true);
		make_bundle_->set_default_epsilon(val1);
		runtime_assert_string_msg( tag->hasOption("epsilon_samples"),
			"When parsing options for the BundleGridSampler mover, found epsilon_min and epsilon_max options, but no epsilon_samples option.  The number of epsilon samples must be specified." );
		core::Size const val3( tag->getOption<core::Size>( "epsilon_samples", 0 ) );
		runtime_assert_string_msg( val3 > 1,
			"The number of epsilon samples must be greater than 1." );
		if ( TR.visible() ) TR << "Setting default epsilon samples to " << val3 << "." << std::endl;
		default_epsilon()->set_samples(val3);
	} else if ( tag->hasOption("epsilon_max") ) {
		runtime_assert_string_msg(tag->hasOption("epsilon_min"),
			"When parsing options for the BundleGridSampler mover, found epsilon_max but no epsilon_min.  This does not make sense.");
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

			//Set crick_params_file, set_bondlengths, set_bondangles, and set_dihedrals options for this helix, based on the tag:
			make_bundle_->set_helix_params_from_tag( helix_index, branch_tag );

			//Set omega1, z1, and delta_omega1 for this helix, based on the tag:
			make_bundle_->set_minor_helix_params_from_tag( helix_index, branch_tag );

			//Set residue_name, invert, and helix_length for this helix, based on the tag:
			make_bundle_->set_other_helix_params_from_tag( helix_index, branch_tag );

			//Parse options that could be sampled for this helix:
			if ( branch_tag->hasOption("r0") ) {
				runtime_assert_string_msg(!branch_tag->hasOption("r0_copies_helix"),
					"When parsing options for the BundleGridSampler mover, found r0 defined alongside an r0_copies_helix statement.  This does not makes sense: either r0 is fixed, or its value copies the r0 value of another helix (not both).");
				runtime_assert_string_msg(!branch_tag->hasOption("r0_min") && !branch_tag->hasOption("r0_max"),
					"When parsing options for the BundleGridSampler mover, found r0 defined alongside r0_min or r0_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
				core::Real const val( branch_tag->getOption<core::Real>("r0", 0.0) );
				runtime_assert_string_msg( val>=(-0.000000001), "When parsing options for the BundleGridSampler mover, found an r0 value less than 0." );
				if ( TR.visible() ) TR << "Setting the r0 value for helix " << helix_index << " to " << val << "." << std::endl;
				r0(helix_index)->set_default_value(val);
				r0(helix_index)->set_perturbable(false);
				r0(helix_index)->set_use_defaults(false);
				make_bundle_->helix(helix_index)->set_r0(val);
			} else if ( branch_tag->hasOption("r0_min") ) {
				runtime_assert_string_msg(!branch_tag->hasOption("r0_copies_helix"),
					"When parsing options for the BundleGridSampler mover, found r0_min defined alongside an r0_copies_helix statement.  This does not makes sense: either r0 is sampled, or its value copies the r0 value of another helix (not both).");
				runtime_assert_string_msg(!branch_tag->hasOption("r0"),
					"When parsing options for the BundleGridSampler mover, found r0 defined alongside r0_min or r0_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
				runtime_assert_string_msg(branch_tag->hasOption("r0_max"),
					"When parsing options for the BundleGridSampler mover, found r0_min but no r0_max.  This does not make sense.");
				core::Real const val1( branch_tag->getOption<core::Real>("r0_min", 0.0) );
				core::Real const val2( branch_tag->getOption<core::Real>("r0_max", 0.0) );
				runtime_assert_string_msg( val1>=(-0.000000001), "When parsing options for the BundleGridSampler mover, found an r0_min value less than 0." );
				runtime_assert_string_msg( val2>=(-0.000000001), "When parsing options for the BundleGridSampler mover, found an r0_max value less than 0." );
				runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an r0_max value less than or equal to the r0_min value.");
				if ( TR.visible() ) TR << "Setting r0_min and r0_max values for helix " << helix_index << " to " << val1 << " and " << val2 << ", respectively." << std::endl;
				r0(helix_index)->set_lower_value(val1);
				r0(helix_index)->set_upper_value(val2);
				r0(helix_index)->set_perturbable(true);
				r0(helix_index)->set_use_defaults(false);
				make_bundle_->helix(helix_index)->set_r0(val1);
				runtime_assert_string_msg( branch_tag->hasOption("r0_samples"),
					"When parsing options for the BundleGridSampler mover, found r0_min and r0_max options, but no r0_samples option.  The number of r0 samples must be specified." );
				core::Size const val3( branch_tag->getOption<core::Size>( "r0_samples", 0 ) );
				runtime_assert_string_msg( val3 > 1,
					"The number of r0 samples must be greater than 1." );
				if ( TR.visible() ) TR << "Setting r0 samples for helix " << helix_index << " to " << val3 << "." << std::endl;
				r0(helix_index)->set_samples(val3);
			} else if ( branch_tag->hasOption("r0_max") ) {
				runtime_assert_string_msg(branch_tag->hasOption("r0_min"),
					"When parsing options for the BundleGridSampler mover, found r0_max but no r0_min.  This does not make sense.");
			} else if ( branch_tag->hasOption("r0_copies_helix") ) {
				core::Size const val( branch_tag->getOption<core::Size>("r0_copies_helix", 0) );
				runtime_assert_string_msg(val > 0 && val < helix_index,
					"When parsing options for the BundleGridSampler mover, found an r0_copies_helix option with a nonsensical value.  Please specify an already-defined helix index.");
				if ( TR.visible() ) TR << "Setting r0 for helix " << helix_index << " to copy the r0 value for helix " << val << "." << std::endl;
				r0(helix_index)->set_perturbable(false);
				r0(helix_index)->set_use_defaults(false);
				r0(helix_index)->set_helix_to_copy(val);
			}

			// Note that omega0 has additional code in it for the special case of copying the pitch angle instead of the omega0 value.
			if ( branch_tag->hasOption("pitch_from_helix") ) {
				runtime_assert_string_msg(
					!branch_tag->hasOption("omega0") &&
					!branch_tag->hasOption("omega0_copies_helix") &&
					!branch_tag->hasOption("omega0_min") &&
					!branch_tag->hasOption("omega0_max") &&
					!branch_tag->hasOption("omega0_samples"),
					"When parsing options for the BundleGridSampler mover, found \"pitch_from_helix\" alongside omega0 options.  This does not make sense.  EITHER a helix copies its pitch angle from another, OR the omega0 value can be set/sampled/copied."
				);
				core::Size const val( branch_tag->getOption<core::Size>("pitch_from_helix", 0) );
				runtime_assert_string_msg(val > 0 && val < helix_index,
					"When parsing options for the BundleGridSampler mover, found a \"pitch_from_helix\" option with a nonsensical value.  Please specify an already-defined helix index.");
				if ( TR.visible() ) TR << "Setting omega0 for helix " << helix_index << " to be set to match the pitch angle for helix " << val << "." << std::endl;
				omega0(helix_index)->set_perturbable(false);
				omega0(helix_index)->set_use_defaults(false);
				omega0(helix_index)->set_helix_to_copy(val);
				omega0(helix_index)->set_omega0_copies_pitch_instead(true); //We're going to copy the pitch angle instead of the omega0 value.
			} else { //All that follows resembles the code for the other parameters.
				if ( branch_tag->hasOption("omega0") ) {
					runtime_assert_string_msg(!branch_tag->hasOption("omega0_copies_helix"),
						"When parsing options for the BundleGridSampler mover, found omega0 defined alongside an omega0_copies_helix statement.  This does not makes sense: either omega0 is fixed, or its value copies the omega0 value of another helix (not both).");
					runtime_assert_string_msg(!branch_tag->hasOption("omega0_min") && !branch_tag->hasOption("omega0_max"),
						"When parsing options for the BundleGridSampler mover, found omega0 defined alongside omega0_min or omega0_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
					core::Real const val( branch_tag->getOption<core::Real>("omega0", 0.0) );
					if ( TR.visible() ) TR << "Setting the omega0 value for helix " << helix_index << " to " << val << "." << std::endl;
					omega0(helix_index)->set_default_value( make_bundle_->convert_angle(val) );
					omega0(helix_index)->set_perturbable(false);
					omega0(helix_index)->set_use_defaults(false);
					make_bundle_->helix(helix_index)->set_omega0( make_bundle_->convert_angle(val) );
				} else if ( branch_tag->hasOption("omega0_min") ) {
					runtime_assert_string_msg(!branch_tag->hasOption("omega0_copies_helix"),
						"When parsing options for the BundleGridSampler mover, found omega0_min defined alongside an omega0_copies_helix statement.  This does not makes sense: either omega0 is sampled, or its value copies the omega0 value of another helix (not both).");
					runtime_assert_string_msg(!branch_tag->hasOption("omega0"),
						"When parsing options for the BundleGridSampler mover, found omega0 defined alongside omega0_min or omega0_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
					runtime_assert_string_msg(branch_tag->hasOption("omega0_max"),
						"When parsing options for the BundleGridSampler mover, found omega0_min but no omega0_max.  This does not make sense.");
					core::Real const val1( branch_tag->getOption<core::Real>("omega0_min", 0.0) );
					core::Real const val2( branch_tag->getOption<core::Real>("omega0_max", 0.0) );
					runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an omega0_max value less than or equal to the omega0_min value.");
					if ( TR.visible() ) TR << "Setting omega0_min and omega0_max values for helix " << helix_index << " to " << val1 << " and " << val2 << ", respectively." << std::endl;
					omega0(helix_index)->set_lower_value( make_bundle_->convert_angle(val1) );
					omega0(helix_index)->set_upper_value( make_bundle_->convert_angle(val2) );
					omega0(helix_index)->set_perturbable(true);
					omega0(helix_index)->set_use_defaults(false);
					make_bundle_->helix(helix_index)->set_omega0( make_bundle_->convert_angle(val1) );
					runtime_assert_string_msg( branch_tag->hasOption("omega0_samples"),
						"When parsing options for the BundleGridSampler mover, found omega0_min and omega0_max options, but no omega0_samples option.  The number of omega0 samples must be specified." );
					core::Size const val3( branch_tag->getOption<core::Size>( "omega0_samples", 0 ) );
					runtime_assert_string_msg( val3 > 1,
						"The number of omega0 samples must be greater than 1." );
					if ( TR.visible() ) TR << "Setting omega0 samples for helix " << helix_index << " to " << val3 << "." << std::endl;
					omega0(helix_index)->set_samples(val3);
				} else if ( branch_tag->hasOption("omega0_max") ) {
					runtime_assert_string_msg(branch_tag->hasOption("omega0_min"),
						"When parsing options for the BundleGridSampler mover, found omega0_max but no omega0_min.  This does not make sense.");
				} else if ( branch_tag->hasOption("omega0_copies_helix") ) {
					core::Size const val( branch_tag->getOption<core::Size>("omega0_copies_helix", 0) );
					runtime_assert_string_msg(val > 0 && val < helix_index,
						"When parsing options for the BundleGridSampler mover, found an omega0_copies_helix option with a nonsensical value.  Please specify an already-defined helix index.");
					if ( TR.visible() ) TR << "Setting omega0 for helix " << helix_index << " to copy the omega0 value for helix " << val << "." << std::endl;
					omega0(helix_index)->set_perturbable(false);
					omega0(helix_index)->set_use_defaults(false);
					omega0(helix_index)->set_helix_to_copy(val);
				}
			}

			if ( branch_tag->hasOption("delta_omega0") ) {
				runtime_assert_string_msg(!branch_tag->hasOption("delta_omega0_copies_helix"),
					"When parsing options for the BundleGridSampler mover, found delta_omega0 defined alongside an delta_omega0_copies_helix statement.  This does not makes sense: either delta_omega0 is fixed, or its value copies the delta_omega0 value of another helix (not both).");
				runtime_assert_string_msg(!branch_tag->hasOption("delta_omega0_min") && !branch_tag->hasOption("delta_omega0_max"),
					"When parsing options for the BundleGridSampler mover, found delta_omega0 defined alongside delta_omega0_min or delta_omega0_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
				core::Real const val( branch_tag->getOption<core::Real>("delta_omega0", 0.0) );
				if ( TR.visible() ) TR << "Setting the delta_omega0 value for helix " << helix_index << " to " << val << "." << std::endl;
				delta_omega0(helix_index)->set_default_value( make_bundle_->convert_angle(val) );
				delta_omega0(helix_index)->set_perturbable(false);
				delta_omega0(helix_index)->set_use_defaults(false);
				make_bundle_->helix(helix_index)->set_delta_omega0( make_bundle_->convert_angle(val) );
			} else if ( branch_tag->hasOption("delta_omega0_min") ) {
				runtime_assert_string_msg(!branch_tag->hasOption("delta_omega0_copies_helix"),
					"When parsing options for the BundleGridSampler mover, found delta_omega0_min defined alongside an delta_omega0_copies_helix statement.  This does not makes sense: either delta_omega0 is sampled, or its value copies the delta_omega0 value of another helix (not both).");
				runtime_assert_string_msg(!branch_tag->hasOption("delta_omega0"),
					"When parsing options for the BundleGridSampler mover, found delta_omega0 defined alongside delta_omega0_min or delta_omega0_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
				runtime_assert_string_msg(branch_tag->hasOption("delta_omega0_max"),
					"When parsing options for the BundleGridSampler mover, found delta_omega0_min but no delta_omega0_max.  This does not make sense.");
				core::Real const val1( branch_tag->getOption<core::Real>("delta_omega0_min", 0.0) );
				core::Real const val2( branch_tag->getOption<core::Real>("delta_omega0_max", 0.0) );
				runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an delta_omega0_max value less than or equal to the delta_omega0_min value.");
				if ( TR.visible() ) TR << "Setting delta_omega0_min and delta_omega0_max values for helix " << helix_index << " to " << val1 << " and " << val2 << ", respectively." << std::endl;
				delta_omega0(helix_index)->set_lower_value( make_bundle_->convert_angle(val1) );
				delta_omega0(helix_index)->set_upper_value( make_bundle_->convert_angle(val2) );
				delta_omega0(helix_index)->set_perturbable(true);
				delta_omega0(helix_index)->set_use_defaults(false);
				make_bundle_->helix(helix_index)->set_delta_omega0( make_bundle_->convert_angle(val1) );
				runtime_assert_string_msg( branch_tag->hasOption("delta_omega0_samples"),
					"When parsing options for the BundleGridSampler mover, found delta_omega0_min and delta_omega0_max options, but no delta_omega0_samples option.  The number of delta_omega0 samples must be specified." );
				core::Size const val3( branch_tag->getOption<core::Size>( "delta_omega0_samples", 0 ) );
				runtime_assert_string_msg( val3 > 1,
					"The number of delta_omega0 samples must be greater than 1." );
				if ( TR.visible() ) TR << "Setting delta_omega0 samples for helix " << helix_index << " to " << val3 << "." << std::endl;
				delta_omega0(helix_index)->set_samples(val3);
			} else if ( branch_tag->hasOption("delta_omega0_max") ) {
				runtime_assert_string_msg(branch_tag->hasOption("delta_omega0_min"),
					"When parsing options for the BundleGridSampler mover, found delta_omega0_max but no delta_omega0_min.  This does not make sense.");
			} else if ( branch_tag->hasOption("delta_omega0_copies_helix") ) {
				core::Size const val( branch_tag->getOption<core::Size>("delta_omega0_copies_helix", 0) );
				runtime_assert_string_msg(val > 0 && val < helix_index,
					"When parsing options for the BundleGridSampler mover, found an delta_omega0_copies_helix option with a nonsensical value.  Please specify an already-defined helix index.");
				if ( TR.visible() ) TR << "Setting delta_omega0 for helix " << helix_index << " to copy the delta_omega0 value for helix " << val << "." << std::endl;
				delta_omega0(helix_index)->set_perturbable(false);
				delta_omega0(helix_index)->set_use_defaults(false);
				delta_omega0(helix_index)->set_helix_to_copy(val);
			}

			if ( branch_tag->hasOption("delta_omega1") ) {
				runtime_assert_string_msg(!branch_tag->hasOption("delta_omega1_copies_helix"),
					"When parsing options for the BundleGridSampler mover, found delta_omega1 defined alongside an delta_omega1_copies_helix statement.  This does not makes sense: either delta_omega1 is fixed, or its value copies the delta_omega1 value of another helix (not both).");
				runtime_assert_string_msg(!branch_tag->hasOption("delta_omega1_min") && !branch_tag->hasOption("delta_omega1_max"),
					"When parsing options for the BundleGridSampler mover, found delta_omega1 defined alongside delta_omega1_min or delta_omega1_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
				core::Real const val( branch_tag->getOption<core::Real>("delta_omega1", 0.0) );
				if ( TR.visible() ) TR << "Setting the delta_omega1 value for helix " << helix_index << " to " << val << "." << std::endl;
				delta_omega1(helix_index)->set_default_value( make_bundle_->convert_angle(val) );
				delta_omega1(helix_index)->set_perturbable(false);
				delta_omega1(helix_index)->set_use_defaults(false);
				make_bundle_->helix(helix_index)->set_delta_omega1_all( make_bundle_->convert_angle(val) );
			} else if ( branch_tag->hasOption("delta_omega1_min") ) {
				runtime_assert_string_msg(!branch_tag->hasOption("delta_omega1_copies_helix"),
					"When parsing options for the BundleGridSampler mover, found delta_omega1_min defined alongside an delta_omega1_copies_helix statement.  This does not makes sense: either delta_omega1 is sampled, or its value copies the delta_omega1 value of another helix (not both).");
				runtime_assert_string_msg(!branch_tag->hasOption("delta_omega1"),
					"When parsing options for the BundleGridSampler mover, found delta_omega1 defined alongside delta_omega1_min or delta_omega1_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
				runtime_assert_string_msg(branch_tag->hasOption("delta_omega1_max"),
					"When parsing options for the BundleGridSampler mover, found delta_omega1_min but no delta_omega1_max.  This does not make sense.");
				core::Real const val1( branch_tag->getOption<core::Real>("delta_omega1_min", 0.0) );
				core::Real const val2( branch_tag->getOption<core::Real>("delta_omega1_max", 0.0) );
				runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an delta_omega1_max value less than or equal to the delta_omega1_min value.");
				if ( TR.visible() ) TR << "Setting delta_omega1_min and delta_omega1_max values for helix " << helix_index << " to " << val1 << " and " << val2 << ", respectively." << std::endl;
				delta_omega1(helix_index)->set_lower_value( make_bundle_->convert_angle(val1) );
				delta_omega1(helix_index)->set_upper_value( make_bundle_->convert_angle(val2) );
				delta_omega1(helix_index)->set_perturbable(true);
				delta_omega1(helix_index)->set_use_defaults(false);
				make_bundle_->helix(helix_index)->set_delta_omega1_all( make_bundle_->convert_angle(val1) );
				runtime_assert_string_msg( branch_tag->hasOption("delta_omega1_samples"),
					"When parsing options for the BundleGridSampler mover, found delta_omega1_min and delta_omega1_max options, but no delta_omega1_samples option.  The number of delta_omega1 samples must be specified." );
				core::Size const val3( branch_tag->getOption<core::Size>( "delta_omega1_samples", 0 ) );
				runtime_assert_string_msg( val3 > 1,
					"The number of delta_omega1 samples must be greater than 1." );
				if ( TR.visible() ) TR << "Setting delta_omega1 samples for helix " << helix_index << " to " << val3 << "." << std::endl;
				delta_omega1(helix_index)->set_samples(val3);
			} else if ( branch_tag->hasOption("delta_omega1_max") ) {
				runtime_assert_string_msg(branch_tag->hasOption("delta_omega1_min"),
					"When parsing options for the BundleGridSampler mover, found delta_omega1_max but no delta_omega1_min.  This does not make sense.");
			} else if ( branch_tag->hasOption("delta_omega1_copies_helix") ) {
				core::Size const val( branch_tag->getOption<core::Size>("delta_omega1_copies_helix", 0) );
				runtime_assert_string_msg(val > 0 && val < helix_index,
					"When parsing options for the BundleGridSampler mover, found an delta_omega1_copies_helix option with a nonsensical value.  Please specify an already-defined helix index.");
				if ( TR.visible() ) TR << "Setting delta_omega1 for helix " << helix_index << " to copy the delta_omega1 value for helix " << val << "." << std::endl;
				delta_omega1(helix_index)->set_perturbable(false);
				delta_omega1(helix_index)->set_use_defaults(false);
				delta_omega1(helix_index)->set_helix_to_copy(val);
			}

			if ( branch_tag->hasOption("delta_t") ) {
				runtime_assert_string_msg(!branch_tag->hasOption("delta_t_copies_helix"),
					"When parsing options for the BundleGridSampler mover, found delta_t defined alongside an delta_t_copies_helix statement.  This does not makes sense: either delta_t is fixed, or its value copies the delta_t value of another helix (not both).");
				runtime_assert_string_msg(!branch_tag->hasOption("delta_t_min") && !branch_tag->hasOption("delta_t_max"),
					"When parsing options for the BundleGridSampler mover, found delta_t defined alongside delta_t_min or delta_t_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
				core::Real const val( branch_tag->getOption<core::Real>("delta_t", 0.0) );
				if ( TR.visible() ) TR << "Setting the delta_t value for helix " << helix_index << " to " << val << "." << std::endl;
				delta_t(helix_index)->set_default_value(val);
				delta_t(helix_index)->set_perturbable(false);
				delta_t(helix_index)->set_use_defaults(false);
				make_bundle_->helix(helix_index)->set_delta_t(val);
			} else if ( branch_tag->hasOption("delta_t_min") ) {
				runtime_assert_string_msg(!branch_tag->hasOption("delta_t_copies_helix"),
					"When parsing options for the BundleGridSampler mover, found delta_t_min defined alongside an delta_t_copies_helix statement.  This does not makes sense: either delta_t is sampled, or its value copies the delta_t value of another helix (not both).");
				runtime_assert_string_msg(!branch_tag->hasOption("delta_t"),
					"When parsing options for the BundleGridSampler mover, found delta_t defined alongside delta_t_min or delta_t_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
				runtime_assert_string_msg(branch_tag->hasOption("delta_t_max"),
					"When parsing options for the BundleGridSampler mover, found delta_t_min but no delta_t_max.  This does not make sense.");
				core::Real const val1( branch_tag->getOption<core::Real>("delta_t_min", 0.0) );
				core::Real const val2( branch_tag->getOption<core::Real>("delta_t_max", 0.0) );
				runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an delta_t_max value less than or equal to the delta_t_min value.");
				if ( TR.visible() ) TR << "Setting delta_t_min and delta_t_max values for helix " << helix_index << " to " << val1 << " and " << val2 << ", respectively." << std::endl;
				delta_t(helix_index)->set_lower_value(val1);
				delta_t(helix_index)->set_upper_value(val2);
				delta_t(helix_index)->set_perturbable(true);
				delta_t(helix_index)->set_use_defaults(false);
				make_bundle_->helix(helix_index)->set_delta_t(val1);
				runtime_assert_string_msg( branch_tag->hasOption("delta_t_samples"),
					"When parsing options for the BundleGridSampler mover, found delta_t_min and delta_t_max options, but no delta_t_samples option.  The number of delta_t samples must be specified." );
				core::Size const val3( branch_tag->getOption<core::Size>( "delta_t_samples", 0 ) );
				runtime_assert_string_msg( val3 > 1,
					"The number of delta_t samples must be greater than 1." );
				if ( TR.visible() ) TR << "Setting delta_t samples for helix " << helix_index << " to " << val3 << "." << std::endl;
				delta_t(helix_index)->set_samples(val3);
			} else if ( branch_tag->hasOption("delta_t_max") ) {
				runtime_assert_string_msg(branch_tag->hasOption("delta_t_min"),
					"When parsing options for the BundleGridSampler mover, found delta_t_max but no delta_t_min.  This does not make sense.");
			} else if ( branch_tag->hasOption("delta_t_copies_helix") ) {
				core::Size const val( branch_tag->getOption<core::Size>("delta_t_copies_helix", 0) );
				runtime_assert_string_msg(val > 0 && val < helix_index,
					"When parsing options for the BundleGridSampler mover, found an delta_t_copies_helix option with a nonsensical value.  Please specify an already-defined helix index.");
				if ( TR.visible() ) TR << "Setting delta_t for helix " << helix_index << " to copy the delta_t value for helix " << val << "." << std::endl;
				delta_t(helix_index)->set_perturbable(false);
				delta_t(helix_index)->set_use_defaults(false);
				delta_t(helix_index)->set_helix_to_copy(val);
			}

			if ( branch_tag->hasOption("z1_offset") ) {
				runtime_assert_string_msg(!branch_tag->hasOption("z1_offset_copies_helix"),
					"When parsing options for the BundleGridSampler mover, found z1_offset defined alongside an z1_offset_copies_helix statement.  This does not makes sense: either z1_offset is fixed, or its value copies the z1_offset value of another helix (not both).");
				runtime_assert_string_msg(!branch_tag->hasOption("z1_offset_min") && !branch_tag->hasOption("z1_offset_max"),
					"When parsing options for the BundleGridSampler mover, found z1_offset defined alongside z1_offset_min or z1_offset_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
				core::Real const val( branch_tag->getOption<core::Real>("z1_offset", 0.0) );
				if ( TR.visible() ) TR << "Setting the z1_offset value for helix " << helix_index << " to " << val << "." << std::endl;
				z1_offset(helix_index)->set_default_value(val);
				z1_offset(helix_index)->set_perturbable(false);
				z1_offset(helix_index)->set_use_defaults(false);
				make_bundle_->helix(helix_index)->set_z1_offset(val);
			} else if ( branch_tag->hasOption("z1_offset_min") ) {
				runtime_assert_string_msg(!branch_tag->hasOption("z1_offset_copies_helix"),
					"When parsing options for the BundleGridSampler mover, found z1_offset_min defined alongside an z1_offset_copies_helix statement.  This does not makes sense: either z1_offset is sampled, or its value copies the z1_offset value of another helix (not both).");
				runtime_assert_string_msg(!branch_tag->hasOption("z1_offset"),
					"When parsing options for the BundleGridSampler mover, found z1_offset defined alongside z1_offset_min or z1_offset_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
				runtime_assert_string_msg(branch_tag->hasOption("z1_offset_max"),
					"When parsing options for the BundleGridSampler mover, found z1_offset_min but no z1_offset_max.  This does not make sense.");
				core::Real const val1( branch_tag->getOption<core::Real>("z1_offset_min", 0.0) );
				core::Real const val2( branch_tag->getOption<core::Real>("z1_offset_max", 0.0) );
				runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an z1_offset_max value less than or equal to the z1_offset_min value.");
				if ( TR.visible() ) TR << "Setting z1_offset_min and z1_offset_max values for helix " << helix_index << " to " << val1 << " and " << val2 << ", respectively." << std::endl;
				z1_offset(helix_index)->set_lower_value(val1);
				z1_offset(helix_index)->set_upper_value(val2);
				z1_offset(helix_index)->set_perturbable(true);
				z1_offset(helix_index)->set_use_defaults(false);
				make_bundle_->helix(helix_index)->set_z1_offset(val1);
				runtime_assert_string_msg( branch_tag->hasOption("z1_offset_samples"),
					"When parsing options for the BundleGridSampler mover, found z1_offset_min and z1_offset_max options, but no z1_offset_samples option.  The number of z1_offset samples must be specified." );
				core::Size const val3( branch_tag->getOption<core::Size>( "z1_offset_samples", 0 ) );
				runtime_assert_string_msg( val3 > 1,
					"The number of z1_offset samples must be greater than 1." );
				if ( TR.visible() ) TR << "Setting z1_offset samples for helix " << helix_index << " to " << val3 << "." << std::endl;
				z1_offset(helix_index)->set_samples(val3);
			} else if ( branch_tag->hasOption("z1_offset_max") ) {
				runtime_assert_string_msg(branch_tag->hasOption("z1_offset_min"),
					"When parsing options for the BundleGridSampler mover, found z1_offset_max but no z1_offset_min.  This does not make sense.");
			} else if ( branch_tag->hasOption("z1_offset_copies_helix") ) {
				core::Size const val( branch_tag->getOption<core::Size>("z1_offset_copies_helix", 0) );
				runtime_assert_string_msg(val > 0 && val < helix_index,
					"When parsing options for the BundleGridSampler mover, found an z1_offset_copies_helix option with a nonsensical value.  Please specify an already-defined helix index.");
				if ( TR.visible() ) TR << "Setting z1_offset for helix " << helix_index << " to copy the z1_offset value for helix " << val << "." << std::endl;
				z1_offset(helix_index)->set_perturbable(false);
				z1_offset(helix_index)->set_use_defaults(false);
				z1_offset(helix_index)->set_helix_to_copy(val);
			}

			if ( branch_tag->hasOption("z0_offset") ) {
				runtime_assert_string_msg(!branch_tag->hasOption("z0_offset_copies_helix"),
					"When parsing options for the BundleGridSampler mover, found z0_offset defined alongside an z0_offset_copies_helix statement.  This does not makes sense: either z0_offset is fixed, or its value copies the z0_offset value of another helix (not both).");
				runtime_assert_string_msg(!branch_tag->hasOption("z0_offset_min") && !branch_tag->hasOption("z0_offset_max"),
					"When parsing options for the BundleGridSampler mover, found z0_offset defined alongside z0_offset_min or z0_offset_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
				core::Real const val( branch_tag->getOption<core::Real>("z0_offset", 0.0) );
				if ( TR.visible() ) TR << "Setting the z0_offset value for helix " << helix_index << " to " << val << "." << std::endl;
				z0_offset(helix_index)->set_default_value(val);
				z0_offset(helix_index)->set_perturbable(false);
				z0_offset(helix_index)->set_use_defaults(false);
				make_bundle_->helix(helix_index)->set_z0_offset(val);
			} else if ( branch_tag->hasOption("z0_offset_min") ) {
				runtime_assert_string_msg(!branch_tag->hasOption("z0_offset_copies_helix"),
					"When parsing options for the BundleGridSampler mover, found z0_offset_min defined alongside an z0_offset_copies_helix statement.  This does not makes sense: either z0_offset is sampled, or its value copies the z0_offset value of another helix (not both).");
				runtime_assert_string_msg(!branch_tag->hasOption("z0_offset"),
					"When parsing options for the BundleGridSampler mover, found z0_offset defined alongside z0_offset_min or z0_offset_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
				runtime_assert_string_msg(branch_tag->hasOption("z0_offset_max"),
					"When parsing options for the BundleGridSampler mover, found z0_offset_min but no z0_offset_max.  This does not make sense.");
				core::Real const val1( branch_tag->getOption<core::Real>("z0_offset_min", 0.0) );
				core::Real const val2( branch_tag->getOption<core::Real>("z0_offset_max", 0.0) );
				runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an z0_offset_max value less than or equal to the z0_offset_min value.");
				if ( TR.visible() ) TR << "Setting z0_offset_min and z0_offset_max values for helix " << helix_index << " to " << val1 << " and " << val2 << ", respectively." << std::endl;
				z0_offset(helix_index)->set_lower_value(val1);
				z0_offset(helix_index)->set_upper_value(val2);
				z0_offset(helix_index)->set_perturbable(true);
				z0_offset(helix_index)->set_use_defaults(false);
				make_bundle_->helix(helix_index)->set_z0_offset(val1);
				runtime_assert_string_msg( branch_tag->hasOption("z0_offset_samples"),
					"When parsing options for the BundleGridSampler mover, found z0_offset_min and z0_offset_max options, but no z0_offset_samples option.  The number of z0_offset samples must be specified." );
				core::Size const val3( branch_tag->getOption<core::Size>( "z0_offset_samples", 0 ) );
				runtime_assert_string_msg( val3 > 1,
					"The number of z0_offset samples must be greater than 1." );
				if ( TR.visible() ) TR << "Setting z0_offset samples for helix " << helix_index << " to " << val3 << "." << std::endl;
				z0_offset(helix_index)->set_samples(val3);
			} else if ( branch_tag->hasOption("z0_offset_max") ) {
				runtime_assert_string_msg(branch_tag->hasOption("z0_offset_min"),
					"When parsing options for the BundleGridSampler mover, found z0_offset_max but no z0_offset_min.  This does not make sense.");
			} else if ( branch_tag->hasOption("z0_offset_copies_helix") ) {
				core::Size const val( branch_tag->getOption<core::Size>("z0_offset_copies_helix", 0) );
				runtime_assert_string_msg(val > 0 && val < helix_index,
					"When parsing options for the BundleGridSampler mover, found an z0_offset_copies_helix option with a nonsensical value.  Please specify an already-defined helix index.");
				if ( TR.visible() ) TR << "Setting z0_offset for helix " << helix_index << " to copy the z0_offset value for helix " << val << "." << std::endl;
				z0_offset(helix_index)->set_perturbable(false);
				z0_offset(helix_index)->set_use_defaults(false);
				z0_offset(helix_index)->set_helix_to_copy(val);
			}

			if ( branch_tag->hasOption("epsilon") ) {
				runtime_assert_string_msg(!branch_tag->hasOption("epsilon_copies_helix"),
					"When parsing options for the BundleGridSampler mover, found epsilon defined alongside an epsilon_copies_helix statement.  This does not makes sense: either epsilon is fixed, or its value copies the epsilon value of another helix (not both).");
				runtime_assert_string_msg(!branch_tag->hasOption("epsilon_min") && !branch_tag->hasOption("epsilon_max"),
					"When parsing options for the BundleGridSampler mover, found epsilon defined alongside epsilon_min or epsilon_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
				core::Real const val( branch_tag->getOption<core::Real>("epsilon", 1.0) );
				if ( TR.visible() ) TR << "Setting the epsilon value for helix " << helix_index << " to " << val << "." << std::endl;
				epsilon(helix_index)->set_default_value(val);
				epsilon(helix_index)->set_perturbable(false);
				epsilon(helix_index)->set_use_defaults(false);
				make_bundle_->helix(helix_index)->set_epsilon(val);
			} else if ( branch_tag->hasOption("epsilon_min") ) {
				runtime_assert_string_msg(!branch_tag->hasOption("epsilon_copies_helix"),
					"When parsing options for the BundleGridSampler mover, found epsilon_min defined alongside an epsilon_copies_helix statement.  This does not makes sense: either epsilon is sampled, or its value copies the epsilon value of another helix (not both).");
				runtime_assert_string_msg(!branch_tag->hasOption("epsilon"),
					"When parsing options for the BundleGridSampler mover, found epsilon defined alongside epsilon_min or epsilon_max.  This does not make sense -- it suggests that the value should both be sampled and not sampled.");
				runtime_assert_string_msg(branch_tag->hasOption("epsilon_max"),
					"When parsing options for the BundleGridSampler mover, found epsilon_min but no epsilon_max.  This does not make sense.");
				core::Real const val1( branch_tag->getOption<core::Real>("epsilon_min", 1.0) );
				core::Real const val2( branch_tag->getOption<core::Real>("epsilon_max", 1.0) );
				runtime_assert_string_msg( val2 > val1, "When parsing options for the BundleGridSampler mover, found an epsilon_max value less than or equal to the epsilon_min value.");
				if ( TR.visible() ) TR << "Setting epsilon_min and epsilon_max values for helix " << helix_index << " to " << val1 << " and " << val2 << ", respectively." << std::endl;
				epsilon(helix_index)->set_lower_value(val1);
				epsilon(helix_index)->set_upper_value(val2);
				epsilon(helix_index)->set_perturbable(true);
				epsilon(helix_index)->set_use_defaults(false);
				make_bundle_->helix(helix_index)->set_epsilon(val1);
				runtime_assert_string_msg( branch_tag->hasOption("epsilon_samples"),
					"When parsing options for the BundleGridSampler mover, found epsilon_min and epsilon_max options, but no epsilon_samples option.  The number of epsilon samples must be specified." );
				core::Size const val3( branch_tag->getOption<core::Size>( "epsilon_samples", 0 ) );
				runtime_assert_string_msg( val3 > 1,
					"The number of epsilon samples must be greater than 1." );
				if ( TR.visible() ) TR << "Setting epsilon samples for helix " << helix_index << " to " << val3 << "." << std::endl;
				epsilon(helix_index)->set_samples(val3);
			} else if ( branch_tag->hasOption("epsilon_max") ) {
				runtime_assert_string_msg(branch_tag->hasOption("epsilon_min"),
					"When parsing options for the BundleGridSampler mover, found epsilon_max but no epsilon_min.  This does not make sense.");
			} else if ( branch_tag->hasOption("epsilon_copies_helix") ) {
				core::Size const val( branch_tag->getOption<core::Size>("epsilon_copies_helix", 0) );
				runtime_assert_string_msg(val > 0 && val < helix_index,
					"When parsing options for the BundleGridSampler mover, found an epsilon_copies_helix option with a nonsensical value.  Please specify an already-defined helix index.");
				if ( TR.visible() ) TR << "Setting epsilon for helix " << helix_index << " to copy the epsilon value for helix " << val << "." << std::endl;
				epsilon(helix_index)->set_perturbable(false);
				epsilon(helix_index)->set_use_defaults(false);
				epsilon(helix_index)->set_helix_to_copy(val);
			}

		} //if ( (*tag_it)->getName() == "Helix" )
	}
	runtime_assert_string_msg(at_least_one_helix, "When parsing options for the BundleGridSampler mover, found no helices!  Define at least one helix with a <Helix...> sub-tag.");

	if ( TR.visible() ) TR.flush();
	return;
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
core::Size BundleGridSampler::add_helix( ) {

	increment_helix_count();
	core::Size const nhelices(n_helices());

	r0_.push_back( PerturbBundleOptionsOP( new PerturbBundleOptions ) );
	r0(nhelices)->set_helix_index(nhelices);
	omega0_.push_back( PerturbBundleOptionsOP( new PerturbBundleOptions ) );
	omega0(nhelices)->set_helix_index(nhelices);
	delta_omega0_.push_back( PerturbBundleOptionsOP( new PerturbBundleOptions ) );
	delta_omega0(nhelices)->set_helix_index(nhelices);
	delta_omega1_.push_back( PerturbBundleOptionsOP( new PerturbBundleOptions ) );
	delta_omega1(nhelices)->set_helix_index(nhelices);
	delta_t_.push_back( PerturbBundleOptionsOP( new PerturbBundleOptions ) );
	delta_t(nhelices)->set_helix_index(nhelices);
	z1_offset_.push_back( PerturbBundleOptionsOP( new PerturbBundleOptions ) );
	z1_offset(nhelices)->set_helix_index(nhelices);
	z0_offset_.push_back( PerturbBundleOptionsOP( new PerturbBundleOptions ) );
	z0_offset(nhelices)->set_helix_index(nhelices);
	epsilon_.push_back( PerturbBundleOptionsOP( new PerturbBundleOptions ) );
	epsilon(nhelices)->set_helix_index(nhelices);

	runtime_assert_string_msg( r0_.size()==nhelices && omega0_.size()==nhelices && delta_omega0_.size()==nhelices && delta_omega1_.size()==nhelices && delta_t_.size()==nhelices && z1_offset_.size()==nhelices && z0_offset_.size()==nhelices && epsilon_.size()==nhelices,
		"In protocols::helical_bundle::BundleGridSampler::add_helix() function: somehow, vector indices are out of sync.  I can't determine how many helices have been defined." );

	//Update the MakeBundle mover, too:
	make_bundle_->add_helix();
	runtime_assert_string_msg( make_bundle_->n_helices()==nhelices,
		"In protocols::helical_bundle::BundleGridSampler::add_helix() function: somehow, the MakeBundle mover is out of sync with the BundleGridSampler mover.  Consult a developer or an exorcist -- this shouldn't happen.");

	return nhelices;
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
core::Size BundleGridSampler::calculate_total_samples() const {
	core::Size total_samples(1);
	core::Size const helixcount( n_helices() );

	//Loop through all helices that are defined for r0 sampling:
	//(Note that all helices are defined.)
	debug_assert( r0_.size()==helixcount );
	debug_assert( omega0_.size()==helixcount );
	debug_assert( delta_omega0_.size()==helixcount );
	debug_assert( delta_omega1_.size()==helixcount );
	debug_assert( delta_t_.size()==helixcount );
	debug_assert( z1_offset_.size()==helixcount );
	debug_assert( z0_offset_.size()==helixcount );
	debug_assert( epsilon_.size()==helixcount );
	for ( core::Size i=1, imax=helixcount; i<=imax; ++i ) {
		if ( r0(i)->is_perturbable() && !r0(i)->use_defaults() && r0(i)->other_helix()==0 && r0(i)->samples()!=0  ) total_samples *= r0(i)->samples();
		else if ( r0(i)->use_defaults() && default_r0()->is_perturbable() && default_r0()->samples()!=0 ) total_samples *= default_r0()->samples();

		if ( omega0(i)->is_perturbable() && !omega0(i)->use_defaults() && omega0(i)->other_helix()==0 && omega0(i)->samples()!=0  ) total_samples *= omega0(i)->samples();
		else if ( omega0(i)->use_defaults() && default_omega0()->is_perturbable() && default_omega0()->samples()!=0 ) total_samples *= default_omega0()->samples();

		if ( delta_omega0(i)->is_perturbable() && !delta_omega0(i)->use_defaults() && delta_omega0(i)->other_helix()==0 && delta_omega0(i)->samples()!=0  ) total_samples *= delta_omega0(i)->samples();
		else if ( delta_omega0(i)->use_defaults() && default_delta_omega0()->is_perturbable() && default_delta_omega0()->samples()!=0 ) total_samples *= default_delta_omega0()->samples();

		if ( delta_omega1(i)->is_perturbable() && !delta_omega1(i)->use_defaults() && delta_omega1(i)->other_helix()==0 && delta_omega1(i)->samples()!=0  ) total_samples *= delta_omega1(i)->samples();
		else if ( delta_omega1(i)->use_defaults() && default_delta_omega1()->is_perturbable() && default_delta_omega1()->samples()!=0 ) total_samples *= default_delta_omega1()->samples();

		if ( delta_t(i)->is_perturbable() && !delta_t(i)->use_defaults() && delta_t(i)->other_helix()==0 && delta_t(i)->samples()!=0  ) total_samples *= delta_t(i)->samples();
		else if ( delta_t(i)->use_defaults() && default_delta_t()->is_perturbable() && default_delta_t()->samples()!=0 ) total_samples *= default_delta_t()->samples();

		if ( z1_offset(i)->is_perturbable() && !z1_offset(i)->use_defaults() && z1_offset(i)->other_helix()==0 && z1_offset(i)->samples()!=0  ) total_samples *= z1_offset(i)->samples();
		else if ( z1_offset(i)->use_defaults() && default_z1_offset()->is_perturbable() && default_z1_offset()->samples()!=0 ) total_samples *= default_z1_offset()->samples();

		if ( z0_offset(i)->is_perturbable() && !z0_offset(i)->use_defaults() && z0_offset(i)->other_helix()==0 && z0_offset(i)->samples()!=0  ) total_samples *= z0_offset(i)->samples();
		else if ( z0_offset(i)->use_defaults() && default_z0_offset()->is_perturbable() && default_z0_offset()->samples()!=0 ) total_samples *= default_z0_offset()->samples();

		if ( epsilon(i)->is_perturbable() && !epsilon(i)->use_defaults() && epsilon(i)->other_helix()==0 && epsilon(i)->samples()!=0  ) total_samples *= epsilon(i)->samples();
		else if ( epsilon(i)->use_defaults() && default_epsilon()->is_perturbable() && default_epsilon()->samples()!=0 ) total_samples *= default_epsilon()->samples();
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
	add_attributes_for_make_bundle_dofs( attlist );
	add_attributes_for_make_bundle_other_defaults( attlist );
	add_attributes_for_make_bundle_minorhelix_defaults( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default("reset", xsct_rosetta_bool, "If reset is set to \"true\" (the default), then input geometry is discarded and the BundleGridSampler builds a pose from scratch.  If \"false\", then parametric geometry is appended to the input geometry.", "true" )
		+ XMLSchemaAttribute::attribute_w_default("nstruct_mode", xsct_rosetta_bool, "If \"true\", sample a different set of mainchain torsions for each RosettaScripts job (with each successive job sampling the next gridpoint in the grid of parameter values to be sampled).  If \"false\" (the default), then each job consists of the whole mainchain sampling effort.", "false" )
		+ XMLSchemaAttribute::attribute_w_default("nstruct_repeats", xsct_non_negative_integer, "In nstruct_mode, this is the number of times each parameter gridpoint will be sampled.  This defaults to 1 (i.e. each successive RosettaScripts job goes on to the next gridpoint), but can be set higher (i.e. successive RosettaScripts jobs repeat gridpoints).", "1" )
		+ XMLSchemaAttribute::attribute_w_default("r0", xsct_real, "A fixed value for r0, the major helix radius.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("r0_min", xsct_real, "Minimum value for r0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("r0_max", xsct_real, "Maximum value for r0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("r0_samples", xsct_non_negative_integer, "Number of samples for r0.", "0" )
		+ XMLSchemaAttribute::attribute_w_default("omega0", xsct_real, "A fixed value for omega0, the twist about the z-axis.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("omega0_min", xsct_real, "Minimum value for omega0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("omega0_max", xsct_real, "Maximum value for omega0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("omega0_samples", xsct_non_negative_integer, "Number of samples for omega0.", "0" )
		+ XMLSchemaAttribute::attribute_w_default("delta_omega0", xsct_real, "A fixed value for delta_omega0, the rigid-body rotation of the minor helix about the z-axis.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("delta_omega0_min", xsct_real, "Minimum value for delta_omega0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("delta_omega0_max", xsct_real, "Maximum value for delta_omega0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("delta_omega0_samples", xsct_non_negative_integer, "Number of samples for delta_omega0.", "0" )
		+ XMLSchemaAttribute::attribute_w_default("delta_omega1", xsct_real, "A fixed value for delta_omega1, the rotation of the minor helix about its own axis.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("delta_omega1_min", xsct_real, "Minimum value for delta_omega1.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("delta_omega1_max", xsct_real, "Maximum value for delta_omega1.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("delta_omega1_samples", xsct_non_negative_integer, "Number of samples for delta_omega1.", "0" )
		+ XMLSchemaAttribute::attribute_w_default("delta_t", xsct_real, "A fixed value for delta_t, the offset along the polypeptide backbone.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("delta_t_min", xsct_real, "Minimum value for delta_t.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("delta_t_max", xsct_real, "Maximum value for delta_t.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("delta_t_samples", xsct_non_negative_integer, "Number of samples for delta_t.", "0" )
		+ XMLSchemaAttribute::attribute_w_default("z1_offset", xsct_real, "A fixed value for z1_offset, the translation of the minor helix along its own axis.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("z1_offset_min", xsct_real, "Maximum value for z1_offset.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("z1_offset_max", xsct_real, "Maximum value for z1_offset.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("z1_offset_samples", xsct_non_negative_integer, "Number of samples for z1_offset.", "0" )
		+ XMLSchemaAttribute::attribute_w_default("z0_offset", xsct_real, "A fixed value for z0_offset, the translation of the minor helix along the z-axis.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("z0_offset_min", xsct_real, "Minimum value for z0_offset.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("z0_offset_max", xsct_real, "Maximum value for z0_offset.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("z0_offset_samples", xsct_non_negative_integer, "Number of samples for z0_offset.", "0" )
		+ XMLSchemaAttribute::attribute_w_default("epsilon", xsct_real, "A fixed value for epsilon, the lateral squash of the bundle.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("epsilon_min", xsct_real, "Minimum value for epsilon.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("epsilon_max", xsct_real, "Maximum value for epsilon.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("epsilon_samples", xsct_non_negative_integer, "Number of samples for epsilon.", "0" );
	if ( ! attribute_w_name_in_attribute_list( "z1", attlist ) ) {
		attlist + XMLSchemaAttribute::attribute_w_default("z1", xsct_real, "A fixed value for z1, the minor helix rise per residue.  Note that this is not normally set by the user!", "0.0" );
	}

	AttributeList subtag_attributes;

	add_attributes_for_helix_params( subtag_attributes );
	add_attributes_for_minor_helix_params( subtag_attributes );
	add_attributes_for_other_helix_params( subtag_attributes );

	subtag_attributes + XMLSchemaAttribute::attribute_w_default( "r0", xsct_real, "Single value for r0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "r0_copies_helix", xsct_non_negative_integer, "Helix index from which r0 should be copied.", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "r0_min", xsct_real, "Minimum value for r0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "r0_max", xsct_real, "Maximum value for r0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "r0_samples", xsct_non_negative_integer, "Number of samples for r0.", "0" )

		+ XMLSchemaAttribute::attribute_w_default( "omega0", xsct_real, "Single value for omega0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "omega0_copies_helix", xsct_non_negative_integer, "Helix index from which omega0 should be copied.", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "omega0_min", xsct_real, "Minimum value for omega0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "omega0_max", xsct_real, "Maximum value for omega0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "omega0_samples", xsct_non_negative_integer, "Number of samples for delta_omega0.", "0" )

		+ XMLSchemaAttribute::attribute_w_default( "delta_omega0", xsct_real, "Single value for delta_omega0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "delta_omega0_copies_helix", xsct_non_negative_integer, "Helix index from which delta_omega0 should be copied.", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "delta_omega0_min", xsct_real, "Minimum value for delta_omega0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "delta_omega0_max", xsct_real, "Maximum value for delta_omega0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "delta_omega0_samples", xsct_non_negative_integer, "Number of samples for delta_omega0.", "0" );
	if ( ! attribute_w_name_in_attribute_list( "delta_omega1", subtag_attributes ) ) {
		subtag_attributes + XMLSchemaAttribute::attribute_w_default( "delta_omega1", xsct_real, "Single value for delta_omega1.", "0.0" );
	}
	subtag_attributes + XMLSchemaAttribute::attribute_w_default( "delta_omega1_copies_helix", xsct_non_negative_integer, "Helix index from which delta_omega1 should be copied.", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "delta_omega1_min", xsct_real, "Minimum value for delta_omega1.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "delta_omega1_max", xsct_real, "Maximum value for delta_omega1.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "delta_omega1_samples", xsct_non_negative_integer, "Number of samples for delta_omega1.", "0" )

		+ XMLSchemaAttribute::attribute_w_default( "delta_t", xsct_real, "Single value for delta_t.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "delta_t_copies_helix", xsct_non_negative_integer, "Helix index from which delta_t should be copied.", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "delta_t_min", xsct_real, "Minimum value for delta_t.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "delta_t_max", xsct_real, "Maximum value for delta_t.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "delta_t_samples", xsct_non_negative_integer, "Number of samples for delta_t.", "0" );
	if ( ! attribute_w_name_in_attribute_list( "z1", attlist ) ) {
		subtag_attributes + XMLSchemaAttribute::attribute_w_default( "z1", xsct_real, "A fixed value for z1, the minor helix rise per residue.  Note that this is not normally set by the user!", "0.0" );
	}

	subtag_attributes
		+ XMLSchemaAttribute::attribute_w_default( "z1_offset", xsct_real, "Single value for z1.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "z1_offset_copies_helix", xsct_non_negative_integer, "Helix index from which z1 should be copied.", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "z1_offset_min", xsct_real, "Minimum value for z1.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "z1_offset_max", xsct_real, "Maximum value for z1.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "z1_offset_samples", xsct_non_negative_integer, "Number of samples for z1.", "0" )

		+ XMLSchemaAttribute::attribute_w_default( "z0_offset", xsct_real, "Single value for z1.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "z0_offset_copies_helix", xsct_non_negative_integer, "Helix index from which z0 should be copied.", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "z0_offset_min", xsct_real, "Minimum value for z0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "z0_offset_max", xsct_real, "Maximum value for z0.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "z0_offset_samples", xsct_non_negative_integer, "Number of samples for z0.", "0" )

		+ XMLSchemaAttribute::attribute_w_default( "epsilon", xsct_real, "Single value for z1.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "epsilon_copies_helix", xsct_non_negative_integer, "Helix index from which epsilon should be copied.", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "epsilon_min", xsct_real, "Minimum value for epsilon.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "epsilon_max", xsct_real, "Maximum value for epsilon.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "epsilon_samples", xsct_non_negative_integer, "Number of samples for epsilon.", "0" )

		+ XMLSchemaAttribute::attribute_w_default( "pitch_from_helix", xsct_non_negative_integer, "Helix index from which pitch should be copied.", "0" );

	utility::tag::XMLSchemaSimpleSubelementList ssl;
	ssl.add_simple_subelement( "Helix", subtag_attributes, "Tags describing individual helices in the bundle"/*, 0 minoccurs*/ )
		.complex_type_naming_func( & subtag_for_bundgrid );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, ssl );
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
