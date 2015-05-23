// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/PerturbBundle.cc
/// @brief  Perturbs a helical bundle by altering the Crick parameters.
/// @details The bundle is centred on the origin, with the outer helix axis pointing along the
/// global z-axis.  This mover calls the PerturbBundleHelix mover.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/helical_bundle/PerturbBundle.hh>
#include <protocols/helical_bundle/PerturbBundleCreator.hh>
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <numeric/crick_equations/BundleParams.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <utility/tag/Tag.hh>

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
#include <protocols/helical_bundle/PerturbBundleHelix.fwd.hh>
#include <protocols/helical_bundle/PerturbBundleHelix.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace helical_bundle {

static thread_local basic::Tracer TR("protocols.helical_bundle.PerturbBundle");

std::string
PerturbBundleCreator::keyname() const
{
	return PerturbBundleCreator::mover_name();
}

protocols::moves::MoverOP
PerturbBundleCreator::create_mover() const {
	return protocols::moves::MoverOP( new PerturbBundle );
}

std::string
PerturbBundleCreator::mover_name()
{
	return "PerturbBundle";
}


/// @brief Creator for PerturbBundle mover.
PerturbBundle::PerturbBundle():
		Mover("PerturbBundle"),
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
		bundleparametersset_index_(1),
		use_degrees_(false)
{}


/// @brief Copy constructor for PerturbBundle mover.
PerturbBundle::PerturbBundle( PerturbBundle const & src ):
	protocols::moves::Mover( src ),
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
	default_z1_offset_( src.default_z1_offset_->clone() ),
	z1_offset_(),
	default_z0_offset_( src.default_z0_offset_->clone() ),
	z0_offset_(),
	bundleparametersset_index_(src.bundleparametersset_index_),
	use_degrees_( src.use_degrees_ )
{
	r0_.clear();
	omega0_.clear();
	delta_omega0_.clear();
	delta_omega1_.clear();
	delta_t_.clear();
	z1_offset_.clear();
	z0_offset_.clear();
	for(core::Size i=1,imax=src.r0_.size(); i<=imax; ++i) r0_.push_back( src.r0_[i]->clone() );
	for(core::Size i=1,imax=src.omega0_.size(); i<=imax; ++i) omega0_.push_back( src.omega0_[i]->clone() );
	for(core::Size i=1,imax=src.delta_omega0_.size(); i<=imax; ++i) delta_omega0_.push_back( src.delta_omega0_[i]->clone() );
	for(core::Size i=1,imax=src.delta_omega1_.size(); i<=imax; ++i) delta_omega1_.push_back( src.delta_omega1_[i]->clone() );
	for(core::Size i=1,imax=src.delta_t_.size(); i<=imax; ++i) delta_t_.push_back( src.delta_t_[i]->clone() );
	for(core::Size i=1,imax=src.z1_offset_.size(); i<=imax; ++i) z1_offset_.push_back( src.z1_offset_[i]->clone() );
	for(core::Size i=1,imax=src.z0_offset_.size(); i<=imax; ++i) z0_offset_.push_back( src.z0_offset_[i]->clone() );
}


/// @brief Destructor for PerturbBundle mover.
PerturbBundle::~PerturbBundle() {}


/// @brief Clone operator to create a pointer to a fresh PerturbBundle object that copies this one.
protocols::moves::MoverOP PerturbBundle::clone() const {
	return protocols::moves::MoverOP( new PerturbBundle( *this ) );
}


/// @brief Fresh_instance operator to create a pointer to a fresh PerturbBundle object that does NOT copy this one.
protocols::moves::MoverOP PerturbBundle::fresh_instance() const {
	return protocols::moves::MoverOP( new PerturbBundle );
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
void PerturbBundle::apply( core::pose::Pose & pose )
{
	bool failed=false;

	if(TR.visible()) TR << "Copying pose." << std::endl;
	core::pose::Pose pose_copy(pose);

	if(TR.visible()) TR << "Finding BundleParametersSet object in pose." << std::endl;
	BundleParametersSetOP params_set;
	core::Size params_set_index(0);
	core::Size n_encountered(0);
	bool breaknow(false);
	for(core::Size i=1, imax=pose_copy.conformation().n_parameters_sets(); i<=imax; ++i) {
		BundleParametersSetOP cur_set( utility::pointer::dynamic_pointer_cast< BundleParametersSet >( pose_copy.conformation().parameters_set(i) ) );
		if(cur_set) {
			++n_encountered; //Increment the number of parameterssets encountered
			if(n_encountered==bundleparametersset_index()) {
				params_set_index=i;
				params_set=cur_set; //Assign the current owning pointer to be the params_set owning pointer.
				breaknow=true;
			}
		}
		if(breaknow) break;
	}
	runtime_assert_string_msg(params_set_index!=0 && params_set, "In protocols::helical_bundle::PerturbBundle::apply() function: BundleparametersSet object with given index not found in pose!");

	write_report( params_set, true); //Write a pre-perturbation report summarizing the initial Crick parameter values.

	if(TR.visible()) TR << "Perturbing Crick equation values." << std::endl;
	if(!perturb_values( params_set )) {
		if(TR.visible()) TR << "Perturbation failed -- senseless values resulted.  Returning input pose." << std::endl;
		failed=true;
	}

	if(!failed) {
		if(TR.visible()) TR << "Rebuilding the helical bundle conformation using the Crick parameters stored in the pose." << std::endl;
		rebuild_conformation(pose_copy, params_set, params_set_index, failed);
	}

	if(!failed) {
		if(TR.visible()) TR << "Perturbation successful.  Copying result to the input pose." << std::endl;
		pose = pose_copy;	
		//write_report( params_set, false); //Write a post-perturbation report summarizing the final Crick parameter values.	
		write_report( utility::pointer::dynamic_pointer_cast< BundleParametersSet >( pose.conformation().parameters_set(params_set_index) ), false); //Write a post-perturbation report summarizing the final Crick parameter values.	
	} else {
		if(TR.visible()) TR << "The current attempt generated Crick parameters that did not permit sensible geometry.  Returning input pose." << std::endl;
	}

	if(TR.Debug.visible()) TR.Debug << "Finished apply function." << std::endl;

	TR.flush(); TR.Debug.flush();

	return;
}

////////////////////////////////////////////////////////////////////////////////


/// @brief Returns the name of this mover ("PerturbBundle").
std::string PerturbBundle::get_name() const{
	return "PerturbBundle";
}

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////

/// @brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
///
void
PerturbBundle::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & /*data_map*/,
		protocols::filters::Filters_map const &/*filters*/,
		protocols::moves::Movers_map const &/*movers*/,
		core::pose::Pose const & /*pose*/
) {

	if ( tag->getName() != "PerturbBundle" ){
		throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible -- the tag name does not match the mover name.");
	}

	if(TR.visible()) TR << "Parsing options for PerturbBundle (\"" << tag->getOption<std::string>("name" ,"") << "\") mover." << std::endl;
	
	//Determine whether input is in degrees or radians.
	set_use_degrees( tag->getOption<bool>( "use_degrees", false ) );
	if(TR.visible()) TR << "Interpreting user-input angles as being in " << (use_degrees() ? "degrees" : "radians") << ".  (Internally, radians are always used, and output will be in radians.)" << std::endl;

	//Set a default perturbation type:
	if( tag->hasOption("default_perturbation_type")) {
		std::string perttype( tag->getOption<std::string>("default_perturbation_type", "") );
		runtime_assert_string_msg( perttype=="gaussian" || perttype=="uniform",
			"In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: allowed perturbation types are \"gaussian\" and \"uniform\"." );
		if(perttype == "gaussian") {
			default_r0()->set_perturbation_type(pt_gaussian);
			default_omega0()->set_perturbation_type(pt_gaussian);
			default_delta_omega0()->set_perturbation_type(pt_gaussian);
			default_delta_omega1()->set_perturbation_type(pt_gaussian);
			default_delta_t()->set_perturbation_type(pt_gaussian);
			default_z1_offset()->set_perturbation_type(pt_gaussian);
			default_z0_offset()->set_perturbation_type(pt_gaussian);
			if(TR.visible()) TR << "Setting default perturbation type to GAUSSIAN." << std::endl;
		} else if (perttype == "uniform") {
			default_r0()->set_perturbation_type(pt_uniform);
			default_omega0()->set_perturbation_type(pt_uniform);
			default_delta_omega0()->set_perturbation_type(pt_uniform);
			default_delta_omega1()->set_perturbation_type(pt_uniform);
			default_delta_t()->set_perturbation_type(pt_uniform);
			default_z1_offset()->set_perturbation_type(pt_uniform);
			default_z0_offset()->set_perturbation_type(pt_uniform);
			if(TR.visible()) TR << "Setting default perturbation type to UNIFORM." << std::endl;
		}
	}

	//Set defaults for the various perturbable degrees of freedom:
	if( tag->hasOption("r0_perturbation") ) {
		core::Real r0pert( tag->getOption<core::Real>("r0_perturbation", 0.0) );
		default_r0()->set_perturbation_magnitude(r0pert);
		if(TR.visible()) TR << "Set r0 perturbation magnitude to " << r0pert << std::endl;
	}
	if( tag->hasOption("r0_perturbation_type") ) {
		if(default_r0()->is_perturbable()) {
			std::string perttype( tag->getOption<std::string>("r0_perturbation_type", "") );
			runtime_assert_string_msg( perttype=="gaussian" || perttype=="uniform",
				"In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: allowed perturbation types are \"gaussian\" and \"uniform\"." );
			if(perttype=="gaussian") default_r0()->set_perturbation_type(pt_gaussian);
			else if (perttype=="uniform") default_r0()->set_perturbation_type(pt_uniform);
			if(TR.visible()) TR << "Set r0 perturbation type to " << perttype << "." << std::endl;
		} else {
			if(TR.Warning.visible()) TR.Warning << "Warning!  The r0_perturbation_type option was specified, but without an r0_perturbation option to set perturbation magnitude, it will be ignored." << std::endl;
		}
	}

	if( tag->hasOption("omega0_perturbation") ) {
		core::Real omega0pert( tag->getOption<core::Real>("omega0_perturbation", 0.0) );
		default_omega0()->set_perturbation_magnitude( convert_angle( omega0pert ) );
		if(TR.visible()) TR << "Set omega0 perturbation magnitude to " << omega0pert << (use_degrees() ? " degrees." : " radians.") << std::endl;
	}
	if( tag->hasOption("omega0_perturbation_type") ) {
		if(default_omega0()->is_perturbable()) {
			std::string perttype( tag->getOption<std::string>("omega0_perturbation_type", "") );
			runtime_assert_string_msg( perttype=="gaussian" || perttype=="uniform",
				"In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: allowed perturbation types are \"gaussian\" and \"uniform\"." );
			if(perttype=="gaussian") default_omega0()->set_perturbation_type(pt_gaussian);
			else if (perttype=="uniform") default_omega0()->set_perturbation_type(pt_uniform);
			if(TR.visible()) TR << "Set omega0 perturbation type to " << perttype << "." << std::endl;
		} else {
			if(TR.Warning.visible()) TR.Warning << "Warning!  The omega0_perturbation_type option was specified, but without an omega0_perturbation option to set perturbation magnitude, it will be ignored." << std::endl;
		}
	}

	if( tag->hasOption("delta_omega0_perturbation") ) {
		core::Real delta_omega0pert( tag->getOption<core::Real>("delta_omega0_perturbation", 0.0) );
		default_delta_omega0()->set_perturbation_magnitude( convert_angle( delta_omega0pert) );
		if(TR.visible()) TR << "Set delta_omega0 perturbation magnitude to " << delta_omega0pert << (use_degrees() ? " degrees." : " radians.") << std::endl;
	}
	if( tag->hasOption("delta_omega0_perturbation_type") ) {
		if(default_delta_omega0()->is_perturbable()) {
			std::string perttype( tag->getOption<std::string>("delta_omega0_perturbation_type", "") );
			runtime_assert_string_msg( perttype=="gaussian" || perttype=="uniform",
				"In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: allowed perturbation types are \"gaussian\" and \"uniform\"." );
			if(perttype=="gaussian") default_delta_omega0()->set_perturbation_type(pt_gaussian);
			else if (perttype=="uniform") default_delta_omega0()->set_perturbation_type(pt_uniform);
			if(TR.visible()) TR << "Set delta_omega0 perturbation type to " << perttype << "." << std::endl;
		} else {
			if(TR.Warning.visible()) TR.Warning << "Warning!  The delta_omega0_perturbation_type option was specified, but without an delta_omega0_perturbation option to set perturbation magnitude, it will be ignored." << std::endl;
		}
	}

	if( tag->hasOption("delta_omega1_perturbation") ) {
		core::Real delta_omega1pert( tag->getOption<core::Real>("delta_omega1_perturbation", 0.0) );
		default_delta_omega1()->set_perturbation_magnitude( convert_angle( delta_omega1pert ) );
		if(TR.visible()) TR << "Set delta_omega1 perturbation magnitude to " << delta_omega1pert << (use_degrees() ? " degrees." : " radians.") << std::endl;
	}
	if( tag->hasOption("delta_omega1_perturbation_type") ) {
		if(default_delta_omega1()->is_perturbable()) {
			std::string perttype( tag->getOption<std::string>("delta_omega1_perturbation_type", "") );
			runtime_assert_string_msg( perttype=="gaussian" || perttype=="uniform",
				"In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: allowed perturbation types are \"gaussian\" and \"uniform\"." );
			if(perttype=="gaussian") default_delta_omega1()->set_perturbation_type(pt_gaussian);
			else if (perttype=="uniform") default_delta_omega1()->set_perturbation_type(pt_uniform);
			if(TR.visible()) TR << "Set delta_omega1 perturbation type to " << perttype << "." << std::endl;
		} else {
			if(TR.Warning.visible()) TR.Warning << "Warning!  The delta_omega1_perturbation_type option was specified, but without an delta_omega1_perturbation option to set perturbation magnitude, it will be ignored." << std::endl;
		}
	}

	if( tag->hasOption("delta_t_perturbation") ) {
		core::Real delta_tpert( tag->getOption<core::Real>("delta_t_perturbation", 0.0) );
		default_delta_t()->set_perturbation_magnitude(delta_tpert);
		if(TR.visible()) TR << "Set delta_t perturbation magnitude to " << delta_tpert << std::endl;
	}
	if( tag->hasOption("delta_t_perturbation_type") ) {
		if(default_delta_t()->is_perturbable()) {
			std::string perttype( tag->getOption<std::string>("delta_t_perturbation_type", "") );
			runtime_assert_string_msg( perttype=="gaussian" || perttype=="uniform",
				"In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: allowed perturbation types are \"gaussian\" and \"uniform\"." );
			if(perttype=="gaussian") default_delta_t()->set_perturbation_type(pt_gaussian);
			else if (perttype=="uniform") default_delta_t()->set_perturbation_type(pt_uniform);
			if(TR.visible()) TR << "Set delta_t perturbation type to " << perttype << "." << std::endl;
		} else {
			if(TR.Warning.visible()) TR.Warning << "Warning!  The delta_t_perturbation_type option was specified, but without an delta_t_perturbation option to set perturbation magnitude, it will be ignored." << std::endl;
		}
	}
	if( tag->hasOption("z1_offset_perturbation") ) {
		core::Real z1_offsetpert( tag->getOption<core::Real>("z1_offset_perturbation", 0.0) );
		default_z1_offset()->set_perturbation_magnitude(z1_offsetpert);
		if(TR.visible()) TR << "Set z1_offset perturbation magnitude to " << z1_offsetpert << std::endl;
	}
	if( tag->hasOption("z1_offset_perturbation_type") ) {
		if(default_z1_offset()->is_perturbable()) {
			std::string perttype( tag->getOption<std::string>("z1_offset_perturbation_type", "") );
			runtime_assert_string_msg( perttype=="gaussian" || perttype=="uniform",
				"In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: allowed perturbation types are \"gaussian\" and \"uniform\"." );
			if(perttype=="gaussian") default_z1_offset()->set_perturbation_type(pt_gaussian);
			else if (perttype=="uniform") default_z1_offset()->set_perturbation_type(pt_uniform);
			if(TR.visible()) TR << "Set z1_offset perturbation type to " << perttype << "." << std::endl;
		} else {
			if(TR.Warning.visible()) TR.Warning << "Warning!  The z1_offset_perturbation_type option was specified, but without an z1_offset_perturbation option to set perturbation magnitude, it will be ignored." << std::endl;
		}
	}
	if( tag->hasOption("z0_offset_perturbation") ) {
		core::Real z0_offsetpert( tag->getOption<core::Real>("z0_offset_perturbation", 0.0) );
		default_z0_offset()->set_perturbation_magnitude(z0_offsetpert);
		if(TR.visible()) TR << "Set z0_offset perturbation magnitude to " << z0_offsetpert << std::endl;
	}
	if( tag->hasOption("z0_offset_perturbation_type") ) {
		if(default_z0_offset()->is_perturbable()) {
			std::string perttype( tag->getOption<std::string>("z0_offset_perturbation_type", "") );
			runtime_assert_string_msg( perttype=="gaussian" || perttype=="uniform",
				"In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: allowed perturbation types are \"gaussian\" and \"uniform\"." );
			if(perttype=="gaussian") default_z0_offset()->set_perturbation_type(pt_gaussian);
			else if (perttype=="uniform") default_z0_offset()->set_perturbation_type(pt_uniform);
			if(TR.visible()) TR << "Set z0_offset perturbation type to " << perttype << "." << std::endl;
		} else {
			if(TR.Warning.visible()) TR.Warning << "Warning!  The z0_offset_perturbation_type option was specified, but without an z0_offset_perturbation option to set perturbation magnitude, it will be ignored." << std::endl;
		}
	}


	//Parse options for specific helices:
  utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	for( utility::vector1< utility::tag::TagCOP >::const_iterator tag_it=branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it) {
		if ( (*tag_it)->getName() == "Helix" ) { //A helix has been added.  Add it, and parse its options.
			runtime_assert_string_msg( (*tag_it)->hasOption("helix_index"), "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was added, but no helix index has been indicated." );
			core::Size helix_index( (*tag_it)->getOption<core::Size>("helix_index", 0) );
			runtime_assert_string_msg(helix_index>0, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was added, but its index was set to 0.  This is not allowed." );
			core::Size const this_helix( add_helix(helix_index) ); //Add and initialize this helix.  By default, no degrees of freedom may be perturbed.  Store the current helix index in this_helix.
			bool has_perturbable_dofs(false);

			if( (*tag_it)->hasOption("r0_copies_helix") ) {
				core::Size copyhelix( (*tag_it)->getOption<core::Size>("r0_copies_helix", 0) );
				runtime_assert_string_msg( copyhelix>0, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was set to 0.  Please specify a sensible helix index for the target that is to be copied." );
				runtime_assert_string_msg( copyhelix!=helix_index, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was the same as the copy helix index.  Please specify a sensible helix index for the target that is to be copied." );
				r0(this_helix)->set_helix_to_copy( copyhelix );
				TR << "Set helix " << helix_index << "'s r0 value to copy that of helix " << copyhelix << std::endl;
				has_perturbable_dofs=true;
			} else {
				if( (*tag_it)->hasOption("r0_perturbation") ) {
					core::Real r0pert( (*tag_it)->getOption<core::Real>("r0_perturbation", 0.0) );
					r0(this_helix)->set_perturbation_magnitude(r0pert);
					if(TR.visible()) TR << "Set r0 perturbation magnitude to " << r0pert << std::endl;
					has_perturbable_dofs=true;
				}
				if( (*tag_it)->hasOption("r0_perturbation_type") ) {
					if(r0(this_helix)->is_perturbable()) {
						std::string perttype( (*tag_it)->getOption<std::string>("r0_perturbation_type", "") );
						runtime_assert_string_msg( perttype=="gaussian" || perttype=="uniform",
							"In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: allowed perturbation types are \"gaussian\" and \"uniform\"." );
						if(perttype=="gaussian") r0(this_helix)->set_perturbation_type(pt_gaussian);
						else if (perttype=="uniform") r0(this_helix)->set_perturbation_type(pt_uniform);
						if(TR.visible()) TR << "Set r0 perturbation type for helix " << helix_index << "to " << perttype << "." << std::endl;
					} else {
						if(TR.Warning.visible())
							TR.Warning << "Warning!  The r0_perturbation_type option was specified for helix " << helix_index << ", but without an r0_perturbation option to set perturbation magnitude, it will be ignored." << std::endl;
					}
				}
			}

			// Note that omega0 has additional code in it for the special case of copying the pitch angle instead of the omega0 value.
			if( (*tag_it)->hasOption("pitch_from_helix") ) {
				runtime_assert_string_msg(
					!(*tag_it)->hasOption("omega0_copies_helix") &&
					!(*tag_it)->hasOption("omega0_perturbation"),
					"When parsing options for the BundleGridSampler mover, found \"pitch_from_helix\" alongside omega0 options.  This does not make sense.  EITHER a helix copies its pitch angle from another, OR the omega0 value can be perturbed/copied."
				);
				core::Size const val( (*tag_it)->getOption<core::Size>("pitch_from_helix", 0) );
				runtime_assert_string_msg( val>0, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was set to 0.  Please specify a sensible helix index for the target that is to be copied." );
				runtime_assert_string_msg( val!=helix_index, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was the same as the copy helix index.  Please specify a sensible helix index for the target that is to be copied." );
				if(TR.visible()) TR << "Setting omega0 for helix " << this_helix << " to be set to match the pitch angle for helix " << val << "." << std::endl;
				omega0(this_helix)->set_helix_to_copy( val );
				omega0(this_helix)->set_omega0_copies_pitch_instead(true); //We're going to copy the pitch angle instead of the omega0 value.
			} else { //All that follows resembles the code for the other parameters.
				if( (*tag_it)->hasOption("omega0_copies_helix") ) {
					core::Size copyhelix( (*tag_it)->getOption<core::Size>("omega0_copies_helix", 0) );
					runtime_assert_string_msg( copyhelix>0, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was set to 0.  Please specify a sensible helix index for the target that is to be copied." );
					runtime_assert_string_msg( copyhelix!=helix_index, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was the same as the copy helix index.  Please specify a sensible helix index for the target that is to be copied." );
					omega0(this_helix)->set_helix_to_copy( copyhelix );
					TR << "Set helix " << helix_index << "'s omega0 value to copy that of helix " << copyhelix << std::endl;
					has_perturbable_dofs=true;
				} else {
					if( (*tag_it)->hasOption("omega0_perturbation") ) {
						core::Real omega0pert( (*tag_it)->getOption<core::Real>("omega0_perturbation", 0.0) );
						omega0(this_helix)->set_perturbation_magnitude( convert_angle( omega0pert ) );
						if(TR.visible()) TR << "Set omega0 perturbation magnitude to " << omega0pert << (use_degrees() ? " degrees." : " radians.") << std::endl;
						has_perturbable_dofs=true;
					}
					if( (*tag_it)->hasOption("omega0_perturbation_type") ) {
						if(omega0(this_helix)->is_perturbable()) {
							std::string perttype( (*tag_it)->getOption<std::string>("omega0_perturbation_type", "") );
							runtime_assert_string_msg( perttype=="gaussian" || perttype=="uniform",
								"In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: allowed perturbation types are \"gaussian\" and \"uniform\"." );
							if(perttype=="gaussian") omega0(this_helix)->set_perturbation_type(pt_gaussian);
							else if (perttype=="uniform") omega0(this_helix)->set_perturbation_type(pt_uniform);
							if(TR.visible()) TR << "Set omega0 perturbation type for helix " << helix_index << "to " << perttype << "." << std::endl;
						} else {
							if(TR.Warning.visible())
								TR.Warning << "Warning!  The omega0_perturbation_type option was specified for helix " << helix_index << ", but without an omega0_perturbation option to set perturbation magnitude, it will be ignored." << std::endl;
						}
					}
				}
			}

			if( (*tag_it)->hasOption("delta_omega0_copies_helix") ) {
				core::Size copyhelix( (*tag_it)->getOption<core::Size>("delta_omega0_copies_helix", 0) );
				runtime_assert_string_msg( copyhelix>0, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was set to 0.  Please specify a sensible helix index for the target that is to be copied." );
				runtime_assert_string_msg( copyhelix!=helix_index, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was the same as the copy helix index.  Please specify a sensible helix index for the target that is to be copied." );
				delta_omega0(this_helix)->set_helix_to_copy( copyhelix );
				TR << "Set helix " << helix_index << "'s delta_omega0 value to copy that of helix " << copyhelix << std::endl;
				has_perturbable_dofs=true;
			} else {
				if( (*tag_it)->hasOption("delta_omega0_perturbation") ) {
					core::Real delta_omega0pert( (*tag_it)->getOption<core::Real>("delta_omega0_perturbation", 0.0) );
					delta_omega0(this_helix)->set_perturbation_magnitude( convert_angle( delta_omega0pert ) );
					if(TR.visible()) TR << "Set delta_omega0 perturbation magnitude to " << delta_omega0pert << (use_degrees() ? " degrees." : " radians.") << std::endl;
					has_perturbable_dofs=true;
				}
				if( (*tag_it)->hasOption("delta_omega0_perturbation_type") ) {
					if(delta_omega0(this_helix)->is_perturbable()) {
						std::string perttype( (*tag_it)->getOption<std::string>("delta_omega0_perturbation_type", "") );
						runtime_assert_string_msg( perttype=="gaussian" || perttype=="uniform",
							"In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: allowed perturbation types are \"gaussian\" and \"uniform\"." );
						if(perttype=="gaussian") delta_omega0(this_helix)->set_perturbation_type(pt_gaussian);
						else if (perttype=="uniform") delta_omega0(this_helix)->set_perturbation_type(pt_uniform);
						if(TR.visible()) TR << "Set delta_omega0 perturbation type for helix " << helix_index << "to " << perttype << "." << std::endl;
					} else {
						if(TR.Warning.visible())
							TR.Warning << "Warning!  The delta_omega0_perturbation_type option was specified for helix " << helix_index << ", but without an delta_omega0_perturbation option to set perturbation magnitude, it will be ignored." << std::endl;
					}
				}
			}

			if( (*tag_it)->hasOption("delta_omega1_copies_helix") ) {
				core::Size copyhelix( (*tag_it)->getOption<core::Size>("delta_omega1_copies_helix", 0) );
				runtime_assert_string_msg( copyhelix>0, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was set to 0.  Please specify a sensible helix index for the target that is to be copied." );
				runtime_assert_string_msg( copyhelix!=helix_index, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was the same as the copy helix index.  Please specify a sensible helix index for the target that is to be copied." );
				delta_omega1(this_helix)->set_helix_to_copy( copyhelix );
				TR << "Set helix " << helix_index << "'s delta_omega1 value to copy that of helix " << copyhelix << std::endl;
				has_perturbable_dofs=true;
			} else {
				if( (*tag_it)->hasOption("delta_omega1_perturbation") ) {
					core::Real delta_omega1pert( (*tag_it)->getOption<core::Real>("delta_omega1_perturbation", 0.0) );
					delta_omega1(this_helix)->set_perturbation_magnitude( convert_angle( delta_omega1pert) );
					if(TR.visible()) TR << "Set delta_omega1 perturbation magnitude to " << delta_omega1pert << (use_degrees() ? " degrees." : " radians.") << std::endl;
					has_perturbable_dofs=true;
				}
				if( (*tag_it)->hasOption("delta_omega1_perturbation_type") ) {
					if(delta_omega1(this_helix)->is_perturbable()) {
						std::string perttype( (*tag_it)->getOption<std::string>("delta_omega1_perturbation_type", "") );
						runtime_assert_string_msg( perttype=="gaussian" || perttype=="uniform",
							"In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: allowed perturbation types are \"gaussian\" and \"uniform\"." );
						if(perttype=="gaussian") delta_omega1(this_helix)->set_perturbation_type(pt_gaussian);
						else if (perttype=="uniform") delta_omega1(this_helix)->set_perturbation_type(pt_uniform);
						if(TR.visible()) TR << "Set delta_omega1 perturbation type for helix " << helix_index << "to " << perttype << "." << std::endl;
					} else {
						if(TR.Warning.visible())
							TR.Warning << "Warning!  The delta_omega1_perturbation_type option was specified for helix " << helix_index << ", but without an delta_omega1_perturbation option to set perturbation magnitude, it will be ignored." << std::endl;
					}
				}
			}

			if( (*tag_it)->hasOption("delta_t_copies_helix") ) {
				core::Size copyhelix( (*tag_it)->getOption<core::Size>("delta_t_copies_helix", 0) );
				runtime_assert_string_msg( copyhelix>0, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was set to 0.  Please specify a sensible helix index for the target that is to be copied." );
				runtime_assert_string_msg( copyhelix!=helix_index, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was the same as the copy helix index.  Please specify a sensible helix index for the target that is to be copied." );
				delta_t(this_helix)->set_helix_to_copy( copyhelix );
				TR << "Set helix " << helix_index << "'s delta_t value to copy that of helix " << copyhelix << std::endl;
				has_perturbable_dofs=true;
			} else {
				if( (*tag_it)->hasOption("delta_t_perturbation") ) {
					core::Real delta_tpert( (*tag_it)->getOption<core::Real>("delta_t_perturbation", 0.0) );
					delta_t(this_helix)->set_perturbation_magnitude(delta_tpert);
					if(TR.visible()) TR << "Set delta_t perturbation magnitude to " << delta_tpert << std::endl;
					has_perturbable_dofs=true;
				}
				if( (*tag_it)->hasOption("delta_t_perturbation_type") ) {
					if(delta_t(this_helix)->is_perturbable()) {
						std::string perttype( (*tag_it)->getOption<std::string>("delta_t_perturbation_type", "") );
						runtime_assert_string_msg( perttype=="gaussian" || perttype=="uniform",
							"In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: allowed perturbation types are \"gaussian\" and \"uniform\"." );
						if(perttype=="gaussian") delta_t(this_helix)->set_perturbation_type(pt_gaussian);
						else if (perttype=="uniform") delta_t(this_helix)->set_perturbation_type(pt_uniform);
						if(TR.visible()) TR << "Set delta_t perturbation type for helix " << helix_index << "to " << perttype << "." << std::endl;
					} else {
						if(TR.Warning.visible())
							TR.Warning << "Warning!  The delta_t_perturbation_type option was specified for helix " << helix_index << ", but without an delta_t_perturbation option to set perturbation magnitude, it will be ignored." << std::endl;
					}
				}
			}
			
			if( (*tag_it)->hasOption("z1_offset_copies_helix") ) {
				core::Size copyhelix( (*tag_it)->getOption<core::Size>("z1_offset_copies_helix", 0) );
				runtime_assert_string_msg( copyhelix>0, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was set to 0.  Please specify a sensible helix index for the target that is to be copied." );
				runtime_assert_string_msg( copyhelix!=helix_index, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was the same as the copy helix index.  Please specify a sensible helix index for the target that is to be copied." );
				z1_offset(this_helix)->set_helix_to_copy( copyhelix );
				TR << "Set helix " << helix_index << "'s z1_offset value to copy that of helix " << copyhelix << std::endl;
				has_perturbable_dofs=true;
			} else {
				if( (*tag_it)->hasOption("z1_offset_perturbation") ) {
					core::Real z1_offsetpert( (*tag_it)->getOption<core::Real>("z1_offset_perturbation", 0.0) );
					z1_offset(this_helix)->set_perturbation_magnitude(z1_offsetpert);
					if(TR.visible()) TR << "Set z1_offset perturbation magnitude to " << z1_offsetpert << std::endl;
					has_perturbable_dofs=true;
				}
				if( (*tag_it)->hasOption("z1_offset_perturbation_type") ) {
					if(z1_offset(this_helix)->is_perturbable()) {
						std::string perttype( (*tag_it)->getOption<std::string>("z1_offset_perturbation_type", "") );
						runtime_assert_string_msg( perttype=="gaussian" || perttype=="uniform",
							"In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: allowed perturbation types are \"gaussian\" and \"uniform\"." );
						if(perttype=="gaussian") z1_offset(this_helix)->set_perturbation_type(pt_gaussian);
						else if (perttype=="uniform") z1_offset(this_helix)->set_perturbation_type(pt_uniform);
						if(TR.visible()) TR << "Set z1_offset perturbation type for helix " << helix_index << "to " << perttype << "." << std::endl;
					} else {
						if(TR.Warning.visible())
							TR.Warning << "Warning!  The z1_offset_perturbation_type option was specified for helix " << helix_index << ", but without an z1_offset_perturbation option to set perturbation magnitude, it will be ignored." << std::endl;
					}
				}
			}
			
			if( (*tag_it)->hasOption("z0_offset_copies_helix") ) {
				core::Size copyhelix( (*tag_it)->getOption<core::Size>("z0_offset_copies_helix", 0) );
				runtime_assert_string_msg( copyhelix>0, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was set to 0.  Please specify a sensible helix index for the target that is to be copied." );
				runtime_assert_string_msg( copyhelix!=helix_index, "In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: a helix was set to copy another, but the target index was the same as the copy helix index.  Please specify a sensible helix index for the target that is to be copied." );
				z0_offset(this_helix)->set_helix_to_copy( copyhelix );
				TR << "Set helix " << helix_index << "'s z0_offset value to copy that of helix " << copyhelix << std::endl;
				has_perturbable_dofs=true;
			} else {
				if( (*tag_it)->hasOption("z0_offset_perturbation") ) {
					core::Real z0_offsetpert( (*tag_it)->getOption<core::Real>("z0_offset_perturbation", 0.0) );
					z0_offset(this_helix)->set_perturbation_magnitude(z0_offsetpert);
					if(TR.visible()) TR << "Set z0_offset perturbation magnitude to " << z0_offsetpert << std::endl;
					has_perturbable_dofs=true;
				}
				if( (*tag_it)->hasOption("z0_offset_perturbation_type") ) {
					if(z0_offset(this_helix)->is_perturbable()) {
						std::string perttype( (*tag_it)->getOption<std::string>("z0_offset_perturbation_type", "") );
						runtime_assert_string_msg( perttype=="gaussian" || perttype=="uniform",
							"In protocols::helical_bundle::PerturbBundle::parse_my_tag() function: allowed perturbation types are \"gaussian\" and \"uniform\"." );
						if(perttype=="gaussian") z0_offset(this_helix)->set_perturbation_type(pt_gaussian);
						else if (perttype=="uniform") z0_offset(this_helix)->set_perturbation_type(pt_uniform);
						if(TR.visible()) TR << "Set z0_offset perturbation type for helix " << helix_index << "to " << perttype << "." << std::endl;
					} else {
						if(TR.Warning.visible())
							TR.Warning << "Warning!  The z0_offset_perturbation_type option was specified for helix " << helix_index << ", but without an z0_offset_perturbation option to set perturbation magnitude, it will be ignored." << std::endl;
					}
				}
			}


			if(!has_perturbable_dofs && TR.visible())
				TR << "No perturbable parameters have been defined for helix " << helix_index << ".  This helix will remain fixed while others are perturbed." << std::endl;

		}
	}

	return;
} //parse_my_tag

////////////////////////////////////////////////////////////////////////////////
//          PUBLIC FUNCTIONS                                                  //
////////////////////////////////////////////////////////////////////////////////


/// @brief Add options for a new helix
///
core::Size PerturbBundle::add_helix( core::Size const helix_index ) {
	runtime_assert_string_msg(helix_not_defined( helix_index ), "In protocols::helical_bundle::PerturbBundle::add_helix() function: could not add helix.  Helix defined multiple times!" );

	r0_.push_back( PerturbBundleOptionsOP( new PerturbBundleOptions ) );
	r0(r0_.size())->set_helix_index(helix_index);
	omega0_.push_back( PerturbBundleOptionsOP( new PerturbBundleOptions ) );
	omega0(omega0_.size())->set_helix_index(helix_index);
	delta_omega0_.push_back( PerturbBundleOptionsOP( new PerturbBundleOptions ) );
	delta_omega0(delta_omega0_.size())->set_helix_index(helix_index);
	delta_omega1_.push_back( PerturbBundleOptionsOP( new PerturbBundleOptions ) );
	delta_omega1(delta_omega1_.size())->set_helix_index(helix_index);
	delta_t_.push_back( PerturbBundleOptionsOP( new PerturbBundleOptions ) );
	delta_t(delta_t_.size())->set_helix_index(helix_index);
	z1_offset_.push_back( PerturbBundleOptionsOP( new PerturbBundleOptions ) );
	z1_offset(z1_offset_.size())->set_helix_index(helix_index);
	z0_offset_.push_back( PerturbBundleOptionsOP( new PerturbBundleOptions ) );
	z0_offset(z0_offset_.size())->set_helix_index(helix_index);

	core::Size const nhelices(r0_.size());

	runtime_assert_string_msg( omega0_.size()==nhelices && delta_omega0_.size()==nhelices && delta_omega1_.size()==nhelices && delta_t_.size()==nhelices && z1_offset_.size()==nhelices && z0_offset_.size()==nhelices,
		"In protocols::helical_bundle::PerturbBundle::add_helix() function: somehow, vector indices are out of sync.  I can't determine how many helices have been defined." );

	return nhelices;
}

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////

/// @brief Is a value in a list?
///
bool PerturbBundle::is_in_list( core::Size const val, utility::vector1 < core::Size> const &list ) const {
	core::Size const listsize(list.size());
	if(listsize==0) return false;
	for(core::Size i=1; i<=listsize; ++i) {
		if(list[i]==val) return true;
	}
	return false;
}


/// @brief Confirms that a helix has not yet been defined.  Returns "true" if the helix
/// has NOT been defined, false otherwise.
bool PerturbBundle::helix_not_defined( core::Size const helix_index) const {
	core::Size r0size=r0_.size();
	core::Size omega0size=omega0_.size();
	core::Size delta_omega0size=delta_omega0_.size();
	core::Size delta_omega1size=delta_omega1_.size();
	core::Size delta_tsize=delta_t_.size();
	core::Size z1_offsetsize=z1_offset_.size();
	core::Size z0_offsetsize=z0_offset_.size();

	if(r0size>0) {
		for(core::Size i=1; i<=r0size; ++i) if(r0(i)->helix_index() == helix_index) return false;
	}
	if(omega0size>0) {
		for(core::Size i=1; i<=omega0size; ++i) if(omega0(i)->helix_index() == helix_index) return false;
	}
	if(delta_omega0size>0) {
		for(core::Size i=1; i<=delta_omega0size; ++i) if(delta_omega0(i)->helix_index() == helix_index) return false;
	}
	if(delta_omega1size>0) {
		for(core::Size i=1; i<=delta_omega1size; ++i) if(delta_omega1(i)->helix_index() == helix_index) return false;
	}
	if(delta_tsize>0) {
		for(core::Size i=1; i<=delta_tsize; ++i) if(delta_t(i)->helix_index() == helix_index) return false;
	}
	if(z1_offsetsize>0) {
		for(core::Size i=1; i<=z1_offsetsize; ++i) if(z1_offset(i)->helix_index() == helix_index) return false;
	}
	if(z0_offsetsize>0) {
		for(core::Size i=1; i<=z0_offsetsize; ++i) if(z0_offset(i)->helix_index() == helix_index) return false;
	}

	return true;
}

/// @brief Perturb the helical bundle parameter values in the pose, subject to the options already set.
/// @details Called by the apply() function.  Returns true for success, false for failure.
bool PerturbBundle::perturb_values( BundleParametersSetOP params_set) const {

	//Get the symmetry of the bundle:
	core::Size const symmetry( params_set->bundle_symmetry() < 2 ? 1 : params_set->bundle_symmetry() );
	//Get the number of symmetry copies in the bundle:
	core::Size const symmetry_copies( params_set->bundle_symmetry_copies()==0 ? symmetry : params_set->bundle_symmetry_copies() );
	//Get the number of helices defined in each symmetry copy in the bundle:
	core::Size const n_helices( params_set->n_helices() );
	runtime_assert_string_msg( n_helices > 0,
		"In protocols::helical_bundle::PerturbBundle::perturb_values() function: the number of helices in the pose is 0.  Unable to proceed -- nothing to perturb." );

	if(TR.Debug.visible()) TR.Debug << "params_set->n_parameters(): " << params_set->n_parameters() << std::endl; //DELETE ME
	if(TR.Debug.visible()) TR.Debug << "n_helices: " << n_helices << std::endl; //DELETE ME
	if(TR.Debug.visible()) TR.Debug << "symmetry_copies: " << symmetry_copies << std::endl; //DELETE ME
	runtime_assert_string_msg( params_set->n_parameters() == n_helices*symmetry_copies,
			"In protocols::helical_bundle::PerturbBundle::perturb_values() function: the pose has corrupted BundleParametersSet data.  Unable to proceed." );

	core::Size helix_index(0); //Counter for the index of the current helix.

	// Lists of helices that have already been processed:
	utility::vector1 < core::Size > helices_processed_r0;
	utility::vector1 < core::Size > helices_processed_omega0;
	utility::vector1 < core::Size > helices_processed_delta_omega0;
	utility::vector1 < core::Size > helices_processed_delta_omega1;
	utility::vector1 < core::Size > helices_processed_delta_t;
	utility::vector1 < core::Size > helices_processed_z1_offset;
	utility::vector1 < core::Size > helices_processed_z0_offset;

	bool loopthrough_failed(false); //Used to break from the nested loops.
	for(core::Size isym=1; isym<=symmetry_copies; ++isym) { //Loop through all of the symmetry copies.
		for(core::Size ihelix=1; ihelix<=n_helices; ++ihelix) { //Loop through all of the helices defined for this symmetry copy
			++helix_index; //Increment the index of the current helix.
			BundleParametersOP params( utility::pointer::dynamic_pointer_cast< BundleParameters >( params_set->parameters(helix_index) ) );
			runtime_assert_string_msg( params,
				"In protocols::helical_bundle::PerturbBundle::perturb_values() function: unable to get an owning pointer for the BundleParameters object." );

			if(helix_index <= n_helices) { //If this is part of the first symmetry repeat, perturb the values
				if(helix_not_defined(helix_index)) { //If custom perturbation properties have NOT been defined for this mover, apply the default perturbation.

					//For each of these, we need to check:
					// -is it perturbable?
					// -has it already been perturbed?
					//Then we perturb it and add it to the list of what's already been perturbed.
					if(default_r0()->is_perturbable() && !is_in_list(helix_index, helices_processed_r0) ) {
						params->set_r0( params->r0() + default_r0()->delta() );
					}
					if(!is_in_list(helix_index, helices_processed_r0)) helices_processed_r0.push_back(helix_index);

					if(default_omega0()->is_perturbable() && !is_in_list(helix_index, helices_processed_omega0) ) {
						params->set_omega0( params->omega0() + default_omega0()->delta() );
					}
					if(!is_in_list(helix_index, helices_processed_omega0)) helices_processed_omega0.push_back(helix_index);

					if(default_delta_omega0()->is_perturbable() && !is_in_list(helix_index, helices_processed_delta_omega0) ) {
						params->set_delta_omega0( params->delta_omega0() + default_delta_omega0()->delta() );
					}
					if(!is_in_list(helix_index, helices_processed_delta_omega0)) helices_processed_delta_omega0.push_back(helix_index);

					if(default_delta_omega1()->is_perturbable() && !is_in_list(helix_index, helices_processed_delta_omega1) ) {
						params->set_delta_omega1_all( params->delta_omega1_all() + default_delta_omega1()->delta() );
					}
					if(!is_in_list(helix_index, helices_processed_delta_omega1)) helices_processed_delta_omega1.push_back(helix_index);
					
					if(default_delta_t()->is_perturbable() && !is_in_list(helix_index, helices_processed_delta_t) )	{
						params->set_delta_t( params->delta_t() + default_delta_t()->delta() );
					}
					if(!is_in_list(helix_index, helices_processed_delta_t)) helices_processed_delta_t.push_back(helix_index);
	
					if(default_z1_offset()->is_perturbable() && !is_in_list(helix_index, helices_processed_z1_offset) )	{
						params->set_z1_offset( params->z1_offset() + default_z1_offset()->delta() );
					}
					if(!is_in_list(helix_index, helices_processed_z1_offset)) helices_processed_z1_offset.push_back(helix_index);

					if(default_z0_offset()->is_perturbable() && !is_in_list(helix_index, helices_processed_z0_offset) )	{
						params->set_z0_offset( params->z0_offset() + default_z0_offset()->delta() );
					}
					if(!is_in_list(helix_index, helices_processed_z0_offset)) helices_processed_z0_offset.push_back(helix_index);

					if(TR.visible()) TR << "Completed default perturbation of helix " << helix_index << "." << std::endl;

				} else { //If custom perturbation properties HAVE been defined for this helix, apply the custom perturbation.

					//Determine the proper index in the r0_, omega_, etc. vectors for this helix:
					core::Size custom_properties_index(0); //The index in the r0_, omega0_, etc. vectors for the custom properties of this helix.
					for(core::Size i=1, imax=r0_.size(); i<=imax; ++i) {
						if(r0(i)->helix_index() == helix_index) {
							custom_properties_index=i;
							break;
						}
					}
					runtime_assert_string_msg( custom_properties_index>0, "Internal error in protocols::helical_bundle::PerturbBundle::perturb_values() function.  Consult a developer or an exorcist -- this shouldn't happen." );

					//Perturbing r0:
					//Case 1: this helix perturbation is from defaults:
					if(r0(custom_properties_index)->use_defaults()) {
						if(default_r0()->is_perturbable() && !is_in_list(helix_index, helices_processed_r0) ) {
							params->set_r0( params->r0() + default_r0()->delta() );
						}
						if(!is_in_list(helix_index, helices_processed_r0)) helices_processed_r0.push_back(helix_index);
					//Case 2: this helix perturbation is based on another helix:
					} else if (r0(custom_properties_index)->is_copy()) {
						//Find the index of the helix that this is a copy of:
						core::Size master_index( r0( custom_properties_index )->other_helix() );
						//For now, require that this helix already be processed:
						runtime_assert_string_msg( is_in_list(master_index, helices_processed_r0 ),
							"In protocols::helical_bundle::PerturbBundle::perturb_values() function: A helix was set to copy another, but the helix that it copies has a higher index. (Helices can only copy lower-index helices, currently.  This might change in a future release.  Note that this message might also result if the helix DoF being copied is never perturbed.)" );
						//Get an owning pointer to the params of the master helix that this one is copying:
						BundleParametersOP master_params( utility::pointer::dynamic_pointer_cast<BundleParameters>( params_set->parameters(master_index) ) );
						params->set_r0( master_params->r0() ); //Copy the r0 parameter.		
						helices_processed_r0.push_back( helix_index );	
					//Case 3: this helix perturbation is independent:
					} else if (r0(custom_properties_index)->is_perturbable()  && !is_in_list(helix_index, helices_processed_r0) ){
						params->set_r0( params->r0() + r0(custom_properties_index)->delta() );
						helices_processed_r0.push_back( helix_index );	
					} else {
						helices_processed_r0.push_back( helix_index );	
					}

					//Perturbing omega0:
					//Case 1: this helix perturbation is from defaults:
					if(omega0(custom_properties_index)->use_defaults()) {
						if(default_omega0()->is_perturbable() && !is_in_list(helix_index, helices_processed_omega0) ) {
							params->set_omega0( params->omega0() + default_omega0()->delta() );
						}
						if(!is_in_list(helix_index, helices_processed_omega0)) helices_processed_omega0.push_back(helix_index);
					//Case 2: this helix perturbation is based on another helix:
					} else if (omega0(custom_properties_index)->is_copy()) {
						//Find the index of the helix that this is a copy of:
						core::Size master_index( omega0( custom_properties_index )->other_helix() );
						//For now, require that this helix already be processed:
						runtime_assert_string_msg( is_in_list(master_index, helices_processed_omega0 ),
							"In protocols::helical_bundle::PerturbBundle::perturb_values() function: A helix was set to copy another, but the helix that it copies has a higher index. (Helices can only copy lower-index helices, currently.  This might change in a future release.  Note that this message might also result if the helix DoF being copied is never perturbed.)" );
						//Get an owning pointer to the params of the master helix that this one is copying:
						BundleParametersOP master_params( utility::pointer::dynamic_pointer_cast<BundleParameters>( params_set->parameters(master_index) ) );

						//Note that there's some extra code here for the special case of copying helical pitch (rise per turn) rather than omega0 (turn per residue).
						if(omega0(custom_properties_index)->omega0_copies_pitch_instead()) {

							core::Real const other_r0( master_params->r0() );
							core::Real const other_omega0( master_params->omega0() );
							core::Real const other_z1( master_params->z1() );
							core::Real const other_sinalpha( other_r0*other_omega0/other_z1 );
							if(other_sinalpha > 1 || other_sinalpha < -1) {
								if(TR.visible()) TR << "Failed to copy pitch angle.  Current parameters do not generate a sensible pitch angle for helix " << master_index << "." << std::endl;
								loopthrough_failed=true;
								break; //Stop looping through the helices.
							}
							//If we've got a good pitch angle, then continue:
							core::Real const other_alpha( asin(other_sinalpha) );
							
							core::Real const this_r0( params->r0() ); //Already set above, if sampled or if copied.
							core::Real const this_z1( params->z1() ); //Cannot be sampled or copied.
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
							params->set_omega0( this_z1/this_r0 * sin(atan(this_r0/other_r0*tan(other_alpha))) );

						} else {
							params->set_omega0( master_params->omega0() ); //Copy the omega0 parameter.
						}
						
						helices_processed_omega0.push_back( helix_index );	
					//Case 3: this helix perturbation is independent:
					} else if (omega0(custom_properties_index)->is_perturbable()  && !is_in_list(helix_index, helices_processed_omega0) ){
						params->set_omega0( params->omega0() + omega0(custom_properties_index)->delta() );
						helices_processed_omega0.push_back( helix_index );	
					} else {
						helices_processed_omega0.push_back( helix_index );	
					}

					//Perturbing delta_omega0:
					//Case 1: this helix perturbation is from defaults:
					if(delta_omega0(custom_properties_index)->use_defaults()) {
						if(default_delta_omega0()->is_perturbable() && !is_in_list(helix_index, helices_processed_delta_omega0) ) {
							params->set_delta_omega0( params->delta_omega0() + default_delta_omega0()->delta() );
						}
						if(!is_in_list(helix_index, helices_processed_delta_omega0)) helices_processed_delta_omega0.push_back(helix_index);
					//Case 2: this helix perturbation is based on another helix:
					} else if (delta_omega0(custom_properties_index)->is_copy()) {
						//Find the index of the helix that this is a copy of:
						core::Size master_index( delta_omega0( custom_properties_index )->other_helix() );
						//For now, require that this helix already be processed:
						runtime_assert_string_msg( is_in_list(master_index, helices_processed_delta_omega0 ),
							"In protocols::helical_bundle::PerturbBundle::perturb_values() function: A helix was set to copy another, but the helix that it copies has a higher index. (Helices can only copy lower-index helices, currently.  This might change in a future release.  Note that this message might also result if the helix DoF being copied is never perturbed.)" );
						//Get an owning pointer to the params of the master helix that this one is copying:
						BundleParametersOP master_params( utility::pointer::dynamic_pointer_cast<BundleParameters>( params_set->parameters(master_index) ) );
						params->set_delta_omega0( master_params->delta_omega0() ); //Copy the delta_omega0 parameter.		
						helices_processed_delta_omega0.push_back( helix_index );	
					//Case 3: this helix perturbation is independent:
					} else if (delta_omega0(custom_properties_index)->is_perturbable()  && !is_in_list(helix_index, helices_processed_delta_omega0) ){
						params->set_delta_omega0( params->delta_omega0() + delta_omega0(custom_properties_index)->delta() );
						helices_processed_delta_omega0.push_back( helix_index );	
					} else {
						helices_processed_delta_omega0.push_back( helix_index );	
					}

					//Perturbing delta_omega1:
					//Case 1: this helix perturbation is from defaults:
					if(delta_omega1(custom_properties_index)->use_defaults()) {
						if(default_delta_omega1()->is_perturbable() && !is_in_list(helix_index, helices_processed_delta_omega1) ) {
							params->set_delta_omega1_all( params->delta_omega1_all() + default_delta_omega1()->delta() );
						}
						if(!is_in_list(helix_index, helices_processed_delta_omega1)) helices_processed_delta_omega1.push_back(helix_index);
					//Case 2: this helix perturbation is based on another helix:
					} else if (delta_omega1(custom_properties_index)->is_copy()) {
						//Find the index of the helix that this is a copy of:
						core::Size master_index( delta_omega1( custom_properties_index )->other_helix() );
						//For now, require that this helix already be processed:
						runtime_assert_string_msg( is_in_list(master_index, helices_processed_delta_omega1 ),
							"In protocols::helical_bundle::PerturbBundle::perturb_values() function: A helix was set to copy another, but the helix that it copies has a higher index. (Helices can only copy lower-index helices, currently.  This might change in a future release.  Note that this message might also result if the helix DoF being copied is never perturbed.)" );
						//Get an owning pointer to the params of the master helix that this one is copying:
						BundleParametersOP master_params( utility::pointer::dynamic_pointer_cast<BundleParameters>( params_set->parameters(master_index) ) );
						params->set_delta_omega1_all( master_params->delta_omega1_all() ); //Copy the delta_omega1 parameter.		
						helices_processed_delta_omega1.push_back( helix_index );	
					//Case 3: this helix perturbation is independent:
					} else if (delta_omega1(custom_properties_index)->is_perturbable()  && !is_in_list(helix_index, helices_processed_delta_omega1) ){
						params->set_delta_omega1_all( params->delta_omega1_all() + delta_omega1(custom_properties_index)->delta() );
						helices_processed_delta_omega1.push_back( helix_index );	
					} else {
						helices_processed_delta_omega1.push_back( helix_index );	
					}

					//Perturbing delta_t:
					//Case 1: this helix perturbation is from defaults:
					if(delta_t(custom_properties_index)->use_defaults()) {
						if(default_delta_t()->is_perturbable() && !is_in_list(helix_index, helices_processed_delta_t) ) {
							params->set_delta_t( params->delta_t() + default_delta_t()->delta() );
						}
						if(!is_in_list(helix_index, helices_processed_delta_t)) helices_processed_delta_t.push_back(helix_index);
					//Case 2: this helix perturbation is based on another helix:
					} else if (delta_t(custom_properties_index)->is_copy()) {
						//Find the index of the helix that this is a copy of:
						core::Size master_index( delta_t( custom_properties_index )->other_helix() );
						//For now, require that this helix already be processed:
						runtime_assert_string_msg( is_in_list(master_index, helices_processed_delta_t ),
							"In protocols::helical_bundle::PerturbBundle::perturb_values() function: A helix was set to copy another, but the helix that it copies has a higher index. (Helices can only copy lower-index helices, currently.  This might change in a future release.  Note that this message might also result if the helix DoF being copied is never perturbed.)" );
						//Get an owning pointer to the params of the master helix that this one is copying:
						BundleParametersOP master_params( utility::pointer::dynamic_pointer_cast<BundleParameters>( params_set->parameters(master_index) ) );
						params->set_delta_t( master_params->delta_t() ); //Copy the delta_t parameter.		
						helices_processed_delta_t.push_back( helix_index );	
					//Case 3: this helix perturbation is independent:
					} else if (delta_t(custom_properties_index)->is_perturbable()  && !is_in_list(helix_index, helices_processed_delta_t) ){
						params->set_delta_t( params->delta_t() + delta_t(custom_properties_index)->delta() );
						helices_processed_delta_t.push_back( helix_index );	
					} else {
						helices_processed_delta_t.push_back( helix_index );	
					}

					//Perturbing z1_offset:
					//Case 1: this helix perturbation is from defaults:
					if(z1_offset(custom_properties_index)->use_defaults()) {
						if(default_z1_offset()->is_perturbable() && !is_in_list(helix_index, helices_processed_z1_offset) ) {
							params->set_z1_offset( params->z1_offset() + default_z1_offset()->delta() );
						}
						if(!is_in_list(helix_index, helices_processed_z1_offset)) helices_processed_z1_offset.push_back(helix_index);
					//Case 2: this helix perturbation is based on another helix:
					} else if (z1_offset(custom_properties_index)->is_copy()) {
						//Find the index of the helix that this is a copy of:
						core::Size master_index( z1_offset( custom_properties_index )->other_helix() );
						//For now, require that this helix already be processed:
						runtime_assert_string_msg( is_in_list(master_index, helices_processed_z1_offset ),
							"In protocols::helical_bundle::PerturbBundle::perturb_values() function: A helix was set to copy another, but the helix that it copies has a higher index. (Helices can only copy lower-index helices, currently.  This might change in a future release.  Note that this message might also result if the helix DoF being copied is never perturbed.)" );
						//Get an owning pointer to the params of the master helix that this one is copying:
						BundleParametersOP master_params( utility::pointer::dynamic_pointer_cast<BundleParameters>( params_set->parameters(master_index) ) );
						params->set_z1_offset( master_params->z1_offset() ); //Copy the z1_offset parameter.		
						helices_processed_z1_offset.push_back( helix_index );	
					//Case 3: this helix perturbation is independent:
					} else if (z1_offset(custom_properties_index)->is_perturbable()  && !is_in_list(helix_index, helices_processed_z1_offset) ){
						params->set_z1_offset( params->z1_offset() + z1_offset(custom_properties_index)->delta() );
						helices_processed_z1_offset.push_back( helix_index );	
					} else {
						helices_processed_z1_offset.push_back( helix_index );	
					}

					//Perturbing z0_offset:
					//Case 1: this helix perturbation is from defaults:
					if(z0_offset(custom_properties_index)->use_defaults()) {
						if(default_z0_offset()->is_perturbable() && !is_in_list(helix_index, helices_processed_z0_offset) ) {
							params->set_z0_offset( params->z0_offset() + default_z0_offset()->delta() );
						}
						if(!is_in_list(helix_index, helices_processed_z0_offset)) helices_processed_z0_offset.push_back(helix_index);
					//Case 2: this helix perturbation is based on another helix:
					} else if (z0_offset(custom_properties_index)->is_copy()) {
						//Find the index of the helix that this is a copy of:
						core::Size master_index( z0_offset( custom_properties_index )->other_helix() );
						//For now, require that this helix already be processed:
						runtime_assert_string_msg( is_in_list(master_index, helices_processed_z0_offset ),
							"In protocols::helical_bundle::PerturbBundle::perturb_values() function: A helix was set to copy another, but the helix that it copies has a higher index. (Helices can only copy lower-index helices, currently.  This might change in a future release.  Note that this message might also result if the helix DoF being copied is never perturbed.)" );
						//Get an owning pointer to the params of the master helix that this one is copying:
						BundleParametersOP master_params( utility::pointer::dynamic_pointer_cast<BundleParameters>( params_set->parameters(master_index) ) );
						params->set_z0_offset( master_params->z0_offset() ); //Copy the z0_offset parameter.		
						helices_processed_z0_offset.push_back( helix_index );	
					//Case 3: this helix perturbation is independent:
					} else if (z0_offset(custom_properties_index)->is_perturbable()  && !is_in_list(helix_index, helices_processed_z0_offset) ){
						params->set_z0_offset( params->z0_offset() + z0_offset(custom_properties_index)->delta() );
						helices_processed_z0_offset.push_back( helix_index );	
					} else {
						helices_processed_z0_offset.push_back( helix_index );	
					}

				}
			} else { //If this is part of a later symmetry repeat, just copy the values from the first symmetry repeat.
				BundleParametersOP ref_params( utility::pointer::dynamic_pointer_cast< BundleParameters >( params_set->parameters(ihelix) ) );
				runtime_assert_string_msg( ref_params,
					"In protocols::helical_bundle::PerturbBundle::perturb_values() function: unable to get an owning pointer for the reference BundleParameters object.  This is odd -- it should not happen." );
				// Copy the r0, omega0, delta_omega0, delta_omega1, delta_t, z1_offset, and z0_offset parameters:
				params->set_r0( ref_params->r0() );
				params->set_omega0( ref_params->omega0() );
				params->set_delta_omega0( ref_params->delta_omega0() + static_cast<core::Real>(isym-1)/static_cast<core::Real>(symmetry)*numeric::constants::d::pi_2 );
				params->set_delta_omega1_all( ref_params->delta_omega1_all() );
				params->set_delta_t( ref_params->delta_t() );
				params->set_z1_offset( ref_params->z1_offset() );
				params->set_z0_offset( ref_params->z0_offset() );
			}

		}
		if(loopthrough_failed) break;
	}

	return (!loopthrough_failed);
}

/// @brief Rebuild the helical bundle conformation using the bundle parameter values in the pose.
///
void PerturbBundle::rebuild_conformation(
	core::pose::Pose &pose,
	BundleParametersSetOP params_set,
	core::Size const params_set_index,
	bool &failed
) const {

	core::Size const n_params( params_set->n_parameters() ); //Get the number of parameters objects (number of helices that we're rebuilding).

	for(core::Size ihelix=1; ihelix<=n_params; ++ihelix) {

		PerturbBundleHelixOP perthelix( new PerturbBundleHelix ); //Construct the mover to perturb this helix.

		perthelix->set_parameters_set_index(params_set_index);
		perthelix->set_parameters_index(ihelix);

		perthelix->apply( pose );

		failed = perthelix->last_apply_failed();
		if(failed) break;

	}

	return;
}

/// @brief Write out the perturbed Crick parameters.
/// @details The "before" parameter determines whether this is a pre-perturbation
/// report or a post-perturbation report.
void PerturbBundle::write_report(
	BundleParametersSetOP params_set,
	bool const before
) const {
	if(!TR.visible()) return; //Do nothing if the tracer isn't visible.

	if(before) {
		TR << "*** Crick parameter values prior to perturbation ***" << std::endl;
	} else {
		TR << "*** Crick parameter values after perturbation ***" << std::endl;
	}

	core::Size const nparams(params_set->n_parameters());
	if(nparams==0) {
		TR << "No Parameters objects in the ParametersSet object!" << std::endl;
	}

	for(core::Size iparams=1; iparams<=nparams; ++iparams) { //Loop through the params objects
		BundleParametersOP params( utility::pointer::dynamic_pointer_cast<BundleParameters>(params_set->parameters(iparams)) );
		if(!params) continue;
		TR << "Helix " << iparams << ":" << std::endl;
		TR << "     r0: " << params->r0() << std::endl;
		TR << "     omega0: " << params->omega0() << std::endl;
		TR << "     delta_omega0: " << params->delta_omega0() << std::endl;
		TR << "     delta_omega1: " << params->delta_omega1_all() << std::endl;
		TR << "     delta_t: " << params->delta_t() << std::endl;
		TR << "     z1_offset: " << params->z1_offset() << std::endl;
		TR << "     z0_offset: " << params->z0_offset() << std::endl;
	}

	TR << "*** END SUMMARY ***" << std::endl;

	TR.flush();

	return;
}


} //namespace helical_bundle
} //namespace protocols
