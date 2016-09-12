// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/MakeBundle.cc
/// @brief  Builds a helical bundle using the Crick parameters.
/// @details The bundle is centred on the origin, with the outer helix axis pointing along the
/// global z-axis.  This mover calls the MakeBundleHelix mover.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/helical_bundle/MakeBundle.hh>
#include <protocols/helical_bundle/MakeBundleCreator.hh>
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
#include <core/pose/PDBInfo.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace helical_bundle {

static THREAD_LOCAL basic::Tracer TR("protocols.helical_bundle.MakeBundle");

std::string
MakeBundleCreator::keyname() const
{
	return MakeBundleCreator::mover_name();
}

protocols::moves::MoverOP
MakeBundleCreator::create_mover() const {
	return protocols::moves::MoverOP( new MakeBundle );
}

std::string
MakeBundleCreator::mover_name()
{
	return "MakeBundle";
}

/// @brief Creator for MakeBundle mover.
MakeBundle::MakeBundle():
	Mover("MakeBundle"),
	reset_pose_(true),
	make_bundle_helix_movers_(),
	bundle_symmetry_(0),
	bundle_symmetry_copies_(0),
	default_r0_(0),
	default_r0_set_(false),
	default_omega0_(0),
	default_omega0_set_(false),
	default_delta_omega0_(0),
	default_delta_omega0_set_(false),
	default_crick_params_file_(""),
	default_crick_params_file_set_(false),
	default_omega1_(0),
	default_omega1_set_(false),
	default_z1_(0),
	default_z1_set_(false),
	default_delta_omega1_all_(0),
	default_delta_omega1_all_set_(false),
	default_residue_name_(),
	default_residue_name_set_(false),
	default_repeating_unit_offset_(0),
	default_repeating_unit_offset_set_(false),
	default_delta_t_(0),
	default_delta_t_set_(false),
	default_z1_offset_(0),
	default_z1_offset_set_(false),
	default_z0_offset_(0),
	default_z0_offset_set_(false),
	default_invert_(false),
	default_invert_set_(false),
	default_helix_length_(0),
	default_helix_length_set_(false),
	default_allow_bondlengths_(false),
	default_allow_bondlengths_set_(false),
	default_allow_bondangles_(false),
	default_allow_bondangles_set_(false),
	default_allow_dihedrals_(true),
	default_allow_dihedrals_set_(false),
	use_degrees_(false),
	last_apply_failed_(false)
{}


/// @brief Copy constructor for MakeBundle mover.
///
/// @brief Creator for MakeBundle mover.
MakeBundle::MakeBundle( MakeBundle const & src ):
	protocols::moves::Mover(src),
	reset_pose_(src.reset_pose_),
	make_bundle_helix_movers_(), //copy this below
	bundle_symmetry_(src.bundle_symmetry_),
	bundle_symmetry_copies_(src.bundle_symmetry_copies_),
	default_r0_(src.default_r0_),
	default_r0_set_(src.default_r0_set_),
	default_omega0_(src.default_omega0_),
	default_omega0_set_(src.default_omega0_set_),
	default_delta_omega0_(src.default_delta_omega0_),
	default_delta_omega0_set_(src.default_delta_omega0_set_),
	default_crick_params_file_(src.default_crick_params_file_),
	default_crick_params_file_set_(src.default_crick_params_file_set_),
	default_omega1_(src.default_omega1_),
	default_omega1_set_(src.default_omega1_set_),
	default_z1_(src.default_z1_),
	default_z1_set_(src.default_z1_set_),
	default_delta_omega1_all_(src.default_delta_omega1_all_),
	default_delta_omega1_all_set_(src.default_delta_omega1_all_set_),
	default_residue_name_(src.default_residue_name_),
	default_residue_name_set_(src.default_residue_name_set_),
	default_repeating_unit_offset_(src.default_repeating_unit_offset_),
	default_repeating_unit_offset_set_(src.default_repeating_unit_offset_set_),
	default_delta_t_(src.default_delta_t_),
	default_delta_t_set_(src.default_delta_t_set_),
	default_z1_offset_(src.default_z1_offset_),
	default_z1_offset_set_(src.default_z1_offset_set_),
	default_z0_offset_(src.default_z0_offset_),
	default_z0_offset_set_(src.default_z0_offset_set_),
	default_invert_(src.default_invert_),
	default_invert_set_(src.default_invert_set_),
	default_helix_length_(src.default_helix_length_),
	default_helix_length_set_(src.default_helix_length_set_),
	default_allow_bondlengths_(src.default_allow_bondlengths_),
	default_allow_bondlengths_set_(src.default_allow_bondlengths_set_),
	default_allow_bondangles_(src.default_allow_bondangles_),
	default_allow_bondangles_set_(src.default_allow_bondangles_set_),
	default_allow_dihedrals_(src.default_allow_dihedrals_),
	default_allow_dihedrals_set_(src.default_allow_dihedrals_set_),
	use_degrees_(src.use_degrees_),
	last_apply_failed_(src.last_apply_failed_)
{
	make_bundle_helix_movers_.clear();
	for ( core::Size i=1, imax=src.make_bundle_helix_movers_.size(); i<=imax; ++i ) {
		make_bundle_helix_movers_.push_back( utility::pointer::dynamic_pointer_cast< MakeBundleHelix >( src.make_bundle_helix_movers_[i]->clone() ) );
	}
}


/// @brief Destructor for MakeBundle mover.
MakeBundle::~MakeBundle() = default;


/// @brief Clone operator to create a pointer to a fresh MakeBundle object that copies this one.
protocols::moves::MoverOP MakeBundle::clone() const {
	return protocols::moves::MoverOP( new MakeBundle( *this ) );
}


/// @brief Fresh_instance operator to create a pointer to a fresh MakeBundle object that does NOT copy this one.
protocols::moves::MoverOP MakeBundle::fresh_instance() const {
	return protocols::moves::MoverOP( new MakeBundle );
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
void MakeBundle::apply (core::pose::Pose & pose)
{
	if ( TR.visible() ) TR << "Building a helical bundle using the Crick equations." << std::endl;

	bool failed=false;

	core::pose::Pose newpose;

	core::Size const total_repeats(symmetry()<2 ? 1 : symmetry());
	core::Size const repeats(symmetry_copies()< 1 ? total_repeats : symmetry_copies());
	core::Real const omega0_offset_increment(numeric::constants::d::pi_2 / static_cast<core::Real>(total_repeats));
	core::Real omega0_offset(0.0);

	//The BundleParametersSet object that will hold the Crick parameters and will be passed into the Conformation object of the pose on which we will be operating.
	BundleParametersSetOP output_parameters_set(BundleParametersSetOP( new BundleParametersSet ));
	output_parameters_set->set_bundle_symmetry( symmetry() ); //Store the symmetry of this bundle.
	output_parameters_set->set_bundle_symmetry_copies( symmetry_copies() ); //Store the number of symmetry copies to actually generate.
	output_parameters_set->set_n_helices( n_helices() ); //Store the number of helices that are defined for each symmetry repeat.

	for ( core::Size irepeat=1; irepeat<=repeats; ++irepeat ) { //Loop through all of the symmetry copies.
		for ( core::Size ihelix=1, ihelixmax=n_helices(); ihelix<=ihelixmax; ++ihelix ) { //Loop through all of the helices defined for each symmetry repeat.
			protocols::helical_bundle::MakeBundleHelixOP make_helix( utility::pointer::dynamic_pointer_cast<protocols::helical_bundle::MakeBundleHelix>(helix(ihelix)->clone()) );
			make_helix->set_reset_pose(true);
			if ( irepeat > 1 ) make_helix->set_delta_omega0( make_helix->delta_omega0() + omega0_offset ); //Offset omega0 for symmetry copies.
			core::pose::Pose helixpose; //A temporary pose.
			make_helix->apply(helixpose);

			failed=make_helix->last_apply_failed(); //Did the apply fail?
			if ( failed ) break;

			//Add the newly-generated helix to the pose for the bundle.
			if ( newpose.empty() ) newpose=helixpose;
			else newpose.append_pose_by_jump(helixpose,1); //This will also append the ParametersSet object -- but this creates a list of ParametersSets, rather than one ParametersSet with a list of Parameters.
		}
		if ( repeats>1 ) omega0_offset += omega0_offset_increment;
	}

	if ( failed ) {
		set_last_apply_failed(true);
		if ( TR.visible() ) TR << "Build of helical bundle failed, likely due to bad Crick parameters resulting in nonsensical geometry.  Returning input pose." << std::endl;
	} else {
		set_last_apply_failed(false);

		//At this point, newpose has a bunch of ParametersSet objects, each with one Parameters object.  We want a single ParametersSet object with all of the Parameters objects.
		//So we do some shuffling around:
		for ( core::Size i=1, imax=newpose.conformation().n_parameters_sets(); i<=imax; ++i ) { //Loop through all ParametersSets
			output_parameters_set->add_parameters( newpose.conformation().parameters_set(i)->parameters(1) ); //Put all of the Parameters objects into this ParametersSet object.
		}
		newpose.conformation().clear_parameters_set_list(); //Delete all of the ParametersSet objects in newpose.
		newpose.conformation().add_parameters_set( utility::pointer::dynamic_pointer_cast<ParametersSet>( output_parameters_set ) );

		if ( reset_pose() ) {
			if ( TR.Debug.visible() ) TR.Debug << "Replacing input pose with helical bundle." << std::endl;
			pose.clear();
			pose=newpose;
		} else {
			if ( TR.Debug.visible() ) TR.Debug << "Appending helical bundle to pose." << std::endl;
			if ( pose.pdb_info() ) {
				pose.pdb_info()->detach_from();
				pose.pdb_info(nullptr);
			}
			pose.append_pose_by_jump(newpose, 1);
		}
	}

	if ( TR.Debug.visible() ) TR.Debug << "Finished apply function." << std::endl;

	return;
}

////////////////////////////////////////////////////////////////////////////////


/// @brief Returns the name of this mover ("MakeBundle").
std::string MakeBundle::get_name() const{
	return "MakeBundle";
}

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////

/// @brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
///
void
MakeBundle::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data_map*/,
	protocols::filters::Filters_map const &/*filters*/,
	protocols::moves::Movers_map const &/*movers*/,
	core::pose::Pose const & /*pose*/
) {

	if ( tag->getName() != "MakeBundle" ) {
		throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible -- the tag name does not match the mover name.");
	}

	if ( TR.visible() ) TR << "Parsing options for MakeBundle (\"" << tag->getOption<std::string>("name" ,"") << "\") mover." << std::endl;

	runtime_assert_string_msg( !tag->hasOption("delta_omega1_all"), "The delta_omega1_all option has been renamed delta_omega1 for simplicity.  Please update your scripts accordingly." );

	//Get whether we're using radians or degrees:
	set_use_degrees( tag->getOption<bool>("use_degrees", false) );
	if ( TR.visible() ) TR << "Reading " << (use_degrees() ? "degrees" : "radians") << " for all input angle measurements.  (Radians are used internally, and output is in radians.)" << std::endl;

	//Set symmetry options for the MakeBundle mover:
	set_symmetry_options_from_tag(tag);

	//Set the reset option:
	if ( tag->hasOption("reset") ) {
		set_reset_pose( tag->getOption<bool>("reset", false) );
		if ( TR.visible() ) TR << "Reset mode set to " << (reset_pose() ? "true" : "false") << "." << std::endl;
	}

	//Set defaults for whether the mover can set bond lengths, bond angles, and dihedrals:
	set_dofs_from_tag(tag);

	//Set defaults for the major helix:
	if ( tag->hasOption("r0") ) {
		set_default_r0( tag->getOption<core::Real>("r0", 0.0) );
		if ( TR.visible() ) TR << "Default r0 (major helix radius) set to " << default_r0() << "." << std::endl;
	}
	if ( tag->hasOption("omega0") ) {
		set_default_omega0( tag->getOption<core::Real>("omega0", 0.0) );
		if ( TR.visible() ) TR << "Default omega0 (major helix turn per residue) set to " << default_omega0() << "." << std::endl;
	}
	if ( tag->hasOption("delta_omega0") ) {
		set_default_delta_omega0( tag->getOption<core::Real>("delta_omega0", 0.0) );
		if ( TR.visible() ) TR << "Default delta_omega0 (major helix rotation) set to " << default_delta_omega0() << "." << std::endl;
	}

	//Set defaults for the minor helix:
	set_minorhelix_defaults_from_tag( tag );
	if ( tag->hasOption("delta_omega1") ) {
		set_default_delta_omega1_all( tag->getOption<core::Real>("delta_omega1", 0) );
		if ( TR.visible() ) TR << "Default delta_omega1 (minor helix rotation) set to " << default_delta_omega1_all() << "." << std::endl;
	}

	//Set defaults for other params:
	if ( tag->hasOption("delta_t") ) {
		set_default_delta_t( tag->getOption<core::Real>("delta_t", 0) );
		if ( TR.visible() ) TR << "Default delta_t (residue offset) set to " << default_delta_t() << "." << std::endl;
	}
	if ( tag->hasOption("z1_offset") ) {
		set_default_z1_offset( tag->getOption<core::Real>("z1_offset", 0) );
		if ( TR.visible() ) TR << "Default z1 offset (helix offset along the minor helix axis) set to " << default_z1_offset() << "." << std::endl;
	}
	if ( tag->hasOption("z0_offset") ) {
		set_default_z0_offset( tag->getOption<core::Real>("z0_offset", 0) );
		if ( TR.visible() ) TR << "Default z0 offset (helix offset along the major helix axis) set to " << default_z0_offset() << "." << std::endl;
	}
	set_other_defaults_from_tag( tag );

	//Check that at least one helix is defined:
	bool at_least_one_helix = false;

	//Parse sub-tags, setting major and minor helix params.  (Note that the add_helix function will automatically set the defaults).
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	for ( auto const & branch_tag : branch_tags ) {
		if ( branch_tag->getName() == "Helix" ) { //A helix has been added.  Add it, and parse its options.
			at_least_one_helix = true;
			runtime_assert_string_msg( !branch_tag->hasOption("delta_omega1_all"), "The delta_omega1_all option has been renamed delta_omega1 for simplicity.  Please update your scripts accordingly." );
			add_helix(); //This will set all of the parameters to defaults
			core::Size const helix_index( n_helices() ); //The index of the current helix
			if ( TR.visible() ) TR << "Added a helix." << std::endl;

			//Set crick_params_file, set_bondlengths, set_bondangles, and set_dihedrals options for this helix, based on the tag.
			set_helix_params_from_tag( helix_index, branch_tag );

			//Major helix params:
			if ( branch_tag->hasOption( "r0" ) ) {
				core::Real const r0val( branch_tag->getOption<core::Real>("r0", 0) );
				helix(helix_index)->set_r0(r0val);
				if ( TR.visible() ) TR << "\tSet r0 value (major helix radius) to " << r0val << "." << std::endl;
			}
			if ( branch_tag->hasOption( "omega0" ) ) {
				core::Real const omega0val( branch_tag->getOption<core::Real>("omega0", 0) );
				helix(helix_index)->set_omega0( convert_angle(omega0val) );
				if ( TR.visible() ) TR << "\tSet omega0 value (major helix turn per residue) to " << omega0val << "." << std::endl;
			}
			if ( branch_tag->hasOption( "delta_omega0" ) ) {
				core::Real const delta_omega0val( branch_tag->getOption<core::Real>("delta_omega0", 0) );
				helix(helix_index)->set_delta_omega0( convert_angle(delta_omega0val) );
				if ( TR.visible() ) TR << "\tSet delta_omega0 value (major helix rotation) to " << delta_omega0val << "." << std::endl;
			}

			//Set minor helix params (omega1, z1, delta_omega1):
			set_minor_helix_params_from_tag(helix_index, branch_tag);

			//Set residue_name, invert, and helix_length for this helix:
			set_other_helix_params_from_tag(helix_index, branch_tag);
			if ( branch_tag->hasOption( "delta_t" ) ) {
				core::Real const delta_tval( branch_tag->getOption<core::Real>("delta_t", 0) );
				helix(helix_index)->set_delta_t(delta_tval);
				if ( TR.visible() ) TR << "\tSet delta_t value (residue offset) to " << delta_tval << "." << std::endl;
			}
			if ( branch_tag->hasOption( "z1_offset" ) ) {
				core::Real const z1_offsetval( branch_tag->getOption<core::Real>("z1_offset", 0) );
				helix(helix_index)->set_z1_offset(z1_offsetval);
				if ( TR.visible() ) TR << "\tSet z1 offset value (helix offset along the minor helix axis) to " << z1_offsetval << "." << std::endl;
			}
			if ( branch_tag->hasOption( "z0_offset" ) ) {
				core::Real const z0_offsetval( branch_tag->getOption<core::Real>("z0_offset", 0) );
				helix(helix_index)->set_z0_offset(z0_offsetval);
				if ( TR.visible() ) TR << "\tSet z0 offset value (helix offset along the major helix axis) to " << z0_offsetval << "." << std::endl;
			}

		}
	}

	runtime_assert_string_msg(at_least_one_helix, "In protocols::helical_bundle::MakeBundle::parse_my_tag(): At least one helix must be defined using a <Helix ...> sub-tag!");

	return;
} //parse_my_tag

/// @brief Set crick_params_file, set_bondlengths, set_bondangles, and set_dihedrals options for a single helix, based on an input tag.
///
void MakeBundle::set_helix_params_from_tag ( core::Size const helix_index, utility::tag::TagCOP tag )
{
	//FIRST, check for a minor helix Crick params file, and set minor helix params accordingly.
	if ( tag->hasOption( "crick_params_file" ) ) {
		std::string const paramsfile = tag->getOption<std::string>( "crick_params_file", "" );
		helix(helix_index)->set_minor_helix_params_from_file( paramsfile );
		if ( TR.visible() ) TR << "\tRead minor helix parameters from " << paramsfile << "." << std::endl;
	}
	//SECOND, check for params that have been set manually, and set accordingly.
	//Params for what DOFs can be set:
	if ( tag->hasOption("set_bondlengths") ) {
		helix(helix_index)->set_allow_bondlengths( tag->getOption<bool>("set_bondlengths", true) );
		if ( TR.visible() ) {
			TR << "\tSet the permission for the mover to set bondlengths for this helix to " << (default_allow_bondlengths() ? "true." : "false.") << std::endl;
		}
	}
	if ( tag->hasOption("set_bondangles") ) {
		helix(helix_index)->set_allow_bondangles( tag->getOption<bool>("set_bondangles", true) );
		if ( TR.visible() ) {
			TR << "\tSet the permission for the mover to set bondangles for this helix to " << (default_allow_bondangles() ? "true." : "false.") << std::endl;
		}
	}
	if ( tag->hasOption("set_dihedrals") ) {
		helix(helix_index)->set_allow_dihedrals( tag->getOption<bool>("set_dihedrals", true) );
		if ( TR.visible() ) {
			TR << "\tSet the permission for the mover to set dihedrals for this helix to " << (default_allow_dihedrals() ? "true." : "false.") << std::endl;
		}
	}
	return;
}

/// @brief Set omega1, z1, and delta_omega1 for a single helix, based on an input tag.
///
void MakeBundle::set_minor_helix_params_from_tag( core::Size const helix_index, utility::tag::TagCOP tag )
{
	if ( tag->hasOption( "omega1" ) ) {
		core::Real const omega1val( tag->getOption<core::Real>("omega1", 0) );
		helix(helix_index)->set_omega1( convert_angle(omega1val) );
		if ( TR.visible() ) TR << "\tSet omega1 value (minor helix turn per residue) to " << omega1val << "." << std::endl;
	}
	if ( tag->hasOption( "z1" ) ) {
		core::Real const z1val( tag->getOption<core::Real>("z1", 0) );
		helix(helix_index)->set_z1(z1val);
		if ( TR.visible() ) TR << "\tSet z1 value (minor helix rise per residue) to " << z1val << "." << std::endl;
	}
	if ( tag->hasOption( "delta_omega1" ) ) {
		core::Real const deltaomega1val( tag->getOption<core::Real>("delta_omega1", 0) );
		helix(helix_index)->set_delta_omega1_all( convert_angle(deltaomega1val) );
		if ( TR.visible() ) TR << "\tSet delta_omega1 value (minor helix rotation) to " << deltaomega1val << "." << std::endl;
	}
	return;
}

/// @brief Set residue_name, invert, and helix_length for a single helix, based on an input tag.
///
void MakeBundle::set_other_helix_params_from_tag( core::Size const helix_index, utility::tag::TagCOP tag )
{
	if ( tag->hasOption( "residue_name" ) ) {
		std::string const resname( tag->getOption<std::string>("residue_name", "ALA") );
		utility::vector1< std::string > resnames;
		parse_resnames(resname, resnames);
		helix(helix_index)->set_residue_name(resnames);
		if ( TR.visible() ) {
			if ( resnames.size()==1 ) {
				TR << "\tSet the residue type for the helix to " << resnames[1] << "." << std::endl;
			} else {
				TR << "\tSet the following residue types for the repeating unit in this helix: ";
				for ( core::Size i=1, imax=resnames.size(); i<=imax; ++i ) {
					TR << resnames[i];
					if ( i<imax ) TR << ", ";
				}
				TR << "." << std::endl;
			}
		}
	}
	if ( tag->hasOption( "repeating_unit_offset" ) ) {
		core::Size const repeat_offset( tag->getOption<core::Size>("repeating_unit_offset", 0 ) );
		helix(helix_index)->set_repeating_unit_offset(repeat_offset);
		if ( TR.visible() ) TR << "\tSet repeat unit offset for this helix to " << repeat_offset << "." << std::endl;
	}
	if ( tag->hasOption( "invert" ) ) {
		bool const invertval( tag->getOption<bool>("invert", false) );
		helix(helix_index)->set_invert_helix(invertval);
		if ( TR.visible() ) TR << "\tSet invert value (should the helix be flipped?) to " << (invertval ? "true" : "false") << "." << std::endl;
	}
	if ( tag->hasOption( "helix_length" ) ) {
		core::Size const helixlength( tag->getOption<core::Size>("helix_length", 0) );
		helix(helix_index)->set_helix_length(helixlength);
		if ( TR.visible() ) TR << "\tSet helix length to " << helixlength << "." << std::endl;
	}
	return;
}

/// @brief Set symmetry and symmetry_copies options based on an input tag.
///
void MakeBundle::set_symmetry_options_from_tag( utility::tag::TagCOP tag )
{
	//Set options for the MakeBundle mover:
	if ( tag->hasOption("symmetry") ) {
		set_symmetry( tag->getOption<core::Size>("symmetry", 0) );
		if ( TR.visible() ) {
			if ( symmetry()<2 ) TR << "Symmetry mode set to false." << std::endl;
			else TR << symmetry() << "-fold symmetry set." << std::endl;
		}
	}
	if ( tag->hasOption("symmetry_copies") ) {
		set_symmetry_copies( tag->getOption<core::Size>("symmetry_copies", 0) );
		if ( TR.visible() ) {
			if ( symmetry_copies()<1 ) TR << "All symmetry copies will be generated." << std::endl;
			if ( symmetry_copies()==1 ) TR << "Only the first symmetry copy will be generated." << std::endl;
			else {
				TR << symmetry_copies() << " symmetry copies will be generated." << std::endl;
				if ( symmetry_copies() > symmetry() ) {
					TR << "Note!  The number of symmetry copies is greater than the bundle symmetry.  This only really makes sense for helical repeat symmetries (i.e. symmetries with z-offset)." << std::endl;
				}
			}
		}
	}

	return;
} //set_global_options_from_tag

/// @brief Set defaults for whether the mover can set bond lengths, bond angles, and dihedrals.
///
void MakeBundle::set_dofs_from_tag( utility::tag::TagCOP tag )
{
	if ( tag->hasOption("set_bondlengths") ) {
		set_default_allow_bondlengths( tag->getOption<bool>("set_bondlengths", true) );
		if ( TR.visible() ) {
			TR << "Set the default permission for the mover to set bondlengths to " << (default_allow_bondlengths() ? "true." : "false.") << std::endl;
		}
	}
	if ( tag->hasOption("set_bondangles") ) {
		set_default_allow_bondangles( tag->getOption<bool>("set_bondangles", true) );
		if ( TR.visible() ) {
			TR << "Set the default permission for the mover to set bondangles to " << (default_allow_bondangles() ? "true." : "false.") << std::endl;
		}
	}
	if ( tag->hasOption("set_dihedrals") ) {
		set_default_allow_dihedrals( tag->getOption<bool>("set_dihedrals", true) );
		if ( TR.visible() ) {
			TR << "Set the default permission for the mover to set dihedrals to " << (default_allow_dihedrals() ? "true." : "false.") << std::endl;
		}
	}
	return;
} //set_dofs_from_tag

/// @brief Set the crick_params_file, omega1, and z1 default values from a tag.
///
void MakeBundle::set_minorhelix_defaults_from_tag( utility::tag::TagCOP tag )
{
	if ( tag->hasOption("crick_params_file") ) {
		set_default_crick_params_file( tag->getOption<std::string>("crick_params_file", "") );
		if ( TR.visible() ) TR << "Default Crick params file set to " << default_crick_params_file() << "." << std::endl;
	}
	if ( tag->hasOption("omega1") ) {
		set_default_omega1( convert_angle( tag->getOption<core::Real>("omega1", 0) ) );
		if ( TR.visible() ) TR << "Default omega1 (minor helix turn per residue) set to " << default_omega1() << "." << std::endl;
	}
	if ( tag->hasOption("z1") ) {
		set_default_z1( tag->getOption<core::Real>("z1", 0) );
		if ( TR.visible() ) TR << "Default z1 (minor helix rise per residue) set to " << default_z1() << "." << std::endl;
	}
	return;
} //set_minorhelix_defaults_from_tag

/// @brief Set the residue_name, invert, and helix_length default values based on an input tag.
///
void MakeBundle::set_other_defaults_from_tag( utility::tag::TagCOP tag )
{
	if ( tag->hasOption("residue_name") ) {
		std::string resname_string( tag->getOption<std::string>("residue_name", "") );
		utility::vector1 < std::string > resname_vect;
		parse_resnames( resname_string, resname_vect );
		set_default_residue_name( resname_vect );
		if ( TR.visible() ) {
			if ( resname_vect.size()==1 ) TR << "Default residue name set to " << default_residue_name(1) << "." << std::endl;
			else {
				TR << "Default residue names set to ";
				for ( core::Size i=1, imax=resname_vect.size(); i<=imax; ++i ) {
					TR << default_residue_name(i);
					if ( i<imax ) TR << ", ";
				}
				TR << "." << std::endl;
			}
		}
	}
	if ( tag->hasOption( "repeating_unit_offset" ) ) {
		core::Size const repeat_offset( tag->getOption<core::Size>("repeating_unit_offset", 0 ) );
		set_default_repeating_unit_offset(repeat_offset);
		if ( TR.visible() ) TR << "Default repeat unit offset set to " << repeat_offset << "." << std::endl;
	}
	if ( tag->hasOption("invert") ) {
		set_default_invert( tag->getOption<bool>("invert", false) );
		if ( TR.visible() ) TR << "Default invert (should the helix be flipped?) set to " << (default_invert() ? "true" : "false") << "." << std::endl;
	}
	if ( tag->hasOption("helix_length") ) {
		set_default_helix_length( tag->getOption<core::Size>("helix_length", 0) );
		if ( TR.visible() ) TR << "Default helix length set to " << default_helix_length() << "." << std::endl;
	}
	return;
} //set_other_defaults_from_tag

/// @brief Function to add a helix.
/// @details This creates a MakeBundleHelix mover that will be called at apply time.
/// The new mover is only initialized if default values are provided by THIS mover.
void MakeBundle::add_helix() {
	make_bundle_helix_movers_.push_back( protocols::helical_bundle::MakeBundleHelixOP(new protocols::helical_bundle::MakeBundleHelix()) );
	core::Size const newindex( make_bundle_helix_movers_.size() );
	if ( default_r0_set() ) helix(newindex)->set_r0( default_r0() );
	if ( default_omega0_set() ) helix(newindex)->set_omega0( default_omega0() );
	if ( default_delta_omega0_set() ) helix(newindex)->set_delta_omega0( default_delta_omega0() );

	if ( default_allow_bondlengths_set() ) helix(newindex)->set_allow_bondlengths( default_allow_bondlengths() );
	if ( default_allow_bondangles_set() ) helix(newindex)->set_allow_bondangles( default_allow_bondangles() );
	if ( default_allow_dihedrals_set() ) helix(newindex)->set_allow_dihedrals( default_allow_dihedrals() );

	if ( default_crick_params_file_set() ) helix(newindex)->set_minor_helix_params_from_file(default_crick_params_file()); //The default Crick params are read in from a file, but...
	//Manual settings can override what was read in:
	if ( default_omega1_set() ) helix(newindex)->set_omega1( default_omega1() );
	if ( default_z1_set() ) helix(newindex)->set_z1( default_z1() );
	if ( default_delta_omega1_all_set() ) helix(newindex)->set_delta_omega1_all( default_delta_omega1_all() );
	if ( default_residue_name_set() ) helix(newindex)->set_residue_name( default_residue_name() );
	if ( default_repeating_unit_offset_set() ) helix(newindex)->set_repeating_unit_offset( default_repeating_unit_offset() );
	if ( default_delta_t_set() ) helix(newindex)->set_delta_t( default_delta_t() );
	if ( default_z1_offset_set() ) helix(newindex)->set_z1_offset( default_z1_offset() );
	if ( default_z0_offset_set() ) helix(newindex)->set_z0_offset( default_z0_offset() );
	if ( default_invert_set() ) helix(newindex)->set_invert_helix( default_invert() );
	if ( default_helix_length_set() ) helix(newindex)->set_helix_length( default_helix_length() );
	return;
}

/// @brief Returns the default r0 value (major helix radius).
///
core::Real MakeBundle::default_r0() const {
	runtime_assert_string_msg( default_r0_set(), "In protocols::helical_bundle::MakeBundle::default_r0() The default r0 value has not been set!" );
	return default_r0_;
}

/// @brief Returns the default omega0 value (major helix rise per residue).
///
core::Real MakeBundle::default_omega0() const {
	runtime_assert_string_msg( default_omega0_set(), "In protocols::helical_bundle::MakeBundle::default_omega0(): The default omega0 value has not been set!" );
	return default_omega0_;
}

/// @brief Returns the default delta_omega0 value (major helix rotation).
///
core::Real MakeBundle::default_delta_omega0() const {
	runtime_assert_string_msg( default_delta_omega0_set(), "In protocols::helical_bundle::MakeBundle::default_delta_omega0(): The default delta_omega0 value has not been set!" );
	return default_delta_omega0_;
}

/// @brief Returns the default Crick params file name.
///
std::string MakeBundle::default_crick_params_file() const {
	runtime_assert_string_msg( default_crick_params_file_set(), "In protocols::helical_bundle::MakeBundle::default_crick_params_file(): The default crick_params_file value has not been set!" );
	return default_crick_params_file_;
}

/// @brief Returns the default omega1 value (minor helix turn per residue).
///
core::Real MakeBundle::default_omega1() const {
	runtime_assert_string_msg( default_omega1_set(), "In protocols::helical_bundle::MakeBundle::default_omega1(): The default omega1 value has not been set!" );
	return default_omega1_;
}

/// @brief Returns the default z1 value (minor helix turn per residue).
///
core::Real MakeBundle::default_z1() const {
	runtime_assert_string_msg( default_z1_set(), "In protocols::helical_bundle::MakeBundle::default_z1(): The default z1 value has not been set!" );
	return default_z1_;
}

/// @brief Returns the default delta_omega1_all value (minor helix rotation).
///
core::Real MakeBundle::default_delta_omega1_all() const {
	runtime_assert_string_msg( default_delta_omega1_all_set(), "In protocols::helical_bundle::MakeBundle::default_delta_omega1_all(): The default delta_omega1_all value has not been set!" );
	return default_delta_omega1_all_;
}

/// @brief Returns the default residue name for a particular index in the repeating unit.
///
std::string MakeBundle::default_residue_name( core::Size const index_in_repeating_unit ) const {
	runtime_assert_string_msg( default_residue_name_set(), "In protocols::helical_bundle::MakeBundle::default_residue_name(): The default residue_name value has not been set!" );
	runtime_assert_string_msg( index_in_repeating_unit <= default_residue_name_.size() && index_in_repeating_unit > 0, "In protocols::helical_bundle::MakeBundle::default_residue_name(): The index provided to this function is out of range." );
	return default_residue_name_[index_in_repeating_unit];
}

/// @brief Returns the default repeating unit offset value.
/// @details An offset of 0 means that the first residue of the helix is the first residue of the repeating unit.  An offset
/// of 1 means that the first residue of the helix is the second residue of the repeating unit, etc.
core::Size MakeBundle::default_repeating_unit_offset() const {
	runtime_assert_string_msg( default_repeating_unit_offset_set(), "In protocols::helical_bundle::MakeBundle::default_repeating_unit_offset():  The default repeating unit offset value has not been set!" );
	return default_repeating_unit_offset_;
}

/// @brief Returns the default residue name vector.
///
utility::vector1 < std::string > MakeBundle::default_residue_name() const {
	runtime_assert_string_msg( default_residue_name_set(), "In protocols::helical_bundle::MakeBundle::default_residue_name(): The default residue_name value has not been set!" );
	return default_residue_name_;
}

/// @brief Returns the default delta_t value (residue offset).
///
core::Real MakeBundle::default_delta_t() const {
	runtime_assert_string_msg( default_delta_t_set(), "In protocols::helical_bundle::MakeBundle::default_delta_t(): The default delta_t value has not been set!" );
	return default_delta_t_;
}

/// @brief Returns the default z1_offset value (helix offset along the minor helix axis).
///
core::Real MakeBundle::default_z1_offset() const {
	runtime_assert_string_msg( default_z1_offset_set(), "In protocols::helical_bundle::MakeBundle::default_z1_offset(): The default z1_offset value has not been set!" );
	return default_z1_offset_;
}

/// @brief Returns the default z0_offset value (helix offset along the major helix axis).
///
core::Real MakeBundle::default_z0_offset() const {
	runtime_assert_string_msg( default_z0_offset_set(), "In protocols::helical_bundle::MakeBundle::default_z0_offset(): The default z0_offset value has not been set!" );
	return default_z0_offset_;
}

/// @brief Returns the default invert value (should the helix be flipped?)
///
bool MakeBundle::default_invert() const {
	runtime_assert_string_msg( default_invert_set(), "In protocols::helical_bundle::MakeBundle::default_invert(): The default invert value has not been set!" );
	return default_invert_;
}

/// @brief Returns the default number of residues in each helix
///
core::Size MakeBundle::default_helix_length() const {
	runtime_assert_string_msg( default_helix_length_set(), "In protocols::helical_bundle::MakeBundle::default_helix_length(): The default number of residues per helix has not been set!" );
	return default_helix_length_;
}

/// @brief Returns the default for whether bond lengths should be set by the mover
///
bool MakeBundle::default_allow_bondlengths() const {
	runtime_assert_string_msg( default_allow_bondlengths_set(), "In protocols::helical_bundle::MakeBundle::default_allow_bondlengths(): The default allow_bondlengths value has not been set!" );
	return default_allow_bondlengths_;
}

/// @brief Returns the default for whether bond angles should be set by the mover
///
bool MakeBundle::default_allow_bondangles() const {
	runtime_assert_string_msg( default_allow_bondangles_set(), "In protocols::helical_bundle::MakeBundle::default_allow_bondangles(): The default allow_bondangles value has not been set!" );
	return default_allow_bondangles_;
}

/// @brief Returns the default for whether dihedral angles should be set by the mover
///
bool MakeBundle::default_allow_dihedrals() const {
	runtime_assert_string_msg( default_allow_dihedrals_set(), "In protocols::helical_bundle::MakeBundle::default_allow_dihedrals(): The default allow_dihedrals value has not been set!" );
	return default_allow_dihedrals_;
}

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////


} //namespace helical_bundle
} //namespace protocols
