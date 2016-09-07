// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/SimpleGlycosylateMover.cc
/// @brief A mover for glycosylation of common biological glycosylations.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/carbohydrates/SimpleGlycosylateMover.hh>
#include <protocols/carbohydrates/SimpleGlycosylateMoverCreator.hh>

#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/selection.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/chemical/AA.hh>
#include <core/kinematics/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/chemical/carbohydrates/database_io.hh>

#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>

#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/io/util.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.carbohydrates.SimpleGlycosylateMover" );


namespace protocols {
namespace carbohydrates {
using namespace core::chemical::carbohydrates;
using namespace core::pose::carbohydrates;
using namespace core::kinematics;

SimpleGlycosylateMover::SimpleGlycosylateMover():
	protocols::moves::Mover( "SimpleGlycosylateMover" ),
	strip_existing_glycans_( true ),
	ref_pose_name_( "" ),
	idealize_glycosylation_( false )

{

}

SimpleGlycosylateMover::~SimpleGlycosylateMover()= default;

SimpleGlycosylateMover::SimpleGlycosylateMover( SimpleGlycosylateMover const & )= default;

void
SimpleGlycosylateMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & pose)
{

	using namespace core::kinematics;

	glycosylations_.clear();
	glycosylation_weights_.clear();
	parsed_positions_.clear();
	positions_.clear();

	if ( tag->hasOption("glycosylation") ) {
		glycosylations_ = utility::string_split_multi_delim( tag->getOption< std::string >("glycosylation"), ",'`~+*&|;. ");
	} else if ( tag->hasOption("glycosylations") ) {
		glycosylations_ = utility::string_split_multi_delim( tag->getOption< std::string >("glycosylations"), ",'`~+*&|;. ");
	} else {
		utility_exit_with_message("Must pass either glycosylation or glycosylations!");
	}

	if ( tag->hasOption("position") ) {
		parsed_positions_ = utility::string_split_multi_delim( tag->getOption< std::string >("position"), ",'`~+*&|;. ");
	} else if ( tag->hasOption("positions") ) {
		parsed_positions_ = utility::string_split_multi_delim( tag->getOption< std::string >("positions"), ",'`~+*&|;. ");
	} else if  ( protocols::rosetta_scripts::has_branch(tag, "MoveMap") ) {
		MoveMapOP mm = MoveMapOP( new MoveMap() );

		//protocols::rosetta_scripts::add_movemaps_to_datamap(tag, pose, data, false);
		protocols::rosetta_scripts::parse_movemap( tag, pose, mm, data, false );
		set_positions_from_movemap(mm);
	} else {
		utility_exit_with_message(" Must pass either position or positions");
	}


	if ( tag->hasOption("weights") ) {
		utility::vector1< std::string > weights = utility::string_split_multi_delim( tag->getOption< std::string >("weights"), ",'`~+*&|;. ");
		for ( core::Size i = 1; i <= weights.size(); ++i ) {
			glycosylation_weights_.push_back( utility::string2Real( weights[ i ]));
		}
	}

	//Convert positions to proper resnums.

	strip_existing_glycans_ = tag->getOption<bool>("strip_existing", strip_existing_glycans_);

	ref_pose_name_ = tag->getOption< std::string >("ref_pose_name", ref_pose_name_);
	idealize_glycosylation_ = tag->getOption< bool >("idealize_glycosylation", idealize_glycosylation_);
}

protocols::moves::MoverOP
SimpleGlycosylateMover::clone() const{
	return protocols::moves::MoverOP( new SimpleGlycosylateMover( *this ) );
}

/*
SimpleGlycosylateMover & SimpleGlycosylateMoveroperator=( SimpleGlycosylateMover const & src){
return SimpleGlycosylateMover( src );
}
*/


moves::MoverOP
SimpleGlycosylateMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new SimpleGlycosylateMover );
}

std::string
SimpleGlycosylateMover::get_name() const {
	return "SimpleGlycosylateMover";
}

void
SimpleGlycosylateMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

std::ostream &operator<< (std::ostream &os, SimpleGlycosylateMover const &mover)
{
	mover.show(os);
	return os;
}

void
SimpleGlycosylateMover::set_position(core::Size position){
	positions_.clear();
	positions_.push_back( position );

}

void
SimpleGlycosylateMover::set_positions(const utility::vector1<core::Size> &positions) {
	positions_ = positions;
}

void
SimpleGlycosylateMover::set_positions(const utility::vector1<bool> &positions){
	positions_.clear();
	if ( positions.size() == 0 ) {
		utility_exit_with_message(" Positions empty! ");
	}

	for ( core::Size i = 1; i <= positions.size(); ++i ) {
		if ( positions[ i ] ) {
			positions_.push_back( i );
		}
	}
}

void
SimpleGlycosylateMover::set_positions_from_movemap(core::kinematics::MoveMapCOP mm){
	positions_ = core::kinematics::get_residues_from_movemap_with_id(core::id::BB, *mm);
}

void
SimpleGlycosylateMover::set_glycosylation(const std::string &iupac_or_common_string){
	glycosylations_.clear();
	glycosylations_.push_back(iupac_or_common_string);
}

void
SimpleGlycosylateMover::set_glycosylations(const utility::vector1<std::string> &iupac_strings){
	glycosylations_ = iupac_strings;
}

void
SimpleGlycosylateMover::set_glycosylation_weights(const utility::vector1<core::Real> &weights){
	glycosylation_weights_ = weights;
}

void
SimpleGlycosylateMover::set_strip_existing_glycans(bool strip_existing){
	strip_existing_glycans_ = strip_existing;
}

void
SimpleGlycosylateMover::remove_index( utility::vector1< core::Size > & current_vector, core::Size resnum) const{

	//Not the fastest way to do this, but it should work for now.
	utility::vector1< core::Size > vector_copy = current_vector;
	current_vector.clear();
	for ( core::Size i = 1; i <= vector_copy.size(); ++i ) {
		core::Size item = vector_copy[ i ];
		if ( item != resnum ) {
			current_vector.push_back( vector_copy[ i ]);
		} else {
			continue;
		}
	}
}

utility::vector1< std::string >
SimpleGlycosylateMover::setup_and_load_iupac_sequences() const {

	using namespace core::chemical::carbohydrates;

	std::map< std::string, std::string > const & short_names = core::chemical::carbohydrates::CarbohydrateInfoManager::get_short_name_to_iupac_strings_map();

	std::string common_names_db = "chemical/carbohydrates/common_glycans/";
	utility::vector1< std::string > full_iupac_glycans;

	for ( core::Size i = 1; i <= glycosylations_.size(); ++i ) {
		std::string set_glycan = glycosylations_[ i ];

		std::map< std::string, std::string >::const_iterator it;
		it  = short_names.find( set_glycan );



		// Check for iupac extension and attempt to load the sequence from the iupac file.
		if ( it != short_names.end() ) {
			full_iupac_glycans.push_back( it->second );
			continue;
		} else if ( set_glycan.find(".iupac") != std::string::npos ) {
			// Check for short name
			std::string iupac_file_path =basic::database::find_database_path( common_names_db, set_glycan );
			std::string full_glycan = read_glycan_sequence_file( iupac_file_path );
			full_iupac_glycans.push_back( full_glycan );

			continue;

		} else {
			// Otherwise, it should be a full iupac glycan name.  Push it back.
			full_iupac_glycans.push_back( set_glycan );
			continue;
		}
	}

	return full_iupac_glycans;

}


void
SimpleGlycosylateMover::apply( core::pose::Pose& pose ){

	using namespace numeric::random;
	using namespace core::pose::carbohydrates;
	using namespace core::chemical::carbohydrates;

	//Since we may be deleting residues, we need to add a reference pose.

	//Convert parsed positions.
	if ( parsed_positions_.size() > 0 ) {
		positions_.clear();
		for ( core::Size i = 1; i <= parsed_positions_.size(); ++i ) {
			core::Size resnum = core::pose::parse_resnum( parsed_positions_[ i ], pose);
			if ( ref_pose_name_ != "" ) {
				resnum = pose.corresponding_residue_in_current( resnum, ref_pose_name_);
			}
			positions_.push_back(resnum);
		}
	}

	if ( glycosylations_.size() == 0 || positions_.size() == 0 ) {
		utility_exit_with_message(" Glycosylation(s) and position(s) need to be set! ");
	}

	if ( glycosylation_weights_.size() > 0 && glycosylation_weights_.size() != glycosylations_.size() ) {
		utility_exit_with_message("Number of weights must equal number of glycosylations!");
	}

	std::string ref_pose_name = "simple_glycosylate_mover";
	pose.reference_pose_from_current( ref_pose_name , true /* override */); //Make a refpose as we may be deleting/adding/etc.

	utility::vector1< std::string > glycosylations =  setup_and_load_iupac_sequences();
	WeightedSampler sampler;


	//Go through positions randomly until all are tried.
	// Each round we remove the residue from local positions.
	utility::vector1< core::Size > local_positions = positions_;
	for ( core::Size round = 1; round <= positions_.size(); ++round ) {


		core::Size local_p = numeric::random::random_range(1, local_positions.size());
		core::Size old_resnum = local_positions[ local_p ];
		core::Size resnum = pose.corresponding_residue_in_current( old_resnum , ref_pose_name);
		if ( local_positions.size() > 1 ) {
			remove_index(local_positions, resnum); //Why can't we have a method that takes an index and deletes it?
		}

		//Does the position already have a glycan attached to it?

		if ( strip_existing_glycans_ == false ) {
			utility_exit_with_message(" Glycan extension not currently implemented!");
		} else {
			delete_carbohydrate_branch( pose, resnum ); //Delete any carbohydrates currently attached.
		}

		std::string glycosylation;

		if ( glycosylations.size() == 1 ) {
			glycosylation = glycosylations[ 1 ];
		} else {
			if ( glycosylation_weights_.size() > 0 ) {
				sampler.weights( glycosylation_weights_ );
				glycosylation = glycosylations[ sampler.random_sample(numeric::random::rg()) ];
			} else {
				glycosylation = glycosylations[ random_range( 1, glycosylations.size() ) ];
			}
		}

		//Glycosylate the pose.
		TR << "Glycosylating at " << old_resnum << " : " << resnum << " " << glycosylation << std::endl;

		glycosylate_pose( pose, resnum, glycosylation, idealize_glycosylation_ /* idealize linkages - Seems to be a bug here!*/);

	}

	/* This doesn't work for O-linked glycosylations
	//Check positions are ASN.  Optionally mutate to ASN and motif?
	simple_moves::MutateResidue mutate = simple_moves::MutateResidue();
	mutate.set_preserve_atom_coords( true );
	for (core::Size i = 1; i <= positions_.size(); ++i ){
	core::Size resnum = positions_[ i ];
	if (! pose.residue( resnum ).aa() == core::chemical::aa_asn ){
	TR << "Mutating "<< resnum << " to ASN for glycosylation. " << std::endl;
	mutate.set_target( resnum );
	mutate.set_res_name( core::chemical::aa_asn );
	mutate.apply( pose );
	}
	}
	*/

}


/////////////// Creator ///////////////

protocols::moves::MoverOP
SimpleGlycosylateMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SimpleGlycosylateMover );
}

std::string
SimpleGlycosylateMoverCreator::keyname() const {
	return SimpleGlycosylateMoverCreator::mover_name();
}

std::string
SimpleGlycosylateMoverCreator::mover_name(){
	return "SimpleGlycosylateMover";
}

} //protocols
} //carbohydrates


