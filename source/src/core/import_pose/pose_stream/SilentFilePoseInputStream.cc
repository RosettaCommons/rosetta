// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/import_pose/pose_stream/SilentFilePoseInputStream.hh
/// @brief
/// @author James Thompson

// libRosetta headers

#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// C++ headers
#include <string>

namespace core {
namespace import_pose {
namespace pose_stream {

static basic::Tracer tr( "core.io.pose_stream.silent" );

typedef std::string string;
typedef utility::file::FileName FileName;
using core::io::silent::SilentFileOptions;

SilentFilePoseInputStream::SilentFilePoseInputStream()
: renumber_decoys_( false ), energy_cut_( 1.0 ), order_by_energy_( false ), record_source_( false )
{
	//  utility::vector1< FileName > empty;
	//  filenames(empty);
	sfd_ = core::io::silent::SilentFileDataOP( new core::io::silent::SilentFileData( SilentFileOptions() ) );
}

SilentFilePoseInputStream::SilentFilePoseInputStream( utility::vector1< FileName > const & fns )
: renumber_decoys_( false ), energy_cut_( 1.0 ), order_by_energy_( false ), record_source_( false )
{
	sfd_ = core::io::silent::SilentFileDataOP( new core::io::silent::SilentFileData( SilentFileOptions() ) );
	filenames(fns);
}

SilentFilePoseInputStream::SilentFilePoseInputStream( std::string const & fn )
: renumber_decoys_( false ), energy_cut_( 1.0 ), order_by_energy_( false ), record_source_( false )
{
	sfd_ = core::io::silent::SilentFileDataOP( new core::io::silent::SilentFileData( SilentFileOptions() ) );
	utility::vector1< FileName > fns;
	fns.push_back( fn );
	filenames(fns);
}

SilentFilePoseInputStream::SilentFilePoseInputStream(
	utility::vector1< FileName > const & fns,
	bool order_by_energy
)
: renumber_decoys_( false ), energy_cut_( 1.0 ), order_by_energy_( order_by_energy ), record_source_( false )
{
	sfd_ = core::io::silent::SilentFileDataOP( new core::io::silent::SilentFileData( SilentFileOptions() ) );
	filenames(fns);
}

SilentFilePoseInputStream::SilentFilePoseInputStream(
	utility::vector1< FileName > const & fns,
	core::Real energy_cut
)
: renumber_decoys_( false ), energy_cut_( energy_cut ), order_by_energy_( false ), record_source_( false )
{
	sfd_ = core::io::silent::SilentFileDataOP( new core::io::silent::SilentFileData( SilentFileOptions() ) );
	filenames(fns);
}

SilentFilePoseInputStream::SilentFilePoseInputStream(
	utility::vector1< FileName > const & fns,
	utility::vector1< string > const & input_tags
) :
	renumber_decoys_( false ), energy_cut_( 1.0 ), order_by_energy_( false ), record_source_( false )
{
	sfd_ = core::io::silent::SilentFileDataOP( new core::io::silent::SilentFileData( SilentFileOptions() ) );
	tags(input_tags);
	filenames(fns);
}

SilentFilePoseInputStream::SilentFilePoseInputStream(
	utility::vector1< FileName > const & fns,
	utility::vector1< string > const & input_tags,
	core::Real energy_cut
) :
	renumber_decoys_( false ), energy_cut_( energy_cut ), order_by_energy_( false ), record_source_( false )
{
	sfd_ = core::io::silent::SilentFileDataOP( new core::io::silent::SilentFileData( SilentFileOptions() ) );
	tags(input_tags);
	filenames(fns);
}

void SilentFilePoseInputStream::tags( utility::vector1< string > const & tags ) {
	tags_ = tags;
}

void SilentFilePoseInputStream::filenames(
	utility::vector1< FileName > const & filenames
) {
	filenames_  = filenames;
	// read in SilentStruct objects from all of the filenames
	read_all_files_();
	current_position_ = sfd_->begin();
}

void SilentFilePoseInputStream::reset()
{
	current_position_ = sfd_->begin();
}

/////////////////////////////////////////////////////
void
SilentFilePoseInputStream::set_silent_file_data(
	core::io::silent::SilentFileDataOP & sfd
) {
	sfd_ = sfd;

	filenames_.clear();

	if ( renumber_decoys() ) {
		tr.Debug << "renumbering decoys." << std::endl;
		sfd_->renumber_all_decoys();
	}

	if ( order_by_energy_ ) {
		tr.Debug << "ordering by energy" << std::endl;
		sfd_->order_by_energy();
	}

	current_position_ = sfd_->begin();
	tags_ = sfd_->tags();
}

void SilentFilePoseInputStream::renumber_decoys( bool const setting ) {
	renumber_decoys_ = setting;
}

utility::vector1< FileName > SilentFilePoseInputStream::filenames() const {
	return filenames_;
}

utility::vector1< std::string > SilentFilePoseInputStream::tags() const {
	return tags_;
}

bool SilentFilePoseInputStream::renumber_decoys() const {
	return renumber_decoys_;
}

core::Real SilentFilePoseInputStream::energy_cut() const {
	return energy_cut_;
}

bool SilentFilePoseInputStream::has_another_pose() {
	if ( current_position_ != sfd_->end() ) {
		return true;
	} else {
		return false;
	}
}

void SilentFilePoseInputStream::fill_pose(
	core::pose::Pose & pose,
	core::chemical::ResidueTypeSet const & residue_set
) {
	// check to make sure that we have more poses!
	if ( !has_another_pose() ) {
		utility_exit_with_message(
			"SilentFilePoseInputStream: called fill_pose, but I have no more Poses!"
		);
	}

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( option[ run::debug ]() ) {
		core::Real debug_rmsd = current_position_->get_debug_rmsd();
		tr.Debug << "RMSD to original coordinates for tag "
			<< current_position_->decoy_tag() << " = " << debug_rmsd << std::endl;
	}

	current_position_->fill_pose( pose, residue_set );

	// set up a tag using decoy_tag from SilentStruct
	using namespace basic::datacache;
	pose.data().set(
		core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG,
		DataCache_CacheableData::DataOP( new basic::datacache::CacheableString( current_position_->decoy_tag() ) )
	);
	tr.Debug << "decoy_tag() == " << current_position_->decoy_tag() << std::endl;

	core::pose::setPoseExtraScore( pose, "silent_score", current_position_->get_energy( "score" ) );

	preprocess_pose( pose );

	++current_position_;
} // fill_pose

void SilentFilePoseInputStream::fill_pose(
	core::pose::Pose & pose
) {
	// check to make sure that we have more poses!
	if ( !has_another_pose() ) {
		utility_exit_with_message(
			"SilentFilePoseInputStream: called fill_pose, but I have no more Poses!"
		);
	}

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( option[ run::debug ]() ) {
		core::Real debug_rmsd = current_position_->get_debug_rmsd();
		tr.Debug << "RMSD to original coordinates for tag "
			<< current_position_->decoy_tag() << " = " << debug_rmsd << std::endl;
	}

	current_position_->fill_pose( pose );

	// set up a tag using decoy_tag from SilentStruct
	using namespace basic::datacache;
	pose.data().set(
		core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG,
		DataCache_CacheableData::DataOP( new basic::datacache::CacheableString( current_position_->decoy_tag() ) )
	);
	tr.Debug << "decoy_tag() == " << current_position_->decoy_tag() << std::endl;

	core::pose::setPoseExtraScore( pose, "silent_score", current_position_->get_energy( "score" ) );

	preprocess_pose( pose );

	++current_position_;
} // fill_pose


core::io::silent::SilentStructOP
SilentFilePoseInputStream::next_struct()
{
	// check to make sure that we have more poses!
	if ( !has_another_pose() ) {
		utility_exit_with_message(
			"SilentFilePoseInputStream: called next_struct, but I have no more structs!"
		);
	}

	core::io::silent::SilentStructOP silent_struct( *(current_position_) );
	++current_position_;
	return silent_struct;
} // fill_pose


void SilentFilePoseInputStream::read_all_files_() {
	using utility::vector1;

	if ( record_source_ ) sfd_->set_record_source( true );

	for ( vector1< FileName >::const_iterator current_fn_ = filenames_.begin();
			current_fn_ != filenames_.end(); ++current_fn_
			) {

		if ( !file_exists( *current_fn_ ) ) {
			tr.Error << "Hey! Could not find " + std::string(*current_fn_)
				<< std::endl;
			continue;
		}

		tr.Debug << "reading " << *current_fn_ << std::endl;
		if ( tags_.size() > 0 ) {
			sfd_->read_file( *current_fn_, tags_ );
		} else {
			sfd_->read_file( *current_fn_ );
		}
	}

	if ( renumber_decoys() ) {
		tr.Debug << "renumbering decoys." << std::endl;
		sfd_->renumber_all_decoys();
	}

	if ( energy_cut() != 1.0 ) {
		sfd_->score_filter( energy_cut() );
	}

	if ( order_by_energy_ ) {
		tr.Debug << "ordering by energy" << std::endl;
		sfd_->order_by_energy();
	}
}

void
SilentFilePoseInputStream::set_order_by_energy( bool const & setting ) {
	order_by_energy_ = setting;
	if ( setting ) sfd_->order_by_energy();
	current_position_ = sfd_->begin();
}

void
SilentFilePoseInputStream::set_record_source( bool const & setting ) {
	record_source_ = setting;
}

} // pose_stream
} // import_pose
} // core
