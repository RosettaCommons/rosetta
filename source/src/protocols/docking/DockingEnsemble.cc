// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DockingEnsemble.cc
/// @brief container class for ensemble docking information to be used with ConformerSwitchMover
/// @author Monica Berrondo

#include <protocols/docking/DockingEnsemble.hh>

// Rosetta Headers
#include <basic/datacache/CacheableStringFloatMap.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>
// AUTO-REMOVED #include <core/id/types.hh>
#include <core/pose/datacache/CacheableDataType.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>

// AUTO-REMOVED #include <protocols/scoring/Interface.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>


// Random number generator
#include <numeric/random/random.hh>
//
#include <string>

#include <basic/Tracer.hh>

#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/Jump.hh>
#include <basic/datacache/BasicDataCache.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.moves.DockingEnsemble" );

namespace protocols {
namespace docking {

/// @details Auto-generated virtual destructor
DockingEnsemble::~DockingEnsemble() {}


//constructor with arguments
DockingEnsemble::DockingEnsemble(
	Size start_res,
	Size end_res,
	Size jump_id,
	std::string ensemble_file_path,
	std::string partner,
	core::scoring::ScoreFunctionCOP scorefxn_low,
	core::scoring::ScoreFunctionCOP scorefxn_high
) :
	start_res_(start_res),
	end_res_(end_res),
	jump_id_(jump_id),
	ensemble_file_path_(ensemble_file_path),
	partner_(partner)
{
	// initialize current conf_num to zero
	conf_num_ = 0;

	scorefxn_low_ = scorefxn_low;
	scorefxn_high_ = scorefxn_high;

	load_ensemble();
	ensemble_size_ = ensemble_list_.size();

	runtime_assert(ensemble_size_ > 0);
	conf_size_ = ensemble_list_[1].total_residue();
	runtime_assert((end_res_ - start_res_ + 1) == conf_size_);

	TR << "ensemble summary: start_res_ " << start_res_ <<
		  " end_res_ " << end_res_ <<
		  " conf_size_ " << conf_size_ <<
		  " ensemble_size_ " << ensemble_size_ << std::endl;
}

void DockingEnsemble::load_ensemble()
{
	utility::vector1< core::Real > data;
	utility::io::izstream file( ensemble_file_path_ );
	std::string line;
	TR << "Loading Ensemble" << std::endl;

	// read file names
	while( getline(file, line) ) {
		// read pdb file names into filename array
		if ( line.find("pdb") != std::string::npos )
			pdb_filenames_.push_back( line );
		// read the rest of the data as Reals into data array
		else data.push_back( atof( line.c_str() ) );
	}

	if ( data.size() > 1 ) {
		core::Size count = 1;
		for ( Size i=count; i<count+pdb_filenames_.size(); ++i )
			lowres_reference_energies_.push_back( data[i] );

		count = count+pdb_filenames_.size();
		for ( Size i=count; i<count+pdb_filenames_.size(); ++i )
			highres_reference_energies_.push_back( data[i] );
	}

	ensemble_list_ = core::import_pose::poses_from_pdbs( pdb_filenames_ );
	ensemble_list_cen_ = core::import_pose::poses_from_pdbs( pdb_filenames_ );	// Add by DK

	//Add by DK
	for (Size i = 1; i <= ensemble_list_cen_.size(); i++){
		protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( core::chemical::CENTROID );
		to_centroid.apply( ensemble_list_cen_[i] );
	}
}

void DockingEnsemble::recover_conformer_sidechains( core::pose::Pose & pose )
{
	using namespace core::pose::datacache;
	core::pose::Pose recover_pose = ensemble_list_[conf_num_];
	protocols::simple_moves::ReturnSidechainMoverOP recover_mover( new protocols::simple_moves::ReturnSidechainMover( recover_pose, start_res_, end_res_ ) );
	recover_mover->apply( pose );

	// make sure that the pose has ARBITRARY_FLOAT_DATA in the DataCache
	if ( !pose.data().has( ( CacheableDataType::ARBITRARY_FLOAT_DATA ) ) ){
		using namespace basic::datacache;
		pose.data().set(
			CacheableDataType::ARBITRARY_FLOAT_DATA,
			DataCache_CacheableData::DataOP( new basic::datacache::CacheableStringFloatMap() )
		);
	}

	basic::datacache::CacheableStringFloatMapOP data
		= utility::pointer::dynamic_pointer_cast< basic::datacache::CacheableStringFloatMap >
		( pose.data().get_ptr(CacheableDataType::ARBITRARY_FLOAT_DATA) );

	data->map()[ partner_ ] = highres_reference_energies_[conf_num_];
	pose.data().set( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA, data );

    //jjg: need to cache pointer to this function in the pose too!
}

void DockingEnsemble::calculate_lowres_ref_energy( core::pose::Pose & pose )
{
	core::Real score_low = ( *scorefxn_low_ )( pose );
	TR << "score_low: " << score_low << std::endl;
	lowres_reference_energies_.push_back( score_low );
}

void DockingEnsemble::calculate_highres_ref_energy( core::Size conf_num )
{
	core::pose::Pose conformer = ensemble_list_[conf_num];
	protocols::simple_moves::SwitchResidueTypeSetMover to_fullatom( core::chemical::FA_STANDARD );
	to_fullatom.apply( conformer );

	pack_operations_->apply( conformer );

	core::Real score_high = ( *scorefxn_high_ )( conformer );
	TR << "score_high: " << score_high << std::endl;
	highres_reference_energies_.push_back( score_high );

	std::string filename( pdb_filenames_[highres_reference_energies_.size()] + ".ppk" );
	TR << "filename: " << filename << std::endl;
	conformer.dump_pdb( filename );
	pdb_filenames_[highres_reference_energies_.size()] = filename;
	TR << "filename in array: " << pdb_filenames_[highres_reference_energies_.size()] << std::endl;
}

void DockingEnsemble::update_pdblist_file()
{
	utility::io::ozstream data;
	data.open( ensemble_file_path_ );

	// write out new filenames
	core::Size conf = 1;
	while ( conf <= pdb_filenames_.size() ) {
		data << pdb_filenames_[conf] << '\n';
		++conf;
	}

	// normalize centroid reference energies
	core::Real adjustment = lowres_reference_energies_[1];
	for ( core::Size i=1; i<=lowres_reference_energies_.size(); ++i ) {
		core::Real ref_energy = lowres_reference_energies_[i];
		adjustment = std::min( ref_energy, adjustment );
	}
	for ( core::Size i=1; i<= lowres_reference_energies_.size(); ++i )
		lowres_reference_energies_[i] = lowres_reference_energies_[i] - adjustment;

	// write out centroid reference energies
	conf = 1;
	while ( conf <= pdb_filenames_.size() ) {
		data << lowres_reference_energies_[conf] << '\n';
		++conf;
	}

	// write out fullatom reference energies
	conf = 1;
	while ( conf <= pdb_filenames_.size() ) {
		data << highres_reference_energies_[conf] << '\n';
		++conf;
	}
}

void DockingEnsemble::set_packer( protocols::moves::SequenceMoverOP packer )
{
	pack_operations_ = packer;
}

}  // namespace docking
}  // namespace protocols
