// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

///
/// @author Oliver Lange
/// Archive class to collect structures such that variances of scores can be computed to determine normalized weights


#ifndef INCLUDED_protocols_abinitio_HedgeArchive_hh
#define INCLUDED_protocols_abinitio_HedgeArchive_hh

// Unit Headers
//#include <protocols/abinitio/IterativeAbrelax.fwd.hh>

// Package Headers
#include <protocols/jd2/archive/EvaluatedArchive.hh>
#include <protocols/jd2/archive/ArchiveManager.fwd.hh>
#include <protocols/abinitio/HedgeArchive.fwd.hh>
// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/silent/silent.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// Third-party Headers

//// C++ headers

#include <string>
#include <set>



namespace protocols {
namespace abinitio {

class HedgeArchive : public jd2::archive::EvaluatedArchive {
  typedef jd2::archive::EvaluatedArchive Parent;
  typedef std::list< std::pair< core::Real, core::io::silent::SilentStructOP > > SilentStructs;
  typedef std::map< core::Size, SilentStructs > BatchStructuresMap;
public:
  HedgeArchive( std::string name );

	//called when new decoy is available
  virtual bool add_evaluated_structure(
		 core::io::silent::SilentStructOP,
		 core::io::silent::SilentStructOP alternative_decoy,
		 jd2::archive::Batch const& batch
	);

  virtual void generate_batch() {};

	virtual void rescore() {}; //do nothing since we don't care about scores

	///@brief save and restore status of archive to file-system
	// will also call save_pending_decoys and restore 'incoming_structures_' from the pending-decoy files
  virtual void save_status( std::ostream& ) const;
  virtual void restore_status( std::istream& );


  ///@brief overloaded to make input decoys appear the same as decoys coming from batches
  virtual void init_from_decoy_set( core::io::silent::SilentFileData const& ) {};
	void collect( jd2::archive::Batch const& batch, core::io::silent::SilentStructOPs& ) const;
protected:

	//a batch is completed ---
	void incorporate_batch( core::Size batch_id );
private:

	//store pending decoys in the 'incoming_structures_' in files ( required for restarts )
  void save_pending_decoys( SilentStructs const& decoys, core::Size batch_id  ) const;

	//remove files with pending decoys
  void remove_pending_decoys( core::Size batch_id  ) const;
  core::Real score_cut_per_batch_;
  core::Real add_fuzzy_; //0 for strictly score based, 1 for totally random

	//pending decoys -- these are decoys from incomplete batches, when batch is completed we incorporate a random subset of them into archive
  BatchStructuresMap incoming_structures_;

  typedef   std::set< core::Size > BatchList;
	BatchList old_batches_;
};


}
}


#endif
