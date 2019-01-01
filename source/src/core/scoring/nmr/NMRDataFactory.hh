// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRDataFactory.hh
/// @brief   singleton class that handles the creation of different types of NMR data
///          through a common interface. At the moment, it supports three paramagnetic
///          NMR interactions (i.e. PCSs, RDCs, PREs) but other NMR data types can
///          be easily added. They should provide a constructor that takes a filename
///          as argument and is called from inside the factory function which returns an
///          owning pointer to the new created NMR data object.
/// @details last Modified: 07/01/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_NMRDataFactory_HH
#define INCLUDED_core_scoring_nmr_NMRDataFactory_HH

// Unit headers
#include <core/scoring/nmr/NMRDataFactory.fwd.hh>

// Package headers
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <core/scoring/nmr/pcs/PCSData.fwd.hh>
#include <core/scoring/nmr/rdc/RDCData.fwd.hh>
#include <core/scoring/nmr/pre/PREData.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/SingletonBase.hh>

// C++ headers
#include <string>
#include <map>

namespace core {
namespace scoring {
namespace nmr {

class NMRDataFactory : public utility::SingletonBase< NMRDataFactory > {
	friend class utility::SingletonBase< NMRDataFactory >;

public: // Methods

	/// @brief Factory function that takes an NMR data type and input filename
	///        as arguments and returns a pointer to the NMR data.
	basic::datacache::CacheableDataOP
	get_nmr_data(
		std::string const & nmr_data_type,
		std::string const & filename,
		pose::Pose const & pose
	);

	/// @brief Factory function that takes an NMR data type and input filename
	///        as arguments, creates NMR data and attaches them directly to the pose.
	void
	create_attach_nmr_data_to_pose(
		std::string const & nmr_data_type,
		std::string const & filename,
		pose::Pose & pose
	);

private: // Methods

	/// @brief Empty default constructor
	NMRDataFactory();

};

} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_NMRDataFactory_HH
