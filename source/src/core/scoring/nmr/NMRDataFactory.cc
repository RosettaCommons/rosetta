// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRDataFactory.cc
/// @brief   Implementation of class NMRDataFactory
/// @details last Modified: 07/01/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/NMRDataFactory.hh>

// Package headers
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/scoring/nmr/pcs/PCSData.hh>
#include <core/scoring/nmr/rdc/RDCData.hh>
#include <core/scoring/nmr/pre/PREData.hh>

// Package headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>

// C++ headers
#include <string>
#include <map>


namespace core {
namespace scoring {
namespace nmr {

static basic::Tracer TR("core.scoring.nmr.NMRDataFactory");

/// @brief Empty constructor
NMRDataFactory::NMRDataFactory() {}

/// @brief Factory function that takes an NMR data type and input filename
///        as arguments and returns a pointer to the NMR data
basic::datacache::CacheableDataOP
NMRDataFactory::get_nmr_data(
	std::string const & nmr_data_type,
	std::string const & filename,
	pose::Pose const & pose)
{
	if ( nmr_data_type == "PCS" || nmr_data_type == "pcs" ) {
		return basic::datacache::CacheableDataOP( new pcs::PCSData(filename, pose) );
	} else if ( nmr_data_type == "RDC" || nmr_data_type == "rdc" ) {
		return basic::datacache::CacheableDataOP( new rdc::RDCData(filename, pose) );
	} else if ( nmr_data_type == "PRE" || nmr_data_type == "pre" ) {
		return basic::datacache::CacheableDataOP( new pre::PREData(filename, pose) );
	} else {
		utility_exit_with_message( "ERROR: Could not create NMR data. Data type was not found by NMRDataFactory." );
	}
}

/// @brief Factory function that takes an NMR data type and input filename
///        as arguments, creates NMR data and attaches them to the pose
void
NMRDataFactory::create_attach_nmr_data_to_pose(
	std::string const & nmr_data_type,
	std::string const & filename,
	pose::Pose & pose)
{
	if ( nmr_data_type == "PCS" || nmr_data_type == "pcs" ) {
		pcs::PCSDataOP pcs_data_ptr = pcs::PCSDataOP( new pcs::PCSData(filename, pose) );
		pose.data().set(core::pose::datacache::CacheableDataType::NMR_PCS_DATA, pcs_data_ptr);
	} else if ( nmr_data_type == "RDC" || nmr_data_type == "rdc" ) {
		rdc::RDCDataOP rdc_data_ptr = rdc::RDCDataOP( new rdc::RDCData(filename, pose) );
		pose.data().set(core::pose::datacache::CacheableDataType::NMR_RDC_DATA, rdc_data_ptr);
	} else if ( nmr_data_type == "PRE" || nmr_data_type == "pre" ) {
		pre::PREDataOP pre_data_ptr = pre::PREDataOP( new pre::PREData(filename, pose) );
		pose.data().set(core::pose::datacache::CacheableDataType::NMR_PRE_DATA, pre_data_ptr);
	} else {
		utility_exit_with_message( "ERROR: Could not create NMR data. Data type was not found by NMRDataFactory." );
	}
}

} // namespace nmr
} // namespace scoring
} // namespace core
