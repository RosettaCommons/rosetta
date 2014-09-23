// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file MinMover.cc
/// @brief
/// @author ashworth

// Unit headers
#include <protocols/simple_moves/symmetry/SymmetricRMSMover.hh>

// Package headers
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DiagnosticData.hh>
#include <core/pose/datacache/CacheableDataType.hh>

// AUTO-REMOVED #include <protocols/jd2/ScoreMap.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
// ObjexxFCL Headers

// C++ Headers

// Utility Headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {
namespace symmetry {

using namespace ObjexxFCL;
static thread_local basic::Tracer TR( "protocols.simple_moves.symmetry.SymmetricRMSMover" );

SymmetricRMSMover::SymmetricRMSMover()
	: protocols::moves::Mover("SymmetricRMSMover") {}

SymmetricRMSMover::~SymmetricRMSMover(){}

void
SymmetricRMSMover::apply( core::pose::Pose & pose )
{
	using namespace basic::datacache;
	//using core::pose::datacache::CacheableDataType::SCORE_MAP;

	std::map < std::string, core::Real > score_map;
	if( !(pose.data().has( core::pose::datacache::CacheableDataType::SCORE_MAP ) ) ) {
		score_map[ "NO_OUTPUT_TAG_CACHED_SORRY" ] = 0.0;
	} else {
		score_map = ( static_cast< DiagnosticData const &>( pose.data().get( core::pose::datacache::CacheableDataType::SCORE_MAP ))).data() ;
	}

	assert( core::pose::symmetry::is_symmetric( pose ) );
	SymmetricConformation const & SymmConf (
		dynamic_cast<SymmetricConformation const & > ( pose.conformation()) );
	core::conformation::symmetry::SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	FArray1D_bool superpos ( pose.total_residue(), false );
	for (Size res=1; res <= symm_info->num_total_residues_without_pseudo(); ++res )
	{
		superpos(res) = true;
	}
	if ( get_native_pose() ) {
		core::Real const rms( core::scoring::rmsd_with_super_subset( *get_native_pose(), pose, superpos, core::scoring::is_protein_CA ) );
		score_map[ "rms" ] = rms;
		using namespace basic::datacache;
		pose.data().set(core::pose::datacache::CacheableDataType::SCORE_MAP, DataCache_CacheableData::DataOP( new basic::datacache::DiagnosticData( score_map ) ));
	}

}

std::string
SymmetricRMSMover::get_name() const {
	return "SymmetricRMSMover";
}


} // symmetry
} // moves
} // protocols
