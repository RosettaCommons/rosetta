// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/RNA_SecStructLegacyInfo.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_protocols_rna_RNA_SecStructLegacyInfo_hh
#define INCLUDED_protocols_rna_RNA_SecStructLegacyInfo_hh


// Project headers
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/CacheableData.hh>

// Utility headers

// Numceric Headers

// C++ headers
#include <string>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace farna {
namespace secstruct {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of RNA centroid information inside the pose.
//// Rhiju move this to its own namespace!
class RNA_SecStructLegacyInfo: public basic::datacache::CacheableData  {

public:

	RNA_SecStructLegacyInfo(){};

	RNA_SecStructLegacyInfo( std::string const & rna_secstruct_legacy_string ):
		rna_secstruct_legacy_( rna_secstruct_legacy_string )
	{}

	RNA_SecStructLegacyInfo( RNA_SecStructLegacyInfo const & src );

	basic::datacache::CacheableDataOP
	clone() const
	{
		return basic::datacache::CacheableDataOP( new RNA_SecStructLegacyInfo( *this ) );
	}

	Size
	size() const {
		return rna_secstruct_legacy_.size();
	}

	void
	set_secstruct( std::string const & secstruct ){ rna_secstruct_legacy_ = secstruct; }

	std::string const &
	get_secstruct() const { return rna_secstruct_legacy_; }

private:

	std::string rna_secstruct_legacy_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

bool
has_rna_secstruct_legacy( core::pose::Pose const & pose );

std::string const &
get_rna_secstruct_legacy( core::pose::Pose & pose );

std::string const &
get_rna_secstruct_legacy_from_const_pose( core::pose::Pose const & pose );

void
set_rna_secstruct_legacy( core::pose::Pose & pose, std::string const & rna_secstruct_legacy_string ); //By default, unknown, actually.

void
clear_rna_secstruct_legacy_info( core::pose::Pose & pose );

} //secstruct
} //farna
} //protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_farna_secstruct_legacy_RNA_SecStructLegacyInfo )
#endif // SERIALIZATION


#endif
