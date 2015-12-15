// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// @file:   core/scoring/facts/FACTSPose.hh
// @brief:  Header-file for class declartions for the FACTS algorithm
//          FACTS: Fast Analytical Continuum Treatment of Solvation by URS HABERTHUR and AMEDEO CAFLISCH
// @author: Hahnbeom Park

#ifndef INCLUDED_core_scoring_facts_FACTSPose_HH
#define INCLUDED_core_scoring_facts_FACTSPose_HH

// Unit headers
#include <core/scoring/facts/FACTSPotential.fwd.hh>
#include <core/scoring/facts/FACTSResidue.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/RotamerSetBase.hh>

#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

////////////////////////////////////////////////////////////////////////////////////////////////////
class FACTSPoseInfo : public basic::datacache::CacheableData {

public:
	typedef conformation::Residue   Residue;
	typedef conformation::ResidueOP ResidueOP;

public:

	FACTSPoseInfo();
	FACTSPoseInfo( FACTSPoseInfo const & src );

	basic::datacache::CacheableDataOP clone() const
	{ return basic::datacache::CacheableDataOP( new FACTSPoseInfo( *this ) );}

	Size size() const { return residue_info_.size(); }

	FACTSResidueInfo & residue_info( Size const i )
	{ return *residue_info_[i]; }

	FACTSResidueInfo const & residue_info( Size const i ) const
	{ return *residue_info_[i]; }

	bool being_packed( Size const seqpos ) const
	{ return being_packed_[ seqpos ]; }

	void set_placeholder( Size const i, ResidueOP rsd, FACTSResidueInfoOP info );

	FACTSResidueInfo const & placeholder_info( Size const seqpos ) const
	{
		debug_assert( placeholder_info_[ seqpos ] );
		return *placeholder_info_[ seqpos ];
	}

	Residue const & placeholder_residue( Size const seqpos ) const
	{
		debug_assert( placeholder_residue_[ seqpos ] );
		return *placeholder_residue_[ seqpos ];
	}

	void initialize( pose::Pose const & pose, FACTSRsdTypeMap &rsdtypemap );

	void set_repack_list( utility::vector1< bool > const & repacking_residues );

	bool is_changed( pose::Pose const &pose );

	void update_enumeration_shell( pose::Pose const &pose,
		bool const enumerate_second_shell = false );

	inline bool context_derivative_empty() { return context_derivative_empty_; }

public:
	utility::vector1< FACTSResidueInfoOP > residue_info_; // these are allocated in initialize
	utility::vector1< ResidueOP > placeholder_residue_; // these may be null pointers
	utility::vector1< FACTSResidueInfoOP > placeholder_info_;
	utility::vector1< bool > being_packed_; // stores info from the packertask when setup_for_packing calls set_repack_list
	bool context_derivative_empty_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // scoring
} // core
#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_facts_FACTSPose )
#endif // SERIALIZATION


#endif
