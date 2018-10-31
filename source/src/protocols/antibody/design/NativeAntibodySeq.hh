// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/design/MutateFrameworkForCluster.hh
/// @brief Mutates Framework regions after insertion of a particular cluster
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_design_NativeAntibodySeq_hh
#define INCLUDED_protocols_antibody_design_NativeAntibodySeq_hh

#include <protocols/antibody/design/NativeAntibodySeq.fwd.hh>

#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/chemical/AA.hh>

// Utility headers
#include <utility/vector1.hh>
#include <core/types.hh>
#include <basic/datacache/CacheableData.hh>
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace antibody {
namespace design {

/// @brief Class that keeps track of the 'native' sequence during design.  Native means that no design has occured.
/// Used for conservative base design to keep it from hopping around residue types.
class NativeAntibodySeq : public basic::datacache::CacheableData {


public:



	NativeAntibodySeq(core::pose::Pose const & pose, AntibodyInfo const & ab_info);

	NativeAntibodySeq( NativeAntibodySeq const & src);

	virtual ~NativeAntibodySeq();

	/// @brief Sets the sequence from the PDB into this class and into the pose.
	void
	set_sequence(core::pose::Pose const & pose, AntibodyInfo const & ab_info );

	/// @brief Set this class to the pose datacache.
	void
	set_to_pose(core::pose::Pose & pose);

	/// @brief Sets the sequence from the CDR of the PDB into this class and into the pose.
	void
	set_from_cdr(core::pose::Pose const & pose, AntibodyInfo const & ab_info, CDRNameEnum cdr);

public:

	///@brief Get a map of PDBInfo residues to AA residues as the sequence
	std::map< std::string, core::chemical::AA> const &
	get_full_sequence() const ;

	///@brief Get a map of each native CDR sequence
	std::map< CDRNameEnum, utility::vector1< core::chemical::AA>> const &
	get_cdr_sequence() const ;

	/// @brief Get the full pose sequence, accounting for length changes to the pose.
	std::string
	get_native_sequence_matching_current_length(core::pose::Pose const & pose, AntibodyInfo const & ab_info ) const;

	virtual basic::datacache::CacheableDataOP
	clone() const;


private:

	std::map< std::string , core::chemical::AA> seq_;
	std::map< CDRNameEnum, utility::vector1< core::chemical::AA > > cdr_seq_;

#ifdef    SERIALIZATION
public:
	///@brief Purely for Serialization.
	NativeAntibodySeq();
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //design
} //antibody
} //protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_antibody_design_NativeAntibodySeq )
#endif // SERIALIZATION

#endif // INCLUDED_protocols_antibody_design_NativeAntibodySeq_hh
