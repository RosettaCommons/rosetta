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



namespace protocols {
namespace antibody {
namespace design {

/// @brief Class that keeps track of the 'native' sequence during design.  Native means that no design has occured.
/// Used for conservative base design to keep it from hopping around residue types.
class NativeAntibodySeq : public basic::datacache::CacheableData {


public:

	NativeAntibodySeq(core::pose::Pose const & pose, AntibodyInfoCOP ab_info);

	NativeAntibodySeq( NativeAntibodySeq const & src);

	virtual ~NativeAntibodySeq();
	
	/// @brief Sets the sequence from the PDB into this class and into the pose.
	void
	set_sequence(core::pose::Pose const & pose);

	/// @brief Set this class to the pose datacache.
	void
	set_to_pose(core::pose::Pose & pose);

	/// @brief Sets the sequence from the CDR of the PDB into this class and into the pose.
	void
	set_from_cdr(core::pose::Pose const & pose, CDRNameEnum cdr);

	/// @brief Get the full pose sequence.
	std::string
	get_sequence(core::pose::Pose const & pose) const;

	virtual basic::datacache::CacheableDataOP
	clone() const;

private:

	NativeAntibodySeq(); // Default constructor is private: Need ab_info to initialize

	AntibodyInfoCOP ab_info_;
	std::map< std::string , core::chemical::AA> seq_;
	std::map< CDRNameEnum, utility::vector1< core::chemical::AA > > cdr_seq_;

};

} //design
} //antibody
} //protocols

#endif // INCLUDED_protocols_antibody_design_NativeAntibodySeq_hh
