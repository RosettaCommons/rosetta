// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/architects/StrandArchitect.hh
/// @brief Architect that creates a beta strand
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_architects_StrandArchitect_hh
#define INCLUDED_protocols_denovo_design_architects_StrandArchitect_hh

// Unit headers
#include <protocols/denovo_design/architects/StrandArchitect.fwd.hh>
#include <protocols/denovo_design/architects/DeNovoArchitect.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <set>

namespace protocols {
namespace denovo_design {
namespace architects {

struct PairedStrandNames {
public:
	PairedStrandNames( std::string const & seg1, std::string const & seg2 ):
		segment1( seg1 ), segment2( seg2 ) {}

	std::string segment1;
	std::string segment2;

	bool
	operator<( PairedStrandNames const & other ) const;

	friend std::ostream &
	operator<<( std::ostream & os, PairedStrandNames const & paired_strands );

private:
	PairedStrandNames() {}
};

/// @brief Individual strands are oriented pointing either "UP" or "DOWN"
///        If two adjacent strands have the same orientation, they are parallel
///        If two adjacent strands have different orientation, they are antiparallel
enum StrandOrientation {
	UP = 1,
	DOWN = 2,
	ORIENTATIONS_END = 3
};
typedef utility::vector1< StrandOrientation > StrandOrientations;

typedef long int RegisterShift;
typedef utility::vector1< RegisterShift > RegisterShifts;

typedef long int StrandBulge;
typedef utility::vector1< StrandBulge > StrandBulges;

///@brief Architect that creates a beta strand
class StrandArchitect : public protocols::denovo_design::architects::DeNovoArchitect {
public:
	typedef std::set< core::Size > LengthSet;
	typedef protocols::denovo_design::architects::DeNovoArchitect DeNovoArchitect;
	typedef protocols::denovo_design::architects::DeNovoArchitectOP DeNovoArchitectOP;
	typedef protocols::denovo_design::components::StructureData StructureData;
	typedef protocols::denovo_design::components::StructureDataOP StructureDataOP;

public:
	StrandArchitect( std::string const & id_value );

	virtual ~StrandArchitect();

	static std::string
	class_name() { return "StrandArchitect"; }

	virtual std::string
	type() const;

	DeNovoArchitectOP
	clone() const;

	virtual StructureDataOP
	design( core::pose::Pose const & pose, core::Real & random ) const;

protected:
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

public:
	void
	set_paired_strands( PairedStrandNames const & strands );

	components::StructureDataCOPs::const_iterator
	motifs_begin() const;

	components::StructureDataCOPs::const_iterator
	motifs_end() const;

	void
	set_length( std::string const & length_str );

	void
	set_length( Lengths const & lengths_val );

	void
	set_orientation( std::string const & orientations_str );

	void
	set_orientation( StrandOrientations const & orientations );

	void
	set_register_shift( std::string const & register_shift_str );

	void
	set_register_shift( RegisterShifts const & register_shifts );

	void
	set_bulge( std::string const & bulges_str );

	void
	set_bulge( StrandBulges const & bulges );

	void
	enumerate_permutations();

public:
	// Data field names
	static std::string const
	paired_strands_keyname();

	static std::string const
	register_shift_keyname();

	static std::string const
	orientation_keyname();

	static std::string const
	bulge_keyname();

	static StrandOrientation
	int_to_orientation( int const integer );

	static PairedStrandNames
	str_to_paired_strands( std::string const & paired_strand_str );

public:
	/// interaction with StructureData
	RegisterShift
	retrieve_register_shift( StructureData const & sd ) const;

	StrandOrientation
	retrieve_orientation( StructureData const & sd ) const;

	StrandBulge
	retrieve_bulge( StructureData const & sd ) const;

	PairedStrandNames
	retrieve_paired_strands( StructureData const & sd ) const;

	void
	store_paired_strands( StructureData & sd, PairedStrandNames const & strands ) const;

private:
	void
	store_register_shift( StructureData & sd, RegisterShift const shift ) const;

	void
	store_orientation( StructureData & sd, StrandOrientation const orient ) const;

	void
	store_bulge( StructureData & sd, StrandBulge const bulge ) const;

private:
	void
	needs_update();

	components::StructureDataCOPs
	compute_permutations() const;

private:
	components::StructureDataCOPs motifs_;
	Lengths lengths_;
	StrandOrientations orientations_;
	RegisterShifts register_shifts_;
	PairedStrandNames paired_strands_;
	StrandBulges bulges_;
	bool updated_;

};

} //protocols
} //denovo_design
} //architects

#endif //INCLUDED_protocols_denovo_design_architects_StrandArchitect_hh

