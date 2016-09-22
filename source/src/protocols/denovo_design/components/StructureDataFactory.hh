// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/components/StructureDataFactory.hh
/// @brief Singleton for creating StructureData objects
/// @author Tom Linsky (tlinsky@uw.edu)
#ifndef INCLUDED_protocols_denovo_design_components_StructureDataFactory_hh
#define INCLUDED_protocols_denovo_design_components_StructureDataFactory_hh

// Unit headers
#include <protocols/denovo_design/components/StructureDataFactory.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/denovo_design/components/StructureDataObserver.fwd.hh>
#include <protocols/denovo_design/types.hh>

// Core headers
#include <core/io/Remarks.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <basic/datacache/WriteableCacheableData.fwd.hh>
#include <utility/SingletonBase.hh>
#include <utility/excn/EXCN_Base.hh>

// C++ headers
#include <map>

namespace protocols {
namespace denovo_design {
namespace components {

///@brief Singleton for creating StructureData objects
class StructureDataFactory : public utility::SingletonBase< StructureDataFactory > {
public:
	friend class utility::SingletonBase< StructureDataFactory >;

	typedef std::map< std::string, std::set< basic::datacache::WriteableCacheableDataOP > > WriteableCacheableDataMap;

public:
	/// @brief stores the data of this permutation into a pose for later retrieval
	///        StructureData stored in a pose this way can be retrieved by calling
	///        get_from_pose(), or get_from_const_pose()
	void
	save_into_pose( core::pose::Pose & pose, StructureData const & sd ) const;

	/// @brief clears all structuredata information (including cached xml) from
	///        this pose.
	void
	clear_from_pose( core::pose::Pose & pose ) const;

	/// @brief retrieves a StructureData object from the pose observable cache
	///        creates one if necessary
	StructureData const &
	get_from_pose( core::pose::Pose & pose ) const;

	/// @brief retrieves a StructureData object from the pose observable cache
	///        Creates a StructureData using the pose (but doesn't attach it) if
	///        the cached StructureData could not be retrieved properly
	StructureData const &
	get_from_const_pose( core::pose::Pose const & pose ) const;


	/// @brief retrieves a StructureData object from the pose data cache
	///        Creates a StructureData using the pose (but doesn't store it) if
	///        the cached StructureData is not present
	/// @param[in]  pose    The input pose
	StructureData
	create_from_pose( core::pose::Pose const & pose ) const;

	/// @brief retrieves a StructureData object from the pose data cache
	///        Creates a StructureData using the pose (but doesn't store it) if
	///        the cached StructureData is not present. Adds id_val as a prefix
	///        for all new segments
	/// @param[in]  pose    The input pose
	/// @param[in]  prefix  A prefix to be added to the segments in the pose. For example, if
	///                     prefix is empty, the first helix would be named 'H01', but if prefix
	///                     is set to 'pose1', the first helix will be named 'pose.H01'
	StructureData
	create_from_pose( core::pose::Pose const & pose, SegmentName const & prefix ) const;

	/// @brief creates from CacheableData stream
	StructureData
	create_from_cacheable_data( std::istream & raw_stream ) const;

	/// @brief creates a StructureData from an xml stringstream
	StructureData
	create_from_xml( std::istream & xmltag ) const;

	/// @brief creates SD from description stored in PDB remarks
	StructureData
	create_from_remarks( core::io::Remarks const & rem ) const;

	/// @brief creates a StructureData from a motif string
	/// @param[in] motifs String of secstruct/abego motifs (e.g. "1LX-5EB-2LG-5EB-1LA-1LB-1LG-10HA-1LX" )
	/// @returns StructureData containing one segment per motif with the specified
	///          secondary structure and abego. Each motif is connected to
	///          the previous and next motifs.  Segments are named by their
	///          secondary structure (e.g. L01, E01, L02, E02, L03, H01, L04 )
	StructureData
	create_from_motifs( std::string const & motif_str ) const;

	/// @brief creates a StructureData from a motif string, optionally prefixing a string to the
	///        name of each segment
	/// @param[in] motifs String of secstruct/abego motifs (e.g. "1LX-5EB-2LG-5EB-1LA-1LB-1LG-10HA-1LX" )
	/// @param[in] prefix Name to be prepended before the segment names. If empty,
	///                   names will be as [ L01, E01, L02, ... ].  If set to "myprefix", names
	///          will be as [ myprefix.L01, myprefix.E01, ... ]
	/// @returns StructureData containing one segment per motif with the specified
	///          secondary structure and abego. Each motif is connected to
	///          the previous and next motifs.  Segments are named by their
	///          secondary structure (e.g. L01, E01, L02, E02, L03, H01, L04 )
	StructureData
	create_from_motifs( std::string const & motif_str, SegmentName const & prefix ) const;

public:
	// TO BE MADE PRIVATE EVENTUALLY

	/// @brief checks for a StructureData in the pose observer cache
	///        Returns true if one is found, false otherwise
	bool
	observer_attached( core::pose::Pose const & pose ) const;

	/// @brief returns observer pointer if the pose has one cached
	///        pointer returned can be null if no StructureData is present in the cache
	StructureDataObserverCOP
	retrieve_observer( core::pose::Pose const & pose ) const;

	/// @brief attaches cacheable observer. Sets up cached SD if there isn't one
	void
	attach_observer( core::pose::Pose & pose ) const;

	/// @brief detaches cacheable observer and removes from observable cache
	void
	detach_observer( core::pose::Pose & pose ) const;

	/// @brief checks whether StructureData is stored in the pose's datacache
	bool
	has_cached_data( core::pose::Pose const & pose ) const;

private:
	// Creation methods

	/// @brief Private constructor for singleton.
	StructureDataFactory();

	/// @brief sets data from a given pose's information, not taking into account PDB remarks
	/// @param[in]  pose    The input pose
	/// @param[in]  prefix  A prefix to be added to the segments in the pose. For example, if
	///                     prefix is empty, the first helix would be named 'H01', but if prefix
	///                     is set to 'pose1', the first helix will be named 'pose.H01'
	StructureData
	infer_from_pose( core::pose::Pose const & pose, SegmentName const & prefix ) const;

private:
	// I/O with datacache

	StructureDataOP
	retrieve_cached_data_ptr( core::pose::Pose & pose ) const;

	/// @brief retrieves a StructureData object from the pose's datacache
	StructureData const &
	retrieve_cached_data( core::pose::Pose const & pose ) const;

	/// @brief stores a StructureData in the pose's datacache
	void
	set_cached_data( core::pose::Pose & pose, StructureData const & sd ) const;

	/// @brief clears stored StructureData from the pose's datacache
	void
	clear_cached_data( core::pose::Pose & pose ) const;

private:
	// I/O with Remarks

	/// @brief adds a remark to remarks object
	void
	add_perm_remark( core::io::Remarks & remarks, std::string const & rem_value ) const;

public:
	// Constants
	static int const REMARK_NUM;
};

/// @brief prepare a string that was stored in string datacache
void
clean_from_storage( std::string & st );

/// @brief prepare a string to be stored in the string datacache
void
clean_for_storage( std::string & ss );

class EXCN_RemarksNotPresent : public utility::excn::EXCN_Base {
public:
	EXCN_RemarksNotPresent( std::string const & msg ):
		utility::excn::EXCN_Base(),
		msg_( msg ) {}

	virtual void
	show( std::ostream & os ) const { os << msg_; }

private:
	std::string msg_;
};

class SegmentCounts {
public:
	enum SegmentType {
		HELIX = 1,
		STRAND = 2,
		LOOP = 3,
		LIGAND = 4,
		NUM_SEGMENT_TYPES = 5
	};

	SegmentCounts();

	SegmentCounts( core::pose::Pose const & pose );

	std::string
	new_segment_name( char const ss_type, core::Size const seqpos );

private:
	SegmentType
	type( char const ss_type, core::Size const seqpos ) const;

private:
	utility::vector1< core::Size > counts_;
	core::select::residue_selector::ResidueSubset ligand_subset_;
};

/// @brief Adds segments to the given SD
/// @param[in]     id_val      Parent ID of the segments.  Can be empty
/// @param[in,out] sd          StructureData to be modified
/// @param[in]     chain_ss    Secondary structure of the full chain
/// @param[in]     chain_abego Abego for the full chain
/// @param[in,out] counts      Secondary structure element counters for naming segments
void
add_segments_for_chain(
	std::string const & id_val,
	StructureData & sd,
	std::string const & chain_ss,
	std::string const & chain_abego,
	SegmentCounts & counts );

} //protocols
} //denovo_design
} //components

#endif //INCLUDED_protocols_denovo_design_components_StructureDataFactory_hh

