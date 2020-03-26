// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/PDBPoseMap.hh
/// @brief  class to allow querying for pose resid with pdb chain/seqpos
/// @author Steven Lewis
/// @author Yih-En Andrew Ban (yab@u.washington.edu)


#ifndef INCLUDED_core_pose_PDBPoseMap_hh
#define INCLUDED_core_pose_PDBPoseMap_hh


// Unit headers
#include <core/pose/PDBPoseMap.fwd.hh>

// Package headers
#include <core/pose/PDBInfo.fwd.hh>

// Project headers

// Utility headers
#include <core/types.hh>
#include <utility/VirtualBase.hh>

// C++ Headers
#include <map>
#include <string>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {




/// @brief PDBPoseMap can be queried with PDB information (chain, sequence position)
///  and returns a pose's resid position.  Useful for handing input/output in terms
///  of PDB positions.  Can be tucked into the pose for repeated access, or generated
///  just-in-time for a single use.  Basically a wrapper class for std::map.
/// @warning danpf@uw.edu  -- In order to avoid breaking workflows I added a secondary map
///  that tracks the map<ResidueKey, Size> without the segmentID field being checked. This is
///  sort of a strange if you use segmentIDs, but if you don't you won't notice any differences.
class PDBPoseMap : public utility::VirtualBase {


public: // typedefs

	typedef core::Size Size;


private: // forward declarations

	struct ResidueKey;


private: // typedefs

	typedef utility::VirtualBase Super;


private: // structs

	/// @brief sortable residue key internal to PDBPoseMap
	struct ResidueKey {
		/// @brief default constructor
		ResidueKey () :
			chainID(' '),
			resSeq(0),
			iCode(' '),
			segmentID("    ")
		{}

		/// @brief value constructor
		ResidueKey( char const c, int const r, char const i, std::string const & s = "    " ) :
			chainID( c ),
			resSeq( r ),
			iCode( i ),
			segmentID(s)
		{}

		/// @brief comparator, lexicographic ordering
		bool
		operator <( ResidueKey const & rval ) const
		{
			return std::tie(chainID, resSeq, iCode, segmentID) < std::tie(rval.chainID, rval.resSeq, rval.iCode, rval.segmentID);
		}

		/// @brief chain id
		char chainID;
		/// @brief residue sequence number
		int resSeq;
		/// @brief insertion code
		char iCode;
		std::string segmentID = "    ";


#ifdef    SERIALIZATION
		template< class Archive > void save( Archive & arc ) const;
		template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION
	};
	struct no_segID_comp {
		bool operator()(ResidueKey const & a, ResidueKey const & b) const {
			return std::tie(a.chainID, a.resSeq, a.iCode) < std::tie(b.chainID, b.resSeq, b.iCode);
		}
	};


public: // constructors

	/// @brief default constructor
	PDBPoseMap();

	/// @brief PDBInfo constructor
	PDBPoseMap( PDBInfo const & info );

	/// @brief copy constructor
	PDBPoseMap( PDBPoseMap const & map );


public: // destructor

	/// @brief default destructor
	~PDBPoseMap() override;


public: // assignment

	/// @brief copy assignment
	PDBPoseMap &
	operator =( PDBPoseMap const & m );


public: // methods

	/// @brief number of mappings
	inline
	Size size() const
	{
		return pdb2pose_.size();
	}

	/// @brief lookup pose numbering
	/// @param[in] chain  chain id
	/// @param[in] pdb_res  pdb residue numbering
	/// @param[in] ins_code  insertion code
	/// @return pose numbering for residue, returns 0 if not found
	/// @warning this will default to return the no-segmentID version if we can't find it!
	Size
	find(
		char const chain,
		int const pdb_res,
		char const ins_code = ' ',
		std::string const & segmentID = "    "
	) const;

	/// @brief insert pdb -> pose number mapping
	/// @param[in] chain  chain id
	/// @param[in] pdb_res  pdb residue numbering
	/// @param[in] ins_code insertion code, use ' ' if no insertion code
	/// @param[in] pose_res  pose numbering for residue
	/// @remarks if the chain is equal to the PDBInfo's empty record character,
	///  the insertion will be skipped
	void
	insert(
		char const chain,
		int const pdb_res,
		char const ins_code,
		std::string const & segmentID,
		Size const pose_res
	);

	/// @brief remove mapping for pdb residue key only if Pose residue matches
	/// @param[in] chain  chain id
	/// @param[in] pdb_res  pdb residue numbering
	/// @param[in] ins_code insertion code, use ' ' if no insertion code
	/// @param[in] pose_res the mapped Pose residue
	/// @return true if key-value pair erase, false otherwise
	bool
	conditional_erase(
		char const chain,
		int const pdb_res,
		char const ins_code,
		std::string const & segmentID,
		Size const pose_res
	);

	/// @brief forcibly remove mapping for pdb residue key
	/// @param[in] chain  chain id
	/// @param[in] pdb_res  pdb residue numbering
	/// @param[in] ins_code insertion code, use ' ' if no insertion code
	void
	erase(
		char const chain,
		int const pdb_res,
		char const ins_code,
		std::string const & segmentID = "    "
	);

	/// @brief clear the current mapping data
	inline
	void
	clear()
	{
		pdb2pose_.clear();
		pdb2pose_noSegmentIDs_.clear();
	}

	/// @brief fill with corresponding pdb -> pose residue mapping
	/// @note does not clear any currently existing mapping data
	void
	fill( PDBInfo const & info );

private: // data

	/// @brief maps ResidueKey -> Pose internal numbering
	std::map< ResidueKey, Size > pdb2pose_;
	std::map< ResidueKey, Size, no_segID_comp> pdb2pose_noSegmentIDs_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //end class PDBPoseMap

} // namespace pose
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pose_PDBPoseMap )
#endif // SERIALIZATION


#endif
