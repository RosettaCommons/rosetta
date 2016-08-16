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

// Core headers
#include <core/io/Remarks.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/SingletonBase.hh>

namespace protocols {
namespace denovo_design {
namespace components {

///@brief Singleton for creating StructureData objects
class StructureDataFactory : public utility::SingletonBase< StructureDataFactory > {
public:
	static StructureDataFactory *
	create_singleton_instance();

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

	/// @brief retrieves a StructureData object from observable cache
	///        creates one if necessary, and sets it in the pose
	StructureData const &
	get_from_pose( core::pose::Pose & pose, std::string const & newid ) const;

	/// @brief retrieves a StructureData object from the pose observable cache
	///        Creates a StructureData using the pose (but doesn't attach it) if
	///        the cached StructureData could not be retrieved properly
	StructureData const &
	get_from_const_pose( core::pose::Pose const & pose ) const;


	/// @brief retrieves a StructureData object from the pose observable cache
	///        Creates a StructureData using the pose (but doesn't attach it) if
	///        the cached StructureData could not be retrieved properly
	StructureDataOP
	create_from_pose( core::pose::Pose const & pose ) const;

	/// @brief retrieves a StructureData object from the pose observable cache
	///        Creates a StructureData using the pose (but doesn't attach it) if
	///        the cached StructureData could not be retrieved properly
	StructureDataOP
	create_from_pose( core::pose::Pose const & pose, std::string const & newid ) const;

	/// @brief creates a StructureData from an xml stringstream
	StructureDataOP
	create_from_xml( std::istream & xmltag ) const;

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

	/// @brief attaches cacheable observer.  Overwrites whatever was there before
	void
	attach_observer( core::pose::Pose & pose, StructureData const & sd ) const;

	/// @brief detaches cacheable observer and removes from observable cache
	void
	detach_observer( core::pose::Pose & pose ) const;

	/// @brief stores a string in the pose's datacache
	bool
	has_cached_string( core::pose::Pose const & pose ) const;

	/// @brief stores a string in the pose's datacache
	std::string
	retrieve_cached_string( core::pose::Pose const & pose ) const;

	/// @brief stores a string in the pose's datacache
	void
	set_cached_string( core::pose::Pose & pose, std::string const & ssorig ) const;

private:
	// Creation methods

	/// @brief sets data from a given pose's pdb remark records. Only looks at remarks which are subcomponents of the given id
	/// if the pdb remarks are not found, the permutation information is inferred from the information that can be gleaned from the pose
	StructureDataOP
	create_new_from_pose( core::pose::Pose const & pose, std::string const & newid ) const;

	/// @brief sets data from a given pose's information, not taking into account PDB remarks
	StructureDataOP
	infer_from_pose( core::pose::Pose const & pose, std::string const & newid ) const;


private:
	// I/O with Remarks

	/// @brief creates SD from description stored in PDB remarks
	StructureDataOP
	create_from_remarks( core::io::Remarks const & rem ) const;

	/// @brief adds a remark to remarks object
	void
	add_perm_remark( core::io::Remarks & remarks, std::string const & rem_value ) const;

private:
	// I/O with cached string

	/// @brief retrieves string from pose's datacache
	std::string
	retrieve_cached_string( core::pose::Pose const & pose, std::string const & data_name ) const;

	/// @brief stores a string in the pose's datacache
	void
	set_cached_string(
		core::pose::Pose & pose,
		std::string const & ssorig,
		std::string const & data_name ) const;

	void
	clear_cached_string( core::pose::Pose & pose ) const;

public:
	// Constants
	static int const REMARK_NUM;
	static std::string const DATA_NAME;

};

/// @brief prepare a string that was stored in string datacache
void
clean_from_storage( std::string & st );

/// @brief prepare a string to be stored in the string datacache
void
clean_for_storage( std::string & ss );

} //protocols
} //denovo_design
} //components

#endif //INCLUDED_protocols_denovo_design_components_StructureDataFactory_hh

