// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//////////////////////////////////////////////////////////////////////
/// @file
///
/// @brief
/// A general mover class for performing BCL mutates on Rosetta small molecule residues
///
/// @details
/// This class makes a mover for BCL small molecule drug design mutates.
/// The design leverages the BCL util::Implementation< template> class, which allows
/// the templated argument type to be constructed from an object data label, which
/// itself can be generated from a string. This is essentially how users interact
/// with the BCL from the command-line: they pass serializable arguments to the BCL
/// that are interpreted as object data labels for construction of objects within the
/// the specified application.
///
/// @author
/// Benjamin P. Brown (benjamin.p.brown17@gmail.com)
////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_protocols_drug_design_bcl_BCLFragmentMutateMover_HH
#define INCLUDED_protocols_drug_design_bcl_BCLFragmentMutateMover_HH

// unit forward includes
#include <protocols/drug_design/bcl/BCLFragmentMutateMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/drug_design/bcl/BCLFragmentBaseMover.hh>

// core includes
#include <core/chemical/bcl/BCLFragmentHandler.hh>
#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/AtomRefMapping.hh>

// package includes
#include <core/chemical/AtomRefMapping.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// utility includes
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// citation manager
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

// BCL includes
#ifdef USEBCL
#include <bcl/include/chemistry/bcl_chemistry_atoms_complete_standardizer.h>
#include <bcl/include/chemistry/bcl_chemistry_bond_isometry_handler.h>
#include <bcl/include/chemistry/bcl_chemistry_fragment_complete.h>
#include <bcl/include/chemistry/bcl_chemistry_fragment_mutate_interface.h>
#include <bcl/include/chemistry/bcl_chemistry_stereocenters_handler.h>
#include <bcl/include/math/bcl_math_mutate_interface.h>
#include <bcl/include/sdf/bcl_sdf.h>
#include <bcl/include/sdf/bcl_sdf_atom_info.h>
#include <bcl/include/sdf/bcl_sdf_bond_info.h>
#include <bcl/include/util/bcl_util_implementation.h>
#include <bcl/include/util/bcl_util_loggers.h>
#endif

namespace protocols
{
namespace drug_design
{
namespace bcl
{
class BCLFragmentMutateMover :
	public protocols::drug_design::bcl::BCLFragmentBaseMover
{
public:

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	/// @brief default constructor
	BCLFragmentMutateMover();

#ifdef USEBCL
	//! @brief default sample confs constructor
	BCLFragmentMutateMover
	(
		core::Size const &n_max_mutate_attempts,
		::bcl::util::Implementation< ::bcl::chemistry::FragmentMutateInterface> const &mutate
	);
#endif

	/// @brief create copy constructor
	protocols::moves::MoverOP clone() const override;

	/// @brief create this type of objectt
	protocols::moves::MoverOP fresh_instance() const override;

	/////////////////
	// data access //
	/////////////////

	static std::string mover_name();

	std::string get_name() const override;

	std::string get_tag( Pose & pose, core::Size resid ) const;

#ifdef USEBCL
	std::string const &get_mutate_label() const;
#else
	void get_mutate_label() const;
#endif

	////////////////
	// operations //
	////////////////

	//! @brief apply BCLFragmentMutateMover
	void apply( core::pose::Pose & pose ) override;

	//! @brief set the maximum number of mutate attempts
	void set_n_max_mutate_attempts( core::Size const &n_max_mutate_attempts);

	//! @brief set the mutate using a BCL object data label string
	void set_mutate( std::string const &object_data_label);

	//////////////////////
	// helper functions //
	//////////////////////

	//! @brief parse xml file
	void parse_my_tag( TagCOP tag, basic::datacache::DataMap &data ) override;

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	//! @brief Provide citations to the passed CitationCollectionList.
	//! @details This mover is unpublished. Hopefully we can publish it soon and I can change
	//! this to a citation for an actual manuscript
	void provide_citation_info( basic::citation_manager::CitationCollectionList &citations) const override;

private:

	//////////
	// data //
	//////////

	//! maximum number of times to retry a particular mutate if it fails.
	//! if the mutate fails for any reason, e.g. the mutated molecule failed the
	//! druglikeness filter or 3D coordinates could not be generated, then the mutate
	//! returns a null pointer, so we will just return the initial fragment
	core::Size max_mutate_retry_;

#ifdef USEBCL
	//! the mutate interface the object will be handling, wrapped in an implementation for generalizability
	::bcl::util::Implementation< ::bcl::chemistry::FragmentMutateInterface > mutate_;
#endif

}; // BCLFragmentMutateMover
} // bcl
} // protocols
} // drug_design

#endif // INCLUDED_protocols_drug_design_bcl_BCLFragmentMutateMover_HH
