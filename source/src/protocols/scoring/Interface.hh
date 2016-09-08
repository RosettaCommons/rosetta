// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Interface - information about the interface between to partners
/// @brief contains the following information:
///  calculate the interface between the two (or other?) partners
///  set the packer task to only pack the interface
/// @author Monica Berrondo


#ifndef INCLUDED_protocols_scoring_Interface_hh
#define INCLUDED_protocols_scoring_Interface_hh

// Unit headers
#include <protocols/scoring/Interface.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh> // replace with .fwd.hh
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Symmetry headers

// Utility headers
#include <ObjexxFCL/FArray1D.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace scoring {

class
	Interface : public utility::pointer::ReferenceCount {
public:

	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::Real Real;


	// default constructor
	Interface() {
		jump_number_ = 1;
		distance_squared_ = 8.0 * 8.0;
		use_input_partners_ = false;
	}

	// constructor with arguments
	Interface(
		Size const jump_number_in
	) : jump_number_( jump_number_in )
	{
		// set the default distance to 8 Angstrom (for cendist)
		// is this the correct number for a default?
		distance_squared_ = 8.0 * 8.0;
		use_input_partners_ = false;
	}

	Interface(
		ObjexxFCL::FArray1D_bool partner
	)
	{
		partner_ = partner;
		// set the default distance to 8 Angstrom (for cendist)
		// is this the correct number for a default?
		distance_squared_ = 8.0 * 8.0;
		use_input_partners_ = true;
	}

	void jump( Size const jump_number);

	core::Size jump_id() const {
		return jump_number_;
	}

	void calculate( core::pose::Pose const & pose );
	void print( core::pose::Pose const & pose );
	void show( core::pose::Pose const & pose );
	void set_pack( core::pose::Pose const & pose, PackerTaskOP task );

	core::Size closest_interface_residue( core::pose::Pose const & pose, core::Size src_rsd, core::Real & distance );

	void distance( Real const distance_in );

	core::Vector center( core::pose::Pose const & pose );

	bool is_interface( core::conformation::Residue const & rsd ) const;// { return is_interface_( rsd.seqpos() ); }

	bool is_interface( Size const position ) const;// { return is_interface_(position); }

	/// @brief returns the total number of residues on both sides of the interface for
	/// a given jump
	Size interface_nres();

	/// @brief returns whether the rsd1 and rsd2 are considered a contact pair
	/// based on contact_list_ array
	bool is_pair(
		core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2
	);
	// symmetric interfaces
	void set_symmetric_pack( core::pose::Pose const & pose, PackerTaskOP task );

	/// @brief Get two lists containing the residues on each side of the interface
	///
	/// You might expect something called pair_list() to output pairs or somehow
	/// relate to the is_pair function, but it does not.
	/// It returns a two element list.
	///  - The first element is a list of all residues on the lower partner of the
	///    interface.
	///  - The second element is a list of all the residues in the interface on the
	///    upper partner.
	/// Both lists are sorted by residue number, so the indices of one list have
	/// no relationship to the other list.
	/// @todo rename this to something more logical & give it a better return type
	/// @todo make a function that returns a vector1<pair<Size,Size> > containing
	/// all interacting pairs. This would be easy to implement.
	utility::vector1 < utility::vector1_int > pair_list() { return pair_list_; }

	/// @brief The contact_list is pose.size() long, each element
	/// (contact_list_[i]) contains a list of residues from the other partner that
	/// interacts with residue i.  This is calculated in protein_calculate, and
	/// used by Interface.is_pair() function as well as the interchain_vdw,
	/// interchain_env, and interchain_pair scoring components.
	utility::vector1 < utility::vector1_size > contact_list() { return contact_list_; }

private:
	Size jump_number_;
	Real distance_squared_;
	ObjexxFCL::FArray1D_bool partner_;
	ObjexxFCL::FArray1D_bool is_interface_;
	utility::vector1 < utility::vector1_int > pair_list_;
	utility::vector1 < utility::vector1_size > contact_list_;
	//std::vector < bool > is_interface_;
	bool use_input_partners_;

	void protein_calculate( core::pose::Pose const & pose );
	void ligand_calculate( core::pose::Pose const & pose );
	void NA_calculate( core::pose::Pose const & ); //currently set to protein_calculate()
	void symmetric_protein_calculate( core::pose::Pose const & pose );

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace scoring
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_scoring_Interface )
#endif // SERIALIZATION


#endif
