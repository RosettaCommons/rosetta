// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/raw_data/DecoyStruct.hh
///
/// @brief silent input file reader for mini

#ifndef INCLUDED_core_io_raw_data_DecoyStruct_hh
#define INCLUDED_core_io_raw_data_DecoyStruct_hh

// mini headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>


#include <core/io/raw_data/RawStruct.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>

// Utility headers

// Numeric headers
#include <numeric/xyzVector.hh>

// C++ Headers
#include <string>
#include <map>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace raw_data {

	class DecoyStruct : public RawStruct {
	public:

		// the constructor:
		DecoyStruct( Size const nres_in )
		{
			nres_        = nres_in;
			fullatom_    = false;
			resize( nres_in );
		}

		DecoyStruct()
		{
			nres_        = 0;
			fullatom_    = false;
			decoy_tag_   = "empty";
		}

		DecoyStruct(
			core::pose::Pose const& pose,
			std::string tag = "empty_tag",
			bool fa = false
		);

		// resize everything appropriately
		void resize(
			Size const nres_in
		);

		// destructor
		~DecoyStruct() {}

		//DecoyStruct & operator= (DecoyStruct const & src);

		void fill_pose(
			core::pose::Pose & pose
		);

		void fill_pose(
			core::pose::Pose & pose,
			core::chemical::ResidueTypeSet const& residue_set
		);

		// fills the res array
		//void set_sequence(std::string const & sequence);

		virtual void print_conformation( std::ostream& output ) const;

		/// @brief calculates the RMSD between the C-alpha atoms of a Pose built from the torsions in this
		/// DecoyStruct and the C-alpha atoms from this DecoyStruct.
		virtual Real get_debug_rmsd();

		/// @brief data getters/setters
		bool fullatom() const {
			return fullatom_;
		}

		Real phi( unsigned int seqpos ) const {
			return phi_[seqpos];
		}

		Real psi( unsigned int seqpos ) const {
			return psi_[seqpos];
		}

		Real omega( unsigned int seqpos ) const {
			return omega_[seqpos];
		}

		char secstruct( unsigned int seqpos ) const {
			return secstruct_[seqpos];
		}

		utility::vector1< Real > chi( unsigned int seqpos ) const {
			return chi_[ seqpos ];
		}

		Vector coords( unsigned int seqpos ) const {
			return coords_[seqpos];
		}


		void fullatom( bool fullatom ) {
			fullatom_ = fullatom;
		}

		void phi( unsigned int seqpos, Real phi ) {
			phi_[seqpos] = phi;
		}

		void psi( unsigned int seqpos, Real psi ) {
			psi_[seqpos] = psi;
		}

		void omega( unsigned int seqpos, Real omega ) {
			omega_[seqpos] = omega;
		}

		void secstruct( unsigned int seqpos, char ss ) {
			secstruct_[seqpos] = ss;
		}

		void chi( unsigned int seqpos, utility::vector1< Real > chis ) {
			chi_[seqpos] = chis;
		}

		void coords( unsigned int seqpos, Vector coords ) {
			coords_[seqpos] = coords;
		}


	protected:
		bool fullatom_;

		utility::vector1< char > secstruct_;
		utility::vector1< Real > phi_;
		utility::vector1< Real > psi_;
		utility::vector1< Real > omega_;
		utility::vector1< utility::vector1< Real > > chi_;
		utility::vector1< Vector > coords_;

	}; // class DecoyStruct
} // namespace silent
} // namespace io
} // namespace core

#endif
