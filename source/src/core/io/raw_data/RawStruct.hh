// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/raw_data/RawStruct.hh
///
/// @brief Struct base class
/// @author James Thompson, Monica Berrondo

#ifndef INCLUDED_core_io_raw_data_RawStruct_hh
#define INCLUDED_core_io_raw_data_RawStruct_hh

// mini headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// AUTO-REMOVED #include <core/io/raw_data/Raw.fwd.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>

#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
// AUTO-REMOVED #include <string>

#include <utility/vector1.hh>
#include <map>
#include <ostream>


namespace core {
namespace io {
namespace raw_data {

	/////////////////////////////////////////////////////////////////////////////
	// holds all the data for a single entry in a silent file
	class RawStruct : public utility::pointer::ReferenceCount {

	public:
		// destructor
		virtual ~RawStruct();

		/// @brief Fill a Pose with the conformation information in this RawStruct and the FA_STANDARD
		/// ResidueTypeSet. This is a virtual method which must be implemented by classes derived from RawStruct.
		virtual void fill_pose(
			core::pose::Pose & pose
		) = 0;

		/// @brief Fill a Pose with the conformation information in this RawStruct and the ResidueTypeSet
		/// provided by the caller. This is a virtual method which must be implemented by classes derived from RawStruct.
		virtual void fill_pose(
			core::pose::Pose & pose,
			core::chemical::ResidueTypeSet const& residue_set
		) = 0;

		/// @brief Do some sort of comparison between the actual RMSD of this silent-struct and
		/// the cached coordinates. Used for RawStruct objects that are rebuild from torsions
		/// or other shortened representations of data.
		virtual Real get_debug_rmsd();

		/// @brief print out a header line to the given ozstream. In a rosetta++ silent-file, this contained the lines
		/// SEQUENCE: <protein sequence>\nSCORE: <list of score-types>.
		virtual void print_header( std::ostream& out, std::map < std::string, core::Real > const & score_map, 
                                 std::map < std::string, std::string > const & string_map = ( std::map < std::string, std::string > () ),
                                 bool print_sequence = true ) const;
		/// @brief print out a SCORE line to the given ozstream.
		virtual void print_scores( std::ostream& out, std::map < std::string, core::Real > const & score_map,
                                 std::map < std::string, std::string > const & string_map = ( std::map < std::string, std::string > () ) ) const;

		virtual void print_conformation( std::ostream& out ) const;

		/// @brief data access methods.
		Size nres() {
			return nres_;
		}

		/// @brief returns the tag associated with this RawStruct
		std::string decoy_tag() {
			return decoy_tag_;
		}

		/// @brief returns the sequence for this RawStruct
		std::string sequence() {
			return sequence_;
		}

		/// @brief returns the number of residues in this RawStruct
		void nres( Size nres ) {
			nres_ = nres;
		}

		/// @brief sets the tag associated with this RawStruct
		void decoy_tag( std::string tag ) {
			decoy_tag_ = tag;
		}

		/// @brief sets the sequence for this RawStruct
		void sequence( std::string sequence ) {
			sequence_ = sequence;
		}

	protected:
		Size nres_;
		std::string decoy_tag_;
		std::string sequence_;

}; // class RawStruct

} // namespace silent
} // namespace io
} // namespace core

#endif
