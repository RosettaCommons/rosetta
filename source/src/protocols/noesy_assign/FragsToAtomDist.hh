// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragstoAtomDist.fwd.hh
/// @brief  simulate fragments and recored atom distances
/// @author Zaiyong Zhang (zaiyong.zhang@tum.de)
/// @date   Thr NOV 22 13:22:31 2012


#ifndef INCLUDED_protocols_noesy_assign_FragsToAtomDist_HH
#define INCLUDED_protocols_noesy_assign_FragsToAtomDist_HH

//unit headers
#include <protocols/noesy_assign/FragsToAtomDist.fwd.hh>

//package headers

//project headers
#include <core/types.hh>
#include <core/fragment/FragSet.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

//C++ headers
#include <map>
#include <cmath>

namespace protocols {
namespace noesy_assign {

class FragsToAtomDist : public utility::pointer::ReferenceCount {
public:
	class DistanceRecord {
	public:
		DistanceRecord()
			: cum_dist6_( 0 ), cum_dist_( 0 ), min_dist_( 1000 ), count_( 0 ) {};
		DistanceRecord( core::Real dist6, core::Real dist, core::Real min_dist ,core::Size count ) :
			cum_dist6_( std::pow( dist6, -6.0 ) ), cum_dist_( dist ), min_dist_( min_dist ), count_( count ) {};

		void record( core::Real );
		void record_inv6( core::Real );
		bool is_valid() const { return count_; }

		core::Real average_dist6() const { return std::pow( cum_dist6_ / count_, -1.0/6.0 ); }
		core::Real average_dist() const { return cum_dist_ / count_; }
		core::Real min_dist() const { return min_dist_; }
		core::Real popular_bin() const ;
		utility::vector1< core::Real > dist_track() const { return dist_track_ ; }
		//		std::map< core::Real, core::Size > hist_dist() const { return hist_dist_ };

	private:
		core::Real cum_dist6_;
		core::Real cum_dist_;
		core::Real min_dist_;
		core::Size count_;
		utility::vector1< core::Real > dist_track_;
		//		std::map< core::Real, core::Size > hist_dist_;
		//		std::map< core::Real, core::Size > hist_dist6_;
	};

	typedef std::map< core::id::NamedAtomID, DistanceRecord > NamedInnerMap;
	typedef std::map< core::id::NamedAtomID, NamedInnerMap > NamedDistanceMap;
	static DistanceRecord NO_CONTACT;
	//( 100.0, 100.0, 100.0, 1 );
public:
	FragsToAtomDist() {};
	FragsToAtomDist( std::string const& filename ) {
		read_from_file( filename );
	}
	// read/write from stream
	void read_from_stream( std::istream& );
	void write_to_stream( std::ostream& ) const;
	void write_hist_to_stream( std::ostream& ) const;
	// wrapper for read/write from file
	void read_from_file( std::string const& filename ); // open stream call read_from_stream
	void write_to_file( std::string const& filename ) const; //open output file, call write_to_stream
	void write_hist_to_file( std::string const& filename ) const;
	// generate distance data from fragments
	void generate_from_fragments(
			 core::fragment::FragSetOP fragments,
			 std::string const& sequence,
			 //	 bool r6_averaging,
			 core::Size cycles,
			 core::Size dump_freq
	);

	// generate distance data from fragments
	void generate_from_frag_file(
			 std::string const& filename,
			 std::string const& sequence,
			 //			 bool r6_averaging,
			 core::Size cycles,
			 core::Size dump_freq
	); //read fragments, call generate_from_fragments

	// query distance between two atoms
	FragsToAtomDist::DistanceRecord const& distance_record( core::id::NamedAtomID atom1, core::id::NamedAtomID atom2 ) const;
	core::Real distance(
		 core::id::NamedAtomID atom1,
		 core::id::NamedAtomID atom2,
		 bool r6_averaged = true
	) const {
		return r6_averaged
			? distance_record( atom1, atom2 ).average_dist6()
			: distance_record( atom1, atom2 ).average_dist();
	}


	// compare sequence of fragment data
	bool check_sequence( std::string const& other_sequence ) const;
	std::string const& sequence() const { return sequence_; };

	// has r6_everaging been used to prudecd
	//	bool r6_averaged() const { return r6_averaged_; }


private:
	// small helper function: swap atoms to have always the smaller residue number/atom number first.
	void swap_atoms( core::id::NamedAtomID& atom1, core::id::NamedAtomID& atom2) const;

	void set_sequence( std::string const& sequence ) {
		sequence_=sequence;
	}
	//	void set_r6_averaged( bool r6_averaged ) {
	//		r6_averaged_=r6_averaged;
	//	}

	//main helper function -- computes the distance map
	void compute_average_distances( core::Size ,core::Size);

	//fill vector with natoms() of each residue
	//	void initialize_atom_number_of_each_residue( Size, utility::vector1< core::Size >&, core::pose::Pose const& );
	//	void initialize_maps(Size, DistanceMap&, DistanceMap&, utility::vector1< core::Size >, core::pose::Pose const& );
	//	void initialize_atomid_all(Size, AtomIDMap&, utility::vector1< core::Size > );

	//	void average_distmap( DistanceMap&,DistanceMap );

	std::string sequence_;
	//	bool r6_averaged_;
	core::fragment::FragSetOP frags_;
	NamedDistanceMap named_distmap_;

};// class FragstoAtomDist

std::ostream& operator<< ( std::ostream&, FragsToAtomDist::DistanceRecord const& );
std::istream& operator>> ( std::istream&, FragsToAtomDist::DistanceRecord const& );

}// namespace noesy_assign
}// namespace protocols
#endif
