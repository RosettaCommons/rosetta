// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loophash/BackboneDB.hh
/// @brief I little in the spirit of fragments but way more memory and speed efficient. Backbones are stored as 16bit fixed comma integers. :)
/// @author Mike Tyka
/// @author Ken Jung


#ifndef INCLUDED_protocols_loophash_BackboneDB_hh
#define INCLUDED_protocols_loophash_BackboneDB_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <boost/unordered_map.hpp>

#include <vector>
#include <map>

#include <protocols/frag_picker/VallChunk.hh>
#include <utility/vector1.hh>
#include <numeric/geometry/BoundingBox.fwd.hh>

//numeric headers
#include <numeric/geometry/hashing/SixDHasher.fwd.hh>


namespace protocols {
namespace loophash {

// surely there's a func like this somewhere ?
template <class T>
T sqr( T a){ return a*a; }

const core::Real MAXIMAL_FLOAT = 100000000.0;

typedef boost::unordered_multimap< boost::uint64_t, core::Size , numeric::geometry::hashing::bin_index_hasher >   BackboneIndexMap;
typedef numeric::geometry::BoundingBox< core::Vector > BoundingBox;
const int HASH_POSITION_GRID_BASE= 5;

short RealAngleToShort( core::Real angle );

core::Real ShortToRealAngle( short angle );


// Structs for storing data
struct BBExtraData{
	std::string pdb_id;
	std::vector < int > rotamer_id;
	std::string sequence;
};
struct BBData {
	std::vector < short > angles;
	core::Size  extra_key;
};

class BackboneSegment {
	public:
		BackboneSegment(){}

		BackboneSegment(
				const std::vector<core::Real> &phi,
				const std::vector<core::Real> &psi,
				const std::vector<core::Real> &omega
				){
			phi_ = phi;
			psi_ = psi;
			omega_ = omega;
		}

		void apply_to_pose( core::pose::Pose &pose, core::Size ir, bool cut  = false ) const;
		void read_from_pose( core::pose::Pose const &pose, core::Size ir, core::Size length );

		void print() const ;

		core::Size length() const { return  phi_.size(); }
		const std::vector<core::Real> &phi()   const { return phi_; }
		const std::vector<core::Real> &psi()   const { return psi_; }
		const std::vector<core::Real> &omega() const { return omega_; }
		
		bool compare(const BackboneSegment &bs1, core::Real tolerance) const;

		bool operator==( const BackboneSegment &bs1 ) const;

		bool operator!=(const BackboneSegment &other) const {
		     return !(*this == other);
		}

		core::Size get_mem_foot_print(){ return phi_.size() * 3 * sizeof( core::Real ) + sizeof( BackboneSegment ); }
	private:
		std::vector<core::Real> phi_;
		std::vector<core::Real> psi_;
		std::vector<core::Real> omega_;

};

core::Real get_rmsd( const BackboneSegment &bs1, const BackboneSegment &bs2 );


class BackboneDB
{
	public:
		BackboneDB(){ extra_ = false; }

		core::Real angle( core::Size index, core::Size offset );

		// Always grab extra data from Vall
		void add_pose( const core::pose::Pose &pose, core::Size nres, core::Size &offset, protocols::frag_picker::VallChunkOP chunk = NULL );

		//obsolete?
		//void add_backbone_segment( const BackboneSegment &bs, core::Size &offset );


		void get_protein( core::Size index, BBData & protein ) const;
		void get_extra_data( core::Size index, BBExtraData & extra ) const;
		void add_protein( BBData new_protein );
		void add_extra_data( BBExtraData extra );

		void get_backbone_segment(
				const core::Size index,
				const core::Size offset,
				const core::Size len,
				BackboneSegment &bs
				) const;

		core::Size size() const { return data_.size(); }
		core::Size extra_size() const { return extra_data_.size(); }

		// Only supports writing to text, extra data is mandatory
		void write_db( std::string filename );

		// Only supports reading from text, extra data must exist in db, loading is optional
		// Returns a range of keys of which protein is read
		// also returns a map of homolog indices
		void read_db( std::string filename, bool load_extra,
				core::Size num_partitions, core::Size assigned_num,
				std::pair< core::Size, core::Size > & loop_range,
				std::map< core::Size, bool > & homolog_index );
		// Hack to give a default parameter to loop_range
		inline void read_db( std::string filename, bool load_extra ) {
			std::pair< core::Size, core::Size > range;
			std::map < core::Size, bool > homolog_index;
			read_db( filename, load_extra, 1, 0, range, homolog_index );
		}

		// For backwards compatability
		void read_legacydb( std::string filename );

		// Reads in homologs into homologs_
		void read_homologs();

		core::Size get_mem_foot_print(){ return data_.size() * sizeof( BBData ) + data_.size() * sizeof( BBExtraData ) + sizeof(BackboneDB); };

	private:

		// Gonna change this to a vector of BBData structs, one for each protein
		// This obviously introduces overhead, but minimal (~2MB), and it makes it easier to store
		// other info as well as sort the database

		// Extra info is stored in a separate struct (BBExtraData) with an index in the BBData struct

		// std::vector < short > data_;
		std::vector < BBData > data_;
		std::vector < BBExtraData > extra_data_;

		// if this flag is true, extra_data_ is populated
		bool extra_;

		/// @brief Homology map, contains all homologs
		std::map< std::string, bool >             homologs_;

};


}
}

#endif

