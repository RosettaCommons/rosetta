// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StrandFragment.hh
///
/// @brief Small helper class that stores the start and end of a strand secondary structure

/// @author Tim jacobs

#ifndef INCLUDED_protocols_features_strand_assembly_STRANDFRAGMENT_HH
#define INCLUDED_protocols_features_strand_assembly_STRANDFRAGMENT_HH

//Core
#include <core/types.hh>
#include <core/pose/Pose.hh>

//Utility
#include <utility/vector1.fwd.hh>

//Devel
//#include <devel/strand_assembly/NativeResidue.hh>

//External
//#include <boost/serialization/access.hpp>
//#include <boost/serialization/map.hpp>
//#include <boost/serialization/vector.hpp>
//#include <boost/serialization/string.hpp>

//C++ Headers
#include <string>
#include <vector>
#include <map>

namespace protocols {
namespace features {
namespace strand_assembly {

class StrandFragment{

public:

	StrandFragment(core::Size start, core::Size end);

	StrandFragment(core::Size start, core::Size end, bool direction);

	//This really shouldn't exist, but boost uses it for serialization (I couldn't figure out the way around this).
	//Either way, start_ and end_ are read-only members, so a Fragment assembled in this way won't be very useful....
	StrandFragment();

	~StrandFragment();

	core::Size get_end() const;
	std::string get_pdb_source() const;
	core::Size get_start() const;
	core::Size get_size() const;
	bool get_direction() const;
	void set_pdb_source(std::string pdb_source_);
	void set_direction(bool direction);

//	void insertResiduesFromPose(const core::pose::Pose & pose,
//			core::Size start, core::Size end, const core::pose::Pose & this_pose);
//
//	std::string print() const;

private:
	/*
	friend class boost::serialization::access;

	// When the class Archive corresponds to an output archive, the
	// & operator is defined similar to <<.	Likewise, when the class Archive
	// is a type of input archive the & operator is defined similar to >>.
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & start_;
		ar & end_;
		ar & pdb_source_;
		ar & direction_;
//		ar & residue_map_;
//		ar & residue_list_;
	} */

	core::Size start_;
	core::Size end_;
	std::string pdb_source_;
	bool direction_;
//	std::map<core::Size, std::vector<NativeResidue> > residue_map_;
//	std::vector<std::vector<NativeResidue> > residue_list_;

};

//namespace boost { namespace serialization {
//template<class Archive>
//inline void save_construct_data(
//		Archive & ar, const STRANDFragment * t, const unsigned int file_version
//){
//		// save data required to construct instance
//		size_t start(t->get_start());
//		size_t end(t->get_end());
//		ar << start;
//		ar << end;
//}
//
//template<class Archive>
//inline void load_construct_data(
//		Archive & ar, STRANDFragment * t, const unsigned int file_version
//){
//		// retrieve data from archive required to construct new instance
//		size_t start;
//		size_t end;
//		ar >> start;
//		ar >> end;
//		// invoke inplace constructor to initialize instance of my_class
//		::new(t)STRANDFragment(start, end);
//}
//}}

} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif /* STRANDFRAGMENT_HH_ */
