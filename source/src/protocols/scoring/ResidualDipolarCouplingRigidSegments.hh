// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ResidualDipolarCouplingRigidSegments.hh
/// @brief  Creates, stores and manages a list of individual RDC objects for Rigid Segment-based RDC scoring
/// @author Nikolas Sgourakis

#ifndef INCLUDED_protocols_scoring_ResidualDipolarCouplingRigidSegments_hh
#define INCLUDED_protocols_scoring_ResidualDipolarCouplingRigidSegments_hh

#include <protocols/scoring/ResidualDipolarCouplingRigidSegments.fwd.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>

#include <basic/datacache/CacheableData.hh>
#include <numeric/numeric.functions.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ResidualDipolarCoupling.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace scoring {

void store_RDC_segments_in_pose(ResidualDipolarCouplingRigidSegmentsOP, core::pose::Pose&);
ResidualDipolarCouplingRigidSegmentsOP retrieve_RDC_segments_from_pose(core::pose::Pose&);
ResidualDipolarCouplingRigidSegmentsCOP retrieve_RDC_segments_from_pose(core::pose::Pose const&);

/// @brief ResidualDipolarCouplingRigidSegmentss are mainly handled by this class
/// @detail related classed: RDC --- a single line in an RDC file - representing a single dipolar coupling
///                         ResidualDipolarCouplingRigidSegmentsEnergy -- an energy method which triggers computations handled by this class.
///
///
class ResidualDipolarCouplingRigidSegments: public basic::datacache::CacheableData {
	// friend class ResidualDipolarCouplingRigidSegmentsEnergy;
public:
	// typedefs
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::scoring::ResidualDipolarCoupling ResidualDipolarCoupling;
	typedef core::scoring::ResidualDipolarCoupling::RDC_lines RDC_lines;
	typedef utility::vector1< core::scoring::ResidualDipolarCouplingOP > RDC_Segments;

public:
	/// @brief standard c'stor -- will access option -in:file:rdc to read RDC data
	//  ResidualDipolarCouplingRigidSegments() :
	//   nex_(0), nrows_(0) {
	//   reserve_buffers();
	//   read_RDC_file();
	//  }

	ResidualDipolarCouplingRigidSegments() {
		//  std::cout << "Calling Constructor" << std::endl;
		read_RDC_segment_file_from_cmdline();
		sort_into_segments( read_RDCs_from_cmdline() );
	}

	/// @brief alternative c'stor if you have a list of RDC lines
	/*  ResidualDipolarCouplingRigidSegments( RDC_lines data_in, Loops segments ) :
	All_RDC_lines_(data_in) {
	//  preprocess_data();
	//  reserve_buffers();
	}*/

	//explicit copy c'stor to initialize buffers
	//  ResidualDipolarCouplingRigidSegments(ResidualDipolarCouplingRigidSegments const& other);

	//explicit assignment operator to initialize buffers
	// ResidualDipolarCouplingRigidSegments& operator=(ResidualDipolarCouplingRigidSegments const & other);

	//explicit destructor because we use raw pointers for buffers
	// virtual ~ResidualDipolarCouplingRigidSegments();

	//this class lives in the PoseCache.... need to provide clone()
	basic::datacache::CacheableDataOP clone() const {
		return basic::datacache::CacheableDataOP( new ResidualDipolarCouplingRigidSegments(*this) );
	}

	/// @brief compute dipolar score for given segment definition
	/// alignment tensor optimization will be performed for each segment individually
	core::Real compute_total_score(core::pose::Pose const& pose)const;
	core::Real compute_pairwise_score() const; ///total score must have been evaluated before calls to this method are made.
	/// @brief read RDC data from file
	//  void read_RDC_file();

	// do you need accessor for individual Tensors... do it like this
	//Tensor tensor_of_segment( Size i );

	void show(std::ostream&) const;

private:
	/// @brief read RDC data from file
	void sort_into_segments(RDC_lines all_rdcs);
	RDC_lines read_RDCs_from_cmdline() const;
	void read_RDC_segment_file_from_cmdline();
	void read_RDC_segment_file(std::string const& );
	Size find_segid_from_RDC_line(core::scoring::RDC const& line) const;
	Size find_effective_plane(core::scoring::RDC const& line) const;

private:
	RDC_Segments rdc_segments_;
	protocols::loops::Loops segment_definitions_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

extern std::ostream& operator<<(std::ostream&, ResidualDipolarCouplingRigidSegments const&);

} //scoring
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_scoring_ResidualDipolarCouplingRigidSegments )
#endif // SERIALIZATION


#endif
