// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
///
/// @file protocols/seeded_abinitio/PlaceFragments.hh
/// @author Eva-Maria Strauch (evas01@u.washington.edu)

#ifndef INCLUDED_protocols_seeded_abinitio_PlaceFragments_hh
#define INCLUDED_protocols_seeded_abinitio_PlaceFragments_hh

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.hh>
#include <utility/vector1.hh>
#include <boost/unordered/unordered_map.hpp>

namespace protocols {
namespace seeded_abinitio {

class PlaceFragments : public protocols::moves::Mover {
 public:
  typedef core::pose::Pose Pose;

  //PlaceFragments(const FragSetOP& fragments);

  PlaceFragments();

  // undefined, commenting out to fix PyRosetta build  bool is_seed  ( protocols::loops::Loops & loops, core::Size & residue );

  void apply( core::pose::Pose & pose );

  virtual std::string get_name() const;

  void parse_my_tag( utility::tag::TagCOP tag,
                     basic::datacache::DataMap &,
                     protocols::filters::Filters_map const &,
                     protocols::moves::Movers_map const &,
                     core::pose::Pose const & );

  protocols::moves::MoverOP clone() const { return( protocols::moves::MoverOP( new PlaceFragments( *this ) ) ); }
  protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new PlaceFragments ); }

  virtual ~PlaceFragments();

	void initialize_fragments( const core::fragment::FragSetOP & fragments );

	//use privat members...and setters/getters
	void apply_frame(  core::pose::Pose & pose, core::fragment::Frame &frame, int aln_len, core::Size seq_start, core::Size max_frag_len  );
	void create_fragments( core::pose::Pose & pose, core::Size insert_start, core::Size insert_stop );

 private:
	void initialize( const core::fragment::FragSetOP  fragments );

  	/// mover that should be applied after placement
    protocols::moves::MoverOP mover_;

    /// filter that should be applied
    protocols::filters::FilterOP filter_;

    //void initialize_library();
		//core::fragment::FragSetOP fragments9_;

		/// storing fragments
		core::fragment::FragSetOP fragments_;
		boost::unordered_map<core::Size, core::fragment::Frame> library_;

		/// how many fragments to be picked
    core::Size nfrags_;
  	/// length of fragment to be picked on the flight
  	core::Size fsize_;
		/// number of residues to align with fragments
		int cartfrag_overlap_;

		/// stub
    std::string input_stubs_;
    /// seeds/segements
    utility::vector1< std::pair < std::string,std::string > > seg_vector_;

		bool frags_onflight_;

		/// secondary structure sequence
		std::string ss_;

		bool use_seq_;

};

}
}

#endif
