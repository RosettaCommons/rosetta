// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/fldsgn/MatchResidues.cc
/// @brief  
/// @author Javier Castellanos ( javiercv@uw.edu )

// Unit Headers
#include <protocols/fldsgn/MatchResidues.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/AtomID_Map.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <utility/exit.hh>

// Parser headers
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/moves/Mover.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>


// Boost Headers
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include "boost/assign.hpp"


static thread_local basic::Tracer TR( "protocols.fldsgn.MatchResidues" );

namespace protocols {
namespace fldsgn {

typedef utility::vector1< core::Size > VecSize;
typedef utility::vector1< VecSize > VecVecSize;


// @brief default constructor
MatchResidues::MatchResidues()
{}


// @brief destructor
MatchResidues::~MatchResidues() {}



core::Real 
MatchResidues::compute_comb( core::pose::Pose const & pose, VecSize const & comb ) const
{
  std::map< core::id::AtomID, core::id::AtomID > atom_id_map;
	for(Size i = 1; i < comb.size(); i++) {
		const core::id::AtomID mod_id(pose.residue_type( comb[i] ).atom_index( "CA" ), comb[i] );
		const core::id::AtomID ref_id(pose.residue_type( reference_residues_indexes_[i] ).atom_index( "CA" ), reference_residues_indexes_[i]);
		atom_id_map.insert( std::make_pair(mod_id, ref_id) );
	}
	return  core::scoring::rms_at_all_corresponding_atoms(pose, reference_pose_, atom_id_map); 
}

core::Real 
MatchResidues::superimpose_comb( core::pose::Pose & pose, VecSize const & comb ) const
{
	core::id::AtomID_Map< core::id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose, core::id::BOGUS_ATOM_ID );
	for ( Size i = 1; i < comb.size(); ++i ) {
		const core::id::AtomID mod_id(pose.residue_type( comb[i] ).atom_index( "CA" ), comb[i] );
		const core::id::AtomID ref_id(pose.residue_type( reference_residues_indexes_[i] ).atom_index( "CA" ), reference_residues_indexes_[i]);
		atom_map.set( mod_id, ref_id);

	}
	return  core::scoring::superimpose_pose(pose, reference_pose_, atom_map); 
}

core::Real 
MatchResidues::compute( core::pose::Pose const & pose, VecSize & best_fit ) const
{
	Real smallest_rms = 99999;
	BOOST_FOREACH(VecSize const & comb, mod_segment_prod_){
		Real rms = compute_comb(pose, comb);
		if( rms < smallest_rms ) {
			smallest_rms = rms;
			best_fit = comb;
		}
	}
	TR << "best fit:" <<std::endl;
	TR <<"ref_pos\tmod_pos" <<  std::endl;
	for(Size i = 1; i <= best_fit.size(); i++)// {
		TR << "\t" << reference_residues_indexes_[i] << "\t" << best_fit[i] << std::endl;
	TR << std::endl << "rmsd matched positions " << smallest_rms << std::endl;
	
	return smallest_rms;	
}


void 
MatchResidues::cart_product(
    VecVecSize& rvvi,  // final result
    VecSize&  rvi,   // current result 
    VecVecSize::const_iterator me, // current input
    VecVecSize::const_iterator end)  const // final input
{
    if(me == end) {
        // terminal condition of the recursion. We no longer have
        // any input vectors to manipulate. 
        rvvi.push_back(rvi);
        return;
    }

    const VecSize& mevi = *me;
    for(VecSize::const_iterator it = mevi.begin();
        it != mevi.end();
        it++) {
        rvi.push_back(*it);  
        cart_product(rvvi, rvi, me+1, end);
        rvi.pop_back(); 
    }
}

VecVecSize 
MatchResidues::cart_product( VecVecSize const & input) const {
  VecVecSize tmp;
  VecVecSize output;
  VecSize outputTemp;
	cart_product(tmp, outputTemp, input.begin(), input.end());
	//remove duplicates in combination
	BOOST_FOREACH(VecSize const & vec, tmp) {
		std::set<Size> comb(vec.begin(), vec.end());
		if(comb.size() == vec.size())
			output.push_back(vec);
	}
	return output;
}


/// @brief parse xml
void
MatchResidues::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	const std::string reference = tag->getOption< std::string >("reference");
  core::import_pose::pose_from_pdb( reference_pose_, reference );
	
	VecVecSize mod_match_segments;

  threshold_ = tag->getOption< core::Real >("threshold", 3.0);
	const std::string blueprint = tag->getOption< std::string >("blueprint", "");
	std::map< std::string, boost::tuple< core::Size, core::Size> > ss_seg;
	if(blueprint != "") {
		// get the segment start and end points from the blueprint	
		protocols::jd2::parser::BluePrint blue;
		blue.read_blueprint( blueprint );
		const std::string ss = blue.secstruct();
		ss_seg = map_ss_segments( ss );
  }

	BOOST_FOREACH( utility::tag::TagCOP pairs_tag, tag->getTags() ) {
		if( pairs_tag->getName() == "match" ) {
			const Size ref_pos = pairs_tag->getOption< Size >("ref_pos");
			Size mod_pos_start= 0;
			Size mod_pos_end = 0;
			if( pairs_tag->hasOption("segment") ) {
				assert( blueprint != "");
				boost::tuple< Size, Size> range = ss_seg[ pairs_tag->getOption< std::string >("segment") ];
				mod_pos_start = boost::get<0>(range);
				mod_pos_end = boost::get<1>(range);

			} else if( pairs_tag->hasOption("mod_pos_start") || pairs_tag->hasOption("mod_pos_end") ) {
				mod_pos_start = pairs_tag->getOption< Size >("mod_pos_start");
				mod_pos_end = pairs_tag->getOption< Size >("mod_pos_end");

			} else if( pairs_tag->hasOption("mod_pos") ) {
				mod_pos_start = mod_pos_end = pairs_tag->getOption< Size >("mod_pos");

			} else {
				utility_exit_with_message("need to specify a segment to be matched in the modified structure");
			}
			VecSize vrange;
			for(Size n = mod_pos_start; n <= mod_pos_end; n++)
				vrange.push_back(n);
			
			TR << "Residue " << ref_pos << " from the reference structure matches against segment [" << mod_pos_start << "," << mod_pos_end <<"]" << std::endl;
			reference_residues_indexes_.push_back( ref_pos );
			mod_match_segments.push_back( vrange );

		} else {
			utility_exit_with_message("wrong tag \"" + pairs_tag->getName() + "\" used in MatchResidues.");
		}
	}
	mod_segment_prod_ = cart_product( mod_match_segments );
}


std::map< std::string, boost::tuple< core::Size, core::Size> > 
MatchResidues::map_ss_segments( std::string const & ss) const
{
	std::map< std::string, boost::tuple<Size, Size> > ss_map;
	std::map< char, Size > ss_seg_count = boost::assign::map_list_of ( 'L', 0 ) ( 'H', 0 ) ( 'E', 0 );
	char last = ss[0];
	ss_seg_count[ last ] += 1;
	Size start = 1;
	Size end = 1;

	for(Size i = 1; i < ss.size(); i++) {
		if( ss[i] == last ) {
			// push if last element
			if(i + 1 == ss.size())
				ss_map.insert( std::make_pair< std::string, boost::tuple<Size, Size> >( last + boost::lexical_cast< std::string >( ss_seg_count[ss[i - 1] ] ), boost::make_tuple(start, i) ) );
			continue;
		} else {
			end = i ;
			ss_map.insert( std::make_pair< std::string, boost::tuple<Size, Size> >( last + boost::lexical_cast< std::string >( ss_seg_count[ss[i - 1] ] ), boost::make_tuple(start, end) ) );
			ss_seg_count[ ss[i] ] += 1;
			start = end + 1;
			last = ss[i +  1];
		}
	}
	return ss_map;
}


} // fldsgn
} // protocols
