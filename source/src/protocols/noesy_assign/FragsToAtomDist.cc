// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
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
///

// Unit Headers
#include <protocols/noesy_assign/FragsToAtomDist.hh>

// Package Headers
#include <protocols/noesy_assign/MethylNames.hh>
#include <protocols/noesy_assign/Exceptions.hh>


// Project Headers
#include <protocols/simple_moves/sidechain_moves/JumpRotamerSidechainMover.hh>

#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/Frame.hh>

#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

#include <core/pose/Pose.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/pose/util.hh>
#include <core/types.hh>

// Utility headers
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>
#include <numeric/xyzVector.io.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <string>
//#include <iostream>

static thread_local basic::Tracer tr( "protocols.noesy_assign.FragsToAtomDist" );

core::Real half_adjust( core::Real in ) {
	return floor(in+0.5);
}

std::map< core::Real,core::Size > simple_histogram ( utility::vector1< core::Real > samples, core::Size bin_number=20 ) {
	core::Real low( *( samples.begin() ) );
	core::Real high( *( samples.begin() ) );
	for ( utility::vector1< core::Real >::const_iterator it=samples.begin(); it!=samples.end(); it++ ) {
		if ( *it < low ) { low= *it; }
		if ( *it > high ) { high= *it; }
	}
	core::Real step( ( high-low )/bin_number );
	std::map< core::Real,core::Size > hist;
	if ( step !=0 ) {
		for (core::Size i=0; i< bin_number; i++) {
			hist[ low+step*( i+0.5 ) ]=0;
		}
	}
	else {
		for (core::Size i=0; i< bin_number; i++) {
			hist[ low+1*( i ) ]=0;
		}
	}

	for ( utility::vector1< core::Real >::const_iterator i=samples.begin(); i !=samples.end(); i++ ) {
		core::Real bin( ( floor( ( *i-low )/step )+0.5)*step+low );
		if ( *i == high ) {	bin=low+(bin_number-0.5)*step; }
		//		if ( bin >= high ) { bin=low+(bin_number-0.5)*step; }
		//		if ( bin <= low ) { bin=low+0.5*step; }
		hist[ bin ]+=1;
	}
	return hist;
}




using namespace std;
using namespace core;
using namespace core::pose;
using namespace core::id;
using namespace fragment;
using namespace core::chemical;

namespace protocols {
namespace noesy_assign {

typedef utility::vector1< core::Size > SizeList;
typedef utility::vector1< std::pair< std::string, SizeList> > AtomGrps;
typedef utility::vector1< AtomGrps > GroupList;
typedef std::map< core::id::AtomID,  FragsToAtomDist::DistanceRecord > InnerMap;
typedef std::map< core::id::AtomID, InnerMap > DistanceMap;



// core::Real FragsToAtomDist::DistanceRecord::popular_bin() const {
// 	core::Size max_bin(0);
// 	core::Real pop_dist(-1.0);
// 	std::map< core::Real, core::Size >::const_iterator it;
// 	tr.Info << "size of hist_dist_  "<< hist_dist_.size() <<std::endl;
// 	for ( it=hist_dist_.begin(); it !=hist_dist_.end(); it++ ) {
// 		tr.Info << " it->second "<< it->second<<std::endl;
// 		tr.Info << " max_bin "<< max_bin<<std::endl;
// 		if ( it->second > max_bin ) {
// 			tr.Info << " it->first "<< it->first<<std::endl;
// 			max_bin=it->second;
// 			pop_dist=it->first;
// 		}
// 	}
// 	tr.Info << " pop_dist "<< pop_dist<<std::endl;
// 	return pop_dist;
// }
FragsToAtomDist::DistanceRecord FragsToAtomDist::NO_CONTACT( 100.0, 100.0, 100.0, 1 );

std::ostream& operator<< ( std::ostream& os, FragsToAtomDist::DistanceRecord const& dr ) {
	return os << " " << setw(10) << setprecision(4) << dr.average_dist6()
						<< " " << setw(10) << setprecision(4) << dr.average_dist()
						<< " " << setw(10) << setprecision(4) << dr.min_dist();
		//						<< " " << setw(10) << setprecision(4) << dr.popular_bin();
}

std::istream& operator>> ( std::istream& is, FragsToAtomDist::DistanceRecord& dr ) {
	core::Real t1, t2, t3;
	is >> t1 >> t2 >> t3;
	dr = FragsToAtomDist::DistanceRecord( pow( t1, -6.0 ), t2, t3, 1 );
	return is;
}

void FragsToAtomDist::DistanceRecord::record_inv6( core::Real inv6 ) {
	core::Real const dist( std::pow( inv6, -1.0/6.0 ) );
	cum_dist6_ += inv6;
	cum_dist_ += dist;
 	min_dist_ = std::min( min_dist_, dist );
	++count_;
	dist_track_.push_back( dist );
	//tr.Info << " dist " << dist << " half adjust " << half_adjust( dist ) <<std::endl;
	//	hist_dist_[ half_adjust( dist ) ]=hist_dist_[ half_adjust( dist ) ]+1;
	//tr.Info << " hist_dist_ " << hist_dist_[ half_adjust( dist ) ] <<std::endl;
	//	hist_dist6_[ half_adjust( inv6 ) ]=hist_dist6_[ half_adjust( inv6 ) ]+1;
}

void FragsToAtomDist::DistanceRecord::record( core::Real dist ) {
	core::Real const inv6( std::pow( dist, -6.0 ) );
	cum_dist6_ += inv6;
	cum_dist_ += dist;
 	min_dist_ = std::min( min_dist_, dist );
	++count_;
	dist_track_.push_back( dist );
}


void initialize_group_list(
  	 GroupList& grp_list,
 		 core::pose::Pose const& pose )
{
	MethylNameLibrary const& methyl_lib( *MethylNameLibrary::get_instance() );
	for ( core::Size pos =1; pos<=pose.total_residue(); pos++) {
		core::chemical::ResidueType const& rsd(pose.residue_type( pos ));
		//		std::string name( utility::trim( pose.residue_type( pos ).atom_name( iatom1 ) ) );
		MethylNames const& methyls( methyl_lib[ rsd.aa() ] );
		AtomGrps new_grps;
		for ( MethylNames::const_iterator it = methyls.begin(); it != methyls.end(); ++it ) {
			SizeList indices;
			tr.Info << "pos " << pos << " " << it->first << " ";
			for ( MethylNames::AtomList::const_iterator ait = it->second.begin(); ait != it->second.end(); ++ait ) {
				std::string atom( *ait );
				if ( pos == 1 && atom == "H" ) {
					atom="1H";
				}
				indices.push_back( rsd.atom_index( atom ) );
				tr.Info << atom << " ";
			}
			new_grps.push_back( std::make_pair( it->first, indices ) );
			//			tr.Info << "CONTROL: " << it->first << " " << new_grps.back().first << std::endl;
			tr.Info << std::endl;
		}
		tr.Info << "generated " << new_grps.size() << " proton-groups for position " << pos << " " << rsd.name3()<< std::endl;
		grp_list.push_back( new_grps );
	}
}

bool FragsToAtomDist::check_sequence(std::string const& other_sequence) const {
	if (sequence_.compare(other_sequence)==0) {
		return true;
	}	else {
		return false;
	}
}

FragsToAtomDist::DistanceRecord const& FragsToAtomDist::distance_record( core::id::NamedAtomID atom1, core::id::NamedAtomID atom2 ) const {
	swap_atoms( atom1, atom2 );
	NamedDistanceMap::const_iterator it1 = named_distmap_.find( atom1 );
	if ( it1 != named_distmap_.end() ) {
		NamedInnerMap::const_iterator it2 = it1->second.find( atom2 );
		if ( it2 != it1->second.end() ) {
			return it2->second;
		}
	}
	//	tr << "my test on store with distance_record 1  "<<atom1 <<std::endl;
	//	tr << "my test on store with distance_record 2  "<<atom2 <<std::endl;
	return NO_CONTACT;
}

void FragsToAtomDist::write_to_stream(std::ostream& output) const {
	Size count =1;
	for (std::string::const_iterator i = sequence_.begin();i != sequence_.end();++i) {
		if ( count%50 == 1 ) {
			output << "DATA SEQUENCE ";
		}
		output << *i;
		if ( count%50 == 0 ) {
			output << endl;
		}
		if ( count%10 == 0 && count%50 != 0 ) {
			output << " ";
		}
		count++;
	}
	output << endl << endl;
	//	if ( FragsToAtomDist::r6_averaged_ ) {
	//		output << "DATA R6_AVERAGED YES" << endl;
	//	}
	//	else {
	//		output << "DATA R6_AVERAGED NO" << endl;
	//	}
	output << "VARS   RESID1 ATOMNAME1 RESID2 ATOMNAME2 DIST_R6_AVERAGED DIST_AVERAGED DIST_MIN" << endl << endl;
	NamedDistanceMap::const_iterator iter1;
	NamedInnerMap::const_iterator iter2;
	for ( iter1 = named_distmap_.begin(); iter1 != named_distmap_.end(); iter1++ )	{
		for (iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++ ) {
			if ( iter2->second.is_valid() ) {
				output << setw(5)<<iter1->first.rsd()
							 << setw(5)<<iter1->first.atom()
							 << setw(5)<<iter2->first.rsd()
							 << setw(8)<<iter2->first.atom()
							 << setw(8)<< iter2->second <<endl;
			}
		}
	}
}

void FragsToAtomDist::write_hist_to_stream(std::ostream& output) const {
	Size count =1;
	for (std::string::const_iterator i = sequence_.begin();i != sequence_.end();++i) {
		if ( count%50 == 1 ) {
			output << "DATA SEQUENCE ";
		}
		output << *i;
		if ( count%50 == 0 ) {
			output << endl;
		}
		if ( count%10 == 0 && count%50 != 0 ) {
			output << " ";
		}
		count++;
	}
	output << endl << endl;
	NamedDistanceMap::const_iterator iter1;
	NamedInnerMap::const_iterator iter2;
	std::map< core::Real, core::Size >::const_iterator iter3;
	output << "VARS   RESID1 ATOMNAME1 RESID2 ATOMNAME2 LOW HIGH" << endl << endl;
	for ( iter1 = named_distmap_.begin(); iter1 != named_distmap_.end(); iter1++ )	{
		for (iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++ ) {
			output << setw(5)<<iter1->first.rsd()
						 << setw(8)<<iter1->first.atom()
						 << setw(5)<<iter2->first.rsd()
						 << setw(8)<<iter2->first.atom();
			utility::vector1< core::Real > dist_track_to_write ( iter2->second.dist_track() );
			std::map< core::Real,core::Size > hist_dist_to_write( simple_histogram( dist_track_to_write ) );
			core::Real low( *( dist_track_to_write.begin() ) );
			core::Real high( *( dist_track_to_write.begin() ) );
			for ( utility::vector1< core::Real >::const_iterator i=dist_track_to_write.begin();i !=dist_track_to_write.end();i++ ) {
				if ( *i < low ) { low= *i; }
				if ( *i > high ) { high= *i; }
			}
			output << setw(10) << setprecision(4) << low
						 << setw(10) << setprecision(4) << high;
			for (iter3 = hist_dist_to_write.begin(); iter3 != hist_dist_to_write.end(); iter3++) {
				output << setw(8) << setprecision(14)<< iter3->second;
			}
			output <<endl;
		}
	}
}



void FragsToAtomDist::write_to_file(std::string const& filename) const {
	utility::io::ozstream output( filename  );
	write_to_stream(output);
}

void FragsToAtomDist::write_hist_to_file(std::string const& filename) const {
	utility::io::ozstream output( filename  );
	write_hist_to_stream(output);
}

void FragsToAtomDist::read_from_stream(std::istream& input) {
	std::string line;
	while ( getline( input, line ) ) {
		if ( line.size() == 0 ) continue;
		std::istringstream tag_stream( line );
		if (line.substr(0,4)=="DATA") {
			if (line.substr(0,13)=="DATA SEQUENCE") {
				line.erase(0,14);
				for (std::string::iterator i = line.begin();i != line.end();++i) {
					if ( *i !=' ' ) {
						sequence_+=*i;
					}
				}
			}
			//			if (line.substr(0,20)=="DATA R6_AVERAGED YES") {
			//r6_averaged_=true;
			//}
			//else if (line.substr(0,19)=="DATA R6_AVERAGED NO") {
			//	r6_averaged_=false;
			//}
		}
    else if ( line.substr(0,4) !="VARS" && line.length()>10 ) {
			int residue1,residue2;
			std::string name1,name2;
			DistanceRecord dist;
			tag_stream >> residue1 >> name1 >> residue2 >> name2 >> dist;
			NamedAtomID atom1( name1,residue1);
			NamedAtomID atom2( name2,residue2);
			named_distmap_[atom1][atom2]=dist;
		}
	}
}

void FragsToAtomDist::read_from_file(std::string const& filename) {
	utility::io::izstream input( filename );
	read_from_stream( input );
}



void FragsToAtomDist::generate_from_fragments( FragSetOP fragments, std::string const& sequence, core::Size cycles, core::Size cycle_mod ) {
	frags_=fragments;
	set_sequence(sequence);
	//set_r6_averaged(r6_averaging);
	compute_average_distances( cycles, cycle_mod );
}

void FragsToAtomDist::generate_from_frag_file( std::string const& filename, std::string const& sequence, core::Size cycles, core::Size cycle_mod ) {
	generate_from_fragments( FragmentIO().read_data( filename ), sequence, cycles, cycle_mod);
}

//re-order atom-ids so that always smaller residue number and/or smaller atom number comes first
void FragsToAtomDist::swap_atoms( NamedAtomID& atom1, NamedAtomID& atom2 ) const {
	bool swap( false );
	if ( atom1.rsd() > atom2.rsd() ) {
		swap = true;
	} else if ( atom1.rsd() == atom2.rsd() ) {
		AA aa( aa_from_oneletter_code( sequence_[ atom1.rsd()-1 ] ) );
		MethylNameLibrary const& methyl_lib( *MethylNameLibrary::get_instance() );
		try {
			if ( methyl_lib[aa].proton_index( atom1.atom() ) > methyl_lib[ aa ].proton_index( atom2.atom() ) ) {
				swap = true;
			}
		} catch( EXCN_UnknownAtomname& excn ) {
			excn.add_msg( "cannot find proton/methyl name at position "+ObjexxFCL::string_of( atom1.rsd() ) );
		}
	}
	if ( swap ) {
		NamedAtomID atom3( atom1 );
		atom1=atom2;
		atom2=atom3;
	}
}

//translate from AtomID into NamedAtomID --- store in named_distmap_
void store_distmap_with_namedatoms(
	 FragsToAtomDist::NamedDistanceMap &named_distmap,
   DistanceMap const& distmap,
	 GroupList const& proton_groups
) {
	DistanceMap::const_iterator iter1;
	InnerMap::const_iterator iter2;
	for ( iter1 = distmap.begin(); iter1 != distmap.end(); iter1++)	{
		for (iter2=iter1->second.begin(); iter2 != iter1->second.end();iter2++) {
			core::id::AtomID const& grp_1(iter1->first);
			core::id::AtomID const& grp_2(iter2->first);
			std::string const& grp_name1( proton_groups[ grp_1.rsd() ][ grp_1.atomno() ].first );
			std::string const& grp_name2( proton_groups[ grp_2.rsd() ][ grp_2.atomno() ].first );
			core::id::NamedAtomID named_atom_1( grp_name1, grp_1.rsd() );
			core::id::NamedAtomID named_atom_2( grp_name2, grp_2.rsd() );
			named_distmap[named_atom_1][named_atom_2]=iter2->second;
		}
	}
}

void store_distance_snapshot(
	 DistanceMap& distmap,
	 //	 DistanceMap& countmap,
	 core::pose::Pose const& short_pose,
	 Size start,
	 Size short_size,
	 GroupList const& grps
	 //	 bool r6_averaging
) {
	for ( core::Size rsd1 =1; rsd1<=short_size; rsd1++) {
		AtomGrps const& group1( grps[ rsd1+start-1 ] );
		for ( core::Size igrp1 =1; igrp1<=group1.size(); igrp1++ ) {
			for ( core::Size rsd2 =rsd1; rsd2<=short_size; rsd2++) {
				AtomGrps const& group2( grps[ rsd2+start-1 ] );
				for ( core::Size igrp2 = ( rsd2==rsd1 ? igrp1+1 : 1 ); igrp2<=group2.size(); igrp2++ ) {
					Real cumdist( 0 );
					for ( SizeList::const_iterator atom1 = group1[igrp1].second.begin(); atom1 != group1[igrp1].second.end(); ++atom1 ) {
						PointPosition const & xyz_1 = short_pose.xyz( id::AtomID( *atom1, rsd1  ) );
						for ( SizeList::const_iterator atom2 = group2[igrp2].second.begin(); atom2 != group2[igrp2].second.end(); ++atom2 ) {
							PointPosition const & xyz_2 = short_pose.xyz( id::AtomID( *atom2, rsd2 ) );
							Real const inv_dist2( 1.0/xyz_1.distance_squared( xyz_2 ) );
							Real const inv_dist6( inv_dist2 * inv_dist2 * inv_dist2 );
							cumdist += inv_dist6;
						}
					}
					id::AtomID grpid1( igrp1, rsd1+start-1 );
					id::AtomID grpid2( igrp2, rsd2+start-1 );
					distmap[ grpid1 ][ grpid2 ].record_inv6( cumdist );
					//tr << "my test on store with snapshot 1  "<<igrp1 << "   "<< rsd1+start-1 <<std::endl;
					//tr << "my test on store with snapshot 2  "<<igrp2 << "   "<< rsd2+start-1 <<std::endl;
				}
			}
		}
	}
}

void FragsToAtomDist::compute_average_distances(core::Size cycles,core::Size dump_freq ) {
	using namespace protocols::simple_moves;

	if (dump_freq>cycles) dump_freq=cycles;

	sidechain_moves::JumpRotamerSidechainMover jrmover;
	core::scoring::ScoreFunctionOP scorefxn(  scoring::get_score_function() );
	scorefxn->set_weight( core::scoring::fa_dun, 0 ); //since we use the JumpRotamer Mover the dunbrack energy is already in the sampling bias
	core::pose::Pose pose;
	core::pose::make_pose_from_sequence(	pose, sequence_, "fa_standard", true	);

	//SizeList natoms( pose.total_residue(), 0);
	//	initialize_natoms( natoms, pose );

	GroupList proton_groups;
	initialize_group_list( proton_groups, pose );

	//	DistanceMap countmap;
	DistanceMap distmap;
	//initialize_maps( distmap, countmap, proton_groups );

	for ( FrameIterator frame = frags_->nonconst_begin(), eframe=frags_->nonconst_end();
				frame != eframe; ++frame ) {
		core::pose::Pose short_pose;
		core::pose::make_pose_from_sequence(
												short_pose,
												sequence_.substr( frame->start()-1, frame->length() ),
												"fa_standard",
												true
		);
		//		tr << " my test   " << short_pose.sequence() << "   "<<std::endl;
		Size const short_size( frame->length() );
		std::string short_sequence = short_pose.sequence();
		Size const short_start( frame->start() );
		tr.Info << "frame position: " << frame->start() << std::endl;

		frame->shift_to(1);

		utility::vector1< core::Real > old_chi(4);
		utility::vector1< core::Real > new_chi(4);
		for ( Size ifrag=1; ifrag<=frame->nr_frags(); ifrag++ ) {
			//	if ( ifrag % 10 == 0 ) tr.Info << "fragment (" << ifrag << "," << frame->nr_frags() << ")" << std::endl;
			frame->apply( ifrag, short_pose );
			core::Real score_old(scorefxn->score(short_pose));

			//start monte-carlo simulation
			for ( core::Size cycle = 1; cycle<=cycles; cycle++ )	{
				//better select rsd randomly
				//for ( Size irsd = 1; irsd<=short_size; irsd++) {
				Size const irsd( numeric::random::rg().random_range(1, short_size) );
				core::conformation::Residue const& rsd( short_pose.residue( irsd ));
				core::Size const n_chi_angles(rsd.nchi());
				old_chi=rsd.chi();

				//make chi-move
				if ( n_chi_angles != 0 ) {
					jrmover.make_chi_move(rsd,old_chi,new_chi);
				}

				//set new chi-angles in pose
				for (core::Size ichi = 1; ichi<=n_chi_angles; ichi++) {
					short_pose.set_chi(ichi,irsd,new_chi[ ichi ]);
				}

				//evaluate energy and accept/reject
				core::Real score_new(scorefxn->score(short_pose));
				if ( score_new > score_old ) {
					core::Real Delta_score(score_new - score_old);
					core::Real p( exp( -1*Delta_score ) );
					if ( numeric::random::rg().uniform() > p) {
						for (core::Size ichi = 1; ichi<=n_chi_angles; ichi++) {
							short_pose.set_chi(ichi,irsd,old_chi[ ichi ]);
						}
					} else {
						score_old = score_new;
						old_chi = new_chi;
					}
				} else {
					score_old = score_new;
					old_chi = new_chi;
				}

				if (cycle%dump_freq == 0 ) {
					tr.Debug << "Store distances " << cycle << " " << dump_freq << std::endl;
					store_distance_snapshot( distmap, short_pose, short_start, short_size, proton_groups );
				}
			}
		}
	}
	//average_distmap( distmap,countmap, r6_averaged_ );
	store_distmap_with_namedatoms(named_distmap_,distmap,proton_groups);
}


} //noesy_assign
} //protocols
