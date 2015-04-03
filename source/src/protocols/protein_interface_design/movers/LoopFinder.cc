// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file protocols/protein_interface_design/LoopFinder.cc
/// @brief Finds loops and adds them to the basic::datacache::DataMap. Parseable options to control how loops are found
/// @author Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/LoopFinder.hh>
#include <protocols/protein_interface_design/movers/LoopFinderCreator.hh>

// Package headers
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <protocols/scoring/Interface.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.fwd.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loopfinder.hh>


// Numeric headers
#include <numeric/xyzVector.hh>


// the following must be included because LoopFinder is derived from DRMover
#include <core/scoring/ScoreFunction.hh>
//#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <basic/datacache/DataMap.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace protocols::moves;

static thread_local basic::Tracer TR( "protocols.protein_interface_design.movers.LoopFinder" );

std::string
LoopFinderCreator::keyname() const
{
	return LoopFinderCreator::mover_name();
}

protocols::moves::MoverOP
LoopFinderCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoopFinder );
}

std::string
LoopFinderCreator::mover_name()
{
	return "LoopFinder";
}

LoopFinder::LoopFinder() :
	simple_moves::DesignRepackMover( LoopFinderCreator::mover_name() )
{}

LoopFinder::LoopFinder(
	bool const interface,
	bool const ch1,
	bool const ch2,
	core::Size const min_length,
	core::Size const max_length,
	core::Size const mingap,
	core::Size const resnum,
	core::Real const ca_ca_distance,
	core::Real const iface_cutoff,
	protocols::loops::LoopsOP loops
) :
	simple_moves::DesignRepackMover( LoopFinderCreator::mover_name() ),
	interface_(interface),
	ch1_(ch1),
	ch2_(ch2),
	min_length_(min_length),
	max_length_(max_length),
	mingap_( mingap ),
	resnum_( resnum ),
	ca_ca_distance_( ca_ca_distance ),
	iface_cutoff_(iface_cutoff)
{
	loops_ = protocols::loops::LoopsOP( new protocols::loops::Loops( *loops ) );
}

LoopFinder::~LoopFinder() {}

protocols::moves::MoverOP LoopFinder::clone() const {
    return( protocols::moves::MoverOP( new LoopFinder( *this ) ) );
}

void
LoopFinder::apply( core::pose::Pose & pose )
{
	using namespace protocols::loops;

	protocols::scoring::Interface iface;
	iface.distance( iface_cutoff_ );
	if( interface_ ) {
		iface.jump( interface_ );
		iface.calculate( pose );
	}

	LoopsOP all_loops( new Loops );
	LoopsOP temp_loops( new Loops );
	LoopsOP ch1ch2_loops( new Loops );
	LoopsOP ch1_loops( new Loops );
	LoopsOP ch2_loops( new Loops );
	loopfinder( pose, *all_loops ); // runs dssp, inserts ss into the pose, and extracts all loops at least 3 residues long

	// all loops, no restrictions
	if( !interface_ && ch1_ && ch2_ && (max_length_ >= 1000) && (min_length_ <= 1) && (mingap_ == 0) ) {
		loops_ =  all_loops ;
		return;
	}

	if( all_loops->size() > 0 ) {
		for( Loops::const_iterator it = all_loops->begin(); it != all_loops->end(); ++it ) {
			LoopCOP loop( LoopOP( new Loop(*it) ) );
			if( pose.residue( loop->start() ).is_upper_terminus() || pose.residue( loop->stop() ).is_lower_terminus() ) continue; // skip if terminal loop
			if( loop->size() < min_length_ || loop->size() > max_length_ ) continue; // skip this loop

			// look at all residues in each loop.
			// do we want interface residues?
			// if so, go through all residues of each loop, checking whether interfacial and which chain it's in.
			// if not, only end up looking at 1st residue in loop to check which chain it's in (break statements prevent multiple addition)

			for( core::Size i = loop->start(); i <= loop->stop(); ++i ) {
				if( interface_ ) {
					// if even one residue is at the interface, count whole loop as interface for purposes of aggressive remodeling
					if( !iface.is_interface(i) ) continue;
				}

				if( resnum_ > 0) {
					TR.Debug <<"residue : " << resnum_ << " was specified" << std::endl;
					if( pose.residue( i ).xyz( "CA" ).distance( pose.residue( resnum_ ).xyz( "CA" ) ) < ca_ca_distance_ ) {
						TR.Debug<< "residue is within " << ca_ca_distance_ << " of the specified target chain(s)" << std::endl;
											}
					else {
						TR.Debug<<"residue is not within " << ca_ca_distance_ << " of this loop residue " << std::endl;
					  break;
					}
				}
				if( ch1_ && ch2_ ) {
					ch1ch2_loops->add_loop( *loop, mingap_ );
					break;
				}
				if( ch1_ ) {
					if( pose.residue(i).chain() == 1 ) ch1_loops->add_loop( *loop, mingap_ );
					break;
				}
				if( ch2_ ) {
					if( pose.residue(i).chain() == 2 ) ch2_loops->add_loop( *loop, mingap_ );
					break;
				}
				else TR << "Neither chain1 nor chain2 specified for loop finding. No loops added!" << std::endl;
			}
		}
		TR.Debug << "ch1ch2 loops " << *ch1ch2_loops << std::endl;
		TR.Debug << "ch1 loops " << *ch1_loops << std::endl;
		TR.Debug << "ch2 loops " << *ch2_loops << std::endl;
		if( ch1_ && ch2_ ) *loops_ = *ch1ch2_loops ; // copy, rather than change pointer address
		else if( ch1_ ) *loops_ = *ch1_loops ;
		else if( ch2_ ) *loops_ = *ch2_loops ;
		//runtime_assert( data_->has( "loops", "loops" ) );
	}
	else
		TR << "No loops found." << std::endl;
}

std::string
LoopFinder::get_name() const {
	return LoopFinderCreator::mover_name();
}


void
LoopFinder::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
	interface_ = tag->getOption<Size>( "interface", 1 );
	ch1_ = tag->getOption<bool>( "ch1", 0 );
	ch2_ = tag->getOption<bool>( "ch2", 1 );
	min_length_ = tag->getOption<core::Size>( "min_length", 3 );
	max_length_ = tag->getOption<core::Size>( "max_length", 1000 );
	iface_cutoff_ = tag->getOption<core::Real>( "iface_cutoff", 8.0 );
	runtime_assert( ch1_ || ch2_ );
	if( interface_ ) runtime_assert( pose.num_jump() >= 1 );
	mingap_ = tag->getOption<core::Size>( "mingap", 1 );
	if( tag->hasOption( "resnum" ) || ( tag->hasOption( "pdb_num" ) ))	{
		resnum_ = core::pose::get_resnum( tag, pose );
		TR<<"user specified residue " << resnum_ << " for distance cutoff" << std::endl;

	}
	else	{
		resnum_ = 0;
		TR<<"no target residue has been specified"<< std::endl;
	}
  ca_ca_distance_ = tag->getOption<core::Real>( "CA_CA_distance", 15 );
	if( resnum_ != 0 )	TR<<"distance cutoff from user defined residue is " << ca_ca_distance_ << std::endl;

	// add loopsOP to the basic::datacache::DataMap
	loops_ = protocols::loops::LoopsOP( new protocols::loops::Loops );
	data.add( "loops", "found_loops", loops_ );

	TR << "LoopFinder mover: interface="<<interface_<<" iface_cutoff="<<iface_cutoff_<<" ch1="<<ch1_<<" ch2="<<ch2_<<" min_length="<<min_length_<<
		" max_length="<<max_length_<< std::endl;

}

} //movers
} //protein_interface_design
} //protocols
