// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file protocols/protein_interface_design/movers/HotspotHasherMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/HotspotHasherMover.hh>
#include <protocols/protein_interface_design/movers/HotspotHasherMoverCreator.hh>

// Package headers

// Project headers
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static thread_local basic::Tracer TR( "protocols.protein_interface_design.movers.HotspotHasherMover" );

std::string
HotspotHasherMoverCreator::keyname() const
{
	return HotspotHasherMoverCreator::mover_name();
}

protocols::moves::MoverOP
HotspotHasherMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new HotspotHasherMover );
}

std::string
HotspotHasherMoverCreator::mover_name()
{
	return "HotspotHasher";
}

HotspotHasherMover::HotspotHasherMover() : protocols::moves::Mover( HotspotHasherMoverCreator::mover_name() ) { }
HotspotHasherMover::HotspotHasherMover(
	std::vector<std::string> const resnames,
	core::scoring::ScoreFunctionCOP scorefxn,
	core::Size const n_stubs,
	core::Size const target_resnum,
	protocols::filters::FilterOP hotspot_filter,
	core::Real const target_distance,
	std::string const hashin_fname,
	std::string const hashout_fname
) :
	protocols::moves::Mover( HotspotHasherMoverCreator::mover_name() ),
	scorefxn_(scorefxn),
	resnames_(resnames),
	n_stubs_(n_stubs),
	target_resnum_(target_resnum),
	target_distance_(target_distance),
	hashin_fname_(hashin_fname),
	hashout_fname_(hashout_fname)
{
	hotspot_filter_ = hotspot_filter;
}

HotspotHasherMover::~HotspotHasherMover() {}

protocols::moves::MoverOP
HotspotHasherMover::clone() const {
	return( protocols::moves::MoverOP( new HotspotHasherMover( *this ) ) );
}

void HotspotHasherMover::apply( core::pose::Pose & pose ) {

	// finding a request for ALL adds the main 18 aa's to
	if ( std::find( resnames_.begin(), resnames_.end(), "ALL" ) != resnames_.end() ) {
		resnames_.push_back( "ALA" );
		resnames_.push_back( "ARG" );
		resnames_.push_back( "ASN" );
		resnames_.push_back( "ASP" );
		resnames_.push_back( "GLU" );
		resnames_.push_back( "GLN" );
		resnames_.push_back( "HIS" );
		resnames_.push_back( "ILE" );
		resnames_.push_back( "LEU" );
		resnames_.push_back( "LYS" );
		resnames_.push_back( "MET" );
		resnames_.push_back( "PHE" );
		resnames_.push_back( "PRO" );
		resnames_.push_back( "SER" );
		resnames_.push_back( "THR" );
		resnames_.push_back( "TRP" );
		resnames_.push_back( "TYR" );
		resnames_.push_back( "VAL" );
	}

	protocols::hotspot_hashing::HotspotStubSet stubset;

	// read existing hashes
	if ( utility::file::file_exists( hashin_fname_ ) ) {
		stubset.read_data( hashin_fname_ );
		TR << "Found hash file " << hashin_fname_ << std::endl;
	}
	if ( utility::file::file_exists( hashout_fname_ ) ) { // useful for interrupted runs
		stubset.read_data( hashout_fname_ );
		TR << "Found hash file " << hashout_fname_ << std::endl;
	}


	// for each residue requested
	for ( std::vector< std::string >::const_iterator it=resnames_.begin() ; it!=resnames_.end(); ++it ) {
		std::string resname = *it;

		TR << "Hash contains " << stubset.size(resname) << " " << resname << " stubs." << std::endl;

		// check to see if we've already finished our hash
		core::Size stubs_left = n_stubs_;
		stubs_left -= stubset.size( resname );
		if ( stubs_left <= 0 ) {
			// perform a scorecut
			if ( basic::options::option[ basic::options::OptionKeys::out::scorecut ].user() ) {
				Real score_cutoff = basic::options::option[ basic::options::OptionKeys::out::scorecut ]();
				std::stringstream i;
				i.str("");
				i << score_cutoff;
				stubset.clear();
				stubset.read_data( hashout_fname_ );
				protocols::hotspot_hashing::HotspotStubSetOP cut_stubs = stubset.subset( score_cutoff );
				std::string newfname = i.str() + "cut_" + hashout_fname_;
				cut_stubs->write_all( newfname );
			}
			return;
		}

		// do hashing in 10-stub cycles to minimize file i/o
		Size n_per(10);

		Size n_cycles = n_stubs_ / n_per;
		// make sure we do at least one cycle
		if ( n_cycles <= 0 ) n_cycles = 1;
		// PERFORM HASHING
		for ( Size i = 1; i <= n_cycles; ++i ) {
			stubset.clear();
			stubset.score_threshold( score_threshold_ );
			TR << "Finding " << n_per*i << "/" << n_stubs_ << " " << resname << " stubs" ;
			if ( target_resnum_ ) {
				TR << " " << target_distance_ << "A from " << target_resnum_ << std::endl;
				stubset.fill( pose, scorefxn_, target_resnum_, target_distance_, resname, n_per );
			} else {
				TR << "." << std::endl;
				stubset.fill( pose, scorefxn_, resname, n_per );
			}
			stubset.write_all( hashout_fname_ );
		}
	} // for each residue

	// perform a scorecut
	if ( basic::options::option[ basic::options::OptionKeys::out::scorecut ].user() ) {
		Real score_cutoff = basic::options::option[ basic::options::OptionKeys::out::scorecut ]();
		std::stringstream i;
		i.str("");
		i << score_cutoff;
		stubset.clear();
		stubset.read_data( hashout_fname_ );
		protocols::hotspot_hashing::HotspotStubSetOP cut_stubs = stubset.subset( score_cutoff );
		std::string newfname = i.str() + "cut_" + hashout_fname_;
		cut_stubs->write_all( newfname );
	}
} // HotspotHasherMover::apply

std::string
HotspotHasherMover::get_name() const {
	return HotspotHasherMoverCreator::mover_name();
}


void
HotspotHasherMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const & filters, Movers_map const &, core::pose::Pose const & pose )
{

	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();

	n_stubs_ = tag->getOption<core::Size>( "nstubs", 1000 );

	// target residue
	// target_resnum gets set below with residues
	target_resnum_ = 0;
	if ( tag->hasOption( "target_residue_pdb_num" ) || tag->hasOption( "target_residue_res_num" ) ) {
		target_resnum_ = core::pose::get_resnum( tag, pose, "target_residue_" );
	}

	target_distance_ = tag->getOption<core::Real>( "target_distance", 15.0 );

	score_threshold_ = tag->getOption<core::Real>( "threshold", -1.0 );

	// hash in/out
	hashin_fname_ = tag->getOption<std::string>( "in", "");
	hashout_fname_ = tag->getOption<std::string>( "out", "hash.stubs");

	// filter
	std::string const hotspot_filter_name( tag->getOption<std::string>( "hotspot_filter", "true_filter" ) );
	protocols::filters::Filters_map::const_iterator find_filter( filters.find( hotspot_filter_name ));
	bool const filter_found( find_filter != filters.end() );
	if ( filter_found ) {
		hotspot_filter_ = find_filter->second->clone();
	} else {
		if ( hotspot_filter_name != "true_filter" ) {
			TR<<"***WARNING WARNING! Filter defined for HotspotHasher not found in filter_list!!!! Defaulting to truefilter***"<<std::endl;
			runtime_assert( filter_found );
		} else {
			hotspot_filter_ = protocols::filters::FilterOP( new protocols::filters::TrueFilter );
		}
	}

	// residues
	utility::vector0< TagCOP > const & hasher_tags( tag->getTags() );
	for ( utility::vector0< TagCOP >::const_iterator hash_it=hasher_tags.begin(); hash_it!=hasher_tags.end(); ++hash_it ) {
		TagCOP const hash_tag_ptr = *hash_it;
		std::string tag_name = hash_tag_ptr->getName();
		if ( tag_name == "residue" ) {
			std::string resname( hash_tag_ptr->getOption< std::string >( "type", "" ) );
			resnames_.push_back( resname );
		}
	}
	runtime_assert( resnames_.size() > 0 );

	TR<<"hashing mover finding residues: ";
	for ( std::vector< std::string >::const_iterator it=resnames_.begin() ; it!=resnames_.end(); ++it ) TR<<*it<<" ";
	if ( target_resnum_ ) TR << target_distance_ << "A away from residue " << target_resnum_ << std::endl;
	TR<<std::endl;
} // HotspotHasherMover::parse_my_tag

} //movers
} //protein_interface_design
} //protocols
