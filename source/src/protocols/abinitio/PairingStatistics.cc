// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Oliver Lange

// Unit Headers
#include <protocols/abinitio/PairingStatistics.hh>

// Package Headers
#include <protocols/abinitio/Template.hh>
#include <protocols/abinitio/Templates.hh>
// AUTO-REMOVED #include <protocols/abinitio/TemplateJumpSetup.hh>

// Project Headers
#include <core/types.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/id/AtomID.hh>

#include <core/fragment/FragSet.fwd.hh>
// AUTO-REMOVED #include <core/fragment/ConstantLengthFragSet.hh>
// AUTO-REMOVED #include <core/fragment/FragData.hh> //to get secondary structure
// AUTO-REMOVED #include <core/fragment/SecstructSRFD.hh> //to get secondary structure
// AUTO-REMOVED #include <core/fragment/FragID_Iterator.hh>

// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>

// AUTO-REMOVED #include <core/sequence/util.hh>
#include <core/id/SequenceMapping.hh>

// AUTO-REMOVED #include <core/fragment/SecondaryStructure.hh>
// AUTO-REMOVED #include <protocols/jumping/JumpSample.hh>


//numeric headers
#include <numeric/random/random.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// Utility headers
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>


// AUTO-REMOVED
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>


// C++ headers
#include <cstdlib>

#ifdef WIN32
#include <iterator>
#endif

#include <string>
#include <vector>


static basic::Tracer tr("protocols.abinitio.PairingStats");
using namespace core;
using namespace basic;
using namespace basic::options;
using namespace ObjexxFCL::format;

void protocols::abinitio::PairingStatistics::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	option.add_relevant(  templates::force_native_topology );
	option.add_relevant(  templates::topology_rank_cutoff  );
}

namespace protocols {
namespace abinitio {

static numeric::random::RandomGenerator RG(5651234);  // <- Magic number, do not change it!

using namespace core;
//using namespace jumping;
using namespace basic::options;
using namespace basic::options::OptionKeys;



PairingStatEntry::PairingStatEntry() {}

PairingStatEntry::PairingStatEntry( core::scoring::dssp::StrandPairing const& strand, Model const& id ) :
  strand_pairing_( strand ),
	weight_( -1.0 )
{
  models_.push_back( id );
}


bool PairingStatEntry::add_pairing( core::scoring::dssp::StrandPairing const& new_strand, Model const& id ) {
  bool success( true );
  if ( models_.size() ) {
		//queries 'mergeable'
    success = strand_pairing_.merge( new_strand, true /*do_merge */);
  } else {
    strand_pairing_ = new_strand;
  }
  if ( success ) models_.push_back( id );
  return success;
}


bool PairingStatEntry::compatible( core::scoring::dssp::StrandPairing const& strand ) const {
  return strand_pairing_.mergeable( strand );
}

bool PairingStatEntry::has_model( std::string const& model ) const {
	for ( ModelList::const_iterator it = models_.begin(), eit = models_.end(); it != eit; ++it ) {
		if ( model == *it ) return true;
	}
	return false;
}

bool PairingStatEntry::operator==( PairingStatEntry const& other) const {
	return strand_pairing_.mergeable( other.pairing() );
}

bool PairingStatEntry::operator!=( PairingStatEntry const& other) const {
	return !strand_pairing_.mergeable( other.pairing() );
}

bool _MergableEntries::operator() (core::scoring::dssp::StrandPairing const& p1, core::scoring::dssp::StrandPairing const& p2) const {
	return p1.mergeable( p2 );
}

void PairingStatistics::add_entry(core::scoring::dssp::StrandPairing const& ps, Model const& id ) {
	StatEntries::iterator itentry = entries_.find( ps );
	bool merged( false );
	if ( itentry != entries_.end() ) {
		//			tr.Trace << "found that it matches with " << itentry->first << std::endl;
		merged = itentry->second.add_pairing( ps, id );
		//			if (!merged) tr.Trace << "strangely couldn't merge the matched one..."<< std::endl;
		//			else tr.Trace << "new merged strand is " << itentry->second.pairing() << std::endl;
		core::scoring::dssp::StrandPairing const& key( itentry->first );
		core::scoring::dssp::StrandPairing const& new_pairing( itentry->second.pairing() );
		if ( merged && ( !entries_.key_eq()( new_pairing, key ) ) ) {
			core::scoring::dssp::StrandPairing new_pairing( itentry->second.pairing() ); //need a copy
			entries_.erase( itentry );
			add_entry( new_pairing, id );
			//			entries_[ itentry->second.pairing() ] = itentry->second;
		}
	}
	if ( !merged ) entries_[ ps ]=PairingStatEntry( ps, id );
		//for ( StatEntries::iterator itentry = entries_.begin(), eitentry = entries_.end();
		//	  itentry != eitentry && !merged; ++itentry ) {
		//      merged = itentry->add_pairing( *it, id );
		//    }
}

void PairingStatistics::add_topology( core::scoring::dssp::StrandPairingSet const& topology, Model const& id ) {
	//	if (! topology.size() ) return; //also add empty sets  -- otherwise the modelname floats around and can't be found in this list
	for ( core::scoring::dssp::StrandPairingSet::const_iterator it = topology.begin(), eit = topology.end();
			it != eit; ++it  ) {
		//bool merged ( false );
		tr.Trace << "adding stand pairing to hash.. " << *it << std::endl;
		add_entry( *it, id );
	}
	topols_[ id ] = topology;
}


void PairingStatistics::compute( Templates const& templates ) {
  ModelFreq model_freq;  //count of the underlying structures ( letters 1-4 of model-name )
  for ( Templates::const_iterator it = templates.begin(), eit = templates.end();
	it != eit; ++it ) {
    Template const& model( *it->second );
    core::scoring::dssp::PairingList template_pairings, target_pairings;
    model.strand_pairings().get_beta_pairs( template_pairings );
    model.map_pairings2target( template_pairings, target_pairings );
    add_topology( core::scoring::dssp::StrandPairingSet( target_pairings ), model.name() );
    model_freq[ model.name().substr(0,4) ] += 1;
  }
  compute_model_weights( model_freq );
}


PairingStatistics::PairingStatistics( core::scoring::dssp::StrandPairingSet const& topology ) {
	 add_topology( topology, "SINGLE_TOP" );
	 ModelFreq model_freq;
	 model_freq[ "SING" ] = 1;
	 compute_model_weights( model_freq );
}

PairingStatistics::~PairingStatistics() {}

core::Real PairingStatistics::strand_weight( core::scoring::dssp::StrandPairing const& pairing ) const {
	StatEntries::const_iterator itentry = entries_.find( pairing );
	if ( itentry != entries_.end() ) {
		return itentry->second.weight();
	}
	//  for ( StatEntries::const_iterator entry= entries_.begin(), eentry = entries_.end();
	//	entry != eentry; ++entry ) {
	//    if ( entry->compatible( pairing ) ) {
	//      return entry->weight();
	//    }
	//  }
  return 0.0;
}


core::Real PairingStatistics::weight( Model id ) const {
	for ( ModelWeight::const_iterator top = model_weight_.begin(); top != model_weight_.end(); ++top ) {
		if ( top->second==id ) return top->first;
	}
	utility_exit_with_message("Model name not known: " + id );
	return 0.0;
}

void PairingStatistics::compute_model_weights(  ModelFreq& model_freq ) {
	Real const contact_order_weight( basic::options::option[ basic::options::OptionKeys::jumps::contact_score ] );
  std::list< std::pair< core::Real, Model > > weight_list;
  for ( Topologies::const_iterator top = topols_.begin(), etop=topols_.end();
	top != etop; ++top ) {
    Real score( 0 );
    Real const norm( sqrt( (Real)  model_freq[ top->first.substr(0,4) ]  ) );
    // loop over pairings and find them in entries_
    for ( core::scoring::dssp::StrandPairingSet::const_iterator pairing = top->second.begin(),
	    epairing = top->second.end(); pairing != epairing; ++pairing ) {
      // find pairing in entries
			StatEntries::iterator itentry = entries_.find( *pairing );
			if ( itentry != entries_.end() ) {
				Real const weight( 1.0/norm * itentry->second.frequency() * ( itentry->second.size()-1.0
						+ contact_order_weight * ( std::max( 0, (int) itentry->second.contact_order() - 20)) ) );
				itentry->second.set_weight( weight );//so far only used for output later on
				score += weight;
			}
		}
//
//       for ( StatEntries::iterator entry= entries_.begin(), eentry = entries_.end();
// 	    entry != eentry; ++entry ) {
// 				if ( entry->compatible( *pairing ) ) {
// 					Real const contact_order_weight( basic::options::option[ basic::options::OptionKeys::jumps::contact_score ] );
// 					Real const weight( 1.0/norm * entry->frequency() * ( entry->size()-1.0
// 							+ contact_order_weight * ( std::max( 0, (int) entry->contact_order() - 20)) ) );
// 					entry->set_weight( weight );//so far only used for output later on
// 					score += weight;
// 					break;
// 				}
//       }
    //} // added score for each pairing
		weight_list.push_back( std::make_pair( score, top->first ) );
    //		weights_[ top->first ] = score; //also want the score in a map NAME --> score
  } // score for each model/topology
  weight_list.sort();
  weight_list.reverse();
  model_weight_.clear();
  copy( weight_list.begin(), weight_list.end(), std::back_inserter( model_weight_ ) );
}

core::scoring::dssp::StrandPairingSet const &
PairingStatistics::suggest_topology( std::string& topol_id ) const {
	runtime_assert( nr_models() );

  if ( !option[ templates::force_native_topology ] ) {
    Real const rel_rank_cutoff( option[ templates::topology_rank_cutoff ] );

	if( 1 > model_weight_.size() ){
		utility_exit_with_message( "The array model_weight_[] has 0 members, yet we need one. The topology is: " + topol_id );
	}
	Real const rank_cutoff( weight( 1 )*rel_rank_cutoff );
    Size top_max;
    Size const nr_models( model_weight_.size() );
    if ( rel_rank_cutoff >= 0.99 ) {
      top_max = std::min( 5, (int) nr_models);
    } else {
      for ( top_max = 1; (top_max < nr_models) && weight( top_max ) > rank_cutoff; top_max++ ) {};
    }
    Size const rg_topol ( static_cast< int >( RG.uniform() * top_max ) + 1 );
    topol_id = ranked_model( rg_topol );
    runtime_assert( topol_id != "BOGUS" );
    tr.Info << "cutoff: " << rank_cutoff << " use topology ** " << rg_topol << ": " << topol_id << " ** to select pairings " << std::endl;
    return topology( topol_id );
  }	else {
    return native_topology_;
  }
}

// void PairingStatistics::remove_stupid_strands() {

// }

std::ostream& operator<< ( std::ostream& out, PairingStatistics const& ps ) {
	out << "PAIRING_STATISTICS " << ps.model_weight_.size();
	//	out << ps.entries_;
//   out << "PAIRSTAT: -----------------------------\n";
//   for ( PairingStatistics::ModelFreq::const_iterator it = ps.model_freq_.begin(), eit = ps.model_freq_.end();
// 	it != eit; ++it ) {
//     out << it->first << " " << it->second << "\nPAIRSTAT: ";
//   }
 //  out << "\nPAIRSTAT: ";
//   for ( PairingStatistics::ModelWeight::const_iterator it = ps.model_weight_.begin(), eit = ps.model_weight_.end();
// 	it != eit; ++it ) {
//     out << it->first << " " << it->second  << "\nPAIRSTAT: ";
//   }
//   out << " -----------------------------\nPAIRSTAT: ";
  for ( PairingStatistics::ModelWeight::const_iterator it = ps.model_weight_.begin(), eit = ps.model_weight_.end();
	it != eit; ++it ) {
		out << "\nSTRAND_TOPOLOGY " << " " << ps.topology( it->second ).size() << " " << it->first << " " << it->second;
    for ( core::scoring::dssp::StrandPairingSet::const_iterator isp = ps.topology( it->second ).begin(),
	    eisp = ps.topology( it->second ).end(); isp != eisp; ++isp ) {
			if ( !isp->range_check() ) {
				tr.Error << "[ERROR] skip inconsistent pairing... " << *isp << std::endl;
				continue;
			}
			out << "\nPAIRSTAT_ENTRY: ";
      out << F(5,2, ps.strand_weight( *isp ) ) << " ";
      if ( ps.is_native_pairing( *isp ) ) out << " * ";
      else out << "  . ";
      out << *isp;
    }
  }

  return out;
}


void PairingStatEntry::show( std::ostream& out ) const {
	if ( !strand_pairing_.range_check() ) {
		tr.Error << "[ERROR] skip inconsistent pairing... " << strand_pairing_ << std::endl;
		return;
	}
  out << "PAIRSTAT_ENTRY: " << F(5,2,weight() ) << " " << frequency() << " " <<strand_pairing_ << " ";
  for ( ModelList::const_iterator it = models_.begin(), eit = models_.end();
	it != eit; ++it ) {
    out << *it << " ";
  }
}

std::istream& operator>> ( std::istream& is, PairingStatEntry& ps ) {
	std::string tag;
	is >> tag;
	if ( tag != "PAIRSTAT_ENTRY:" ) {
		tr.Trace << "failed reading PAIRSTAT_ENTRY: --- found instead: " << tag << std::endl;
		is.setstate( std::ios_base::failbit );
		return is;
	}
	Size nr_models;
	is >> ps.weight_ >> nr_models >> ps.strand_pairing_;
	for ( Size ct = 1; ct <= nr_models; ct ++) {
		is >> tag;
		ps.models_.push_back( tag );
	}
	return is;
}

std::ostream& operator<< (std::ostream& out, StatEntries const& ps ) {
  for ( StatEntries::const_iterator it = ps.begin(), eit = ps.end();
				it != eit; ++it ) {
		out << it->second << "\n";
  }
	return out;
}

std::istream& operator>> (std::istream& is, StatEntries& pslist) {
	pslist.clear();
	PairingStatEntry ps;
	while ( is >> ps ) { //not perfect because it fucks up the stream and the next line...
		pslist[ ps.pairing() ] = ps;
	}
	return is;
}

std::istream & operator>>( std::istream &is, PairingStatistics &ps) {
	std::string tag;
	Size ntops;
	is >> tag >> ntops;
	if ( tag != "PAIRING_STATISTICS" ) {
		tr.Trace << "failed reading PAIRING_STATISTIC --- found instead: " << tag << std::endl;
		is.setstate( std::ios_base::failbit );
		return is;
	}
	tr.Trace << " read " << ntops << " topologies from file... " << std::endl;
	ps.model_weight_.reserve( ntops+10 );

	//cause the hash-container to have at least ntops*3+10 buckets
	ps.entries_.rehash( ntops*3+10 );

	for ( Size ct_top = 1; ct_top <= ntops; ct_top++ ) {
		Size nstrand;
		Real model_weight;
		std::string model_ID;
		is >> tag >> nstrand >> model_weight >> model_ID;
		if ( tag != "STRAND_TOPOLOGY" ) {
			tr.Trace << "failed reading STRAND_TOPOLOGY --- found instead: " << tag << std::endl;
			is.setstate( std::ios_base::failbit );
			return is;
		}
		tr.Debug << "reading strand-topology " << model_ID << " with " << nstrand << " strands... "<<std::endl;
		core::scoring::dssp::StrandPairingSet sps;
		for ( Size ct = 1; ct <= nstrand; ct++ ) {
			core::scoring::dssp::StrandPairing pairing;
			Real weight;
			char ntag;
			std::string entry_tag;
			is >> entry_tag >> weight >> ntag >> pairing;
			if ( is.fail() ) return is;
			if ( !pairing.range_check() ) {
				tr.Error << "[ERROR] read inconsistent pairing in " << tag << " " << model_weight << " " << model_ID << std::endl;
				tr.Error << "offending pairing " << pairing << std::endl;
				is.setstate( std::ios_base::failbit );
				return is;
			}
			if ( entry_tag != "PAIRSTAT_ENTRY:" ) {
				tr.Trace << "failed reading PAIRSTAT_ENTRY: --- found instead: " << entry_tag << std::endl;
				is.setstate( std::ios_base::failbit );
				return is;
			}
			runtime_assert( pairing.range_check() );
			sps.push_back( pairing );

			//maintain also an extra list of all individual strand-pairings found... "PairingStatEntry"
			//too slow: option one, skip this condensing test -> check memory
			// otpion two,
			//	bool found( false );
			// for ( StatEntries::iterator try_entry= ps.entries_.begin(), eentry = ps.entries_.end();
// 						try_entry != eentry; ++try_entry ) {
// 				if ( try_entry->compatible( pairing ) ) {
// 					found = true;
// 					if ( try_entry->weight() != weight ) {
// 						tr.Warning << "inconsistent weights in topology " << *try_entry << std::endl;
// 						tr.Warning << "new weight is ignored: " << weight << " for strand " << pairing << "which had weight " << try_entry->weight() << std::endl;
// 					}
// 					try_entry->models().push_back( model_ID );
// 					break;
// 				} // if
// 			}
			StatEntries::iterator try_entry = ps.entries_.find( pairing );
			if ( try_entry != ps.entries_.end() ) {
				if ( try_entry->second.weight() != weight ) {
					tr.Warning << "inconsistent weights in topology " << try_entry->second << std::endl;
					tr.Warning << "new weight is ignored: " << weight << " for strand " << pairing << "which had weight " << try_entry->second.weight() << std::endl;
				}
				try_entry->second.models().push_back( model_ID );
			} else {
				PairingStatEntry entry( pairing, model_ID );
				entry.models().reserve( ntops );
				entry.set_weight( weight );
				ps.entries_[ pairing ] = entry;
			}

		} // finished reading this topology
		ps.topols_[ model_ID ]=sps;
		ps.model_weight_.push_back( std::make_pair( model_weight, model_ID ) );
	} //for all expected topologies
	return is;
}

}
}
