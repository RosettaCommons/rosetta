// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/antibody/grafting/scs_functor.cc
/// @brief Structural Component Selector (SCS): implemetation of predicates for filtering and sorting
/// @author Sergey Lyskov

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/scs_functor.hh>

#include <protocols/antibody/grafting/scs_blast.hh>
#include <protocols/antibody/grafting/exception.hh>
#include <protocols/antibody/grafting/scs_helper.hh>

#include <basic/Tracer.hh>

//#include <unordered_map>

#include <iomanip>

// Options
#include <basic/options/option.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>


namespace protocols {
namespace antibody {
namespace grafting {

using std::string;

using utility::CSI_Reset;  using utility::CSI_Bold;  using utility::CSI_Red;


static THREAD_LOCAL basic::Tracer TR("protocols.antibody.grafting");

// /// \details return true if s is in array
// inline bool is_in(string const &s, std::vector<string> const &array)
// {
// 	return std::find(array.begin(), array.end(), s) != array.end();
// }



void SCS_Comparator::apply(AntibodySequence const &antibody_sequence, SCS_ResultsOP results) const
{
	struct {
		string name;
		SCS_ResultVector & result;
	} J[] {
		{"frh", results->frh}, {"h1", results->h1}, {"h2", results->h2}, {"h3", results->h3},
		{"frl", results->frl}, {"l1", results->l1}, {"l2", results->l2}, {"l3", results->l3},
	};

	for(auto &j : J) {
		// std::function<bool(SCS_ResultOP const &, SCS_ResultOP const &)> f = [](SCS_ResultOP const &a, SCS_ResultOP const &b) {
		// 	SCS_BlastResult const *aa = dynamic_cast< SCS_BlastResult const *>( a.get() );
		// 	SCS_BlastResult const *bb = dynamic_cast< SCS_BlastResult const *>( b.get() );
		// 	if( aa and bb ) return compare(*aa, *bb);
		// 	else throw _AE_scs_failed_("SCS_BlastComparator::compare: Error! Could not cast SCS_Results to SCS_BlastResult!");
		// };

		std::sort(j.result.begin(), j.result.end(), [this, &antibody_sequence](SCS_ResultOP const &a, SCS_ResultOP const &b) { return compare(antibody_sequence, *a, *b); } );
	}
}


bool SCS_BlastComparator::compare(AntibodySequence const &antibody_sequence, SCS_Result const &a, SCS_Result const &b) const
{
	SCS_BlastResult const *aa = dynamic_cast< SCS_BlastResult const *>( &a );
	SCS_BlastResult const *bb = dynamic_cast< SCS_BlastResult const *>( &b );
	if( aa and bb ) return compare(antibody_sequence, *aa, *bb);
	else throw _AE_scs_failed_("SCS_BlastComparator::compare: Error! Could not cast SCS_Results to SCS_BlastResult!");
}




/// @details filter helper function: generate string with results sizes
string result_sizes(SCS_ResultsOP r, int width)
{
	std::stringstream s;

	s << std::left;
	s << "h1:" << std::setw(width) << r->h1.size() << " h2:" << std::setw(width) << r->h2.size() << " h3:" << std::setw(width) << r->h3.size() << ' ';
	s << "l1:" << std::setw(width) << r->l1.size() << " l2:" << std::setw(width) << r->l2.size() << " l3:" << std::setw(width) << r->l3.size() << ' ';

	s << "frh:" << std::setw(width) << r->frh.size() << " frl:" << std::setw(width) << r->frl.size() << ' ';

	s << "orientation:" << std::setw(width) << r->orientation.size();

	return s.str();
}


// // def filter_by_sequence_length(k, results, cdr_query, cdr_info):
// //     """ Template filter by sequence length
// //     """
// //     if Options.verbose: print 'filtering by sequence length...'
// //     for r in results[:]:
// //         pdb = r['subject-id']
// //         if   k in ['heavy','light_heavy']: return
// //         elif k in ['L1','L2','L3','H1','H2','H3'] and  len(cdr_query[k]) != len( cdr_info[pdb][k] ):
// //             results.remove(r)
// //             if Options.verbose and r not in results: print 'Filter sequence_length, removing:%s %s_query:%s %s_info:%s' % (pdb, k, len(cdr_query[k]), k, len( cdr_info[pdb][k] ))
// //         elif k == 'light' and 'light_length' in cdr_info[pdb]  and  len(cdr_query[k]) != int( cdr_info[pdb]['light_length'] ): results.remove(r)
// //         elif k == 'FRH':
// //             template_length = 61 if pdb == 'pdb2x7l_chothia.pdb' else 63
// //             if not len(cdr_query[k]) == template_length: results.remove(r)
// //             if Options.verbose and r not in results: print 'Filter sequence_length, removing:%s %s_query:%s %s_info:%s' % (pdb, k, len(cdr_query[k]), k, template_length)
// //             #elif k == 'FRH'  and  (not  len(cdr_query[k])-8 <= 67 <= len(cdr_query[k])+8): results.remove(r)
// //         elif k == 'FRL':
// //             template_length = 60 if pdb == 'pdb3h0t_chothia.pdb' else 58
// //             if not len(cdr_query[k]) == template_length: results.remove(r)
// //             if Options.verbose and r not in results: print 'Filter sequence_length, removing:%s %s_query:%s %s_info:%s' % (pdb, k, len(cdr_query[k]), k, template_length)
// //             #if not  len(cdr_query[k])-8 <= template_length <= len(cdr_query[k])+8: results.remove(r)


void SCS_BlastFilter_by_sequence_length::apply(AntibodySequence const &A, SCS_ResultsOP results) const
{
	TR.Debug << "SCS_BlastFilter_by_sequence_length: Results count before filtering  " << CSI_Red() << result_sizes(results) << CSI_Reset() << std::endl;

	struct {
		SCS_ResultVector &r;
		string query_sequence;
		string const SCS_BlastResult::*result_sequence;
	} h1_h2_h3_l1_l2_l3[] {
		{ results->h1, A.h1_sequence(), &SCS_BlastResult::h1 },
		{ results->h2, A.h2_sequence(), &SCS_BlastResult::h2 },
		{ results->h3, A.h3_sequence(), &SCS_BlastResult::h3 },
		{ results->l1, A.l1_sequence(), &SCS_BlastResult::l1 },
		{ results->l2, A.l2_sequence(), &SCS_BlastResult::l2 },
		{ results->l3, A.l3_sequence(), &SCS_BlastResult::l3 },
	};
	for(auto &region : h1_h2_h3_l1_l2_l3) {
		for(auto p = region.r.rbegin(); p != region.r.rend(); ) {
			SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
			if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_alignment_length::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

			if( region.query_sequence.size() != (br->*region.result_sequence).size() ) {
				TR.Trace << CSI_Red() << "SCS_BlastFilter_by_sequence_length: Filtering " << br->pdb << ": " << region.query_sequence.size() << "!=" << (br->*region.result_sequence).size() << "..." << CSI_Reset() << std::endl;
				region.r.erase( std::next(p++).base() );
			}
			else {
				// TR.Debug << "SCS_BlastFilter_by_sequence_length:" << br->pdb << ' ' << region.query_sequence << ' ' << br->*region.result_sequence << std::endl;
				++p;
			}
		}
	}

	for(auto p = results->frh.rbegin(); p != results->frh.rend(); ) {
		SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
		if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_alignment_length::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

		uint template_length = (br->pdb == "2x7l") ? 61 : 63;

		if( br->frh.size() != template_length ) {
			TR.Trace << CSI_Red() << "SCS_BlastFilter_by_sequence_length: Filtering " << br->pdb << "..." << CSI_Reset() << std::endl;
			results->frh.erase( std::next(p++).base() );
		}
		else ++p;
	}

	for(auto p = results->frl.rbegin(); p != results->frl.rend(); ) {
		SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
		if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_alignment_length::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

		uint template_length = (br->pdb == "3h0t") ? 60 : 58;

		if( br->frl.size() != template_length ) {
			TR.Trace << CSI_Red() << "SCS_BlastFilter_by_sequence_length: Filtering " << br->pdb << "..." << CSI_Reset() << std::endl;
			results->frl.erase( std::next(p++).base() );
		}
		else ++p;
	}

	TR.Debug << "SCS_BlastFilter_by_sequence_length: Results count after filtering   " << CSI_Red() << result_sizes(results) << CSI_Reset() << std::endl;
}




// // def filter_by_alignment_length(k, results, cdr_query, cdr_info):
// //     """ Template filter by alignment length
// //     """
// //     for r in results[:]:
// //         pdb = r['subject-id']
// //         if k == 'H3' and int(r['alignment-length']) < 0.10 * len( cdr_info[pdb]['H3'] ): results.remove(r)
// //         if k == 'H2' and int(r['alignment-length']) < 0.55 * len( cdr_info[pdb]['H2'] ): results.remove(r)
// //         if k in ['L1', 'L2', 'L3', 'H1']  and  int(r['alignment-length']) < 0.70 *  len( cdr_info[pdb][k] ): results.remove(r)
// //         if Options.verbose and r not in results: print 'Filter alignment_length removing:%s %s_query:%s alignment-length:%s ' % (pdb, k, len(cdr_info[pdb][k]), r['alignment-length'])
void SCS_BlastFilter_by_alignment_length::apply(AntibodySequence const &/*antibody_sequence*/, SCS_ResultsOP results) const
{
	TR.Debug << "SCS_BlastFilter_by_alignment_length: Results count before filtering " << CSI_Red() << result_sizes(results) << CSI_Reset() << std::endl;

	struct {
		SCS_ResultVector &r;
		double v;
		string const SCS_BlastResult::*sequence;  //std::function< string const &(SCS_BlastResult const&) > sequence; //[](SCS_BlastResult const &r) -> string const & { return r.h1; }
	} R[] {
		{ results->h1, 0.70, &SCS_BlastResult::h1 },
		{ results->h2, 0.55, &SCS_BlastResult::h2 },
		{ results->h3, 0.10, &SCS_BlastResult::h3 },
		{ results->l1, 0.70, &SCS_BlastResult::l1 },
		{ results->l2, 0.70, &SCS_BlastResult::l2 },
		{ results->l3, 0.70, &SCS_BlastResult::l3 },
	};

	for(auto &region : R) {
		for(auto p = region.r.rbegin(); p != region.r.rend(); ) {
			SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
			if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_alignment_length::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

			//TR << "alignment_length: " << br->alignment_length << "   " << region.v * (br->*region.sequence).size() << std::endl;

			if( br->alignment_length < region.v * (br->*region.sequence).size() ) {
				TR.Trace << CSI_Red() << "SCS_BlastFilter_by_alignment_length: Filtering " << br->pdb << "..." << CSI_Reset() << std::endl;
				region.r.erase( std::next(p++).base() );
			}
			else ++p;
		}
	}
	TR.Debug << "SCS_BlastFilter_by_alignment_length: Results count after filtering  " << CSI_Red() << result_sizes(results) << CSI_Reset() << std::endl;
}

//  def filter_by_template_resolution(k, results, cdr_query, cdr_info):
//  """ Template filter by resolution
//  """
//  for r in results[:]:
//    pdb = r['subject-id']
//    if float(cdr_info[pdb]['resolution']) > 2.8: results.remove(r)
//    if Options.verbose and r not in results: print 'Filter template_resolution, removing:%s resolution:%s' % (pdb, cdr_info[pdb]['resolution'])

SCS_BlastFilter_by_template_resolution::SCS_BlastFilter_by_template_resolution() {
	set_resolution_cutoff( 2.8 );
}

core::Real SCS_BlastFilter_by_template_resolution::get_resolution_cutoff() const {
	return resolution_cutoff_;
}

void SCS_BlastFilter_by_template_resolution::set_resolution_cutoff(core::Real cutoff) {
	resolution_cutoff_=cutoff;
}

void SCS_BlastFilter_by_template_resolution::apply(AntibodySequence const &/*antibody_sequence*/, SCS_ResultsOP results) const
{
	TR.Debug << "SCS_BlastFilter_by_template_resolution: Results count before filtering " << CSI_Red() << result_sizes(results) << CSI_Reset() << std::endl;

	struct {
		SCS_ResultVector &r;
	} R[] {
		{ results->h1 },
		{ results->h2 },
		{ results->h3 },
		{ results->l1 },
		{ results->l2 },
		{ results->l3 },
	};

	for(auto &region : R) { // loop over all members of array of structs, R?
		for(auto p = region.r.rbegin(); p != region.r.rend(); ) { // for each region, h1-l3, loop over all alignments
			SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() ); // use iterator to "get" alignment result and cast from Result to BlastResult
			if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_template_resolution::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

			//TR << "resolution: " << br->resolution << "   " << region.v * (br->*region.sequence).size() << std::endl;

			if( br->resolution > get_resolution_cutoff() ) {
				TR.Trace << CSI_Red() << "SCS_BlastFilter_by_template_resolution: Filtering " << br->pdb << "..." << CSI_Reset() << std::endl;
				region.r.erase( std::next(p++).base() );
			}
			else ++p;
		}
	}

	// Filter heavy framework region by resolution
	for(auto p = results->frh.rbegin(); p != results->frh.rend(); ) {
		SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
		if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_sequence_identity::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

		if( br->resolution > get_resolution_cutoff() ) {
			TR.Trace << CSI_Red() << "SCS_BlastFilter_by_template_resolution: Filtering " << br->pdb << "..." << CSI_Reset() << std::endl;
			results->frh.erase( std::next(p++).base() );
		}
		else ++p;
	}

	// Filter light framework region by resolution
	for(auto p = results->frl.rbegin(); p != results->frl.rend(); ) {
		SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
		if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_sequence_identiy::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

		if( br->resolution > get_resolution_cutoff() ) {
			TR.Trace << CSI_Red() << "SCS_BlastFilter_by_template_resolution: Filtering " << br->pdb << "..," << CSI_Reset() << std::endl;
			results->frl.erase( std::next(p++).base() );
		}
		else ++p;
	}

	// Filter orientation by resolution
	for(auto p = results->orientation.rbegin(); p != results->orientation.rend(); ) {
		SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
		if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_sequence_identiy::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

		if( br->resolution > get_resolution_cutoff() ) {
			TR.Trace << CSI_Red() << "SCS_BlastFilter_by_template_resolution: Filtering " << br->pdb << "..," << CSI_Reset() << std::endl;
			results->orientation.erase( std::next(p++).base() );
		}
		else ++p;
	}

	TR.Debug << "SCS_BlastFilter_by_template_resolution: Results count after filtering  " << CSI_Red() << result_sizes(results) << CSI_Reset() << std::endl;
}
  //  def filter_by_sequence_homolog(k, results, cdr_query, cdr_info):
  //  """ Template filter for sequence identity
  //  """
  //
  //  if Options.verbose: print 'filtering by sequence identity...'
  //#print results
  //
  //    for r in results[:]:
  //      pdb = r['subject-id']
  //
  //      if k in ['L1','L2','L3','H1','H2','H3'] and sid_checker(cdr_query[k], cdr_info[pdb][k]) >= sid_cutoff_cdr:
  //      results.remove(r)
  //      if Options.verbose and r not in results: print 'Filter sequence_identity (%s%% cut-off), removing:%s %s_query:%s %s_info:%s %s%% identity' % (sid_cutoff_cdr, pdb, k, len(cdr_query[k]), k, len( cdr_info[pdb][k] ), round(sid_checker(cdr_query[k], cdr_info[pdb][k]), 2) )
  //
  //      elif  k in ['FRL'] and sid_checker(cdr_query[k], frl_info[pdb][k]) >= sid_cutoff_fr:
  //      results.remove(r)
  //      if Options.verbose and r not in results: print 'Filter sequence_identity (%s%% cut-off), removing:%s %s_query:%s %s_info:%s %s%% identity' % (sid_cutoff_fr, pdb, k, len(cdr_query[k]), k, len( frl_info[pdb][k] ), round(sid_checker(cdr_query[k], frl_info[pdb][k]), 2) )
  //
  //      elif  k in ['FRH'] and sid_checker(cdr_query[k], frh_info[pdb][k]) >= sid_cutoff_fr:
  //      results.remove(r)
  //      if Options.verbose and r not in results: print 'Filter sequence_identity (%s%% cut-off), removing:%s %s_query:%s %s_info:%s %s%% identity' % (sid_cutoff_fr, pdb, k, len(cdr_query[k]), k, len( frh_info[pdb][k] ), round(sid_checker(cdr_query[k], frh_info[pdb][k]), 2) )
  //      elif  k in ['light_heavy'] and float(r['%-identity'])  >= sid_cutoff_fr - 5:
  //      results.remove(r)
  //      if Options.verbose and r not in results: print 'Filter sequence_identity (%s%% cut-off), removing:%s %s %s%% identity' % (sid_cutoff_fr, pdb, k, r['%-identity'])
  //

/// @details Helper Function: returns fraction of identical residues in sequences of identical length
///          otherwise, the function returns 0. Modeled on Python equivalent.
core::Real sid_checker(std::string seq_q, std::string seq_t)
{
  core::Real ratio = 0; // default, if lengths are unequal
  core::Size id_counter = 0;

  if( seq_q.size() == seq_t.size() and seq_q.size()) {
    for( core::Size i=0; i!=seq_q.size(); ++i ) {
      if( seq_q[i] == seq_t[i] ){
        id_counter++;
      }
    }
    ratio = 100.*id_counter/seq_q.size();
  }
  return ratio;
}

SCS_BlastFilter_by_sequence_identity::SCS_BlastFilter_by_sequence_identity() {
	init_from_options();
}

core::Real SCS_BlastFilter_by_sequence_identity::get_sid_cutoff_cdr() const {
	return sid_cutoff_cdr_;
}

core::Real SCS_BlastFilter_by_sequence_identity::get_sid_cutoff_fr() const {
	return sid_cutoff_fr_;
}

void SCS_BlastFilter_by_sequence_identity::set_sid_cutoff_cdr(core::Real cutoff) {
	sid_cutoff_cdr_=cutoff;
}

void SCS_BlastFilter_by_sequence_identity::set_sid_cutoff_fr(core::Real cutoff) {
	sid_cutoff_fr_=cutoff;
}

void SCS_BlastFilter_by_sequence_identity::init_from_options() {

	using namespace basic::options;

	set_sid_cutoff_cdr( 100.0 ); // default if flag is not given, i.e. do not filter
	set_sid_cutoff_fr( 100.0 ); // default if flag is not given, i.e. do not filter

	if ( option[ basic::options::OptionKeys::antibody::exclude_homologs ] ) {
		// note the default for the above flag is false and the below flags is 80.0 %
		set_sid_cutoff_cdr( option[ basic::options::OptionKeys::antibody::exclude_homologs_cdr_cutoff]() );
		set_sid_cutoff_fr( option[ basic::options::OptionKeys::antibody::exclude_homologs_fr_cutoff]() );
	}

}

void SCS_BlastFilter_by_sequence_identity::apply(AntibodySequence const& A,
  SCS_ResultsOP results) const
{

	core::Real sid_ratio;

  TR.Debug << "SCS_BlastFilter_by_sequence_identity: Results count before filtering " << CSI_Red() << result_sizes(results) << CSI_Reset() << std::endl;

  struct {
    SCS_ResultVector &r;
    string query_sequence; // maybe do sid_cutoff as a struct in scs_blast, storing by region?
    string const SCS_BlastResult::*result_sequence;
  } h1_h2_h3_l1_l2_l3[] {
    { results->h1, A.h1_sequence(), &SCS_BlastResult::h1 },
    { results->h2, A.h2_sequence(), &SCS_BlastResult::h2 },
    { results->h3, A.h3_sequence(), &SCS_BlastResult::h3 },
    { results->l1, A.l1_sequence(), &SCS_BlastResult::l1 },
    { results->l2, A.l2_sequence(), &SCS_BlastResult::l2 },
    { results->l3, A.l3_sequence(), &SCS_BlastResult::l3 },
  };

  // Filter CDRs by sequence identity
  for(auto &region : h1_h2_h3_l1_l2_l3) {
    for(auto p = region.r.rbegin(); p != region.r.rend(); ) {
      SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
      if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_sequence_identity::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

			sid_ratio = sid_checker( region.query_sequence, br->*region.result_sequence );
      if( sid_ratio > get_sid_cutoff_cdr() ) {
        TR.Trace << CSI_Red() << "SCS_BlastFilter_by_sequence_identity: Filtering " << br->pdb << "... with SID ratio of " << sid_ratio << CSI_Reset() << std::endl;
        region.r.erase( std::next(p++).base() );
      }
	  else ++p;

      // else {
      // 	TR.Debug << "SCS_BlastFilter_by_sequence_length:" << br->pdb << ' ' << region.query_sequence << ' ' << br->*region.result_sequence << std::endl;
      // }
    }
  }

	// Get FRH, FRL, orientation from input sequence
	FRH_FRL query_fr( calculate_frh_frl(A) );

	string query_frh = query_fr.frh1 + query_fr.frh2 + query_fr.frh3 + query_fr.frh4;
	string query_frl = query_fr.frl1 + query_fr.frl2 + query_fr.frl3 + query_fr.frl4;

	//string query_heavy_orientation = query_fr.frh1 + A.h1_sequence() + query_fr.frh2 + A.h2_sequence() + query_fr.frh3 + A.h3_sequence() + query_fr.frh4;
	//string query_light_orientation = query_fr.frl1 + A.l1_sequence() + query_fr.frl2 + A.l2_sequence() + query_fr.frl3 + A.l3_sequence() + query_fr.frl4;

	//string query_orientation = query_light_orientation + query_heavy_orientation; // Note: must be in this order

	// Filter heavy framework region by sequence identity
  for(auto p = results->frh.rbegin(); p != results->frh.rend(); ) {
    SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
    if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_sequence_identity::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

		sid_ratio = sid_checker( query_frh, br->frh);
    if( sid_ratio > get_sid_cutoff_fr() ) {
        TR.Trace << CSI_Red() << "SCS_BlastFilter_by_sequence_identity: Filtering " << br->pdb << "... with SID ratio of " << sid_ratio << CSI_Reset() << std::endl;
      results->frh.erase( std::next(p++).base() );
    }
	else ++p;
  }

	// Filter light framework region by sequence identity
  for(auto p = results->frl.rbegin(); p != results->frl.rend(); ) {
    SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
    if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_sequence_identiy::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

		sid_ratio = sid_checker( query_frl, br->frl);
    if( sid_ratio > get_sid_cutoff_fr() ) {
        TR.Trace << CSI_Red() << "SCS_BlastFilter_by_sequence_identity: Filtering " << br->pdb << "... with SID ratio of " << sid_ratio << CSI_Reset() << std::endl;
      results->frl.erase( std::next(p++).base() );
    }
	else ++p;
  }

	// Filter orientation by sequence identity?
	// How? Cannot reassemble orientation sequence from result.
  TR.Debug << "SCS_BlastFilter_by_sequence_identity: Results count after filtering   " << CSI_Red() << result_sizes(results) << CSI_Reset() << std::endl;

}

  //    def filter_by_outlier(k, results, cdr_query, cdr_info):
  //    """ Template filter by outlier
  //
  //    FIXME: Some better description may help
  //    """
  //
  //    outlier = {}
  //    for line in file( _script_path_ + '/info/outlier_list' ): outlier[tuple(line.split()[:2])] = line.split()[2] == 'true'
  //    for r in results[:]:
  //    pdb = r['subject-id']
  //
  //    if outlier.get( (pdb, k), False ): results.remove(r)
  //    if Options.verbose and r not in results: print 'Filter outlier, removing:%s' % pdb
  //
  //    def filter_by_orientational_distance(k, results, cdr_query, cdr_info):
  //    """ Template filter by operational distance
  //
  //    This filter is trajectory dependent; should be kept as last filter implemented
  //
  //    """

void SCS_BlastFilter_by_outlier::apply(AntibodySequence const& /* A */,
                                                 SCS_ResultsOP results) const
{

  TR.Debug << "SCS_BlastFilter_by_outlier: Results count before filtering " << CSI_Red() << result_sizes(results) << CSI_Reset() << std::endl;

  // Get outlier data from file
  // regions: "FRL", "FRH", "L1", "L2", "L3", "H1", "H2"
  std::map< std::string, std::map<std::string, bool> > outlier_map;
  outlier_map = SCS_Helper::get_ab_region_outliers();

  struct {
    SCS_ResultVector &r;
    string ab_region;
  } frh_h1_h2_frl_l1_l2_l3[] {
    { results->frh, "FRH" },
    { results->h1, "H1" },
    { results->h2, "H2" },
    { results->frl, "FRL" },
    { results->l1, "L1" },
    { results->l2, "L2" },
    { results->l3, "L3" },
  };

  for(auto &region : frh_h1_h2_frl_l1_l2_l3) {
    for(auto p = region.r.rbegin(); p != region.r.rend(); ) {
      SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
      if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_outlier::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

			// check if pdb is contained in outlier list... somewhat worrying that list doesn't contain all pdbs
			if( outlier_map.find(br->pdb) == outlier_map.end() ) {
				TR.Trace << "SCS_BlastFilter_by_outlier: Could not find " << br->pdb << " in outlier list." << std::endl;
				++p;
			}
			else if( outlier_map.at(br->pdb).at(region.ab_region) ) {
				TR.Trace << CSI_Red() << "SCS_BlastFilter_by_outlier: Filtering " << br->pdb << CSI_Reset() << std::endl;
				region.r.erase( std::next(p++).base() );
			}
			else ++p;
    }
  }

  TR.Debug << "SCS_BlastFilter_by_outlier: Results count after filtering   " << CSI_Red() << result_sizes(results) << CSI_Reset() << std::endl;

}

  //    def filter_by_template_bfactor(k, results, cdr_query, cdr_info):
  //    """ Template filter by bfactor
  //    """
  //    bfactor = {}
  //
  //    if k in ['L1', 'L2', 'L3', 'H1', 'H2', 'H3']:
  //    for line in file( _script_path_ + '/info/list_bfactor50' ):
  //    for i, e in enumerate(['L1', 'L2', 'L3', 'H1', 'H2', 'H3']):
  //    if k == e: bfactor[ 'pdb'+line.split()[0]+'_chothia.pdb' ] = line.split()[i+1] == 'True'
  //
  //    for r in results[:]:
  //    pdb = r['subject-id']
  //    if bfactor.get( pdb, False ): results.remove(r)
  //    if Options.verbose and r not in results: print 'Filter B-factor50, removing:%s' % pdb
  //
  //


void SCS_BlastFilter_by_template_bfactor::apply(AntibodySequence const& /* A */,
                                                 SCS_ResultsOP results) const
{
	TR.Debug << "SCS_BlastFilter_by_template_bfactor: Results count before filtering " << CSI_Red() << result_sizes(results) << CSI_Reset() << std::endl;

	// Get bfactor data from file
	// regions: "l1", "l2", "l3", "h1", "h2", "h3"
	std::map< std::string, std::map<std::string, bool> > bfactor_map;
	bfactor_map = SCS_Helper::get_ab_cdr_bfactors();

	struct {
		SCS_ResultVector &r;
		string ab_region;
	} h1_h2_h3_l1_l2_l3[] {
		{ results->h1, "h1" },
		{ results->h2, "h2" },
		{ results->h3, "h3" },
		{ results->l1, "l1" },
		{ results->l2, "l2" },
		{ results->l3, "l3" },
	};

	for(auto &region : h1_h2_h3_l1_l2_l3) {
		for(auto p = region.r.rbegin(); p != region.r.rend(); ) {
			SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
			if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_template_bfactor::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

			// check if pdb is contained in bfactor list... somewhat worrying that list doesn't contain all pdbs
			if( bfactor_map.find(br->pdb) == bfactor_map.end() ) {
				TR.Trace << "SCS_BlastFilter_by_template_bfactor: Could not find " << br->pdb << " in outlier list." << std::endl;
				++p;
			}
			else if( bfactor_map.at(br->pdb).at(region.ab_region) ) {
				TR.Trace << CSI_Red() << "SCS_BlastFilter_by_template_bfactor: Filtering " << br->pdb << CSI_Reset() << std::endl;
				region.r.erase( std::next(p++).base() );
			}
			else ++p;
			// else {
			// 	TR.Debug << "SCS_BlastFilter_by_sequence_length:" << br->pdb << ' ' << region.query_sequence << ' ' << br->*region.result_sequence << std::endl;
			// }
		}
	}

	TR.Debug << "SCS_BlastFilter_by_template_mbfactor: Results count after filtering   " << CSI_Red() << result_sizes(results) << CSI_Reset() << std::endl;
}

  //    if Options.verbose: print 'filtering by orientational distance...'
  //      if k in ['light_heavy']:
  //        for i in range(0, Options.number_of_templates):
  //          if (i+1) == len(results):
  //            print "Warning: may not be enough distinct light_heavy orientations; some may be used repeatedly"
  //            break
  //            top_pdb = results[i]['subject-id']
  //            (grep_status, grep_output) = commands.getstatusoutput("grep " + top_pdb + " " + script_dir + "/comparisons.txt")
  //            orientational_dictionary = {}
  //            (keystring, valstring) = grep_output.split("\n")
  //            orient_keys = keystring.split(" ")
  //            orient_vals = valstring.split(" ")
  //            for j in range(1, len(orient_keys)):
  //            orientational_dictionary[orient_keys[j]] = orient_vals[j]
  //            for r in results[(i+1):]:
  //            pdb = r['subject-id']
  //            if float(orientational_dictionary[pdb]) < float(Options.orientational_distance_cutoff):
  //            results.remove(r)
  //            if (i+1) == len(results):
  //            print "Warning: may not be enough distinct light_heavy orientations; some may be used repeatedly"
  //            break
  //            else: return
  //

SCS_BlastFilter_by_OCD::SCS_BlastFilter_by_OCD() {

	init_from_options();
}

core::Size SCS_BlastFilter_by_OCD::get_n_orientational_templates() const{
	return n_orientational_templates_;
}

core::Real SCS_BlastFilter_by_OCD::get_ocd_cutoff() const{
	return ocd_cutoff_;
}

void SCS_BlastFilter_by_OCD::set_n_orientational_templates(core::Size n) {
	n_orientational_templates_=n;
}

void SCS_BlastFilter_by_OCD::set_ocd_cutoff(core::Real cutoff) {
	ocd_cutoff_=cutoff;
}

void SCS_BlastFilter_by_OCD::init_from_options(){

	using namespace basic::options;

	// note the defaults are 10 and 0.5 respectively
	set_n_orientational_templates( option[ basic::options::OptionKeys::antibody::n_multi_templates ]() );
	set_ocd_cutoff( option[ basic::options::OptionKeys::antibody::ocd_cutoff ]() );

}

///@details Filter by OCD will remove all templates of lower bit score that are within X.X OCD of the "current" template.
///         This may result in less than the number of templates requested, so templates will repeat!
void SCS_BlastFilter_by_OCD::apply(AntibodySequence const& /* A */,
                                                 SCS_ResultsOP results) const
{
	// The idea here is to iterate until the number of multiple templates is found
	// if not, repeat templates.

	SCS_ResultVector &r = results->orientation;

	TR.Debug << "SCS_BlastFilter_by_OCD: Results count before filtering " << CSI_Red() << result_sizes(results) << CSI_Reset() << std::endl;

	// Get OCD data from file
	std::map< std::string, std::map<std::string, core::Real> > ocd_map;
	ocd_map = SCS_Helper::get_ab_OCDs();

	// may work better as a for loop as in the python code
	for (core::Size j=0; j < get_n_orientational_templates(); ++j) { // we want to prune templates within x OCD of the top N

		if ( j+1 == r.size() ) {
			TR.Warning << TR.Red << "SCS_BlastFilter_by_OCD: there may not be enough distinct light_heavy orientations! Some may be used repeatedly!" << std::endl;
			break;
		}

		auto p = r.begin()+j; //get best aligned model, then second best, then third best, then ...
		SCS_BlastResult const *top_br = dynamic_cast< SCS_BlastResult const *>( p->get() );
		if( !top_br ) throw _AE_scs_failed_("SCS_BlastFilter_by_OCD::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

		for ( core::Size i=j+1; i < r.size(); ++i ) { // loop over result vector, but update size each iteration as we delete things?

			auto p = r.begin()+i;
			SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
			if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_OCD::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

			if ( ocd_map.at(top_br->pdb).at(br->pdb) < get_ocd_cutoff() ) {
				TR.Trace << CSI_Red() << "SCS_BlastFilter_by_OCD: Filtering " << br->pdb << " with OCD: " << ocd_map.at(top_br->pdb).at(br->pdb) << CSI_Reset() << std::endl;
				r.erase( r.begin()+i );
				--i; //check index position again after removing template
			}
			else {
				TR.Trace << CSI_Red() << "SCS_BlastFilter_by_OCD: Not filtering " << br->pdb << " with OCD: " << ocd_map.at(top_br->pdb).at(br->pdb) << CSI_Reset() << std::endl;
			}

			// check to make sure we're not out of templates, break if we are
			// this gets trigged once number of results is equal to number of templates
			if ( i==r.size() ) {
				if ( i+1>=get_n_orientational_templates() ) break; // in the case we do have enough
				TR.Warning << TR.Red << "SCS_BlastFilter_by_OCD: there may not be enough distinct light_heavy orientations! Some may be used repeatedly!" << std::endl;
				break;
			}

		}

	}

	TR.Debug << "SCS_BlastFilter_by_OCD: Results count after filtering   " << CSI_Red() << result_sizes(results) << CSI_Reset() << std::endl;
}

} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__
