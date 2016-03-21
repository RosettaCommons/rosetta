// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/antibody/grafting/scs_functor.cc
/// @brief Structural Component Selector (SCS): implemetation of predicates for filtering and sorting
/// @author Sergey Lyskov

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/scs_functor.hh>

#include <protocols/antibody/grafting/scs_blast.hh>
#include <protocols/antibody/grafting/exception.hh>

#include <basic/Tracer.hh>

//#include <unordered_map>

#include <iomanip>

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
		SCS_ResultsVector & result;
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
	TR.Debug << "SCS_BlastFilter_by_sequence_length: Results count before filtering  " << CSI_Red << result_sizes(results) << CSI_Reset << std::endl;

	struct {
		SCS_ResultsVector &r;
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
		for(auto p = region.r.rbegin(); p != region.r.rend(); ++p) {
			SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
			if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_alignment_length::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

			if( region.query_sequence.size() != (br->*region.result_sequence).size() ) {
				TR.Trace << CSI_Red << "SCS_BlastFilter_by_sequence_length: Filtering " << br->pdb << "..." << CSI_Reset << std::endl;
				region.r.erase( std::next(p).base() );
			}
			// else {
			// 	TR.Debug << "SCS_BlastFilter_by_sequence_length:" << br->pdb << ' ' << region.query_sequence << ' ' << br->*region.result_sequence << std::endl;
			// }
		}
	}

	for(auto p = results->frh.rbegin(); p != results->frh.rend(); ++p) {
		SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
		if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_alignment_length::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

		uint template_length = (br->pdb == "2x7l") ? 61 : 63;

		if( br->frh.size() != template_length ) {
			TR.Trace << CSI_Red << "SCS_BlastFilter_by_sequence_length: Filtering " << br->pdb << "..." << CSI_Reset << std::endl;
			results->frh.erase( std::next(p).base() );
		}
	}

	for(auto p = results->frl.rbegin(); p != results->frl.rend(); ++p) {
		SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
		if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_alignment_length::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

		uint template_length = (br->pdb == "3h0t") ? 60 : 58;

		if( br->frl.size() != template_length ) {
			TR.Trace << CSI_Red << "SCS_BlastFilter_by_sequence_length: Filtering " << br->pdb << "..." << CSI_Reset << std::endl;
			results->frl.erase( std::next(p).base() );
		}
	}

	TR.Debug << "SCS_BlastFilter_by_sequence_length: Results count after filtering   " << CSI_Red << result_sizes(results) << CSI_Reset << std::endl;
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
	TR.Debug << "SCS_BlastFilter_by_alignment_length: Results count before filtering " << CSI_Red << result_sizes(results) << CSI_Reset << std::endl;

	struct {
		SCS_ResultsVector &r;
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
		for(auto p = region.r.rbegin(); p != region.r.rend(); ++p) {
			SCS_BlastResult const *br = dynamic_cast< SCS_BlastResult const *>( p->get() );
			if( !br ) throw _AE_scs_failed_("SCS_BlastFilter_by_alignment_length::apply: Error! Could not cast SCS_Results to SCS_BlastResult!");

			//TR << "alignment_length: " << br->alignment_length << "   " << region.v * (br->*region.sequence).size() << std::endl;

			if( br->alignment_length < region.v * (br->*region.sequence).size() ) {
				TR.Trace << CSI_Red << "SCS_BlastFilter_by_alignment_length: Filtering " << br->pdb << "..." << CSI_Reset << std::endl;
				region.r.erase( std::next(p).base() );
			}
		}
	}
	TR.Debug << "SCS_BlastFilter_by_alignment_length: Results count after filtering  " << CSI_Red << result_sizes(results) << CSI_Reset << std::endl;
}


} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__
