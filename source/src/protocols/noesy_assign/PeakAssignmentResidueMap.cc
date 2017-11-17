// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/CrossPeakList.hh>
#include <protocols/noesy_assign/PeakAssignmentResidueMap.hh>
#include <protocols/noesy_assign/ResonanceList.hh>

// Package Headers
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>
#include <protocols/noesy_assign/Exceptions.hh>
#include <protocols/noesy_assign/util.hh>
// Project Headers


// Utility headers

// #include <utility/excn/Exceptions.hh>
// #include <utility/vector1.fwd.hh>
// #include <utility/pointer/ReferenceCount.hh>
// #include <numeric/numeric.functions.hh>
// #include <core/util/prof.hh>
#include <basic/Tracer.hh>
// #include <core/options/option.hh>
// #include <core/options/keys/abinitio.OptionKeys.gen.hh>
// #include <core/options/keys/run.OptionKeys.gen.hh>
//#include <core/options/keys/templates.OptionKeys.gen.hh>
#include <basic/prof.hh>
#include <ObjexxFCL/FArray2D.hh>

//// C++ headers
#include <cstdlib>
#include <string>

#ifdef WIN32
#include <iterator>
#endif

#include <protocols/noesy_assign/CrossPeak.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <cmath>
#include <set>


static basic::Tracer tr( "protocols.noesy_assign.assignments" );

using core::Real;
using namespace core;
using namespace basic;
//using namespace core::options;
//using namespace OptionKeys;

namespace protocols {
namespace noesy_assign {

PeakAssignmentResidueMap::PeakAssignmentResidueMap() : residues_(), BOGUS_ASSIGNMENTS(), atoms_() {}
PeakAssignmentResidueMap::~PeakAssignmentResidueMap() = default;

void PeakAssignmentResidueMap::add( PeakAssignmentOP const& assignment ) {
	Size const res1( assignment->resid( 1 ) );
	Size const res2( assignment->resid( 2 ) );
	while ( residues_.size() < res1 ) {
		residues_.push_back( PeakAssignmentMap() );
	}
	auto res2_entry( residues_[ res1 ].insert( PeakAssignmentMap::value_type( res2, PeakAssignments() ) ).first );
	//res2_entry points now to an element in the res2 indexed map
	res2_entry->second.push_back( assignment ); //currently no check for redundancy!

	while ( atoms_.size() < res1 ) {
		atoms_.push_back( AtomList() );
	}
	while ( atoms_.size() < res2 ) {
		atoms_.push_back( AtomList() );
	}
	atoms_[ res1 ].insert( assignment->atom( 1 ) );
	atoms_[ res2 ].insert( assignment->atom( 2 ) );

}

void PeakAssignmentResidueMap::add_all_atoms( ResonanceList const& rslist ) {
	for ( auto const & it : rslist ) {
		Size res1( it.second->resid() );
		while ( atoms_.size() < res1 ) {
			atoms_.push_back( AtomList() );
		}
		atoms_[ res1 ].insert( it.second->atom() );
	}
}
/// @brief remove assignment...
void PeakAssignmentResidueMap::remove( PeakAssignment const& assignment ) {
	Size const res1( assignment.resid( 1 ) );
	Size const res2( assignment.resid( 2 ) );

	PeakAssignments& res2_entry( assignments( res1, res2 ) );
	bool success( false );
	for ( auto it = res2_entry.begin(); !success && it != res2_entry.end(); ++it ) {
		PeakAssignment const& assi( **it );
		if ( assi == assignment ) {
			res2_entry.erase( it );
			success = true;
		}
	}
	if ( !success ) {
		throw CREATE_EXCEPTION(EXCN_AssignmentNotFound,  assignment, "remove: PeakAssignment not found -- no entry with exact assignment" );
	}
}

PeakAssignmentResidueMap::PeakAssignments const& PeakAssignmentResidueMap::assignments( core::Size res1, core::Size res2 ) const {
	basic::ProfileThis doit( basic::NOESY_ASSIGN_NETWORK_FIND_RAW_ASSIGN );
	if ( residues_.size() < res1 || res1 <= 0 ) {
		throw CREATE_EXCEPTION(EXCN_AssignmentNotFound,  BOGUS_ASSIGNMENT, "PeakAssignment not found -- no entry with res1 " + ObjexxFCL::string_of( res1 ) );
	}
	auto res2_entry( residues_[ res1 ].find( res2 ) );
	if ( res2_entry ==  residues_[ res1 ].end() ) {
		throw CREATE_EXCEPTION(EXCN_AssignmentNotFound,  BOGUS_ASSIGNMENT, "PeakAssignment not found -- no entry with res2 " + ObjexxFCL::string_of( res2 ) );
	}
	return res2_entry->second;
}

PeakAssignmentResidueMap::PeakAssignments& PeakAssignmentResidueMap::assignments( core::Size res1, core::Size res2 ) {
	basic::ProfileThis doit( basic::NOESY_ASSIGN_NETWORK_FIND_RAW_ASSIGN );
	if ( residues_.size() < res1 || res1 <= 0 ) {
		throw CREATE_EXCEPTION(EXCN_AssignmentNotFound,  BOGUS_ASSIGNMENT, "PeakAssignment not found -- no entry with res1 " + ObjexxFCL::string_of( res1 ) );
	}
	auto res2_entry( residues_[ res1 ].find( res2 ) );
	if ( res2_entry ==  residues_[ res1 ].end() ) {
		throw CREATE_EXCEPTION(EXCN_AssignmentNotFound,  BOGUS_ASSIGNMENT, "PeakAssignment not found -- no entry with res2 " + ObjexxFCL::string_of( res2 ) );
	}
	return res2_entry->second;
}

PeakAssignmentResidueMap::PeakAssignments& PeakAssignmentResidueMap::_assignments( core::Size res1, core::Size res2 ) {
	basic::ProfileThis doit( basic::NOESY_ASSIGN_NETWORK_FIND_RAW_ASSIGN );
	if ( residues_.size() < res1 || res1 <= 0 ) {
		return BOGUS_ASSIGNMENTS;
	}
	auto res2_entry( residues_[ res1 ].find( res2 ) );
	if ( res2_entry ==  residues_[ res1 ].end() ) {
		return BOGUS_ASSIGNMENTS;
	}
	return res2_entry->second;
}

bool PeakAssignmentResidueMap::has( core::Size res1, core::Size res2 ) {
	return _assignments( res1, res2 ).size();
}

void PeakAssignmentResidueMap::add( CrossPeakList const& cpl ) {
#ifndef WIN32
	for ( auto const & it : cpl.peaks() ) {
		for ( const auto & ait : it->assignments() ) {
			add( ait );
		}
	}
#endif
}

/*void PeakAssignmentResidueMap::invalidate_competitors_to_sequential_NOE( CrossPeakList& ) {
basic::ProfileThis doit( basic::NOESY_ASSIGN_NETWORK_INVALIDATE_SEQ_NOE );
for ( Size i=1; i<total_residue(); i++ ) {
tr.Info << "focus on sequential " << i << " " << i+1 << std::endl;
PeakAssignments sequential;
if ( has( i, i+1 )) {
PeakAssignments& sequential_1( assignments( i, i+1 ) );
copy( sequential_1.begin(), sequential_1.end(), back_inserter( sequential ) );
}
if ( has( i+1, i )) {
PeakAssignments& sequential_2( assignments( i+1, i ) );
copy( sequential_2.begin(), sequential_2.end(), back_inserter( sequential ) );
}
if ( has( i, i ) ) {
PeakAssignments& sequential_3( assignments( i, i ) );
copy( sequential_3.begin(), sequential_3.end(), back_inserter( sequential ) );
}
if ( sequential.size() > 1 ) tr.Info << "multiple sequential assignments found" << std::endl;
for ( PeakAssignments::iterator it = sequential.begin(); it!=sequential.end(); ++it ) {
CrossPeak& cp( (*it)->crosspeak() );
tr.Debug << "keep "
<< (*it)->atom(  1 ) << "--"
<< (*it)->atom(  2 ) << std::endl;
for ( CrossPeak::PeakAssignments::const_iterator cit = cp.assignments().begin(); cit!=cp.assignments().end(); ++cit ) {
if ( (**cit) == (**it) ) continue;
if ( tr.Debug.visible() ) {
tr.Debug << "remove "
<< (*cit)->atom(  1 ) << "--"
<< (*cit)->atom(  2 ) << std::endl;
}
//cp.invalidate_assignment( ai );
}
}
}
}
*/
void PeakAssignmentResidueMap::check_for_symmetric_peaks( CrossPeakList& cpl, bool accumulate_symmetry ) {
#ifndef WIN32
	for ( auto const & it : cpl.peaks() ) {
		CrossPeak::PeakAssignments const& assignments( it->assignments() );
		for ( auto const & assignment : assignments ) {
			PeakAssignment& current( *assignment );
			// tr.Debug << "sym-check: current " << current << " "
			//        << current.crosspeak().peak_id() << " " << current.crosspeak().filename() << std::endl;
			//check if a symmetric partner for current exists
			Size const res1( current.resid(  1 ) );
			Size const res2( current.resid(  2 ) );
			Real found( 0 );
			PeakAssignments symmetric_candidates_inv( _assignments( res1, res2 ) );
			if ( symmetric_candidates_inv.size() ) { //BOGUS_ASSIGNMENTS would be empty
				for ( auto & sym_it : symmetric_candidates_inv ) {
					//  tr.Debug << "sym-check: other " << **sym_it;
					if ( sym_it->is_symmetric_partner_of( current ) ) {
						//       tr.Debug << " MATCH 1 " << (*sym_it)->crosspeak().peak_id() << " " << (*sym_it)->crosspeak().filename() << std::endl;
						if ( !accumulate_symmetry ) {
							//        tr.Debug << std::endl;
							found = 1;
							break;
						}
						found += sym_it->normalized_peak_volume();
					}
					//       tr.Debug << std::endl;
				}
			}
			PeakAssignments symmetric_candidates( _assignments( res2, res1 ) );
			if ( symmetric_candidates.size() ) { //BOGUS_ASSIGNMENTS would be empty
				for ( auto & symmetric_candidate : symmetric_candidates ) {
					//      tr.Debug << "sym-check: other " << **sym_it;
					if ( symmetric_candidate->is_symmetric_partner_of( current ) ) {
						//       tr.Debug << " MATCH 2 " << (*sym_it)->crosspeak().peak_id() << " " << (*sym_it)->crosspeak().filename() << std::endl;
						if ( !accumulate_symmetry ) {
							found = 1;
							//          tr.Debug << std::endl;
							break;
						}
						found += symmetric_candidate->normalized_peak_volume();
					}
					//     tr.Debug << std::endl;
				}
			}
			current.set_symmetry( found );
		}
	}
#endif
}

typedef PeakAssignmentResidueMap::PeakAssignments PeakAssignments;
void retrieve_assignment( PeakAssignments const& list, Size resonance_id1, Size resonance_id2, PeakAssignments& intra_res_NOEs ) {
#ifndef WIN32
	basic::ProfileThis doit( basic::NOESY_ASSIGN_NETWORK_RETRIEVE_ASSIGN );
	for ( auto const & it : list ) {
		if ( ( it->resonance_id( 1 ) == resonance_id1 && it->resonance_id( 2 ) == resonance_id2 )
				&& ( it->resonance_id( 1 ) == resonance_id2 && it->resonance_id( 2 ) == resonance_id1 ) ) {
			intra_res_NOEs.push_back( it );
		}
	}
#endif
}

Real sum_IntraNOE( PeakAssignments const& list, id::NamedAtomID const& atom1, id::NamedAtomID const& atom2 ) {
	Real intra_sumV( 0.0 );
#ifndef WIN32
	for ( auto const & it : list ) {
		if ( ( it->atom( 1 ) == atom1 && it->atom( 2 ) == atom2 )
				&& ( it->atom( 1 ) == atom2 && it->atom( 2 ) == atom1 ) ) {
			intra_sumV+=it->normalized_peak_volume();
		}
	}
#endif
	return intra_sumV;
}


void PeakAssignmentResidueMap::assignments( core::Size res1, core::Size res2, PeakAssignments& collector ) const {
	basic::ProfileThis doit( basic::NOESY_ASSIGN_NETWORK_FIND_RAW_ASSIGN );
	if ( res1 <= 0 || res2 <= 0 ) return;
	if ( residues_.size() < res1 ) {
		return;
	}
	auto res2_entry( residues_[ res1 ].find( res2 ) );
	if ( res2_entry ==  residues_[ res1 ].end() ) {
		return;
	}
	copy( res2_entry->second.begin(), res2_entry->second.end(), back_inserter( collector ) );

}

void PeakAssignmentResidueMap::fill_covalent_gammas( Size alpha_resid, std::map< core::id::NamedAtomID, bool >& collector ) const {
	basic::ProfileThis doit( basic::NOESY_ASSIGN_NETWORK_FILL_COV_GAMMA );
	//collector.get_allocator().allocate( 50 ); //no aminoacid with more than 16 atoms...
	for ( Size resid = (alpha_resid > 1 ? alpha_resid - 1 : alpha_resid);
			resid <= ( alpha_resid < atoms_.size() ? alpha_resid + 1 : atoms_.size() ); ++resid ) {
		//  tr.Trace << "fill_covalent_gammas: probe residue: " << resid << std::endl;
		for ( auto const & it : atoms_[ resid ] ) {
			//   tr.Trace << "fill_covalent_gammas: check atom: " << *it << std::endl;
			//to get integration test back: take only resonance ides within -30+30 of alpha
			//   Size new_reso_id( resonances[ *it ].label() );
			//if ( new_reso_id<= alpha-30 ) continue;
			//   if ( new_reso_id>= alpha+30 ) continue;
			collector.insert( std::make_pair( it, false ) );
		}
	}
	//assume that residues are sequential in resonance list, as they usually are...
	//  core::Size const start( alpha >= ( 30 + resonances().start_key() ) ? alpha-30 : resonances().start_key() );
	//  core::Size const end( alpha + 30 < resonances().last_key() ? alpha+30 : resonances().last_key() );
	//  for ( Size ind = start; ind <= end; ind++ ) {
	//   try {
	//    Size const new_resid( resonances()[ ind ].resid() );
	//    if ( new_resid == alpha_resid || new_resid == alpha_resid - 1 || new_resid == alpha_resid + 1 ) {
	//     collector[ ind ] = false;
	//    }
	//   } catch (EXCN_UnknownResonance&  ) {
	//    continue;
	//   }
	//  }
}

Real PeakAssignmentResidueMap::compute_Nk(
	PeakAssignment const& alpha_beta,
	id::NamedAtomID const& gamma_atom,
	// Size gamma_ind, //is 1 or 2  relation to previous code: gamma_old = ag_or_bg.resonance_id( gamma_new );
	bool connect_in_i,
	bool connect_in_j,
	bool sequential,
	PeakAssignments const& close_to_i_assignments,
	PeakAssignments const& close_to_j_assignments,
	Real longrange_peak_volume
) const {
	basic::ProfileThis doit( basic::NOESY_ASSIGN_NETWORK_COMPUTE_NK );
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	Real const vmax( params.vmax_ );
	Real const vmin( params.vmin_ );

	//  Size const alpha( alpha_beta.resonance_id( 1 ) );
	//  Size const beta( alpha_beta.resonance_id( 2 ) );
	// Size const gamma( ag_or_bg.resonance_id( gamma_ind ) );
	//find out if we have an initial assignment to gamma
	///WHAT IF I FIND MULTIPLE alpha - gamma assignments

	// PeakAssignments intra_residue_NOEs_ag;
	// intra_residue_NOEs_ag.reserve( 1000 );
	Real intra_sumV_ag( 0.0 );
	if ( connect_in_i ) intra_sumV_ag = sum_IntraNOE( close_to_i_assignments, alpha_beta.atom( 1 ), gamma_atom ); //find all assignments with alpha and gamma (both directions)
	//  for ( PeakAssignments::const_iterator it = intra_residue_NOEs_ag.begin(); it != intra_residue_NOEs_ag.end(); ++it ) {
	//   intra_sumV_ag+=(*it)->normalized_peak_volume();
	//  }

	PeakAssignments intra_residue_NOEs_bg;
	// intra_residue_NOEs_bg.reserve( 1000 );
	Real intra_sumV_bg( 0.0 );
	if ( connect_in_j ) intra_sumV_bg = sum_IntraNOE( close_to_j_assignments, alpha_beta.atom( 2 ), gamma_atom ); //find all assignments with alpha and gamma
	//  for ( PeakAssignments::const_iterator it = intra_residue_NOEs_bg.begin(); it != intra_residue_NOEs_bg.end(); ++it ) {
	//   intra_sumV_bg+=(*it)->normalized_peak_volume();
	//  }

	//( ag_or_bg.atom( gamma_ind ) );
	/// WARNING WARNING WARNING THREAD UNSAFE covalent_compliance RELIES ON THREAD-UNSAFE SINGLETON CovalentCompliance
	bool const is_covalent_ag( covalent_compliance( alpha_beta.atom( 1 ), gamma_atom ) );
	bool const is_covalent_bg( covalent_compliance( gamma_atom, alpha_beta.atom( 2 ) ) );
	Real const Vcc_ag( is_covalent_ag ? vmax : ( ( connect_in_i || sequential ) ? vmin : 0 ) );
	Real const Vcc_bg( is_covalent_bg ? vmax : ( ( connect_in_j || sequential ) ? vmin : 0 ) );
	Real const vag( std::max( Vcc_ag, connect_in_i ? intra_sumV_ag : longrange_peak_volume ) );
	Real const vbg( std::max( Vcc_bg, connect_in_j ? intra_sumV_bg : longrange_peak_volume ) );
	//  tr.Trace << alpha_beta.atom( 1 ) << "-" << gamma_atom << "-" << alpha_beta.atom( 2 )
	//       << ": " << Vcc_ag << " " << Vcc_bg << " " << vag << " " << vbg << " / (vmin: " << vmin << ") "
	//       << intra_sumV_ag << " " << intra_sumV_bg << " " << longrange_peak_volume << " ==> " << sqrt( (vag<vmin || vbg < vmin )? 0 : vag*vbg ) << std::endl;

	if ( vag < vmin || vbg < vmin ) return 0.0;//Heavyside function //verbatim would mean that in case of vag == vmin , 0.5*vmin is returned... but that seems stupid.
	return sqrt( vag*vbg );
}


Real sum_peak_volumes( PeakAssignments const& list, core::Size i, core::Size j ) {
	Real intra_sumV( 0.0 );
	for ( auto const & it : list ) {
		if ( it->resonance_id( 1 ) == i && it->resonance_id( 2 ) == j ) {
			intra_sumV+=it->normalized_peak_volume();
		}
	}
	return intra_sumV;
}

void PeakAssignmentResidueMap::network_analysis2() { // ResonanceList const& resonances ) {
	tr.Info << "start network analysis (type2)..." << std::endl;
	typedef std::pair< core::Size, core::Size > ResonancePair;
	typedef std::map< ResonancePair, core::Real > AnchorMap;
	// typedef std::map< ResonancePair, bool > LazyMap;
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	AnchorMap anchor_weights; //Farray would also be fine, but we need to map resonance keys to 1...N instead of using keys from input file
	ObjexxFCL::FArray2D_double residue_sum( atoms_.size(), atoms_.size(), 0.0 );
	// for ( ResonanceList::const_iterator rit=resonances.begin(); rit!=resonance.end(); ++rit ) {
	//  for ( ResonanceList::const_iterator rit=resonances.begin(); rit!=resonance.end(); ++rit ) {
	//   anchor_weights
	//}
	for ( ResidueList::const_iterator it = residues_.begin(); it != residues_.end(); ++it ) {
		for ( auto const & mit : *it ) {
			//list of assignments between a residue i (given by it) and another residue (given by mit)
			PeakAssignments const& assignments_ij( mit.second );

			//determine which residues we are talking about:
			PeakAssignment const& first_assignment( **assignments_ij.begin() );
			Size const resi( first_assignment.resid(  1 ) );
			Size const resj( first_assignment.resid(  2 ) );

			//get the list of inverse assignments
			PeakAssignments const& assignments_ji( _assignments( resj, resi ) );

			//compute anchor_weight \nu(i,j)
			bool sequential( std::abs( (int) resi- (int) resj ) <= 1 );
			for ( auto ait = assignments_ij.begin(); ait != assignments_ij.end(); ++ait ) {
				//ait - is now the alpha, beta assignment whose Nk we want to get --- referring to guntert paper.
				//    PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
				core::Size const i( (*ait)->resonance_id( 1 ) );
				core::Size const j( (*ait)->resonance_id( 2 ) );
				bool covalent( false );
				if ( sequential ) {
					/// WARNING WARNING WARNING THREAD UNSAFE covalent_compliance RELIES ON THREAD-UNSAFE SINGLETON CovalentCompliance
					covalent = covalent_compliance( (*ait)->atom(1), (*ait)->atom(2) );
				}
				ResonancePair ij( i,j );
				ResonancePair ji( j,i );
				if ( anchor_weights.find(ij) == anchor_weights.end() ) {
					core::Real sum_pv( sum_peak_volumes( assignments_ij, i, j ) );
					sum_pv += sum_peak_volumes( assignments_ji, j, i );
					core::Real cov_contrib = covalent ? params.vmax_ : ( sequential ? params.vmin_ : 0 );
					core::Real n_ij = std::max( sum_pv, cov_contrib );
					if ( n_ij < params.vmin_ ) {
						n_ij = 0;
					}
					anchor_weights[ij]= n_ij;
					anchor_weights[ji]= n_ij;
				} //not in AnchorMap yet.
			} //all assignments i,j between resi, resj
		} //iterate over all residue resj
	}//iterate over all residues resi


	//now reap the benefits...
	for ( ResidueList::const_iterator it = residues_.begin(); it != residues_.end(); ++it ) {
		for ( auto const & mit : *it ) {
			//list of assignments between a residue i (given by it) and another residue (given by mit)
			PeakAssignments const& assignments_ij( mit.second );

			//determine which residues we are talking about:
			PeakAssignment const& first_assignment( **assignments_ij.begin() );
			Size const resi( first_assignment.resid(  1 ) );
			Size const resj( first_assignment.resid(  2 ) );

			//are resi and resj sequential
			//bool sequential( std::abs( (int) resi- (int) resj ) <= 1 );

			//get all resonances K that could potentiall give us a triangle i, j, k
			ResonanceList::Resonances resK;
			std::set< core::Size > neighbor_residues;
			for ( Size i = ( resi > 1 ? resi -1 : resi ); i<= ( resi < atoms_.size() ? resi + 1 : resi ); ++i ) {
				neighbor_residues.insert( i );
			}
			for ( Size i = ( resj > 1 ? resj -1 : resj ); i<= ( resj < atoms_.size() ? resj + 1 : resj ); ++i ) {
				neighbor_residues.insert( i );
			}
			for ( core::Size neighbor_residue : neighbor_residues ) {
				try {
					ResonanceList::Resonances const& retrieved( first_assignment.resonances().resonances_at_residue( neighbor_residue ) );
					//     tr.Debug << resi << " " << resj << " resonance k @ " << *sit << std::endl;
					copy( retrieved.begin(), retrieved.end(), back_inserter( resK ) );
				} catch (EXCN_UnknownResonance const & excn ) {
				}
			}
			//   for ( ResonanceList::Resonances::const_iterator itK = resK.begin(); itK != resK.end(); ++itK ) {
			//    if ( (*itK)->is_proton() ) tr.Debug << (*itK)->atom() << std::endl;
			//   }

			core::Real sumNK_resij( 0.0 );
			//now compute network anchoring for each assignment between residue-i and residue-j
			for ( auto const & ait : assignments_ij ) {
				core::Size const i( ait->resonance_id( 1 ) );
				core::Size const j( ait->resonance_id( 2 ) );

				core::Real sumNK( 0 );
				for ( ResonanceList::Resonances::const_iterator itK = resK.begin(); itK != resK.end(); ++itK ) {
					if ( !(*itK)->is_proton() ) continue;
					if ( (*itK)->label() == i || (*itK)->label() == j ) continue;
					core::Size const k( (*itK)->label() );
					ResonancePair ik( i,k );
					ResonancePair kj( k,j );
					core::Real wik, wkj;

					wik=0.0;
					AnchorMap::const_iterator nik( anchor_weights.find(ik) );
					if ( nik == anchor_weights.end() ) { //not computed yet --- definitely no assignments for ik, just add covalent crap
						bool sequential = std::abs( (int) resi - (int) (*itK)->resid() ) <= 1;
						bool covalent( false );
						if ( sequential ) {
							/// WARNING WARNING WARNING THREAD UNSAFE covalent_compliance RELIES ON THREAD-UNSAFE SINGLETON CovalentCompliance
							covalent = covalent_compliance( ait->atom(1), (*itK)->atom() );
						}
						wik = covalent ? params.vmax_ : ( sequential ? params.vmin_ : 0 );
						anchor_weights[ik] = wik;
						anchor_weights[ResonancePair(k,i)] = wik;
					} else {
						wik=std::max( wik, nik->second );
					}

					wkj=0.0;
					AnchorMap::const_iterator nkj( anchor_weights.find(kj) );
					if ( nkj == anchor_weights.end() ) { //not computed yet --- definitely no assignments for kj, just add covalent crap
						bool sequential = std::abs( (int) resj - (int) (*itK)->resid() ) <= 1;
						bool covalent( false );
						if ( sequential ) {
							/// WARNING WARNING WARNING THREAD UNSAFE covalent_compliance RELIES ON THREAD-UNSAFE SINGLETON CovalentCompliance
							covalent = covalent_compliance( ait->atom(2), (*itK)->atom() );
						}
						wkj = covalent ? params.vmax_ : ( sequential ? params.vmin_ : 0 );
						anchor_weights[kj] = wkj;
						anchor_weights[ResonancePair(j,k)] = wkj;
					} else {
						wkj = std::max( wkj, nkj->second );
					}

					sumNK += sqrt( wik*wkj );

				}
				ait->set_network_anchoring( sumNK );
				//    tr.Debug << sumNK << " " << **ait << std::endl;
				sumNK_resij+=sumNK;
			} // for ait
			residue_sum( resi, resj )+=sumNK_resij;
			residue_sum( resj, resi )+=sumNK_resij;
		} //resj
	} //resi
	for ( ResidueList::const_iterator it = residues_.begin(); it != residues_.end(); ++it ) {
		for ( auto const & mit : *it ) {
			PeakAssignments const& assignments_ij( mit.second );
			//now cycle through all assignments between resi and resj
			PeakAssignment const& first_assignment( **assignments_ij.begin() );
			Size const resi( first_assignment.resid(  1 ) );
			Size const resj( first_assignment.resid(  2 ) );
			for ( auto const & ait : assignments_ij ) {
				ait->set_network_anchoring_per_residue( residue_sum( resi, resj ) );
			}
		}
	}
}

void PeakAssignmentResidueMap::network_analysis( Size nr_assignments ) {
	tr.Info << "network_analysis..." << std::endl;
	PROF_START( NOESY_ASSIGN_NETWORK );


	utility::vector1< core::Real > Nk_buf( nr_assignments );
	utility::vector1< core::Real > reswise_Nk_buf( nr_assignments, 0.0 ); //this is longer than required...

	Size ct_buf( 1 );
	Size ct_respair_buf( 0 );
	// Real const vmax( 1.0 );
	// Real const vmin( 0.1 );
	//cycle through assigned pairs of residues
	//compute Nk but store it temporarily in Nk_buf, since we don't want to change peak volumina mid-computation.
	//store new Nks in assignments at end

	// this is more convenient because we need to retrieve i-i, i-j, and j-j information and can share this between
	// all assignments that have the same pair i-j
	for ( ResidueList::const_iterator it = residues_.begin(); it != residues_.end(); ++it ) {
		for ( auto const & mit : *it ) {
			PeakAssignments const& assignments_ij( mit.second );
			runtime_assert( assignments_ij.size() > 0 ); // should not be in list if empty

			PeakAssignment const& first_assignment( **assignments_ij.begin() );
			Size const resi( first_assignment.resid(  1 ) );
			Size const resj( first_assignment.resid(  2 ) );

			//   tr.Trace << "network ana of " << resi << " " << resj << std::endl;
			++ct_respair_buf;
			//   Real residue_wise_Nk( 0.0 );

			PeakAssignments assignments_around_ij;
			assignments( resi-1, resj-1, assignments_around_ij );//why these ?
			assignments( resi-1, resj, assignments_around_ij );
			assignments( resi-1, resj+1, assignments_around_ij );//why these ?
			assignments( resi  , resj-1, assignments_around_ij );
			assignments( resi  , resj, assignments_around_ij );
			assignments( resi  , resj+1, assignments_around_ij );
			assignments( resi+1, resj-1, assignments_around_ij );//why these ?
			assignments( resi+1, resj, assignments_around_ij );
			assignments( resi+1, resj+1, assignments_around_ij );//why these ?

			//makes results worse! rescoring tests on 4 targets (sr10, sgr145, rhodopsin, mbp) -- taken out on Nov 8th
			//    PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
			//    if ( params.network_include_reverse_dir_ ) {
			//     //20th Oct 2010: why not the opposite direction?
			//     assignments( resj-1, resi-1, assignments_around_ij );//why these ?
			//     assignments( resj-1, resi, assignments_around_ij );
			//     assignments( resj-1, resi+1, assignments_around_ij );//why these ?
			//     assignments( resj  , resi-1, assignments_around_ij );
			//     assignments( resj  , resi, assignments_around_ij );
			//     assignments( resj  , resi+1, assignments_around_ij );
			//     assignments( resj+1, resi-1, assignments_around_ij );//why these ?
			//     assignments( resj+1, resi, assignments_around_ij );
			//     assignments( resj+1, resi+1, assignments_around_ij );//why these ?

			//    }

			//get i-i and j-j assignments - use the method _X that returns BOGUS_ASSIGNMENTS and doesn't throw exception
			//    PeakAssignments const& intra_i_assignments( _assignments( resi, resi ) );
			//    PeakAssignments const& intra_j_assignments( _assignments( resj, resj ) );

			//we want to look for gamma atoms that help us to connect i -- > j
			//these can come from same or neighbouring residue, hence we require assignments i -> i-1/i/i+1 or j-1/j/j+1 -> j
			//to start in i or end in j.

			//   std::cout << "START NEW RESIDUE PAIR " << resi << " " << resj << std::endl;

			/*   for ( PeakAssignments::const_iterator yit = assignments_ij.begin(); yit != assignments_ij.end(); ++yit ) {
			std::cout << **yit << std::endl;
			}
			std::cout << "--------------------------------------------------" << std::endl;
			for ( PeakAssignments::const_iterator yit = assignments_around_ij.begin(); yit != assignments_around_ij.end(); ++yit ) {
			std::cout << **yit << std::endl;
			}
			*/
			PeakAssignments close_to_i_assignments;
			assignments( resi, resi-1, close_to_i_assignments );
			assignments( resi, resi, close_to_i_assignments );
			assignments( resi, resi+1, close_to_i_assignments );

			PeakAssignments close_to_j_assignments;
			assignments( resj-1, resj, close_to_j_assignments );
			assignments( resj, resj, close_to_j_assignments );
			assignments( resj+1, resj, close_to_j_assignments );

			//now cycle through all assignments between resi and resj
			for ( auto const & ait : assignments_ij ) {
				//ait - is now the alpha, beta assignment whose Nk we want to get --- referring to guntert paper.

				//    PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );

				bool is4D( ait->crosspeak().is4D() );
				if ( is4D ) {
					Real sum_Nk = 1.0;
					Nk_buf[ ct_buf++ ]=sum_Nk;
					reswise_Nk_buf[ ct_respair_buf ] += sum_Nk;
					continue;
				}

				//Size const alpha( (*ait)->resonance_id( 1 ) );
				//Size const beta( (*ait)->resonance_id( 2 ) );
				Size const alpha_resid( ait->resid( 1 ) );
				Size const beta_resid( ait->resid( 2 ) );
				core::id::NamedAtomID const& alpha_atom( ait->atom( 1 ) );
				core::id::NamedAtomID const& beta_atom( ait->atom( 2 ) );
				//     if ( tr.Trace.visible() ) {
				//      PeakAssignment const& assignment( **ait );
				//      Size const resi( assignment.resid(  1 ) );
				//      Size const resj( assignment.resid(  2 ) );
				//      tr.Trace << alpha << " " << resi << " " << beta << " " << resj << std::endl;
				//      tr.Trace << "human readable: " << resonances()[ alpha ].atom() << " -- " << resonances()[ beta ].atom() << std::endl;
				//     }

				core::Size const seq_dist( resi < resj ? resj - resi : resi - resj );
				std::map< id::NamedAtomID, bool > covalent_gammas; //visited this gamma atom already ?
				if ( seq_dist <= 2 ) {
					fill_covalent_gammas( alpha_resid, covalent_gammas );
					fill_covalent_gammas( beta_resid, covalent_gammas );
				}
				//make list of covalent structure gamma atoms

				Real sum_Nk( 0.0 );
				///cycle through all "gamma" atoms

				/// first count how often each peak is in the gamma list... multiple ambiguous assignments of a single peak should not just add up to provide a strong network
				///   thus we count how often each peak appears and then normalize by that number.
				std::map< std::string, core::Size > offset_map; //to provide a quick peak_id multiplicator.
				std::map< core::Size, core::Size > peak_count_map;
				core::Size next_offset( 100000000 ); //expect to have peak_id smaller than that value.
				core::Size last_offset( next_offset );
				{  basic::ProfileThis doit( basic::NOESY_ASSIGN_NETWORK_PEAK_COUNT );
					offset_map[ ait->crosspeak().filename() ] = 0;
					peak_count_map[ 0+ait->crosspeak().peak_id() ] = 1;
					for ( PeakAssignments::const_iterator yit = assignments_around_ij.begin(); yit != assignments_around_ij.end(); ++yit ) {
						auto file_it = offset_map.find( (*yit)->crosspeak().filename() );
						core::Size offset;
						if ( file_it != offset_map.end() ) {
							offset = file_it->second;
						} else { //new filename
							offset_map[ (*yit)->crosspeak().filename() ] = last_offset;
							offset = last_offset;
							last_offset += next_offset;
						}
						peak_count_map[ offset+(*yit)->crosspeak().peak_id() ] += 1;
					}
				} //scope

				for ( PeakAssignments::const_iterator yit = assignments_around_ij.begin(); yit != assignments_around_ij.end(); ++yit ) {
					//ignore assignments that are from the same peak... how can these possibly make the case stronger:
					//offset of current crosspeak (*ait)->crosspeak() is 0 by construction
					if ( 0+ait->crosspeak().peak_id() == (offset_map[ (*yit)->crosspeak().filename() ]+(*yit)->crosspeak().peak_id()) ) {
						//       if ( tr.Trace.visible() ) {      //        tr.Trace << "ignore assignment of same peak: " << resonances()[ (*yit)->resonance_id( 1 ) ].atom() << " " << resonances()[ (*yit)->resonance_id( 2 ) ].atom() << std::endl;
						//       }
						continue;
					}
					//the i-j assignment must be either gamma-alpha or gamma-beta, find out which
					Size gamma_sel( 0 ); //which resonance id for atom gamma?
					bool connect_in_i( false );
					bool connect_in_j( false );
					bool sequential( false );
					Size my_res( 0 );


					//// NOTE this comparison has to be via NamedAtomID ...
					//// actaually we need somthing like a NMRAtomID to take care of QD1 ambiguities?
					////

					// if yit is gamma-beta in (i-j), gamma is residue i, ie., connect_in_i
					id::NamedAtomID const& g_atom1( (*yit)->atom( 1 ) );
					id::NamedAtomID const& g_atom2( (*yit)->atom( 2 ) );
					if ( g_atom1 != alpha_atom && g_atom2 == beta_atom ) {
						connect_in_i = true;
						connect_in_j = false;
						gamma_sel = 1;
						my_res = resj;
					}

					if ( g_atom2 != alpha_atom && g_atom1 == beta_atom ) {
						connect_in_i = true;
						connect_in_j = false;
						gamma_sel = 2;
						my_res = resj;
						//  // makes results worse : rescoring test on 4 targets --- removed Nov 8th.
						//     std::cout << "working on: " << **ait << std::endl;
						//     std::cout << "gamma: " << **yit << std::endl;
						//     std::cout << gamma_sel << " " << my_res << " " << connect_in_i << std::endl;
					}

					// if yit is alpha-gamma in (i-j), gamma is residue j, ie., !connect_in_i
					if ( g_atom1 == alpha_atom && g_atom2 != beta_atom ) {
						connect_in_i = false;
						connect_in_j = true;
						gamma_sel = 2;
						my_res = resi;
					}

					if ( g_atom2 == alpha_atom && g_atom1 != beta_atom ) {
						connect_in_i = false;
						connect_in_j = true;
						gamma_sel = 1;
						my_res = resi;
						//  // makes results worse : rescoring test on 4 targets --- removed Nov 8th.
						//     std::cout << "working on: " << **ait << std::endl;
						//     std::cout << "gamma: " << **yit << std::endl;
						//     std::cout << gamma_sel << " " << my_res << " " << connect_in_i << std::endl;
					}
					if ( !gamma_sel ) continue; //this is a gamma, delta peak or another alpha, beta peak... no use
					int const gamma_resid( (*yit)->resid( gamma_sel ) );
					sequential = ( std::abs( (int) my_res - gamma_resid) <= 1 );

					// PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
					//      if ( params.network_allow_same_residue_connect_ ) {

					//       //Extension Oct 2010
					//       // add also connections that go to same residue as alpha beta:
					//       if ( (*yit)->resonance_id( 1 )!=alpha && resonances()[ (*yit)->resonance_id( 2 ) ].resid() == resj ) {
					//        connect_in_i = true;
					//        connect_in_j = false;
					//        gamma = (*yit)->resonance_id( 1 );
					//        int const gamma_resid( (*yit)->resid(  1 ) );
					//        sequential = ( std::abs( (int) resj - gamma_resid) <= 1 );
					//       }
					//       if ( (*yit)->resonance_id( 2 )!=alpha && resonances()[ (*yit)->resonance_id( 1 ) ].resid() == resj ) {
					//        connect_in_i = true;
					//        connect_in_j = false;
					//        gamma = (*yit)->resonance_id( 2 );
					//        int const gamma_resid( (*yit)->resid(  2 ) );
					//        sequential = ( std::abs( (int) resj - gamma_resid) <= 1 );
					//       }


					//       if (  resonances()[ (*yit)->resonance_id( 1 ) ].resid() == resi && (*yit)->resonance_id( 2 )!=beta ) {
					//        connect_in_i = false;
					//        connect_in_j = true;
					//        gamma = (*yit)->resonance_id( 2 );
					//        int const gamma_resid( (*yit)->resid(  2 ) );
					//        sequential = ( std::abs( (int) resi - gamma_resid) <= 1 );
					//       }
					//       if (  resonances()[ (*yit)->resonance_id( 2 ) ].resid() == resi && (*yit)->resonance_id( 1 )!=beta ) {
					//        connect_in_i = false;
					//        connect_in_j = true;
					//        gamma = (*yit)->resonance_id( 1 );
					//        int const gamma_resid( (*yit)->resid(  1 ) );
					//        sequential = ( std::abs( (int) resi - gamma_resid) <= 1 );
					//       }
					//      } //same_residue_connect ---


					// if ( tr.Trace.visible() && !gamma ) {
					//       tr.Trace << "no gamma found for " << resonances()[ (*yit)->resonance_id( 1 ) ].atom() << " " << resonances()[ (*yit)->resonance_id( 2 ) ].atom() << std::endl;
					// }

					id::NamedAtomID const& gamma_atom( (*yit)->atom( gamma_sel ) );
					if ( seq_dist <= 2 ) covalent_gammas[ gamma_atom ] = true; //this one has been visited

					Size ambiguity_factor( peak_count_map[ (*yit)->crosspeak().peak_id() + offset_map[ (*yit)->crosspeak().filename() ]] );

					//      if ( tr.Trace.visible() ) {
					//       tr.Trace  << "add: " << (*yit)->atom( 1 ) << " "
					//            <<  (*yit)->atom( 2 )
					//             << " from peak: " << (*yit)->crosspeak().peak_id() << " " << (*yit)->crosspeak().filename() << " ambiguity: " << ambiguity_factor
					//             << std::endl;
					//       }
					sum_Nk += 1.0/ambiguity_factor* compute_Nk( *ait, gamma_atom, connect_in_i, connect_in_j,
						sequential, close_to_i_assignments, close_to_j_assignments, (*yit)->normalized_peak_volume() );
					// tr.Trace << gamma << " " << Vcc_ag << " " << vag << " " << vbg << " "
					//    << intra_sumV << " " << (*yit)->normalized_peak_volume() << std::endl;
				}
				//     tr.Trace << "add covalent gammmas" << std::endl;
				//iterate over missing covalent gammas
				for ( std::map< core::id::NamedAtomID, bool >::const_iterator cgit = covalent_gammas.begin(), cgend = covalent_gammas.end(); cgit != cgend; ++cgit ) {
					if ( cgit->second ) continue; // it has been visited
					id::NamedAtomID const& gamma_atom( cgit->first );
					int gamma_resid( gamma_atom.rsd() );
					//int res1( resi < resj ? resi : resj );
					//int res2( resi >= resj ? resi : resj );
					bool connect_in_i( true );
					bool connect_in_j( true );
					//     tr.Trace << "try gamma atom " << gamma_atom << std::endl;
					//if i, i+2 gamma needs to be exactly between two residues to be effective
					if ( seq_dist == 2 ) {
						int res1( resi < resj ? resi : resj );
						int res2( resi >= resj ? resi : resj );
						if ( !( gamma_resid-res1 == 1 || res2 - gamma_resid == 1 ) ) continue;
						//       connect_in_i = true;
						//       connect_in_j = true;
					}

					if ( seq_dist == 1 ) { //gamma needs to be either in i or in j.
						if ( gamma_resid != (int) resi && gamma_resid != (int) resj ) continue;
						//        if ( gamma_resid == resi ) connect_in_i = true;
						//        if ( gamma_resid == resj ) connect_in_j = true;
						//        if ( !( connect_in_i || connect_in_j ) ) continue;
					}
					sum_Nk += compute_Nk( *ait, gamma_atom, connect_in_i, connect_in_j, true, close_to_i_assignments, close_to_j_assignments, 0.0 /*no longrange NOE*/);
				} //covalent gammas
				//    tr.Trace << "sum: " << sum_Nk << std::endl;
				Nk_buf[ ct_buf++ ]=sum_Nk;
				reswise_Nk_buf[ ct_respair_buf ] += sum_Nk;
			}// all assignments between i-j

		} //all i
	}// all assignments

	//now fill Nks into the PeakAssignments. use same sequence all i, all j (in i-j), all in i-j
	ct_buf = 1;
	ct_respair_buf = 1;
	for ( ResidueList::const_iterator it = residues_.begin(); it != residues_.end(); ++it ) {
		for ( auto const & mit : *it ) {
			PeakAssignments const& assignments_ij( mit.second );

			//now cycle through all assignments between resi and resj
			for ( auto const & ait : assignments_ij ) {
				ait->set_network_anchoring( Nk_buf[ ct_buf++ ] );
				ait->set_network_anchoring_per_residue( reswise_Nk_buf[ ct_respair_buf ] );
			}
			++ct_respair_buf;
		}
	}

	PROF_STOP( NOESY_ASSIGN_NETWORK );

}


} //noesy_assign
} //protocols
