// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ResidueNeighborList.hh
/// @brief  A container class for use by the Etable and FA_Elec classes for storing
///         lists of atom neighbors
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/scoring/ResidueNeighborList.hh>

// Package headers
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

// Project headers
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>

/// Utility headers
#include <utility/vector0.hh>

#include <utility/vector1.hh>


/// #define APL_TEMP_DEBUG

namespace core {
namespace scoring {

/*
const Real ResiduePairNeighborList::narrow_reach = { 0.55 };
const Real ResiduePairNeighborList::wide_reach = { 2.0 };

Real const nr2 = 2 * ResiduePairNeighborList::narrow_reach;
Real const wr2 = 2 * ResiduePairNeighborList::wide_reach;

void
initialize_cut2(
Real heavy_heavy_dist_cutoff,
Real heavy_hydrogen_dist_cutoff,
Real hydrogen_hydrogen_dist_cutoff,
Real * ncut2,
Real * wcut2
)
{

Real ncut[ 3 ] = { heavy_heavy_dist_cutoff + nr2, heavy_hydrogen_dist_cutoff + nr2, hydrogen_hydrogen_dist_cutoff + nr2 };
Real wcut[ 3 ] = { heavy_heavy_dist_cutoff + wr2, heavy_hydrogen_dist_cutoff + wr2, hydrogen_hydrogen_dist_cutoff + wr2 };

ncut2[ 0 ] = ncut[ 0 ] * ncut[ 0 ];
ncut2[ 1 ] = ncut[ 1 ] * ncut[ 1 ];
ncut2[ 2 ] = ncut[ 2 ] * ncut[ 2 ];

wcut2[ 0 ] = wcut[ 0 ] * wcut[ 0 ];
wcut2[ 1 ] = wcut[ 1 ] * wcut[ 1 ];
wcut2[ 2 ] = wcut[ 2 ] * wcut[ 2 ];

}


#ifdef APL_TEMP_DEBUG
Size profile_natoms = 0;
Size profile_nwide_update = 0;
Size profile_nnarrow_update = 0;
Size profile_nwide_update_total = 0;
Size profile_nnarrow_update_total = 0;
bool profile_nodes_processed = false;
Size profile_update_calls = 0;
Size n_narrow_dist_evals = 0;
Size n_wide_dist_evals = 0;
Size n_update_dist_evals = 0;
Size n_update_dist_evals_total = 0;
Size n_empty_edges = 0;
Size n_zero_neighb_atoms = 0;
Size n_some_neighb_atoms = 0;
#endif


ResidueNblistData::ResidueNblistData() : nnarrow_ood_( 0 ), nwide_ood_( 0 ) {}
ResidueNblistData::~ResidueNblistData()
{
#ifdef APL_TEMP_DEBUG
/// APL DEBUG TEMP
if ( profile_update_calls != 0 || n_zero_neighb_atoms != 0 ) {
std::cout << "#updates = " << profile_update_calls << " #narrow "
<< profile_nnarrow_update_total << " #wide " << profile_nwide_update_total << " %narrow "
<< ( (double) profile_nnarrow_update_total ) / ( profile_natoms * profile_update_calls )
<< " %wide " <<  ( (double) profile_nwide_update_total ) / (profile_natoms * profile_update_calls )
<< " dis/update: " << (double) n_update_dist_evals_total / (profile_natoms * profile_update_calls)
<< " dis/wideood: " << (double) n_wide_dist_evals / (profile_nwide_update_total)
<< " dis/narrowood: " << (double) n_narrow_dist_evals / (profile_nnarrow_update_total) << std::endl;
std::cout << "#zero neighb atoms: " << n_zero_neighb_atoms << " #some neighb atoms: " << n_some_neighb_atoms << std::endl;
profile_natoms = 0;
profile_nwide_update = 0;
profile_nnarrow_update = 0;
profile_nwide_update_total = 0;
profile_nnarrow_update_total = 0;
profile_nodes_processed = false;
profile_update_calls = 0;
n_narrow_dist_evals = 0;
n_wide_dist_evals = 0;
n_update_dist_evals = 0;
n_update_dist_evals_total = 0;
n_zero_neighb_atoms = 0;
n_some_neighb_atoms = 0;
}
#endif
}

ResidueNblistData::CacheableDataOP ResidueNblistData::clone() const
{
return new ResidueNblistData( *this );
}

void
ResidueNblistData::setup_for_intrares_nblist(
CountPairFunctionCOP cpfxn,
Real heavy_heavy_dist_cutoff,
Real heavy_hydrogen_dist_cutoff,
Real hydrogen_hydrogen_dist_cutoff
)
{
cpfxn_ = cpfxn;
initialize_cut2( heavy_heavy_dist_cutoff, heavy_hydrogen_dist_cutoff, hydrogen_hydrogen_dist_cutoff, ncut2_, wcut2_ );
}

void ResidueNblistData::initialize( conformation::Residue const & res )
{

Size const natoms( res.natoms() );
is_H_.resize( natoms );
ood_status_.resize( natoms );
narrow_ood_list_.resize( natoms );
wide_ood_list_.resize( natoms );
narrow_coord_.resize( natoms );
wide_coord_.resize( natoms );
at_neighbors_.resize( natoms );

#ifdef APL_TEMP_DEBUG
/// APL DEBUG TEMP
profile_natoms += natoms;
#endif

Real maxd2_to_nbrat;
Vector nbrxyz = res.xyz( res.nbr_atom() );
for ( Size ii = 1; ii <= is_H_.size(); ++ii ) {
Real d2_to_nbrat = res.xyz( ii ).distance_squared( nbrxyz );
if ( ii == 1 || d2_to_nbrat > maxd2_to_nbrat ) {
maxd2_to_nbrat = d2_to_nbrat;
}
is_H_[ ii ] = res.atom_is_hydrogen( ii ) ? 1 : 0;
ood_status_[ ii ] = uptodate;
wide_coord_[ ii ] = narrow_coord_[ ii ] = res.xyz( ii );
}

if ( ! cpfxn_ ) return;

if ( narrow_.size() < natoms ) {
narrow_.resize( natoms );
wide_.resize( natoms );
}
for ( Size ii = 1; ii <= natoms; ++ii ) {
narrow_[ ii ].resize( natoms );
wide_[ ii ].resize( natoms );
}

for ( Size ii = 1; ii <= natoms; ++ii ) {
for ( Size jj = ii + 1; jj <= natoms; ++jj ) {
Real weight( 1.0 );
Size path_dist( 0 );
if ( cpfxn_->count( ii, jj, weight, path_dist ) ) {
Real d2 = res.xyz( ii ).distance_squared( res.xyz( jj ) );
if ( d2 < wcut2_[ is_H_[ ii ] + is_H_[ jj ] ] ) {

wide_[ ii ][ jj ].data() = AtomNeighbor( res.seqpos(), jj, path_dist, weight );
wide_[ jj ][ ii ].data() = AtomNeighbor( res.seqpos(), ii, path_dist, weight );
wide_[ ii ].move_to_back( jj );
wide_[ jj ].move_to_back( ii );
if ( d2 < ncut2_[ is_H_[ ii ] + is_H_[ jj ] ] ) {
narrow_[ ii ][ jj ].data() = wide_[ ii ][ jj ].data();
narrow_[ jj ][ ii ].data() = wide_[ jj ][ ii ].data();
narrow_[ ii ].move_to_back( jj );
narrow_[ jj ].move_to_back( ii );
at_neighbors_[ ii ].push_back( wide_[ ii ][ jj ].data() );
//at_neighbors_[ jj ].push_back( wide_[ jj ][ ii ].data() );

}
}
}
}
}
}


void ResidueNblistData::update( conformation::Residue const & res )
{
//return; /// TEMP!!!! Short circuit the whole thing
#ifdef APL_TEMP_DEBUG
/// APL TEMP
if ( ! profile_nodes_processed ) {
std::cout << "Last update: " << n_update_dist_evals << " distance evaluations" << std::endl;
n_update_dist_evals = 0;
++profile_update_calls;
profile_nodes_processed = true;
}
#endif
Real wide_movement_cutoff( ResiduePairNeighborList::wide_reach - ResiduePairNeighborList::narrow_reach );
Real wide_movement_cutoff2( wide_movement_cutoff * wide_movement_cutoff );

nnarrow_ood_ = 0;
nwide_ood_ = 0;
narrow_ood_list_.clear();
wide_ood_list_.clear();

for ( Size ii = 1; ii <= res.natoms(); ++ii ) {
Real d2narrow = narrow_coord_[ ii ].distance_squared( res.xyz( ii ) );
if ( d2narrow > ResiduePairNeighborList::narrow_reach ) {
Real d2wide = wide_coord_[ ii ].distance_squared( res.xyz( ii ) );
if ( d2wide > wide_movement_cutoff2 ) {
wide_ood_list_.move_to_back( ii );
++nwide_ood_;

#ifdef APL_TEMP_DEBUG
++profile_nwide_update;
++profile_nwide_update_total;
#endif

ood_status_[ ii ] = wide_ood;
wide_coord_[ ii ] = res.xyz( ii );
} else {
#ifdef APL_TEMP_DEBUG
/// APL TEMP
++profile_nnarrow_update;
++profile_nnarrow_update_total;
#endif
++nnarrow_ood_;
narrow_ood_list_.move_to_back( ii );
ood_status_[ ii ] = narrow_ood;
}
narrow_coord_[ ii ] = res.xyz( ii );
} else {
ood_status_[ ii ] = uptodate;
}
}

if ( ! cpfxn_ ) return;

for ( Size ii = 1; ii <= res.natoms(); ++ii ) {
if ( ood_status_[ ii ] == uptodate ) continue;
if ( ood_status_[ ii ] == narrow_ood ) {
for ( Size jj = narrow_[ ii ].head(), jjend = narrow_[ ii ].end(); jj != jjend; jj = narrow_[ ii ][ jj ].next() ) {
Size jjatno = narrow_[ ii ][ jj ].data().atomno();
if ( narrow_[ jjatno ][ ii ].in_list() ) { narrow_[ jjatno ].remove( ii ); }
}
narrow_[ ii ].clear();
} else {
debug_assert( ood_status_[ ii ] == wide_ood );
for ( Size jj = wide_[ ii ].head(), jjend = wide_[ ii ].end(); jj != jjend; jj = wide_[ ii ][ jj ].next() ) {
Size jjatno = wide_[ ii ][ jj ].data().atomno();
if ( wide_[ jjatno ][ ii ].in_list() ) { wide_[ jjatno ].remove( ii ); }
}
wide_[ ii ].clear();
for ( Size jj = narrow_[ ii ].head(), jjend = narrow_[ ii ].end(); jj != jjend; jj = narrow_[ ii ][ jj ].next() ) {
Size jjatno = narrow_[ ii ][ jj ].data().atomno();
if ( narrow_[ jjatno ][ ii ].in_list() ) { narrow_[ jjatno ].remove( ii ); }
}
narrow_[ ii ].clear();
}
}
for ( Size ii = 1; ii <= res.natoms(); ++ii ) {
if ( ood_status_[ ii ] == uptodate ) continue;
int ii_isH = is_H_[ ii ] ? 1 : 0;
if ( ood_status_[ ii ] == narrow_ood ) {
AtomNeighbors const & iiwide( wide_[ ii ] );
for ( Size jj = iiwide.head(), jjend = iiwide.end(); jj != jjend; jj = iiwide[ jj ].next() ) {
Size jjatno = iiwide[ jj ].data().atomno();
if ( ii > jjatno && ood_status_[ jjatno ] != uptodate ) continue; // we've already visited this atpair
if ( ood_status_[ jjatno ] == wide_ood ) continue; // wait until the wide list is up to date
debug_assert( ! narrow_[ ii ][ jjatno ].in_list() );
debug_assert( ! narrow_[ jjatno ][ ii ].in_list() );
debug_assert( wide_[ jjatno ][ ii ].in_list() );
int jj_isH = is_H_[ jjatno ];
Real d2 = narrow_coord_[ ii ].distance_squared( narrow_coord_[ jj ] );
if ( d2 < ncut2_[ ii_isH + jj_isH ] ) {
narrow_[ ii ][ jjatno ].data() = iiwide[ jjatno ].data();
narrow_[ ii ].move_to_back( jjatno );
narrow_[ jjatno ][ ii ].data() = wide_[ jjatno ][ ii ].data();
narrow_[ jjatno ].move_to_back( ii );
}
}
} else {
int ii_isH = is_H_[ ii ];
for ( Size jj = 1; jj <= res.natoms(); ++jj ) {
if ( ii > jj && ood_status_[ jj ] == wide_ood ) continue; // already visited
debug_assert( ! narrow_[ ii ][ jj ].in_list() );
debug_assert( ! narrow_[ jj ][ ii ].in_list() );
debug_assert( ! wide_[ ii ][ jj ].in_list() );
debug_assert( ! wide_[ jj ][ ii ].in_list() );
Real weight( 1.0 ); Size path_dist( 0 );
if ( cpfxn_->count( ii, jj, weight, path_dist ) ) {
Real d2w = wide_coord_[ ii ].distance_squared( wide_coord_[ jj ] );
int jj_isH = is_H_[ jj ];
if ( d2w < wcut2_[ ii_isH + jj_isH ] ) {
wide_[ ii ][ jj ].data() = AtomNeighbor( res.seqpos(), jj, path_dist, weight );
wide_[ jj ][ ii ].data() = AtomNeighbor( res.seqpos(), ii, path_dist, weight );
wide_[ ii ].move_to_back( jj );
wide_[ jj ].move_to_back( ii );
Real d2n = narrow_coord_[ ii ].distance_squared( narrow_coord_[ jj ] );
if ( d2n < ncut2_[ ii_isH + jj_isH ] ) {
narrow_[ ii ][ jj ].data() = wide_[ ii ][ jj ].data();
narrow_[ jj ][ ii ].data() = wide_[ jj ][ ii ].data();
narrow_[ ii ].move_to_back( jj );
narrow_[ jj ].move_to_back( ii );
}
}
}
}
}
}
}

//////////////////////////////////////////////////////

ResiduePairNeighborList::ResiduePairNeighborList() {}
ResiduePairNeighborList::~ResiduePairNeighborList() {}

ResiduePairNeighborList::CacheableDataOP
ResiduePairNeighborList::clone() const {
return new ResiduePairNeighborList( *this );
}

void
ResiduePairNeighborList::initialize_from_residues(
Real heavy_heavy_dist_cutoff,
Real heavy_hydrogen_dist_cutoff,
Real hydrogen_hydrogen_dist_cutoff,
conformation::Residue const & r1,
conformation::Residue const & r2,
ResidueNblistData const & r1dat,
ResidueNblistData const & r2dat,
etable::count_pair::CountPairFunctionCOP cpfxn
)
{
cpfxn_ = cpfxn; // save this for use in update();
initialize_cut2( heavy_heavy_dist_cutoff, heavy_hydrogen_dist_cutoff, hydrogen_hydrogen_dist_cutoff, ncut2_, wcut2_ );

r1_at_neighbors_.resize( r1.natoms() );
for ( Size ii = 1; ii <= r1.natoms(); ++ii ) {
r1_at_neighbors_[ ii ].reserve( r2.natoms() );
}
r2_at_neighbors_.resize( r2.natoms() );
r1_w_non_empty_narrow_.resize( r1.natoms() );

for ( Size ii = 0; ii < 2; ++ii ) {
Size ii_natoms    = ii == 0 ? r1.natoms() : r2.natoms();
Size other_natoms = ii == 1 ? r1.natoms() : r2.natoms();
if ( narrow(ii).size() < ii_natoms ) {
/// only grow the arrays
narrow(ii).resize( ii_natoms );
wide( ii ).resize( ii_natoms );
for ( Size jj = 1; jj <= ii_natoms; ++jj ) {
if ( narrow(ii)[ jj ].size() < other_natoms ) {
narrow(ii)[ jj ].resize( other_natoms );
wide(ii)[   jj ].resize( other_natoms );
} else {
narrow(ii)[ jj ].clear();
wide(ii)[   jj ].clear();
}
}
} else {
for ( Size jj = 1; jj <= ii_natoms; ++jj ) {
if ( narrow(ii)[ jj ].size() < other_natoms ) {
narrow(ii)[jj].resize( other_natoms );
wide(ii)[  jj].resize( other_natoms );
} else {
narrow(ii)[jj].clear();
wide(ii)[  jj].clear();
}
}
}
}

for ( Size ii = 1; ii <= r1.natoms(); ++ii ) {
int ii_isH = r1dat.is_H( ii );
for ( Size jj = 1; jj <= r2.natoms(); ++jj ) {
int jj_isH = r2dat.is_H( jj );
Real weight( 1.0 );
Size pathdist( 0 );
if ( cpfxn_->count( ii, jj, weight, pathdist ) ) {
/// Set the atom-neighbor data even if these atoms never get scored together;
/// this saves us numerous copy operations later
wide(0)[ ii ][ jj ].data() = AtomNeighbor( r2.seqpos(), jj, pathdist, weight );
wide(1)[ jj ][ ii ].data() = AtomNeighbor( r1.seqpos(), ii, pathdist, weight );
narrow(0)[ ii ][ jj ].data() = wide(0)[ ii ][ jj ].data();
narrow(1)[ jj ][ ii ].data() = wide(1)[ jj ][ ii ].data();

Real d2 = r1.xyz( ii ).distance_squared( r2.xyz( jj ) );
if ( d2 < wcut2_[ ii_isH + jj_isH ] ) {
wide(0)[ ii ].move_to_back( jj );
wide(1)[ jj ].move_to_back( ii );
if ( d2 < ncut2_[ ii_isH + jj_isH ] ) {
narrow(0)[ ii ].move_to_back( jj );
narrow(1)[ jj ].move_to_back( ii );
r1_at_neighbors_[ ii ].push_back( wide(0)[ ii ][ jj ].data() );
r2_at_neighbors_[ jj ].push_back( wide(1)[ jj ][ ii ].data() );
}
}
}
}
}
for ( Size ii = 1; ii <= r1.natoms(); ++ii ) {
if ( narrow(0)[ii].head() != 0 ) {
r1_w_non_empty_narrow_.move_to_back( ii );
}
}
#ifdef APL_TEMP_DEBUG
for ( Size ii = 1; ii <= r1.natoms(); ++ii ) {
if ( narrow(0)[ii].head() == 0 ) {
++n_zero_neighb_atoms;
} else {
++n_some_neighb_atoms;
}
}
#endif

//debug_assert( accurate(
// heavy_heavy_dist_cutoff,
// heavy_hydrogen_dist_cutoff,
// hydrogen_hydrogen_dist_cutoff,
// r1, r2, r1dat, r2dat ) );
}

void
ResiduePairNeighborList::update(
Real heavy_heavy_dist_cutoff,
Real heavy_hydrogen_dist_cutoff,
Real hydrogen_hydrogen_dist_cutoff,
conformation::Residue const & r1,
conformation::Residue const & r2,
ResidueNblistData const & r1dat,
ResidueNblistData const & r2dat
)
{
//return; // TEMP!!! short circuit the whole thing
#ifdef APL_TEMP_DEBUG
//std::cout << "update for pair " << r1.seqpos() << " " << r2.seqpos() << std::endl;
if ( profile_nodes_processed ) {
profile_nodes_processed = false;
std::cout << "Previous update: #atoms " << profile_natoms << " narrow " << profile_nnarrow_update << " wide " << profile_nwide_update << std::endl;
profile_nwide_update = 0;
profile_nnarrow_update = 0;
}
#endif

if ( r1dat.nnarrow_ood() == 0 && r1dat.nwide_ood() == 0 && r2dat.nnarrow_ood() == 0 && r2dat.nwide_ood() == 0 ) return;

/// 1. Remove all atom-neighbors in the narrow list for narrow_ood atoms
/// and all atom-neighbors from both the narrow and wide lists for wide_ood atoms

/// 2. Add in atom neighbors for narrow_ood and wide_ood atoms

/// 1a clear the neighbors for res1
//for ( Size ii = 1; ii <= r1.natoms(); ++ii ) {
for ( Size ii = r1dat.narrow_ood_list().head(); ii != 0; ii = r1dat.narrow_ood_list()[ ii ].next() ) {
//if ( r1dat.atom_status(ii) == uptodate ) continue;
debug_assert( r1dat.atom_status(ii) == narrow_ood );
AtomNeighbors & iinarrow( narrow(0)[ ii ] );
for ( Size jj = iinarrow.head(), jjend = 0;
jj != jjend; jj = iinarrow[ jj ].next() ) {
if ( r2dat.atom_status( jj ) == uptodate ) {
narrow(1)[ jj ].remove( ii );
}
}
iinarrow.clear(); // O(k)
}
for ( Size ii = r1dat.wide_ood_list().head(); ii != 0; ii = r1dat.wide_ood_list()[ ii ].next() ) {
debug_assert( r1dat.atom_status(ii) == wide_ood );
AtomNeighbors & iiwide( wide(0)[ ii ] );
for ( Size jj = iiwide.head(), jjend = 0;
jj != jjend; jj = iiwide[ jj ].next() ) {
if ( r2dat.atom_status( jj ) < wide_ood ) {
wide(1)[ jj ].remove( ii );
}
if ( r2dat.atom_status( jj ) == uptodate && narrow(1)[ jj ][ ii ].in_list() ) {
narrow(1)[ jj ].remove( ii );
}
}
iiwide.clear(); // O(k)
narrow(0)[ii].clear(); // O(k)
}

/// 1b Clear the out of date lists for the res2
//for ( Size ii = 1; ii <= r2.natoms(); ++ii ) {
for ( Size ii = r2dat.narrow_ood_list().head(); ii != 0; ii = r2dat.narrow_ood_list()[ ii ].next() ) {
//if ( r2dat.atom_status(ii) == uptodate ) continue;
debug_assert( r2dat.atom_status(ii) == narrow_ood );
AtomNeighbors & iinarrow( narrow(1)[ ii ] );
for ( Size jj = iinarrow.head(), jjend = 0;
jj != jjend; jj = iinarrow[ jj ].next() ) {
if ( r1dat.atom_status( jj ) == uptodate ) {
narrow(0)[ jj ].remove( ii );
}
}
iinarrow.clear(); // O(k)
}
for ( Size ii = r2dat.wide_ood_list().head(); ii != 0; ii = r2dat.wide_ood_list()[ ii ].next() ) {
debug_assert( r2dat.atom_status(ii) == wide_ood );
AtomNeighbors & iiwide( wide(1)[ ii ] );
for ( Size jj = iiwide.head(), jjend = 0;
jj != jjend; jj = iiwide[ jj ].next() ) {
if ( r1dat.atom_status( jj ) < wide_ood ) {
wide(0)[ jj ].remove( ii );
}
if ( r1dat.atom_status( jj ) == uptodate && narrow(0)[ jj ][ ii ].in_list() ) {
narrow(0)[ jj ].remove( ii );
}
}
iiwide.clear(); // O(k)
narrow(1)[ii].clear(); // O(k)
}

// 2a: r1 neighbor list updates
for ( Size ii = r1dat.narrow_ood_list().head(); ii != 0; ii = r1dat.narrow_ood_list()[ ii ].next() ) {
//if ( r1dat.atom_status(ii) == uptodate ) continue;
int ii_isH = r1dat.is_H( ii );
debug_assert( r1dat.atom_status(ii) == narrow_ood );
AtomNeighbors const & iiwide( wide(0)[ ii ] );
for ( Size jj = iiwide.head(), jjend = 0;
jj != jjend; jj = iiwide[ jj ].next() ) {
int jj_isH = r2dat.is_H( jj );
Real d2 = r1dat.narrow_coord(ii).distance_squared( r2dat.narrow_coord(jj) );

#ifdef APL_TEMP_DEBUG
/// APL TEMP
++n_narrow_dist_evals;
++n_update_dist_evals;
++n_update_dist_evals_total;
#endif

if ( d2 < ncut2_[ ii_isH + jj_isH ] ) {
debug_assert( ! narrow(0)[ ii ][ jj ].in_list() );
debug_assert( ! narrow(1)[ jj ][ ii ].in_list() );
//narrow(0)[ ii ][ jjatno ].data() = iiwide[ jjatno ].data();
narrow(0)[ ii ].move_to_back( jj );
//narrow(1)[ jjatno ][ ii ].data() = wide(1)[ jjatno ][ ii ].data();
narrow(1)[ jj ].move_to_back( ii );
}
}
}
for ( Size ii = r1dat.wide_ood_list().head(); ii != 0; ii = r1dat.wide_ood_list()[ ii ].next() ) {
//if ( r1dat.atom_status(ii) == uptodate ) continue;
int ii_isH = r1dat.is_H( ii );
debug_assert( r1dat.atom_status(ii) == wide_ood );
for ( Size jj = 1; jj <= r2.natoms(); ++jj ) {
int jj_isH = r2dat.is_H( jj );
Real weight( 1.0 );
Size path_dist( 0 );
if ( cpfxn_->count( ii, jj, weight, path_dist ) ) {
#ifdef APL_TEMP_DEBUG
/// APL TEMP
++n_wide_dist_evals;
++n_update_dist_evals;
++n_update_dist_evals_total;
#endif
Real d2w = r1dat.wide_coord( ii ).distance_squared( r2dat.wide_coord( jj ) );
if ( d2w < wcut2_[ ii_isH + jj_isH ] ) {
debug_assert( ! wide(0)[ ii ][ jj ].in_list() ); // should have already been removed
debug_assert( ! wide(1)[ jj ][ ii ].in_list() ); // should have already been removed
//wide(0)[ ii ][ jj ].data() = AtomNeighbor( r2.seqpos(), jj, path_dist, weight );
//wide(1)[ jj ][ ii ].data() = AtomNeighbor( r1.seqpos(), ii, path_dist, weight );
wide(0)[ ii ].move_to_back( jj );
wide(1)[ jj ].move_to_back( ii );
Real d2n = r1dat.narrow_coord( ii ).distance_squared( r2dat.narrow_coord( jj ) );
#ifdef APL_TEMP_DEBUG
/// APL TEMP
++n_wide_dist_evals;
++n_update_dist_evals;
++n_update_dist_evals_total;
#endif

if ( d2n < ncut2_[ ii_isH + jj_isH ] ) {
debug_assert( ! narrow(0)[ ii ][ jj ].in_list() ); // should have already been removed
debug_assert( ! narrow(1)[ jj ][ ii ].in_list() ); // should have already been removed
//narrow(0)[ ii ][ jj ].data() = wide(0)[ ii ][ jj ].data();
//narrow(1)[ jj ][ ii ].data() = wide(1)[ jj ][ ii ].data();
narrow(0)[ ii ].move_to_back( jj );
narrow(1)[ jj ].move_to_back( ii );
}
}
}
}
}

// 2b res2 neighborlist updates -- some neighbors have already been detected from part 2a above
// (i.e. when an neighbor from res1 had status wide_ood)
for ( Size ii = r2dat.narrow_ood_list().head(); ii != 0; ii = r2dat.narrow_ood_list()[ ii ].next() ) {
int ii_isH = r2dat.is_H( ii );
debug_assert ( r2dat.atom_status( ii ) == narrow_ood );
AtomNeighbors const & iiwide( wide(1)[ ii ] );
for ( Size jj = iiwide.head(), jjend = 0;
jj != jjend; jj = iiwide[ jj ].next() ) {
if ( r1dat.atom_status( jj ) != uptodate ) continue;
int jj_isH = r1dat.is_H( jj );
Real d2 = r2dat.narrow_coord(ii).distance_squared( r1dat.narrow_coord(jj) );
/// APL TEMP
#ifdef APL_TEMP_DEBUG
++n_narrow_dist_evals;
++n_update_dist_evals;
++n_update_dist_evals_total;
#endif
if ( d2 < ncut2_[ ii_isH + jj_isH ] ) {
debug_assert( ! narrow(1)[ ii ][ jj ].in_list() );
debug_assert( ! narrow(0)[ jj ][ ii ].in_list() );
//narrow(1)[ ii ][ jjatno ].data() = iiwide[ jjatno ].data();
narrow(1)[ ii ].move_to_back( jj );
//narrow(0)[ jjatno ][ ii ].data() = wide(0)[ jjatno ][ ii ].data();
narrow(0)[ jj ].move_to_back( ii );
}
}
}
for ( Size ii = r2dat.wide_ood_list().head(); ii != 0; ii = r2dat.wide_ood_list()[ ii ].next() ) {
debug_assert( r2dat.atom_status(ii) == wide_ood );
int ii_isH = r2dat.is_H( ii );
for ( Size jj = 1; jj <= r1.natoms(); ++jj ) {
if ( r1dat.atom_status(jj) == wide_ood ) continue; // already updated!
int jj_isH = r1dat.is_H( jj );
Real weight( 1.0 );
Size path_dist( 0 );
if ( cpfxn_->count( jj, ii, weight, path_dist ) ) {
Real d2w = r2dat.wide_coord( ii ).distance_squared( r1dat.wide_coord( jj ) );
#ifdef APL_TEMP_DEBUG
/// APL TEMP
++n_wide_dist_evals;
++n_update_dist_evals;
++n_update_dist_evals_total;
#endif
if ( d2w < wcut2_[ ii_isH + jj_isH ] ) {
debug_assert( ! wide(1)[ ii ][ jj ].in_list() ); // should have already been removed
debug_assert( ! wide(0)[ jj ][ ii ].in_list() ); // should have already been removed
//wide(1)[ ii ][ jj ].data() = AtomNeighbor( r1.seqpos(), jj, path_dist, weight );
//wide(0)[ jj ][ ii ].data() = AtomNeighbor( r2.seqpos(), ii, path_dist, weight );
wide(1)[ ii ].move_to_back( jj );
wide(0)[ jj ].move_to_back( ii );
if ( r1dat.atom_status( jj ) < wide_ood ) {
Real d2n = r2dat.narrow_coord( ii ).distance_squared( r1dat.narrow_coord( jj ) );
#ifdef APL_TEMP_DEBUG
/// APL TEMP
++n_wide_dist_evals;
++n_update_dist_evals;
++n_update_dist_evals_total;
#endif
if ( d2n < ncut2_[ ii_isH + jj_isH ] ) {
debug_assert( ! narrow(1)[ ii ][ jj ].in_list() ); // should have already been removed
debug_assert( ! narrow(0)[ jj ][ ii ].in_list() ); // should have already been removed
//narrow(1)[ ii ][ jj ].data() = wide(1)[ ii ][ jj ].data();
//narrow(0)[ jj ][ ii ].data() = wide(0)[ jj ][ ii ].data();
narrow(1)[ ii ].move_to_back( jj );
narrow(0)[ jj ].move_to_back( ii );
}
}
}
}
}
}

//debug_assert( accurate(
// heavy_heavy_dist_cutoff,
// heavy_hydrogen_dist_cutoff,
// hydrogen_hydrogen_dist_cutoff,
// r1, r2, r1dat, r2dat ) );
}


bool
ResiduePairNeighborList::accurate(
Real heavy_heavy_dist_cutoff,
Real heavy_hydrogen_dist_cutoff,
Real hydrogen_hydrogen_dist_cutoff,
conformation::Residue const & r1,
conformation::Residue const & r2,
ResidueNblistData const & r1dat,
ResidueNblistData const & r2dat
) const
{

for ( Size ii = 1; ii <= r1.natoms(); ++ii ) {
for ( Size jj = 1; jj <= r2.natoms(); ++jj ) {
Real weight( 1.0 );
Size path_dist( 0 );
if ( cpfxn_->count( ii, jj, weight, path_dist )) {
Real d2w = r1dat.wide_coord( ii ).distance_squared( r2dat.wide_coord( jj ) );
if ( d2w < wcut2_[ r1dat.is_H( ii ) + r2dat.is_H( jj ) ] ) {
debug_assert( wide(0)[ ii ][ jj ].in_list() );
debug_assert( wide(1)[ jj ][ ii ].in_list() );
} else {
debug_assert( ! wide(0)[ ii ][ jj ].in_list() );
debug_assert( ! wide(1)[ jj ][ ii ].in_list() );
}
Real d2n = r1dat.narrow_coord( ii ).distance_squared( r2dat.narrow_coord( jj ) );
if ( d2n < ncut2_[ r1dat.is_H( ii ) + r2dat.is_H( jj ) ] ) {
debug_assert( narrow(0)[ ii ][ jj ].in_list() );
debug_assert( narrow(1)[ jj ][ ii ].in_list() );
debug_assert( wide(0)[ ii ][ jj ].in_list() );
debug_assert( wide(1)[ jj ][ ii ].in_list() );
} else {
debug_assert( ! narrow(0)[ ii ][ jj ].in_list() );
debug_assert( ! narrow(1)[ jj ][ ii ].in_list() );
}
} else {
debug_assert( ! narrow(0)[ ii ][ jj ].in_list() );
debug_assert( ! narrow(1)[ jj ][ ii ].in_list() );
debug_assert( ! wide(0)[ ii ][ jj ].in_list() );
debug_assert( ! wide(1)[ jj ][ ii ].in_list() );
}
}
}
return true;
}


/// @details Do not deallocate the space representing the atom-neighbors of either residue.
/// This is a time-saving technique to avoid the use of new and delete where possible.
void ResiduePairNeighborList::clear()
{
for ( Size ii = 1; ii <= res1_neighbors_.size(); ++ii ) res1_neighbors_[ ii ].clear();
for ( Size ii = 1; ii <= res2_neighbors_.size(); ++ii ) res2_neighbors_[ ii ].clear();
}

/// @details  This class does not hold the residue id as part of the AtomNeighbor for either residue.
/// That information must be available to the class that's using this data through some other means
void ResiduePairNeighborList::add_atom_neighbor(
Size ind1,
Size ind2,
Size path_distance,
Real weight,
Real weight_func
)
{
debug_assert( ind1 <= res1_neighbors_.size() );
debug_assert( ind2 <= res2_neighbors_.size() );
AtomNeighbor at2( 1, ind2, path_distance, weight, weight_func );
AtomNeighbor at1( 1, ind1, path_distance, weight, weight_func );
res1_neighbors_[ ind1 ].push_back( at2 );
res2_neighbors_[ ind2 ].push_back( at1 );
}
*/

ResidueNblistData::ResidueNblistData() {}
ResidueNblistData::~ResidueNblistData() {}

ResiduePairNeighborList::CacheableDataOP
ResidueNblistData::clone() const
{
	return ResiduePairNeighborList::CacheableDataOP( new ResidueNblistData( *this ) );
}

void ResidueNblistData::initialize(
	conformation::Residue const & res,
	CountPairFunctionCOP cpfxn,
	Real vvd2,
	Real hvd2,
	Real hhd2
)
{
	atom_neighbors_.clear();

	if ( ! cpfxn ) return;

	utility::vector0< Real > cutoffs( 3 );
	cutoffs[ 0 ] = vvd2;
	cutoffs[ 1 ] = hvd2;
	cutoffs[ 2 ] = hhd2;

	for ( Size ii = 1; ii < res.natoms(); ++ii ) {
		int ii_isH = res.atom_is_hydrogen( ii );
		for ( Size jj = ii+1; jj <= res.natoms(); ++jj ) {
			int jj_isH = res.atom_is_hydrogen( jj );
			Real weight( 1.0 );
			Size path_dist( 0 );
			if ( cpfxn->count( ii, jj, weight, path_dist ) ) {
				Real d2 = res.xyz( ii ).distance_squared( res.xyz( jj ) );
				if ( d2 < cutoffs[ ii_isH + jj_isH ] ) {
					atom_neighbors_.push_back( SmallAtNb( ii, jj, weight ) );
				}
			}
		}
	}
}

ResiduePairNeighborList::ResiduePairNeighborList() {}
ResiduePairNeighborList::~ResiduePairNeighborList() {}

ResiduePairNeighborList::CacheableDataOP ResiduePairNeighborList::clone() const { return ResiduePairNeighborList::CacheableDataOP( new ResiduePairNeighborList( *this ) ); }

void
ResiduePairNeighborList::initialize_from_residues(
	Real vvd2,
	Real hvd2,
	Real hhd2,
	conformation::Residue const & r1,
	conformation::Residue const & r2,
	etable::count_pair::CountPairFunctionCOP cpfxn
)
{
	atom_neighbors_.clear();

	//std::cout << "ResiduePairNeighborList::initialize_from_residues " << r1.seqpos() << " " << r2.seqpos() << std::endl;

	utility::vector0< Real > cutoffs( 3 );
	cutoffs[ 0 ] = vvd2;
	cutoffs[ 1 ] = hvd2;
	cutoffs[ 2 ] = hhd2;

	for ( Size ii = 1; ii <= r1.natoms(); ++ii ) {
		int ii_isH = r1.atom_is_hydrogen( ii );
		for ( Size jj = 1; jj <= r2.natoms(); ++jj ) {
			int jj_isH = r2.atom_is_hydrogen( jj );
			Real weight( 1.0 );
			Size path_dist( 0 );
			if ( cpfxn->count( ii, jj, weight, path_dist ) ) {
				Real d2 = r1.xyz( ii ).distance_squared( r2.xyz( jj ) );
				///std::cout << "  atoms "  << ii << " " << jj << " " << r1.atom_name( ii ) << " " << r2.atom_name( jj )
				/// << " d2: " << d2 << " d " << std::sqrt( d2 ) <<  " cutoff2:" <<  cutoffs[ ii_isH + jj_isH ]
				/// << " cutoff " << std::sqrt( cutoffs[ ii_isH + jj_isH ] ) << std::endl;
				if ( d2 < cutoffs[ ii_isH + jj_isH ] ) {
					atom_neighbors_.push_back( SmallAtNb( ii, jj, weight ) );
				}
			}
		}
	}
}


}
}
