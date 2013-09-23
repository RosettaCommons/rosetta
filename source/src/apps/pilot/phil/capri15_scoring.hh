// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <devel/dna/protocols.hh>
#include <devel/dna/util.hh>//_public.hh>
//#include <devel/dna/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/frags/TorsionFragment.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/rigid_body_moves.hh>

#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/dna/DNA_BasePotential.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/etable/Etable.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/AtomVDW.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
//#include <core/scoring/elec/FA_ElecEnergy.hh>
//#include <core/scoring/methods/WholeStructureEnergy.hh>
//#include <core/scoring/etable/EtableEnergy.hh>
//#include <core/scoring/etable/count_pair/CountPairAll.hh>
//#include <core/scoring/etable/count_pair/CountPairFunction.hh>
//#include <core/scoring/etable/count_pair/CountPairFactory.hh>
//#include <core/scoring/etable/count_pair/CountPair1BC4.hh>

#include <core/types.hh>

//#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/VariantType.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/AA.hh>

#include <core/conformation/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/visualize.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/id/AtomID_Map.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/util.hh>
//#include <basic/options/after_opts.hh>

#include <basic/prof.hh> // profiling
#include <basic/basic.hh>
#include <core/id/SequenceMapping.hh>

#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <set>
#include <cstdlib>
#include <sstream>
#include <math.h>

//silly using/typedef


#include <basic/Tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;

//static numeric::random::RandomGenerator RG(93530); // <- Magic number, do not change it!!!

using namespace core;
using namespace protocols;

using utility::vector1;
using std::string;
using std::cout;
using std::endl;
using io::pdb::dump_pdb;


//static basic::Tracer TR( "apps.pilot.phil.loop_model"  );

////////////////////////////////////////////////
// danger USING ////////////////////////////////
using namespace core;
using namespace protocols;
using namespace pose;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace basic::options;
using namespace id;
namespace OK = OptionKeys;
using utility::vector1;
using std::string;


/// @details  Helper for the constraints setup below
id::AtomID
pid( std::string const & name, Size const seqpos, pose::Pose const & pose )
{
	return id::AtomID( pose.residue(seqpos).atom_index( name ), seqpos );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CapriData
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum CapriScoreType {
	BLUE = 1,
	CYAN,
	PINK,
	RED,
	n_experimental_constraint_scores = RED
};




class CapriData : public basic::datacache::CacheableData {

public:

	basic::datacache::CacheableDataOP
	clone() const
	{
		return new CapriData( *this );
	}

public:

	/// bonus for contacting conserved residues
	typedef utility::vector1< std::pair< Size, char > > ConservedResidues;
	ConservedResidues conserved_residues;

};

typedef utility::pointer::owning_ptr< CapriData > CapriDataOP;

///
inline
CapriData const &
retrieve_capri_data_from_pose( pose::Pose const & pose )
{
	assert( pose.data().has( basic::CAPRI_DATA ) );
	assert( dynamic_cast< CapriData const *>( &( pose.data().get( basic::CAPRI_DATA ))));
	return ( static_cast< CapriData const &>(    pose.data().get( basic::CAPRI_DATA )));
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CapriTwoBodyEnergy
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/**

	 centroid score for docking / loop building:

	 "cyan" experimental data
	 reward RNA contacts to conserved residues in protein

	 for the latter, seems like we need an alignment to t033_.fasta ? In case we've done some trimming??

	 K/R     centroid contacts to phosphate backbone OP2/OP1
	 S/T/N/Q centroid contacts to phosphate backbone OP2/OP1
	 D/E/H   centroid contacts to O2'

	 vdw/hybrid vdw

	 backbone O to O2'
	 backbone N to OP2, OP1

	 distance between SAM CE and rGU N1

**/
///
///

class CapriTwoBodyEnergy : public scoring::methods::ContextIndependentTwoBodyEnergy {
public:

	/// @brief  C-tor
	CapriTwoBodyEnergy()
	{
		add_score_type( capri_cen  );
		add_score_type( capri_bb   );
	};


	/// clone
	virtual
	scoring::methods::EnergyMethodOP
	clone() const
	{
		return new CapriTwoBodyEnergy();
	}

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;


	virtual
	void
	eval_atom_derivative(
											 id::AtomID const &,
											 pose::Pose const &,
											 kinematics::DomainMap const &,
											 ScoreFunction const &,
											 EnergyMap const &,
											 Vector &,
											 Vector &
											 ) const {}

	virtual
	void
	eval_intrares_energy(
											 conformation::Residue const &,
											 pose::Pose const &,
											 ScoreFunction const &,
											 EnergyMap &
											 ) const {}


	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	virtual
	Distance
	atomic_interaction_cutoff() const { return 5.0; } // guess

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const { return; }
};


/// @note  This is a short-ranged energy, only scores neighboring residues, only scores protein-rna intxns !!

void
CapriTwoBodyEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &, // pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	// total guesses
	Real const min_centroid_dis2( 2.0 * 2.0 );
	Real const lys_arg_dis2_threshold( 5.0 * 5.0 );
	Real const lys_arg_bonus( -2.0 );
	Real const polar_dis2_threshold( 5.0 * 5.0 );
	Real const polar_bonus( -1.0 );
	Real const asp_glu_dis2_threshold( 5.0 * 5.0 );
	Real const asp_glu_bonus( -1.0 );

	Real const backbone_dis2_threshold( 4.0 * 4.0 );
	Real const min_backbone_dis2( 1.75 * 1.75 );
	Real const backbone_O_bonus( -2.0 );
	Real const backbone_N_bonus( -2.0 );

	if ( rsd1.is_protein() && rsd2.is_protein() ) {
		// protein-protein
		return;

	} else if ( rsd1.is_RNA() && rsd2.is_RNA() ) {
		// rna-rna
		return;

	} else if ( ( rsd1.is_RNA() && rsd2.is_protein() ) ||
							( rsd2.is_RNA() && rsd1.is_protein() ) ) {
		// protein-rna

		Residue const & protein_rsd( rsd1.is_protein() ? rsd1 : rsd2 );
		Residue const &     rna_rsd( rsd1.is_protein() ? rsd2 : rsd1 );

		Vector const & xyz_o1p   ( rna_rsd.xyz("OP2" ) );
		Vector const & xyz_o2p   ( rna_rsd.xyz("OP1" ) );
		Vector const & xyz_o2prime( rna_rsd.xyz("O2'") );

		// protein centroid (protein nbr-atom) contacts to rna
		Real centroid_score( 0.0 );
		{
			AA const & aa( protein_rsd.aa() );
			if ( aa == aa_lys || aa == aa_arg ) {
				// pos charged
				Real const dis2_a( xyz_o1p.distance_squared( protein_rsd.nbr_atom_xyz() ) );
				Real const dis2_b( xyz_o2p.distance_squared( protein_rsd.nbr_atom_xyz() ) );
				if ( dis2_a > min_centroid_dis2 && dis2_a < lys_arg_dis2_threshold ) centroid_score += lys_arg_bonus;
				if ( dis2_b > min_centroid_dis2 && dis2_b < lys_arg_dis2_threshold ) centroid_score += lys_arg_bonus;

			} else if ( aa == aa_ser || aa == aa_thr || aa == aa_asn || aa == aa_gln ) {
				// polar
				Real const dis2_a(    xyz_o1p.distance_squared( protein_rsd.nbr_atom_xyz() ) );
				Real const dis2_b(    xyz_o2p.distance_squared( protein_rsd.nbr_atom_xyz() ) );
				Real const dis2_c( xyz_o2prime.distance_squared( protein_rsd.nbr_atom_xyz() ) );
				if ( dis2_a > min_centroid_dis2 && dis2_a < polar_dis2_threshold ) centroid_score += polar_bonus;
				if ( dis2_b > min_centroid_dis2 && dis2_b < polar_dis2_threshold ) centroid_score += polar_bonus;
				if ( dis2_c > min_centroid_dis2 && dis2_c < polar_dis2_threshold ) centroid_score += polar_bonus;

			} else if ( aa == aa_asp || aa == aa_glu ) {
				// neg charged
				Real const dis2_a( xyz_o2prime.distance_squared( protein_rsd.nbr_atom_xyz() ) );
				if ( dis2_a > min_centroid_dis2 && dis2_a < asp_glu_dis2_threshold ) centroid_score += asp_glu_bonus;
			}
		}

		// protein backbone contacts to rna
		Real backbone_score( 0.0 );
		{
			Real const dis2_a( protein_rsd.xyz( "O" ).distance_squared( xyz_o2prime ) );
			Real const dis2_b( protein_rsd.xyz( "N" ).distance_squared( xyz_o1p ) );
			Real const dis2_c( protein_rsd.xyz( "N" ).distance_squared( xyz_o2p ) );
			if ( dis2_a > min_backbone_dis2 && dis2_a < backbone_dis2_threshold ) backbone_score += backbone_O_bonus;
			if ( dis2_b > min_backbone_dis2 && dis2_b < backbone_dis2_threshold ) backbone_score += backbone_N_bonus;
			if ( dis2_c > min_backbone_dis2 && dis2_c < backbone_dis2_threshold ) backbone_score += backbone_N_bonus;
		}

		//
		if ( false ) { // true ) {//verbose ) {
			using namespace std; Size const pos1( protein_rsd.seqpos() ); Size const pos2( rna_rsd.seqpos() );
			if (     centroid_score != Real(0.0) ) cout << "cen: "  << pos1 << ' ' << pos2 << ' ' << centroid_score << endl;
			if (     backbone_score != Real(0.0) ) cout << "bb: "   << pos1 << ' ' << pos2 << ' ' << backbone_score << endl;
		}

		emap[ capri_cen  ] += centroid_score;
		emap[ capri_bb   ] += backbone_score;


	} else {
		return;
	}


}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CapriTotalEnergy
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/**

	 Updated version of experimental constraints scoring, plus distance score between SAM CE and rGU N1 and conservation
	 scoring.

	 Color code:

	 Blue == "must" make an atomic contact with protein
	 Cyan == would be nice to hace protein nearby (say within 8-10 Angstroms)
	 Pink == no protein within 8-10 A
	 Red == no protein in atomic contact.

**/
///
///

class CapriTotalEnergy : public scoring::methods::WholeStructureEnergy {
public:

	/// @brief  C-tor
	CapriTotalEnergy()
	{
		read_datafile();
		add_score_types();
	};


	/// clone
	virtual
	scoring::methods::EnergyMethodOP
	clone() const
	{
		return new CapriTotalEnergy( *this );
	}

	CapriTotalEnergy( CapriTotalEnergy const & src ):
		WholeStructureEnergy(),
		rna_atoms_( src.rna_atoms_ )
	{
		add_score_types();
	}

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	finalize_total_energy(
												pose::Pose & pose,
												ScoreFunction const &,
												EnergyMap & emap
												) const;

	virtual
	Distance
	atomic_interaction_cutoff() const { return 5.0; } // guess

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const { return; }

private:
	void
	add_score_types()
	{
		add_score_type( capri_blue );
		add_score_type( capri_cyan );
		add_score_type( capri_pink );
		add_score_type( capri_red  );
		add_score_type( capri_cons );
		add_score_type( capri_dist );
	}


	void
	read_datafile();


	void
	score_experimental_constraints( pose::Pose const & pose, EnergyMap & emap ) const;


	void
	setup_rna_atoms_for_pose( pose::Pose const & pose, utility::vector1< utility::vector1< AtomID > > & rna_atoms ) const;

private:
	utility::vector1< utility::vector1< std::pair< std::string/*atomname*/, Size /*seqpos*/> > > rna_atoms_;

};

///
void
CapriTotalEnergy::read_datafile()
{
	std::map< std::string, CapriScoreType > capri_score_type_map;
	capri_score_type_map[ "blue" ] = BLUE;
	capri_score_type_map[ "cyan" ] = CYAN;
	capri_score_type_map[ "pink" ] = PINK;
	capri_score_type_map[ "red"  ] = RED;
	std::ifstream data( "input/new_constraints.txt" );
	if ( !data.good() ) utility_exit_with_message( "missing new_constraints.txt" );
	std::string line, color, atom_name;
	rna_atoms_.clear();
	rna_atoms_.resize( n_experimental_constraint_scores );
	while ( getline( data, line ) ) {
		std::istringstream l(line);
		Size seqpos;
		l >> color >> seqpos;
		l >> atom_name;
		if ( l.fail() || !capri_score_type_map.count( color ) ) {
			utility_exit_with_message( "bad format in new_constraints.txt" );
		}
		CapriScoreType const capri_score_type( capri_score_type_map[ color ] );
		while ( !l.fail() ) {
			rna_atoms_[ capri_score_type ].push_back( std::make_pair( atom_name, seqpos ) );
			l >> atom_name;
		}
	}
	data.close();
}

///
void
CapriTotalEnergy::setup_rna_atoms_for_pose(
																					 pose::Pose const & pose,
																					 utility::vector1< utility::vector1< AtomID > > & rna_atoms
																					 ) const
{
	rna_atoms.clear();
	rna_atoms.resize( n_experimental_constraint_scores );

	Size sam_pos( 0 );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).name3() == "SAM" ) {
			sam_pos = i;
			break;
		}
	}
	assert( sam_pos );

	for ( Size t=1; t<= n_experimental_constraint_scores; ++t ) {
		for ( Size i=1; i<= rna_atoms_[t].size(); ++i ) {
			rna_atoms[t].push_back( pid( rna_atoms_[t][i].first, rna_atoms_[t][i].second + sam_pos, pose ) );
		}
	}
}

/// @details  Calculates the cyan score and the conservation score
void
CapriTotalEnergy::finalize_total_energy(
																				pose::Pose & pose,
																				ScoreFunction const &,
																				EnergyMap & emap
																				) const
{
	// total guesses
	Real const conservation_bonus( -1.0 );
	Size const max_rna_nbr_count( 4 );

	CapriData const & capri_data( retrieve_capri_data_from_pose( pose ) );

	/// score the experimental constraints blue->red
	score_experimental_constraints( pose, emap );

	/// figure out where SAM and rna root position are
	Size sam_pos( 0 ), rna_root_pos( 0 ), nres_rna( 0 );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose.residue(i).name3() == "SAM" ) {
			sam_pos = i;
		} else if ( pose.residue(i).is_RNA() ) {
			++nres_rna;
			assert( sam_pos );
			if ( pose.residue(i).is_bonded( pose.residue( sam_pos ) ) ) {
				assert( !rna_root_pos );
				rna_root_pos = i;
				//break;
			}
		}
	}

	if ( !sam_pos || !rna_root_pos ) utility_exit_with_message("no sam in pose!" );
	assert( nres_rna == ( 74 /*rna length*/ - 8 /* rsds missing in bound rna strx*/ ) );
	assert( sam_pos && sam_pos == pose.total_residue() - nres_rna );


	// now the conservation bonus stuff
	Real conservation_score( 0.0 );

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	for ( CapriData::ConservedResidues::const_iterator it= capri_data.conserved_residues.begin();
				it != capri_data.conserved_residues.end(); ++it ) {
		Size const seqpos( it->first );
		if ( pose.residue( seqpos ).name1() != it->second ) {
			utility_exit_with_message( "sequence mismath in capri15 conservation score" );
		}
		Size rna_nbr_count(0);
		for ( graph::Graph::EdgeListConstIter
						iru  = energy_graph.get_node( seqpos )->const_edge_list_begin(),
						irue = energy_graph.get_node( seqpos )->const_edge_list_end();
					iru != irue; ++iru ) {
			EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
			Size const other_pos( edge->get_other_ind( seqpos ) );
			if ( pose.residue( other_pos ).is_RNA() ) ++rna_nbr_count;
		}
		conservation_score += conservation_bonus * Real( std::min( rna_nbr_count, max_rna_nbr_count ) ) / max_rna_nbr_count;
	}

	// distance constraint
	Real const constraint_distance( pose.residue( sam_pos ).xyz( "CE" ).distance( pose.residue(rna_root_pos).xyz("N1")));
	Real const dist_score( constraint_distance < 6.0 ? 0.0 : numeric::square( constraint_distance - 6.0 ) );

	emap[ capri_cons ] = conservation_score;
	emap[ capri_dist ] = dist_score;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
CapriTotalEnergy::score_experimental_constraints(
																								 pose::Pose const & pose,
																								 //bool const protein_is_fullatom,
																								 EnergyMap & emap
																								 ) const
{
	/**
		 Blue == "must" make an atomic contact with protein
		 Cyan == would be nice to hace protein nearby (say within 8-10 Angstroms)
		 Pink == no protein within 8-10 A
		 Red == no protein in atomic contact.
	**/

	utility::vector1< Real > constraint_distance_by_capri_score_type( n_experimental_constraint_scores, 0.0 );
	utility::vector1< Real > penalty( n_experimental_constraint_scores, 0.0 );
	utility::vector1< int > n_to_exclude( n_experimental_constraint_scores, 0 );

	constraint_distance_by_capri_score_type[ BLUE ] = 4.5 + 1.2; // 1.2 ~ hydrogen-heavyatom bond length?
	constraint_distance_by_capri_score_type[ CYAN ] = 9.0 + 1.2;
	constraint_distance_by_capri_score_type[ PINK ] = 10.0;
	constraint_distance_by_capri_score_type[ RED  ] = 4.5;
	n_to_exclude[ BLUE ] = 1;
	n_to_exclude[ CYAN ] = 3;
	n_to_exclude[ PINK ] = 2;
	n_to_exclude[ RED  ] = 2;
	penalty[ PINK ] = 1.0;
	penalty[ RED  ] = 1.0;

	// this is a nuisance
	utility::vector1< ScoreType > score_type_from_capri_score_type( n_experimental_constraint_scores );
	score_type_from_capri_score_type[ BLUE ] = capri_blue;
	score_type_from_capri_score_type[ CYAN ] = capri_cyan;
	score_type_from_capri_score_type[ PINK ] = capri_pink;
	score_type_from_capri_score_type[ RED  ] = capri_red;


	utility::vector1< utility::vector1< AtomID > > rna_atoms;
	setup_rna_atoms_for_pose( pose, rna_atoms );

	// this could be rewritten to be 4 times faster, duh
	for ( Size t=1; t<= n_experimental_constraint_scores; ++t ) {
		CapriScoreType const capri_score_type = CapriScoreType( t );
		Real const constraint_distance( constraint_distance_by_capri_score_type[ capri_score_type ] );
		utility::vector1< Real > alldevs, alldis2;

		for ( vector1< AtomID >::const_iterator it= rna_atoms[t].begin(), ite=rna_atoms[t].end(); it != ite; ++it ) {
			Vector const & rna_atom_xyz( pose.xyz( *it ) );
			// measures the possibility that there's an atom within constraint_distance[ capri_score_type ]
			// is 0.0 if there could be a heavyatom within that distance, based on rsd.nbr_radius()
			// increasingly positive as the closest possible atom is further away than that distance
			Real mindev ( 999.9 );
			// the smallest dis2 between the constrained rna atom and the nbr atom of a protein residue
			Real mindis2( 999.9 );
			for ( Size i=1; i<= pose.total_residue(); ++i ) {
				Residue const & rsd( pose.residue( i ) );
				if ( rsd.is_protein() ) {
					Real const nbr_dis2( rna_atom_xyz.distance_squared( rsd.nbr_atom_xyz() ) );
					Real const nbr_dis2_threshold( ( rsd.nbr_radius() + constraint_distance ) *
																				 ( rsd.nbr_radius() + constraint_distance ) );
					Real const nbr_dis2_dev( std::max( 0.0, nbr_dis2 - nbr_dis2_threshold ) );
					mindev  = std::min( mindev , nbr_dis2_dev );
					mindis2 = std::min( mindis2, nbr_dis2 );
				}
			}
			alldevs.push_back( mindev );
			alldis2.push_back( mindis2 );
		} // loop over rna_atoms[ capri_score_type ]

		std::sort( alldevs.begin(), alldevs.end() );
		std::sort( alldis2.begin(), alldis2.end() );
		assert( alldevs[1] <= alldevs[2] && alldevs[2] <= alldevs.back() );

		// allow a few mismatches
		Size const n_to_count( alldevs.size() - n_to_exclude[ capri_score_type ] );

		Real score(0.0);
		if ( capri_score_type == PINK || capri_score_type == RED ) {
			// only use the nbr atom distances, look for ones that are smaller than constraint_distance ** 2
			std::reverse( alldis2.begin(), alldis2.end() ); // now in decreasing order
			Real const dis2_threshold( numeric::square( constraint_distance ) );
			for ( Size i=1; i<= n_to_count; ++i ) {
				if ( alldis2[i] < dis2_threshold ) {
					score += std::min( penalty[ capri_score_type ], dis2_threshold - alldis2[i] );
				}
			}
		} else {
			Real totaldev( 0.0 );
			for ( Size i=1; i<= n_to_count; ++i ) {
				totaldev += alldevs[ i ];
			}
			score = std::sqrt( totaldev / n_to_count );
		}
		emap[ score_type_from_capri_score_type[ capri_score_type ] ] = score;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void
setup_capri_data( pose::Pose & pose, std::string const & pssm_file, id::SequenceMapping const & fasta2pose )
{

	CapriDataOP capri_data( new CapriData() );



	// now fill in the conservation data, read a PSI-blast pssm file
	capri_data->conserved_residues.clear();
	utility::vector1< Real > pose_IC( pose.total_residue(), 0.0 );
	utility::vector1< Real > all_ics;
	{ // read IC values for pose residues from the file, using the fasta2pose mapping
		std::ifstream data( pssm_file.c_str() );
		if ( !data.good() ) utility_exit_with_message( "cant open pssm file! "+pssm_file );
		std::string line, tag;
		while ( getline( data,line  ) ) {
			std::istringstream l( line );
			Size seqpos;
			Real information_content;
			char seq;
			l >> seqpos >> seq;
			for ( Size i=1; i<= 40; ++i ) l >> tag; // read the pssm and counts
			l >> information_content;
			if ( !l.fail() ) {
				Size const pose_pos( fasta2pose[ seqpos ] );
				if ( pose_pos ) {
					assert( pose.residue( pose_pos ).name1() == seq );
					pose_IC[ pose_pos ] = information_content;
					all_ics.push_back( information_content );
				}
			}
		}
	} // scope

	std::sort( all_ics.begin(), all_ics.end() );
	std::reverse( all_ics.begin(), all_ics.end() );
	assert( all_ics[1] >= all_ics[2] && all_ics[2] >= all_ics.back() );
	Size const n_conserved( all_ics.size() / 10 );
	Real const ic_threshold( all_ics[ n_conserved ] + 1e-3 );
	//std::cout << "IC-threshold: " << n_conserved << ' ' << ic_threshold << std::endl;
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( pose_IC[ i ] >= ic_threshold ) {
			//std::cout << "Conserved residue: " << i << ' ' << pose.residue(i).name() << ' ' << pose_IC[i] << std::endl;
			capri_data->conserved_residues.push_back( std::make_pair( i, pose.residue(i).name1() ) );
		}
	}


	// store the data in the pose
	pose.data().set( basic::CAPRI_DATA, capri_data );

}


// 	/// first fill in the cyan constraint data
// 	{
// 		Size sam_pos( 0 );
// 		for ( Size i=1; i<= pose.total_residue(); ++i ) {
// 			if ( pose.residue(i).name3() == "SAM" ) {
// 				sam_pos = i;
// 				break;
// 			}
// 		}
// 		if ( !sam_pos ) utility_exit_with_message("no sam in pose!" );

// 		capri_data->cyan_atom_map[ sam_pos +  8 ] = "N1";
// 		capri_data->cyan_atom_map[ sam_pos +  9 ] = "N1";
// 		capri_data->cyan_atom_map[ sam_pos + 10 ] = "N3";
// 		capri_data->cyan_atom_map[ sam_pos + 15 ] = "N1";
// 		capri_data->cyan_atom_map[ sam_pos + 17 ] = "N3";
// 		capri_data->cyan_atom_map[ sam_pos + 18 ] = "N1";
// 		capri_data->cyan_atom_map[ sam_pos + 30 ] = "N3";
// 		capri_data->cyan_atom_map[ sam_pos + 45 ] = "N1";
// 		capri_data->cyan_atom_map[ sam_pos + 46 ] = "N1";
// 		capri_data->cyan_atom_map[ sam_pos + 59 ] = "N1";
// 		capri_data->cyan_atom_map[ sam_pos + 60 ] = "N3";
// 		capri_data->cyan_atom_map[ sam_pos + 61 ] = "N1";
// 		capri_data->cyan_atom_map[ sam_pos + 71 ] = "N3";
// 	}
