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
/// @author James Thompson

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/NamedAtomID.hh>
#include <basic/database/open.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/import_pose/import_pose.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <protocols/jumping/util.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <cmath>

#include <core/id/AtomID.hh>

#ifndef JAMES_UTIL
#define JAMES_UTIL

class DistancePrediction : public utility::pointer::ReferenceCount {
public:
	DistancePrediction(
		core::id::AtomID atom_i,
		core::id::AtomID atom_j,
		core::Real predicted_distance,
		utility::vector1< core::Real > vars) :
		atom_i_( atom_i ),
		atom_j_( atom_j ),
		predicted_distance_( predicted_distance ),
		vars_( vars )
	{}

	core::Real predicted_distance() const {
		return predicted_distance_;
	}

	core::id::AtomID const & atom_i() const {
		return atom_i_;
	}

	core::id::AtomID const & atom_j() const {
		return atom_j_;
	}

	utility::vector1< core::Real > vars() const {
		return vars_;
	}

private:
	core::id::AtomID atom_i_, atom_j_;
	core::Real predicted_distance_;
	utility::vector1< core::Real > vars_;
};

class DistanceMatrix : public utility::pointer::ReferenceCount {

public:
	DistanceMatrix(
		core::pose::Pose const & pose,
		std::string const & atom_name)
	: atom_name_( atom_name )
	{
		calculate_distances( pose );
	}

	void calculate_distances( core::pose::Pose const & pose ) {
		using core::Real;
		using utility::vector1;
		distances_.resize(
			pose.total_residue(), vector1< Real >( pose.total_residue(), 0.0 ));

		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			for ( Size jj = ii + 1; jj <= pose.total_residue(); ++jj ) {
				core::id::NamedAtomID atom_ii( atom_name(), ii );
				core::id::NamedAtomID atom_jj( atom_name(), jj );

				core::Real const dist(
					pose.xyz(atom_ii).distance( pose.xyz(atom_jj) )
				);

				distances_[ii][jj] = dist;
				distances_[jj][ii] = dist;
			} // jj
		} // ii
	} // calculate_distances

	std::string atom_name() const {
		return atom_name_;
	}

	core::Real distance( core::Size const ii, core::Size const jj ) const {
		return distances_[ii][jj];
	}

private:
	std::string atom_name_;
	utility::vector1< utility::vector1< core::Real > > distances_;
};

class TorsionList : public utility::pointer::ReferenceCount {

public:
	TorsionList(
		core::pose::Pose const & pose) {
		calculate_torsions( pose );
	}

	void calculate_torsions( core::pose::Pose const & pose ) {
		using core::Real;
		using utility::vector1;
		torsions_.resize( pose.total_residue() );

		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			torsions_[ii] = pose.residue(ii).mainchain_torsions();
		}
	} // calculate_torsions

	core::Real torsion(
		core::Size const residue,
		core::Size const torsion_idx) const {
		return torsions_[residue][torsion_idx];
	}

	Size n_torsions() const {
		return torsions_.front().size();
	}

private:
	// mainchain torsions from Pose residues
	utility::vector1< utility::vector1< core::Real > > torsions_;
};

utility::vector1< utility::vector1< core::Real > > get_ca_distances(
	core::pose::Pose pose,
	std::string const & atom_name) {
	using core::Size;
	using core::Real;
	using utility::vector1;

	vector1< vector1< Real > > distances(
		pose.total_residue(), vector1< Real >( pose.total_residue(), 0.0 ));
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		for ( Size j = i + 1; j <= pose.total_residue(); ++j ) {
			core::conformation::Residue resi = pose.residue(i);
			core::conformation::Residue resj = pose.residue(j);

			distances[i][j]
				= resi.xyz( atom_name ).distance( resj.xyz( atom_name ) );
			distances[j][i]
				= resi.xyz( atom_name ).distance( resj.xyz( atom_name ) );
		}
	}
	return distances;
}

utility::vector1< int > calculate_burial(
	core::pose::Pose & mypose,
	core::Real const dist_cutoff) {
	utility::vector1< int > burial;
	burial.resize( mypose.total_residue() );

	using core::Size;
	for ( Size i = 1; i <= mypose.total_residue(); ++i ) {
		for ( Size j = i + 1; j <= mypose.total_residue(); ++j ) {
			core::conformation::Residue resi = mypose.residue(i);
			core::conformation::Residue resj = mypose.residue(j);

			core::Real const dist(
				resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) )
			);
			if ( dist < dist_cutoff ) {
				burial[ i ]++;
				burial[ j ]++;
			}
		}
	}
	return burial;
}

// burial is represented as number of number of neighbor atoms within 8 angstroms
utility::vector1< int > calculate_burial(
	core::pose::Pose & mypose) {

	return calculate_burial( mypose, 8 );
}

utility::vector1< core::Real > calc_blosum_scores(
	core::sequence::SequenceAlignment const & aln) {
	using core::Size;
	using core::Real;
	using utility::vector1;
	using utility::file::FileName;
	using namespace core::sequence;

	utility::io::izstream stream;
	basic::database::open( stream, "sequence/BLOSUM62" );
	ScoringSchemeOP blosum_score( new MatrixScoringScheme( 0, 0, FileName("/work/tex/minirosetta_database/sequence/BLOSUM62")) );

	vector1< Real > scores( aln.length(), 0.0 );
	for ( Size ii = 1; ii <= aln.length(); ++ii ) {
		core::Real score(0.0);
		if ( !aln.sequence(1)->is_gap(ii) && !aln.sequence(2)->is_gap(ii) ) {
			score = blosum_score->score(
				aln.sequence(1), aln.sequence(2), ii, ii
			);
		}
		scores[ii] = score;
	}
	return scores;
}

utility::vector1< char > get_ss( core::pose::Pose & pose ) {
	using core::Size;
	protocols::jumping::assign_ss_dssp( pose );

	utility::vector1< char > ss( pose.total_residue(), 'L' );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		ss[ii] = pose.secstruct(ii);
	}

	return ss;
}

core::Real mean( utility::vector1< core::Real > const values ) {
	core::Real total = 0;
	for ( utility::vector1< core::Real >::const_iterator it = values.begin(), end = values.end(); it != end; ++it ) {
		total += *it;
	}
	core::Real mean = total / values.size();
	return mean;
}

core::Real sd( utility::vector1< core::Real > const values ) {
	core::Real mymean = mean( values );
	core::Real total_diff_squared = 0;

	for ( utility::vector1< core::Real >::const_iterator it = values.begin(), end = values.end(); it != end; ++it ) {
		core::Real diff = mymean - *it;
		total_diff_squared += diff * diff;
	}

	return std::sqrt( total_diff_squared );
}

core::pose::PoseOP get_pose_by_id(
	std::string const & id) {
	using std::map;
	using std::string;
	using core::Size;
	using core::pose::Pose;
	using core::pose::PoseOP;
	using namespace core::chemical;

	static Size max_entries( 500 );
	static map< string, PoseOP > pose_map;

	map< string, PoseOP >::iterator iter = pose_map.find( id );
	if ( iter != pose_map.end() ) {
		return iter->second;
	}

	if ( pose_map.size() > max_entries ) pose_map.clear();

	string const prefix( "/work/tex/projects/nn/pdbs/" );
	string const sub_dir( id.substr(1,2) );
	string full_fn( prefix + sub_dir + "/" + id + ".pdb" );

	// hack for pdb id's that have the wrong chain - try guessing chain A!
	if ( !utility::file::file_exists( full_fn ) ) {
		string id_copy( id.substr(0,4) + 'A' );
		full_fn = prefix + sub_dir + "/" + id_copy + ".pdb";
		if ( !utility::file::file_exists( full_fn ) ) {
			utility_exit_with_message( string( "Error: " + full_fn + " doesn't exist!" ) );
		}
	}

	ResidueTypeSetCOP rsd_set(
		ChemicalManager::get_instance()->residue_type_set( "fa_standard" ) );
	PoseOP pose_op( new Pose );
	core::import_pose::pose_from_pdb( *pose_op, *rsd_set, full_fn );
	protocols::jumping::assign_ss_dssp( *pose_op );
	pose_map[ id ] = pose_op;
	return pose_op;
}

utility::vector1< core::Real > calc_dgaps(
	core::sequence::SequenceAlignment const & aln) {
	using core::Size;
	using core::Real;
	using utility::vector1;

	// calculate gap positions
	vector1< Size > gap_positions;
	for ( Size ii = 1; ii <= aln.length(); ++ii ) {
		if ( aln.is_gapped( ii ) ) gap_positions.push_back( ii );
	}

	vector1< Real > dgap_avgs( aln.length(), 0 );
	for ( Size ii = 1; ii <= aln.length(); ++ii ) {
		// find minimum distance to a position in gap_positions
		Real closest_dgap( static_cast< Real > ( aln.length() ) );
		for ( vector1< Size >::const_iterator it = gap_positions.begin(),
				end = gap_positions.end(); it != end; ++it
				) {
			Real const this_dgap(
				std::abs( static_cast< Real > (*it - ii) )
			);
			if ( this_dgap < closest_dgap ) {
				closest_dgap = this_dgap;
			}
		}
		dgap_avgs[ii] = ( closest_dgap );
	}
	return dgap_avgs;
}

#endif
