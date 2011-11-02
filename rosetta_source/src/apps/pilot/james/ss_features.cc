// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ss_features.cc
/// @brief
/// @author James Thompson

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>

#include <core/io/pdb/pose_io.hh>
#include <basic/database/open.hh>
#include <devel/init.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/sequence/util.hh>

#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>

#include <protocols/comparative_modeling/features/ResidueFeature.hh>
#include <protocols/comparative_modeling/features/TorsionFeature.hh>
#include <protocols/comparative_modeling/features/SSFeature.hh>

#include <apps/pilot/james/james_util.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

// option includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

static basic::Tracer tr( "ss_features" );

// featurelist:
// ln_evalue
// blosum score
// distance from a gap
// burial
// predicted distance (or other feature value in template)

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char* argv [] ) {
	// options, random initialization
	devel::init( argc, argv );

	using core::Size;
	using core::Real;
	using std::string;
	using utility::vector1;
	using utility::file::FileName;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::sequence;
	using namespace core::import_pose::pose_stream;

	Real const max_ln_evalue( 0.0 );

	string outfile( option[ out::file::silent ] () );
	std::ofstream output( outfile.c_str() );
	if ( !output.is_open() ) {
		utility_exit_with_message( "Unable to open file: " + outfile + '\n' );
	}

	Size const width( 10 );
	Size const precision( 3 );

	using namespace ObjexxFCL::fmt;

	output
		<< A( width, "ln_ev" )
		<< A( width, "q_resn" )
		<< A( width, "t_resn" )
		<< A( width, "seq_sc" )
		<< A( width, "dgap" )
		<< A( width, "t_bur" )
		<< A( width, "ss_pred" )
		<< A( width, "ss_nat" )
		<< A( width, "q_id" )
		<< A( width, "t_id" )
		<< std::endl;

	vector1< FileName > aln_files( option[ in::file::alignment ]() );
	for ( vector1< FileName >::const_iterator file_it = aln_files.begin(),
				file_end = aln_files.end(); file_it != file_end; ++file_it
	) {

		vector1< SequenceAlignment > alns(
			read_aln( option[ cm::aln_format ](), *file_it )
		);
		tr.Debug << "read " << alns.size() << " alignments from " << *file_it
			<< std::endl;

		typedef vector1< SequenceAlignment >::iterator align_iter;
		for ( align_iter it = alns.begin(), end = alns.end(); it != end; ++it ) {
			Real ln_e_value = it->score();
			if ( ln_e_value > max_ln_evalue ) continue;
			if ( it->identities() == it->length() ) continue;
			tr.Debug << "processing alignment between " << it->sequence(1)->id()
				<< " and " << it->sequence(2)->id() << std::endl;

			// get the query and template poses
			core::pose::PoseOP query_pose = get_pose_by_id( it->sequence(1)->id() );
			core::pose::PoseOP templ_pose = get_pose_by_id( it->sequence(2)->id() );

			// intermediate alignments from alignment sequences to pdb sequences
			SequenceMapping q_aln_to_pdb = map_seq1_seq2(
				it->sequence(1),
				new Sequence(
					query_pose->sequence(), it->sequence(1)->id() + ".pdb", 1
				)
			);
			SequenceMapping t_aln_to_pdb = map_seq1_seq2(
				it->sequence(2),
				new Sequence(
					templ_pose->sequence(), it->sequence(2)->id() + ".pdb", 1
				)
			);

			vector1< Real > scores  ( calc_blosum_scores( *it ) );
			vector1< char > query_ss( get_ss( *query_pose ) );
			vector1< char > templ_ss( get_ss( *templ_pose ) );
			vector1< int >  t_burial( calculate_burial( *templ_pose ) );
			vector1< Real > dgaps( calc_dgaps( *it ) );

			// iterate over the aligned residue pairs,
			for ( Size ii = 1; ii <= it->length(); ++ii ) {
				if ( it->is_gapped( ii ) ) continue;
				Size const q_resi_aln( it->sequence(1)->resnum(ii) );
				Size const t_resi_aln( it->sequence(2)->resnum(ii) );

				Size const t_resi_pdb( t_aln_to_pdb[ t_resi_aln ] );
				if ( t_resi_pdb == 0 ) continue; // skip residues not in PDB

				Size const q_resi_pdb( q_aln_to_pdb[ q_resi_aln ] );
				if ( q_resi_pdb == 0 ) continue; // skip residues not in PDB

				output
					<< F( width, precision, ln_e_value )
					<< A( width, query_pose->residue(q_resi_pdb).name1() )
					<< A( width, templ_pose->residue(t_resi_pdb).name1() )
					<< F( width, precision, scores[ ii ] )
					<< F( width, precision, dgaps[ ii ] )
					<< I( width, t_burial[ t_resi_pdb ] )
					<< A( width, query_ss[ q_resi_pdb ] )
					<< A( width, templ_ss[ t_resi_pdb ] )
					<< A( width, it->sequence(1)->id() )
					<< A( width, it->sequence(2)->id() )
					<< std::endl;
			} // for sequence positions
		} // for alns
	} // for aln_files
} // int main( int argc, char * argv [] )
