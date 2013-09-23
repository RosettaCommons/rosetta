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
// AUTO-REMOVED #include <utility/string_util.hh>
// AUTO-REMOVED #include <utility/file/file_sys_util.hh>

#include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
#include <devel/init.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/sequence/util.hh>

// AUTO-REMOVED #include <core/import_pose/pose_stream/util.hh>
// AUTO-REMOVED #include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>

#include <protocols/comparative_modeling/util.hh>

#include <core/sequence/SequenceAlignment.hh>
#include <core/conformation/Residue.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

// option includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>

#include <utility/excn/Exceptions.hh>

utility::vector1< int > calculate_burial(
	core::pose::Pose & mypose,
	core::Real const dist_cutoff
) {
	utility::vector1< int > burial;
	burial.resize( mypose.total_residue() );

	using core::Size;
	for ( Size i = 1; i <= mypose.total_residue(); ++i ) {
		for ( Size j = i + 1; j <= mypose.total_residue(); ++j ) {
			core::conformation::Residue const& resi = mypose.residue(i);
			core::conformation::Residue const& resj = mypose.residue(j);
			core::Real const dist(resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ));

			if ( dist < dist_cutoff ) {
				burial[i]++;
				burial[j]++;
			}
		}
	}
	return burial;
}

int
main( int argc, char* argv [] ) {
	try {

	// options, random initialization
	devel::init( argc, argv );

	using std::map;
	using core::Size;
	using core::Real;
	using std::string;
	using utility::vector1;
	using core::pose::Pose;
	using utility::file::FileName;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::sequence;
	using namespace core::import_pose::pose_stream;

	//basic::Tracer tr( "ss_features" );

	std::ostream& output( std::cout );
	Size const width( 10 );
	//Size const precision( 3 );  // unused ~Labonte
	Real dist_cutoff( std::numeric_limits< Real >::max() );
	if ( option[ james::dist_thresholds ].user() ) {
		dist_cutoff = option[ james::dist_thresholds ]().front();
	}

	//Size min_seqsep = option[ james::min_seqsep ]();

	using namespace ObjexxFCL::format;

	map< string, Pose > poses = protocols::comparative_modeling::poses_from_cmd_line( option[ in::file::template_pdb ]() );

	output
		<< A( width, "q_resi" )
		<< A( width, "q_resj" )
		<< A( 6, "aa_i" )
		<< A( 6, "aa_j" )
		<< A( 8, "atomi" )
		<< A( 8, "atomj" )
		<< A( width, "buriali" )
		<< A( width, "burialj" )
		<< A( width, "distance" )
		<< A( width, "aln_id" )
		<< std::endl;

	vector1< FileName > aln_files( option[ in::file::alignment ]() );
	for ( vector1< FileName >::const_iterator file_it = aln_files.begin(),
				file_end = aln_files.end(); file_it != file_end; ++file_it
	) {

		vector1< SequenceAlignment > alns(
			read_aln( option[ cm::aln_format ](), *file_it )
		);
		//tr.Debug << "read " << alns.size() << " alignments from " << *file_it
		//	<< std::endl;

		typedef vector1< SequenceAlignment >::iterator align_iter;
		for ( align_iter it = alns.begin(), end = alns.end(); it != end; ++it ) {
			//tr.Debug << "processing alignment between " << it->sequence(1)->id()
			//	<< " and " << it->sequence(2)->id() << std::endl;

			std::string const aln_id = it->sequence(2)->id();
			std::string const template_id = aln_id.substr(0,5);

			map< string, Pose >::iterator pose_it = poses.find( template_id );
			if ( pose_it == poses.end() ) {
				//print_seq_map( std::cerr, seqs );
				string msg( "Error: can't find seq (id = " + template_id + ")" );
				//utility_exit_with_message(msg);
				//tr.Error << msg << std::endl;
				//continue;
			} else {
				core::pose::Pose pose = pose_it->second;;
				utility::vector1< int > burial = calculate_burial( pose, 8 );

				for ( Size ii = 1; ii <= it->length(); ++ii ) {
						for ( Size jj = ii+1; jj <= it->length(); ++jj ) {
							if ( it->is_gapped( ii ) || it->is_gapped(jj) ) continue;

							Size const q_resi( it->sequence(1)->resnum(ii) );
							Size const q_resj( it->sequence(1)->resnum(jj) );
							Size const t_resi( it->sequence(2)->resnum(ii) );
							Size const t_resj( it->sequence(2)->resnum(jj) );

							bool const identity(
								it->sequence(1)->at(ii) == it->sequence(2)->at(ii) &&
								it->sequence(1)->at(jj) == it->sequence(2)->at(jj)
							);
							const core::conformation::Residue& resi = pose.residue(t_resi);
							const core::conformation::Residue& resj = pose.residue(t_resj);

							for ( Size m = 1; m <= resi.natoms(); ++m ) {
									for ( Size n = 1; n <= resj.natoms(); ++n ) {

										// skip hydrogen atoms
										if ( resi.atom_type(m).is_hydrogen() || resj.atom_type(n).is_hydrogen() )
											continue;

										// only take sc-sc and bb-bb pairs
										if ( resi.atom_is_backbone(m) != resj.atom_is_backbone(n) )
											continue;

										// skip non-identical side-chain atoms
										if ( (!resi.atom_is_backbone(m) || !resj.atom_is_backbone(n)) && !identity )
											continue;

										// only keep CA backbone atoms for now
										std::string const & atom_m( resi.atom_type(m).name() );
										std::string const & atom_n( resj.atom_type(n).name() );
										if ( resi.atom_is_backbone(m) && resj.atom_is_backbone(n) && (atom_m != "CAbb" || atom_n != "CAbb" ) )
											continue;

										// determine whether distance is below threshold
										core::Real const distance(resi.xyz(m).distance(resj.xyz(n)));
										if ( distance > dist_cutoff ) continue;

										output
												<< I( 10, q_resi)
												<< I( 10, q_resj)
												<< A(  6, resi.name1() )
												<< A(  6, resj.name1() )
												<< A(  8, resi.atom_name(m) )
												<< A(  8, resj.atom_name(n) )
												<< I( 10, burial[t_resi] )
												<< I( 10, burial[t_resj] )
												<< F( 10, 4, distance )
												<< A( 14, aln_id )
												<< std::endl;
							} // n
						} // m
					} // jj
				} // ii
			} // found a pose
		} // for alns
	} // for aln_files
	
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
} // int main( int argc, char * argv [] )
