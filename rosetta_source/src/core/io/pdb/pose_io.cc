// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
///
/// @brief
/// @author

// Unit headers
#include <core/io/pdb/pose_io.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/orbitals/OrbitalType.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/chemical/residue_io.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/scoring/Energies.hh>

// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>
// AUTO-REMOVED #include <core/id/AtomID_Mask.hh>

#include <ObjexxFCL/format.hh>
// AUTO-REMOVED #include <ObjexxFCL/ObjexxFCL.hh>


//#include <core/io/pdb/pdb_dynamic_reader.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
// Auto-header: duplicate removed #include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/chemical/AtomTypeSet.hh>
#include <core/id/AtomID_Map.hh>
//#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

//Auto using namespaces
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end




//#include <fstream>

/// A temporary copy of the pose_from_pdb code from the demo directory.
/// Will be phased out in favor of file_data routines soon.
///

namespace core {
namespace io {
namespace pdb {


/// special Tracer instance acting as special param for all traced_dump_pdb functions
basic::Tracer TR_dump_pdb_dummy( "core.io.pdb.pose_io.dump_pdb_dummy" );

basic::Tracer TR("core.io.pose_io");

using utility::vector1;

void
dump_pdb(
	pose::Pose const & pose,
	std::ostream & out,
	id::AtomID_Mask const & mask,
	std::string const & tag
) {
	Size const nres( pose.total_residue() );

	Size number(0);

	static std::string const chains( " ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890" );

	out << "MODEL     " << tag << "\n";
	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );

 			if ( ! mask[ id::AtomID( j,i ) ] ) continue;

			//skip outputing virtual atom unless specified.
 			//fixed so that the last atom in atom type set can be something other than a virtual atom --steven combs
			if ( !basic::options::option[ basic::options::OptionKeys::out::file::output_virtual ]() &&
				rsd.atom_type(j).is_virtual() ) continue;

			++number;
			runtime_assert( rsd.chain() < chains.size() ); // silly restriction
			char const chain( chains[ rsd.chain() ] );
			out << "ATOM  " << I(5,number) << ' ' << rsd.atom_name(j) << ' ' <<
				rsd.name3() << ' ' << chain << I(4,rsd.seqpos() ) << "    " <<
				F(8,3,atom.xyz()(1)) <<
				F(8,3,atom.xyz()(2)) <<
				F(8,3,atom.xyz()(3)) <<
				F(6,2,1.0) << F(6,2,1.0) << '\n';

			//now add orbitals if the atom type has orbitals
			if(basic::options::option[ basic::options::OptionKeys::out::file::output_orbitals] &&
					rsd.atom_type(j).atom_has_orbital()){
				utility::vector1<core::Size> orbital_indices(rsd.bonded_orbitals(j));
				foreach(core::Size orbital_index, orbital_indices){
					++number;
					Vector orbital_xyz(rsd.orbital_xyz(orbital_index));
					out << "ATOM  " << I(5,number) << ' ' << rsd.orbital_name(orbital_index) << ' ' <<
						rsd.name3() << ' ' << chain << I(4,rsd.seqpos() ) << "    " <<
						F(8,3,orbital_xyz.x()) <<
						F(8,3,orbital_xyz.y()) <<
						F(8,3,orbital_xyz.z()) <<
						F(6,2,1.0) << F(6,2,1.0) << '\n';
				}


			}



		}
	}
	out << "ENDMDL\n";
}


///////////////////////////////////////////////////////////////////////////////
void
dump_bfactor_pdb(
	pose::Pose const & pose,
	id::AtomID_Map< Real > const & bfactor,
	std::ostream & out,
	std::string const & tag
)
{
	int const nres( pose.total_residue() );

	int number(0);

	static std::string const chains( " ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890" );

	if(tag!="NO_MODEL_LINE_IN_OUTPUT") out << "MODEL     " << tag << "\n";
	for ( int i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );

			++number;
			runtime_assert( rsd.chain() < chains.size() ); // silly restriction
			char const chain( chains[ rsd.chain() ] );
			out << "ATOM  " << I(5,number) << ' ' << rsd.atom_name(j) << ' ' <<
				rsd.name3() << ' ' << chain << I(4,rsd.seqpos() ) << "    " <<
				F(8,3,atom.xyz()(1)) <<
				F(8,3,atom.xyz()(2)) <<
				F(8,3,atom.xyz()(3)) <<
				F(6,2,1.0) << F(6,2, bfactor[ id::AtomID(j,i) ] ) << '\n';

			//now add orbitals if the atom type has orbitals
			if(basic::options::option[ basic::options::OptionKeys::out::file::output_orbitals] &&
					rsd.atom_type(j).atom_has_orbital()){
				utility::vector1<core::Size> orbital_indices(rsd.bonded_orbitals(j));
				foreach(core::Size orbital_index, orbital_indices){
					++number;
					Vector orbital_xyz(rsd.orbital_xyz(orbital_index));
					out << "ATOM  " << I(5,number) << ' ' << rsd.orbital_name(orbital_index) << ' ' <<
						rsd.name3() << ' ' << chain << I(4,rsd.seqpos() ) << "    " <<
						F(8,3,orbital_xyz.x()) <<
						F(8,3,orbital_xyz.y()) <<
						F(8,3,orbital_xyz.z()) <<
						F(6,2,1.0) << F(6,2,1.0) << '\n';
				}


			}


		} // 	for ( int i=1; i<= nres; ++i )
	} // 	for ( Size j=1; j<= rsd.natoms(); ++j )
	if(tag!="NO_MODEL_LINE_IN_OUTPUT") out << "ENDMDL\n";
} // dump_bfactor_pdb


///////////////////////////////////////////////////////////////////////////////
void
dump_pdb_residue(
	conformation::Residue const & rsd,
	Size & atom_number,
	std::ostream & out
) {

	static std::string const chains( " ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890" );

	for ( Size j=1; j<= rsd.natoms(); ++j ) {
		conformation::Atom const & atom( rsd.atom(j) );

		++atom_number;
		runtime_assert( rsd.chain() < chains.size() ); // silly restriction
		char const chain( chains[ rsd.chain() ] );
		out << "ATOM  " << I(5,atom_number) << ' ' << rsd.atom_name(j) << ' ' <<
			rsd.name3() << ' ' << chain << I(4,rsd.seqpos() ) << "    " <<
			F(8,3,atom.xyz()(1)) <<
			F(8,3,atom.xyz()(2)) <<
			F(8,3,atom.xyz()(3)) <<
			F(6,2,1.0) << F(6,2,0.0) << '\n';
		if(basic::options::option[ basic::options::OptionKeys::out::file::output_orbitals] &&
				rsd.atom_type(j).atom_has_orbital()){
			utility::vector1<core::Size> orbital_indices(rsd.bonded_orbitals(j));

			foreach(core::Size orbital_index, orbital_indices){
				Vector orbital_xyz(rsd.orbital_xyz(orbital_index));
				out << "ATOM  " << I(5,atom_number) << ' ' << rsd.orbital_name(orbital_index) << ' ' <<
					rsd.name3() << ' ' << chain << I(4,rsd.seqpos() ) << "    " <<
					F(8,3,orbital_xyz.x()) <<
					F(8,3,orbital_xyz.y()) <<
					F(8,3,orbital_xyz.z()) <<
					F(6,2,1.0) << F(6,2,1.0) << '\n';
			}


		}

	}
}



///////////////////////////////////////////////////////////////////////////////
void
dump_pdb(
	pose::Pose const & pose,
	std::ostream & out,
	std::string const & tag
) {
	pose.dump_pdb(out, tag);
}


void
dump_pdb(
	pose::Pose const & pose,
	std::string const & filename,
	std::string const & tag
) {
	pose.dump_pdb(filename, tag);
	/*
	std::ofstream out( filename.c_str() );
	dump_pdb( pose, out );
	out.close();
	*/
}


/// @brief dump_pdb depending on visibility of tracer
/// @param[in] tr   output performed if tracer is visible or if passed dummy
///  tracer core::io::pdb::TR_dump_pdb_dummy
void
traced_dump_pdb(
	basic::Tracer const & tr,
	pose::Pose const & pose,
	std::ostream & out,
	std::string const & tag
)
{
	if ( ( &tr ) != ( &core::io::pdb::TR_dump_pdb_dummy ) ) {
		if ( !tr.visible() ) {
			return;
		}
	}

	pose.dump_pdb( out, tag );
}


/// @brief dump_pdb depending on visibility of tracer
/// @param[in] tr   output performed if tracer is visible or if passed dummy
///   tracer core::io::pdb::TR_dump_pdb_dummy
void
traced_dump_pdb(
	basic::Tracer const & tr,
	pose::Pose const & pose,
	std::string const & filename,
	std::string const & tag
)
{
	if ( ( &tr ) != ( &core::io::pdb::TR_dump_pdb_dummy ) ) {
		if ( !tr.visible() ) {
			return;
		}
	}

	pose.dump_pdb( filename, tag );
}


// the old way: build termini variants after appending non-terminus residues onto the pose in from the pdb...

// void
// core::import_pose::pose_from_pdb(
// 	pose::Pose & pose,
// 	chemical::ResidueTypeSet const & residue_set,
// 	std::string const & filename
// )
// {
// 	//using namespace core;
// 	using namespace conformation;

// 	typedef numeric::xyzVector< Real > Vector;

// 	// reset current data
// 	pose.clear();

// 	Coords coords;
// 	Strings resids, sequence, pose_resids;

// 	read_pdb( filename, resids, sequence, coords );

// 	int const nres_pdb( resids.size() );

// 	char prev_chain('?');


// 	for ( int i=1; i<= nres_pdb; ++i ) {
// 		std::string const pdb_name( sequence[i] );
// 		std::string const resid( resids[i] );
// 		runtime_assert( resid.size() == 6 );
// 		char const chain( resid[5] );

// 		ResidueCoords const & xyz( coords.find( resid )->second );

// 		ResidueTypeCAPs const & rsd_type_list( residue_set.name3_map( pdb_name ) );
// 		if ( rsd_type_list.empty() ) {
// 			std::cout << "Unrecognized aa: " << pdb_name << '\n';
// 			continue;
// 		}

// 		// look for perfect match:
// 		bool matched( false );
// 		for ( Size j=1; j<= rsd_type_list.size(); ++j ) {
// 			ResidueType const & rsd_type( *(rsd_type_list[j]) );

// 			if ( rsd_type.is_terminus() ) continue; // no termini at this stage

// 			int rsd_missing(0), xyz_missing(0);

// 			for ( Size k=1; k<= rsd_type.natoms(); ++k ) {
// 				if ( xyz.count( rsd_type.atom_name(k) ) == 0 ) ++xyz_missing;
// 			}

// 			for ( ResidueCoords::const_iterator iter=xyz.begin(), iter_end=xyz.end(); iter!= iter_end; ++iter ) {
// 				if ( !rsd_type.has( iter->first ) ) ++rsd_missing;
// 			}

// 			if ( rsd_missing ) continue;

// 			std::cout << "match: " << i << ' ' << rsd_type.name() << ' ' << xyz_missing << std::endl;

// 			matched = true;

// 			// found a perfect match! fill in the coords
// 			ResidueOP new_rsd( ResidueFactory::create_residue( rsd_type ) );

// 			for ( ResidueCoords::const_iterator iter=xyz.begin(), iter_end=xyz.end(); iter!= iter_end; ++iter ) {
// 				new_rsd->atom( iter->first ).xyz( iter->second );
// 			}


// 			if ( chain != prev_chain && pose.total_residue() ) {
// 				pose.append_residue( new_rsd, true, pose.total_residue() );
// 			} else {
// 				pose.append_residue( new_rsd );
// 			}
// 			pose_resids.push_back( resid );

// 			break;
// 		} // j=1,rsd_type_list.size()


// 		if ( !matched ) {
// 			// unforgiving for testing purposes
// 			std::cout << "Unrecognized residue: " << pdb_name << std::endl;
// 			utility_exit();
// 		}

// 		// handle termini
// 		if ( chain != prev_chain ) {
// 			prev_chain = chain;
// 			int const seqpos( pose.total_residue() );
// 			if ( seqpos > 1 ) {
// 				pose.conformation().insert_chain_ending( seqpos - 1 );
// 				// make previous residue a terminus
// 				make_upper_terminus( pose, residue_set, seqpos-1 );
// 			}
// 			make_lower_terminus( pose, residue_set, seqpos );

// 		}

// 	} // i=1,nres_pdb

// 	make_upper_terminus( pose, residue_set, pose.total_residue() );

// 	// now handle missing atoms
// 	id::AtomID_Mask missing( false );

// 	id::initialize( missing, pose ); // dimension the missing-atom mask

// 	for ( Size i=1; i<= pose.total_residue(); ++i ) {
// 		ResidueCoords const & xyz( coords.find( pose_resids[i] )->second );
// 		Residue const & rsd( pose.residue(i) );
// 		for ( Size j=1; j<= rsd.natoms(); ++j ) {
// 			if ( xyz.count( rsd.atom_name(j) ) == 0 ) missing[ id::AtomID( j, i ) ] = true;
// 		}
// 	}

// 	pose.conformation().fill_missing_atoms( missing );
// }

// @brief Write energies information into an output stream (e.g. the tail of a pdb file)
void extract_scores(
	core::pose::Pose const & pose,
	utility::io::ozstream & out
)
{
// 	if(!pose.energies().energies_updated()){
// 		out << "Pose's energies were not current, PDBJobOutputter will not force update" << std::endl;
// 		return;
// 	}

	//This is shamelessly refactored from the older JobDistributor; Jobdistributors.hh:1018; SVN 25940
	// APL: Moving this job-independent code into a central location.
	// Which score terms to use
	core::scoring::EnergyMap weights = pose.energies().weights();
	typedef utility::vector1<core::scoring::ScoreType> ScoreTypeVec;
	ScoreTypeVec score_types;
	for(int i = 1; i <= core::scoring::n_score_types; ++i) {
		core::scoring::ScoreType ii = core::scoring::ScoreType(i);
		if ( weights[ii] != 0 ) score_types.push_back(ii);
	}
	// This version is formatted for easy parsing by R, Excel, etc.
	out << "# All scores below are weighted scores, not raw scores.\n";
	out << "#BEGIN_POSE_ENERGIES_TABLE " << out.filename() << "\n";
	out << "label";
	foreach(core::scoring::ScoreType score_type, score_types){
		out << " " << name_from_score_type(score_type);
	}
	out << " total\n";
	out << "weights";
	foreach(core::scoring::ScoreType score_type, score_types){
		out << " " << weights[score_type];
	}
	out << " NA\n";
	out << "pose";
	core::Real pose_total = 0.0;
	if ( pose.energies().energies_updated() ) {
		foreach(core::scoring::ScoreType score_type, score_types){
			core::Real score = (weights[score_type] * pose.energies().total_energies()[ score_type ]);
			out << " " << score;
			pose_total += score;
		}
		out << " " << pose_total << "\n";
		for(core::Size j = 1, end_j = pose.total_residue(); j <= end_j; ++j) {
			core::Real rsd_total = 0.0;
			out << pose.residue(j).name() << "_" << j;
			foreach(core::scoring::ScoreType score_type, score_types){
				core::Real score = (weights[score_type] * pose.energies().residue_total_energies(j)[ score_type ]);
				out << " " << score;
				rsd_total += score;
			}
			out << " " << rsd_total << "\n";
		}
	}
	out << "#END_POSE_ENERGIES_TABLE " << out.filename() << "\n";
}


} // namespace pdb
} // namespace io
} // namespace core
