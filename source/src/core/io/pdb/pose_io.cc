// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/pose_io.cc
/// @brief  method definitions for input/output functions for use with Pose
/// @author

// Unit headers
#include <core/io/pdb/pose_io.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>

// External headers
#include <ObjexxFCL/format.hh>
#include <boost/foreach.hpp>


using namespace ObjexxFCL::format; // AUTO USING NS

// A temporary copy of the pose_from_pdb code from the demo directory.
// Will be phased out in favor of file_data routines soon.

namespace core {
namespace io {
namespace pdb {


// special Tracer instance acting as special param for all traced_dump_pdb functions
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

			//skip outputting virtual atom unless specified.
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
				utility::vector1<core::Size> const & orbital_indices(rsd.bonded_orbitals(j));
				BOOST_FOREACH(core::Size orbital_index, orbital_indices){
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
				utility::vector1<core::Size> const & orbital_indices(rsd.bonded_orbitals(j));
				BOOST_FOREACH(core::Size orbital_index, orbital_indices){
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
			utility::vector1<core::Size> const & orbital_indices(rsd.bonded_orbitals(j));

			BOOST_FOREACH(core::Size orbital_index, orbital_indices){
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

/// @note Python compatible wrapper avoiding reference parameter
void
dump_pdb_residue(
	conformation::Residue const & rsd,
	std::ostream & out,
	Size start_atom_number)
{
	dump_pdb_residue(rsd, start_atom_number, out);
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


// dump_pdb depending on visibility of tracer
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


// dump_pdb depending on visibility of tracer
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

/// @brief Utility function to round a real value to the given precisions (number of digits after the decimal place) for output.
/// For use solely by extract_scores()
/// @details Apparently, there isn't an easy way to do this in C++, or even the general goal
/// of limiting the precision in output streams. (setprecision() with default formatting doesn't
/// correctly handle very small numbers, and with fixed precision outputs superfluous zeros.)

core::Real restrict_prec( core::Real inval )
{
	if( inval >= 1 || inval <= -1 ) { // Don't alter value, as the default precision of 6 works fine, and we avoid rounding artifacts
		return inval;
	}
	core::Real outval;
	std::stringstream temp;
	temp << std::fixed << std::setprecision(5) << inval;
	temp >> outval;
	return outval;
}

// @brief Write energies information into an output stream (e.g. the tail of a pdb file)
void extract_scores(
	core::pose::Pose const & pose,
	utility::io::ozstream & out
)
{
	extract_scores( pose, out.filename(), out );
}


void extract_scores(
	core::pose::Pose const & pose,
	std::string const & filename,
	std::ostream & out
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
	out << "#BEGIN_POSE_ENERGIES_TABLE " << filename << "\n";
	out << "label";
	BOOST_FOREACH(core::scoring::ScoreType score_type, score_types){
		out << " " << name_from_score_type(score_type);
	}
	out << " total\n";
	out << "weights";
	BOOST_FOREACH(core::scoring::ScoreType score_type, score_types){
		out << " " << weights[score_type];
	}
	out << " NA\n";
	out << "pose";
	core::Real pose_total = 0.0;
	if ( pose.energies().energies_updated() ) {
		BOOST_FOREACH(core::scoring::ScoreType score_type, score_types){
			core::Real score = (weights[score_type] * pose.energies().total_energies()[ score_type ]);
			out << " " << restrict_prec(score);
			pose_total += score;
		}
		out << " " << restrict_prec(pose_total) << "\n";
		for(core::Size j = 1, end_j = pose.total_residue(); j <= end_j; ++j) {
			core::Real rsd_total = 0.0;
			out << pose.residue(j).name() << "_" << j;
			BOOST_FOREACH(core::scoring::ScoreType score_type, score_types){
				core::Real score = (weights[score_type] * pose.energies().residue_total_energies(j)[ score_type ]);
				out << " " << restrict_prec(score);
				rsd_total += score;
			}
			out << " " << restrict_prec(rsd_total) << "\n";
		}
	}
	out << "#END_POSE_ENERGIES_TABLE " << filename << "\n";
}


/////////////////////////////////////////////////////////////////////
void
dump_connect_info(
	pose::Pose const & pose,
	std::ostream & out,
	std::map< id::AtomID, Size > & atom_id_output ){

	// It might be better to put this in file_data.cc, and to make use of the file_data to get the
	// numbering exactly right.
	Real const CUTOFF( basic::options::option[ basic::options::OptionKeys::inout::connect_info_cutoff ]() ); //3.0 by default, but can be set by user
	Real const CUTOFFSQ ( CUTOFF*CUTOFF );

	for ( Size i=1, imax=pose.total_residue(); i<=imax; ++i ) {

		conformation::Residue const & rsd( pose.residue(i) );

		for ( Size j=1, jmax=rsd.natoms(); j<=jmax; ++j ) {

			id::AtomID const atom_id( j, i );

			if ( !atom_id_output[ atom_id ] ) continue;

			utility::vector1<core::id::AtomID>	const & nbr_ids(  pose.conformation().bonded_neighbor_all_res( atom_id, true /*virt*/) );
			for ( Size n = 1, nmax=nbr_ids.size(); n <=nmax; n++ ) {

				id::AtomID const & nbr_id = nbr_ids[ n ];

				if ( !atom_id_output[ nbr_id ] ) continue;

				if ( atom_id.rsd() > nbr_id.rsd()   ) continue;
				if ( atom_id.rsd() == nbr_id.rsd()  && atom_id.atomno() > nbr_id.atomno()  ) continue;

				if ( ( pose.xyz( atom_id ) - pose.xyz( nbr_id ) ).length_squared() < CUTOFFSQ ) continue;

				// Final check: actually look for a connection in the atom tree.
				// V. Mulligan, 2 March 2014: This is a silly check!  It guarantees that the only connections we write out are between atoms and their parents,
				// which means that we NEVER write inter-residue connections between side-chains, or even connections that close rings.  This means that most
				// of the business above, checking residue numbers and whatnot, is for naught.  Commented out for now.
				//core::kinematics::tree::AtomCOP atom( & pose.atom_tree().atom_dont_do_update( atom_id ) );
				//core::kinematics::tree::AtomCOP  nbr( & pose.atom_tree().atom_dont_do_update( nbr_id ) );
				//if( ( nbr->parent() != atom ) && ( atom->parent() != nbr ) ) continue;

				out << "CONECT" << I(5,atom_id_output[ atom_id ]) << I(5,atom_id_output[  nbr_id  ]) << std::endl;

			}
		}
	}
}

/////////////////////////////////////////////////////////////////////
// useful for centroid poses -- spit out "CONECT" information on bonded atoms
//  that might otherwise be ignored by pymol/rasmol.
void
dump_connect_info(
	pose::Pose const & pose,
	std::ostream & out ){

	std::map< id::AtomID, Size  > atom_id_output;

	// It might be better to put this in file_data.cc, and to make use of the file_data to get the
	// numbering exactly right.
	Size count( 0 );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		 for ( Size j=1; j<= pose.residue(i).natoms(); ++j ) {

			 //skip outputting virtual atom unless specified
			 atom_id_output[ id::AtomID(j,i) ] = 0;

			 if ( !basic::options::option[ basic::options::OptionKeys::out::file::output_virtual ]() &&	pose.residue(i).is_virtual( j ) ) continue;

			 count++;
			 atom_id_output[ id::AtomID(j,i) ] = count;
		 }
	 }

	 dump_connect_info( pose, out, atom_id_output );
}

} // namespace pdb
} // namespace io
} // namespace core
