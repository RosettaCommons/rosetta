// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/andrew/find_buns.cc
/// @brief  Application that finds buried unsatisfied hydrogen bonding groups in proteins using the variable-
///         distance to solvent logic.  Outputs a kinemage that highlights the buried-unsatisfied polar groups
///         and also displays all the exposed dots.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <devel/init.hh>
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/ResidueKinWriter.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/id/AtomID.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>

// Utility headers
#include <utility/file/FileName.hh>

// C++ headers
#include <fstream>


utility::vector1< devel::vardist_solaccess::VarSolDRotamerDotsOP >
varsoldist_rotamer_dots_for_pose(
	devel::vardist_solaccess::VarSolDistSasaCalculatorCOP calc,
	core::pose::Pose const & pose
)
{
	using namespace core;
	using namespace core::graph;
	using namespace core::scoring;
	using namespace devel::vardist_solaccess;

	utility::vector1< VarSolDRotamerDotsOP > rotamer_dots( pose.total_residue() );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		rotamer_dots[ ii ] = VarSolDRotamerDotsOP( new VarSolDRotamerDots( core::conformation::ResidueCOP( core::conformation::ResidueOP( new core::conformation::Residue( pose.residue( ii ) ) ) ), calc ) );
		rotamer_dots[ ii ]->increment_self_overlap();
	}


	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		for ( Graph::EdgeListConstIter
				iru  = energy_graph.get_node(ii)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(ii)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			rotamer_dots[ ii ]->intersect_residues( *rotamer_dots[ (*iru)->get_second_node_ind() ] );
		}
	}
	return rotamer_dots;
}

std::list< core::id::AtomID >
buns_for_pose(
	core::pose::Pose const & pose,
	utility::vector1< devel::vardist_solaccess::VarSolDRotamerDotsOP > const & rotamer_dots,
	core::scoring::hbonds::HBondSet const & hbset
)
{
	using namespace core;
	using namespace core::id;
	using namespace core::graph;
	using namespace core::scoring;
	using namespace core::scoring::hbonds;

	std::list< AtomID > buns;
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		for ( Size jj = 1; jj <= pose.residue(ii).nheavyatoms(); ++jj ) {
			AtomID atid( jj, ii );

			if ( pose.residue(ii).atom_type( jj ).is_acceptor() || pose.residue(ii).atom_type( jj ).is_donor() ) {

				if ( pose.residue(ii).atom_type( jj ).is_donor() &&
						pose.residue(ii).type().attached_H_begin( jj ) > pose.residue(ii).type().attached_H_end( jj ) ) {
					// i.e. proline backbone N; not really a donor
					continue;
				}

				if ( rotamer_dots[ ii ]->any_exposed_dots( jj ) ) continue;
				bool exposed_hydrogens = false;
				for ( Size kk = pose.residue(ii).type().attached_H_begin( jj ); kk<= pose.residue(ii).type().attached_H_end( jj ); kk++ ) {
					if ( rotamer_dots[ ii ]->any_exposed_dots( kk ) ) { exposed_hydrogens = true; break; }
				}
				if ( exposed_hydrogens ) continue;

				utility::vector1< HBondCOP > const & jjhbonds( hbset.atom_hbonds( atid ));
				bool attached_h_form_hbonds = false;
				if ( jjhbonds.empty() ) {
					for ( Size kk = pose.residue(ii).type().attached_H_begin( jj ); kk<= pose.residue(ii).type().attached_H_end( jj ); ++kk ) {
						utility::vector1< HBondCOP > const & kkhbonds( hbset.atom_hbonds( AtomID( kk, ii ) ));
						if ( ! kkhbonds.empty() ) {
							attached_h_form_hbonds = true;
							break;
						}
					}
					if ( attached_h_form_hbonds ) continue;
					buns.push_back( atid );
				}
			}
		}
	}
	return buns;
}

core::Vector
pose_nbratom_center_of_mass(
	core::pose::Pose const & pose
)
{
	core::Vector com( 0.0 ); // center of mass
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		com += pose.residue( ii ).xyz( pose.residue( ii ).nbr_atom() );
	}
	com /= pose.total_residue();
	return com;
}

void
write_buns_and_dots_kinemage(
	std::string const & outfile_name,
	std::string const & structure_tag,
	core::pose::Pose const & pose,
	utility::vector1< devel::vardist_solaccess::VarSolDRotamerDotsOP > const & rotamer_dots,
	std::list< core::id::AtomID > const & buns
)
{

	core::conformation::ConformationKinWriter writer;
	std::ofstream fout( outfile_name.c_str() );
	core::conformation::write_kinemage_header( fout, 1, structure_tag, pose_nbratom_center_of_mass( pose ) );
	writer.write_coords( fout, pose.conformation() );
	fout << "@group { Dots } dominant off\n";
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) rotamer_dots[ ii ]->write_dot_kinemage( fout );
	fout << "@group { BUns } dominant\n";
	fout << "@balllist radius= 0.6 color= white\n";
	for ( std::list< core::id::AtomID >::const_iterator iter = buns.begin(); iter != buns.end(); ++iter ) {
		fout << "{" << iter->rsd() << " " << pose.residue( iter->rsd() ).atom_name( iter->atomno() ) << "} P " << pose.xyz( *iter ).x() << " " << pose.xyz( *iter ).y() << " " << pose.xyz( *iter ).z() << "\n";
	}


}


int main( int argc, char * argv [] )
{
	using namespace core;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::id;
	using namespace core::io;
	using namespace core::graph;
	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::scoring::hbonds;

	using namespace devel::vardist_solaccess;

	using namespace basic::options;

	try {
		devel::init( argc, argv );

		utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();

		for ( core::Size ii = 1; ii <= input_jobs.size(); ++ii ) {
			pose::Pose pose;
			core::import_pose::pose_from_file( pose, input_jobs[ ii ]->input_tag() , core::import_pose::PDB_file);
			ScoreFunctionOP sfxn = get_score_function();
			(*sfxn)(pose);

			HBondSet hbset;
			fill_hbond_set( pose, false, hbset );

			VarSolDistSasaCalculatorOP calc( new VarSolDistSasaCalculator );
			utility::vector1< VarSolDRotamerDotsOP > rotamer_dots = varsoldist_rotamer_dots_for_pose( calc, pose );
			std::list< core::id::AtomID > buns = buns_for_pose( pose, rotamer_dots, hbset );

			utility::file::FileName inputname = input_jobs[ ii ]->input_tag();
			std::string outname = inputname.base();
			write_buns_and_dots_kinemage( outname + "_buried_unsats.kin", outname, pose, rotamer_dots, buns );
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;

}

