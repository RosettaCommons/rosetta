// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///
/// @brief Given two residue sets, or "faces", computes their interface energy as the sum of pairwise
///  residue energies over all residue pairs (R1, R2), R1 belonging to face #1 and R2 belonging to
///  face #2.
///
/// @author Andrea Bazzoli (ndrbzz@gmail.com)
///
/// ARGUMENTS
///  ---> -s <POSE_PDB>
///  ---> -face1 <FACE1>, where <FACE1> is the path to a file describing face #1 as specified in function
///                       load_set() (see below).
///  ---> -face2 <FACE2>, where <FACE2> is the path to a file describing face #2; the format is the same.
///
///  ---> -score:hbond_bb_per_residue_energy, to include backbone-backbone Hbond energies.
///
/// OUTPUT:
///  On screen, the total interface energy together with the contributions from all pairs of residues.
///

#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>
#include <string>
#include <iostream>
#include <fstream>

using namespace basic::options::OptionKeys;
using basic::options::option;
using core::Size;


static basic::Tracer TR( "interface_energy.main" );


OPT_KEY( String, face1 )
OPT_KEY( String, face2 )


///
/// @brief loads a set of residues from file. The file format is as follows:
///
///  C1 R1 I1\n
///  ...
///  CN RN IN\n
///
///   Here, Ci, Ri, and Ii indicate the chain identifier, residue index, and
///  insertion code (as specified in the pose's input PDB file) of the ith
///  residue in the set (i=1,...,N; N>=1).
///
/// @param[in] path to the input file.
/// @param[out] rset vector to hold the residues. The vector must be passed
///  empty.
/// @param[in] ps the pose.
///
/// @details: after this function has been called, rset[i] is the pose index
///  of the residue specified by the ith input line (i=1,...,N).
///
/// @details: blank chain identifiers and insertion codes must be specified
///  with the '_' character.
///
void load_set(std::string setf, utility::vector1<Size>& rset,
	core::pose::Pose &ps) {

	std::ifstream setfs(setf.c_str());
	if ( !setfs ) {
		TR << "can't open " << setf << std::endl;
		exit(0);
	}

	char cid;
	int idx;
	char ico;
	while ( setfs >> cid >> idx >> ico ) {

		if ( cid == '_' ) {
			cid = ' ';
		}

		if ( ico == '_' ) {
			ico = ' ';
		}

		rset.push_back( ps.pdb_info()->pdb2pose( cid, idx, ico ) );
	}
}


///
/// @brief given a residue in a pose, prints its identifier in the PDB file
///  from which the pose was loaded.
///
/// @param[in] ridx index of one residue in the pose
/// @param[in] ps the pose
///
/// @details blank chain identifiers and insertion codes are printed as '_'.
///
void print_pdb_info(Size ridx, core::pose::Pose const& ps) {

	// one residue
	char ico = ps.pdb_info()->icode(ridx);
	if ( ico == ' ' ) {
		ico = '_';
	}

	char cid = ps.pdb_info()->chain(ridx);
	if ( cid == ' ' ) {
		cid = '_';
	}

	std::cout <<
		core::chemical::name_from_aa(ps.aa(ridx)) <<
		ps.pdb_info()->number(ridx) <<
		ico <<
		cid <<
		"(" << ridx << ")";
}


////////////////////////////////////////////////////////////////////////////////
//                                    MAIN                                    //
////////////////////////////////////////////////////////////////////////////////

int main( int argc, char * argv [] )
{

	try {

		NEW_OPT( face1, "set of residues in one face of the interface", "" );
		NEW_OPT( face2, "set of residues in the other face of the interface", "" );

		devel::init(argc, argv);

		// load pose from pdb file
		core::pose::Pose ps;
		std::string const input_pdb_name( basic::options::start_file() );
		core::import_pose::pose_from_file( ps, input_pdb_name, core::import_pose::PDB_file );

		utility::vector1<Size> fv1;
		load_set(option[face1], fv1, ps);

		utility::vector1<Size> fv2;
		load_set(option[face2], fv2, ps);

		// score pose
		core::scoring::ScoreFunctionOP scorefxn(
			core::scoring::get_score_function());

		(*scorefxn)(ps);

		core::Real tot_wnrg = 0;
		core::scoring::EnergyMap weights = ps.energies().weights();;

		// SHORT RANGE: the output can be seen as a sequence of N1 blocks, where N1 is
		// number of residues in face 1 of the interface. The ith block contains the
		// pairwise energies of the ith residue in face 1's input file. Within the ith
		// block, line j contains the pairwise energy of the jth residue in face 2's
		// input file (i=1,...,N1; j=1,...,N2, where N2 is the number of residues in
		// face 2.
		std::cout << "##### PAIRWISE SHORT-RANGE ENERGIES #####" << std::endl;

		core::scoring::EnergyGraph const & energy_graph( ps.energies().energy_graph() );

		Size const N1 = fv1.size();
		Size const N2 = fv2.size();

		for ( Size i=1; i<=N1; ++i ) {

			Size ri1 = fv1[i];

			for ( Size j=1; j<=N2; ++j ) {

				Size ri2 = fv2[j];

				core::Real wnrg = 0;
				core::scoring::EnergyEdge const* edge = energy_graph.find_energy_edge(ri1, ri2);
				if ( edge ) {
					core::scoring::EnergyMap nrgs = edge->fill_energy_map();
					wnrg = nrgs.dot(weights);
				}

				print_pdb_info(ri1, ps);
				std::cout << " --- ";
				print_pdb_info(ri2, ps);
				std::cout << ": " << wnrg << std::endl;

				tot_wnrg += wnrg;
			}
		}

		// CONTEXT-INDEPENDENT LONG RANGE: the output can be seen as a sequence of N
		// super-blocks, where N is the number of active energy methods, namely,
		// methods such that (!lrec || lrec->empty()) (see line below) evaluates to
		// false. The ith super-block contains the energies for the ith such active
		// method (i=1,...,N; the ci_lr_2b_methods order is assumed). Each super-block
		// is composed of blocks in the same way as the output of short-range energies
		// (described above).
		for ( core::scoring::ScoreFunction::CI_LR_2B_Methods::const_iterator
				iter = scorefxn->ci_lr_2b_methods_begin(),
				iter_end = scorefxn->ci_lr_2b_methods_end();
				iter != iter_end; ++iter ) {

			core::scoring::LREnergyContainerOP lrec =
				ps.energies().nonconst_long_range_container( (*iter)->long_range_type() );
			if ( !lrec || lrec->empty() ) continue;

			std::cout << std::endl <<
				"##### PAIRWISE ENERGIES FOR CI_LR_2B METHOD IN CHARGE OF " <<
				(*iter)->score_types()[1] << " #####" <<
				std::endl;

			for ( Size i=1; i<=N1; ++i ) {

				Size ri1 = fv1[i];

				for ( Size j=1; j<=N2; ++j ) {

					Size ri2 = fv2[j];

					core::Size rmin = (ri1 < ri2) ? ri1 : ri2;
					core::Size rmax = (ri1 < ri2) ? ri2 : ri1;

					core::Real wnrg = 0;
					core::scoring::EnergyMap emap;

					// are the two residues interacting for the current method?
					for ( core::scoring::ResidueNeighborIteratorOP
							rni = lrec->upper_neighbor_iterator_begin( rmin ),
							rniend = lrec->upper_neighbor_iterator_end( rmin );
							(*rni) != (*rniend); ++(*rni) ) {
						if ( rni->upper_neighbor_id() == rmax ) {
							rni->retrieve_energy(emap);
							wnrg = emap.dot(weights);
							break;
						}
					}

					print_pdb_info(ri1, ps);
					std::cout << " --- ";
					print_pdb_info(ri2, ps);
					std::cout << ": " << wnrg << std::endl;

					tot_wnrg += wnrg;
				}
			}
		}

		// CONTEXT-DEPENDENT LONG RANGE: the output can be seen as a sequence of N
		// super-blocks, where N is the number of context-dependent long-range energy
		// methods. The ith super-block contains the energies for the ith such
		// method (i=1,...,N; the cd_lr_2b_methods order is assumed). Each super-block
		// is composed of blocks in the same way as the output of short-range energies
		// (described above).
		for ( core::scoring::ScoreFunction::CD_LR_2B_Methods::const_iterator
				iter = scorefxn->cd_lr_2b_methods_begin(),
				iter_end = scorefxn->cd_lr_2b_methods_end();
				iter != iter_end; ++iter ) {

			core::scoring::LREnergyContainerOP lrec =
				ps.energies().nonconst_long_range_container( (*iter)->long_range_type() );

			std::cout << std::endl <<
				"##### PAIRWISE ENERGIES FOR CD_LR_2B METHOD IN CHARGE OF " <<
				(*iter)->score_types()[1] << " #####" <<
				std::endl;

			for ( Size i=1; i<=N1; ++i ) {

				Size ri1 = fv1[i];

				for ( Size j=1; j<=N2; ++j ) {

					Size ri2 = fv2[j];

					core::Size rmin = (ri1 < ri2) ? ri1 : ri2;
					core::Size rmax = (ri1 < ri2) ? ri2 : ri1;

					core::Real wnrg = 0;
					core::scoring::EnergyMap emap;

					// are the two residues interacting for the current method?
					for ( core::scoring::ResidueNeighborIteratorOP
							rni = lrec->upper_neighbor_iterator_begin( rmin ),
							rniend = lrec->upper_neighbor_iterator_end( rmin );
							(*rni) != (*rniend); ++(*rni) ) {
						if ( rni->upper_neighbor_id() == rmax ) {
							(*iter)->residue_pair_energy(
								ps.residue(rmin), ps.residue(rmax), ps, *scorefxn, emap);
							wnrg = emap.dot(weights);
							break;
						}
					}

					print_pdb_info(ri1, ps);
					std::cout << " --- ";
					print_pdb_info(ri2, ps);
					std::cout << ": " << wnrg << std::endl;

					tot_wnrg += wnrg;
				}
			}
		}


		std::cout << std::endl;
		std::cout << "##### TOTAL INTERFACE ENERGY: " <<  tot_wnrg << std::endl;

		return 0;

	} // try
catch ( utility::excn::EXCN_Base const & e ) {
	std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
}

} // main

