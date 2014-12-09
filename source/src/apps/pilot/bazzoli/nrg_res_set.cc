// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Prints the energies of a set of residues in a pose
///
/// @param[in] -s <PDB>, where <PDB> is the path to the PDB file containing
/// 	the pose
///
/// @param[in] -nrg_resfile <TGT>, where <TGT> is the path to a file listing the
/// 	residues in the set. Residues are listed one per line using the following
///   format:
///
/// 	I1 C1\n
/// 	...
/// 	IN CN\n
///
///   Here, Ii and Ci indicate the residue index and chain identifier of the ith
/// 	residue in the set (i=1,...,N; N>=1).
///
/// @param[in] -extra_res_fa <PARAMS>. This parameter is required only when the
/// 	residue set contains a ligand for which a .params file is needed. In such a
/// 	 case, <PARAMS> provides the path to the .params file.
///
/// @author Andrea Bazzoli

#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/graph/Graph.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>
#include <string>
#include <iostream>
#include <fstream>

using namespace basic::options::OptionKeys;
using basic::options::option;
using core::Size;


static thread_local basic::Tracer TR( "apps.pilot.nrg_res_set" );


///
/// @brief Returns the residue number of a residue in a pose.
///
/// @parm[in] pdbnum residue number of the residue in its PDB file.
/// @parm[in] pdbchn chain identifier of the residue in the PDB file.
/// @parm[in] ps pose that the residue has been loaded into.
///
core::Size get_pose_resnum(int const pdbnum, char const pdbchn, core::pose::Pose& ps) {

	for ( Size j = 1; j <= ps.total_residue(); ++j )
		if ( ( ps.pdb_info()->chain(j) == pdbchn ) && (ps.pdb_info()->number(j) == pdbnum) )
			return j;

	// residue not found
	TR << "ERROR!! Could not find residue " << pdbnum << " and chain " << pdbchn << std::endl;
	exit(1);
}


///
/// @brief Returns true if a vector of Size contains a given element; returns false
/// 	otherwise
/// @param[in] vec: the vector
/// @param[in] tgt: the target element
///
///
bool in_set(utility::vector1<Size> const& vec, Size const tgt) {

	Size const N = vec.size();
	for(Size i=1; i<=N; i++)
		if(vec[i] == tgt)
			return true;

	return false;
}


OPT_KEY( String, nrg_resfile )

////////////////////////////////////////////////////////////////////////////////
//                                    MAIN                                    //
////////////////////////////////////////////////////////////////////////////////

int main( int argc, char * argv [] )
{

try {

	NEW_OPT( nrg_resfile, "set of residues whose energy is to be computed", "nrg_resfile.txt" );

	devel::init(argc, argv);

	// load pose from pdb file
	core::pose::Pose ps;
	std::string const input_pdb_name( basic::options::start_file() );
	core::import_pose::pose_from_pdb( ps, input_pdb_name );

	// load set of residues whose energy is to be computed. rset[i] contains the
	// ith residue in the input file (i=1,...,N, where N is the number of
	// residues in the file)
	std::string const setf = basic::options::option[nrg_resfile];
	std::ifstream setfs(setf.c_str());
	utility::vector1<Size> rset;
	int ri;
	char rc;
	while(setfs >> ri >> rc) {
		if(rc == '_')
			rc = ' ';
		rset.push_back(get_pose_resnum(ri, rc, ps));
	}

	// score pose
	core::scoring::ScoreFunctionOP scorefxn(
		core::scoring::get_score_function());

	(*scorefxn)(ps);

	//// compute energy of the set (do not count twice intra-set interactions)
	core::scoring::EnergyMap tot_u_nrgs;
	tot_u_nrgs.clear();
	core::scoring::EnergyMap weights = ps.energies().weights();;
	core::Real tot_w_nrg;

	// short-range
	TR << "##### RESIDUE SHORT-RANGE ENERGIES #####\n" << std::endl;

	core::scoring::EnergyGraph const & energy_graph( ps.energies().energy_graph() );

	Size const NSET = rset.size();
	for(Size i=1; i<=NSET; ++i) {

		Size ri = rset[i];
		TR << "### " << ps.pdb_info()->chain(ri) << "-" << ps.pdb_info()->number(ri) << " ###" << std::endl;

		for ( core::graph::Graph::EdgeListConstIter
			nit = energy_graph.get_node(ri)->const_edge_list_begin(),
			nite = energy_graph.get_node(ri)->const_edge_list_end();
			nit != nite; ++nit ) {

			Size ni = (*nit)->get_other_ind(ri);

			core::scoring::EnergyEdge const* edge( static_cast< core::scoring::EnergyEdge const*> (*nit) );
			core::scoring::EnergyMap dnrg = edge->fill_energy_map();

			if(in_set(rset, ni))
				dnrg *= 0.5;

			tot_u_nrgs += dnrg;
			TR << ps.pdb_info()->chain(ni) << "-" << ps.pdb_info()->number(ni) << ":" << dnrg.dot(weights) << std::endl;
		}

		TR << std::endl;
	}

	// long-range (doesn't seem to be active)
	using core::scoring::methods::LongRangeEnergyType;
	TR << std::endl;
	TR << "##### RESIDUE LONG-RANGE ENERGIES #####\n" << std::endl;
	for(Size lr = 1; lr <= core::scoring::methods::n_long_range_types; lr++) {

		LongRangeEnergyType lr_type = LongRangeEnergyType( lr );
		core::scoring::LREnergyContainerCOP lrec = ps.energies().long_range_container( lr_type );
		TR << "#### LR energy " << lr << " ####\n" << std::endl;

		if( !lrec || lrec->empty() )
			continue;

		for(Size i=1; i<=NSET; ++i) {

			Size ri = rset[i];
			TR << "### " << ps.pdb_info()->chain(ri) << "-" << ps.pdb_info()->number(ri) << " ###" << std::endl;

			for ( core::scoring::ResidueNeighborConstIteratorOP
				rni = lrec->const_neighbor_iterator_begin( ri );
				*rni != *( lrec->const_neighbor_iterator_end( ri ) );
				++(*rni) ) {

				Size ni = rni->neighbor_id();

				core::scoring::EnergyMap dnrg;
				rni->retrieve_energy( dnrg );

				if(in_set(rset, ni))
					dnrg *= 0.5;

				tot_u_nrgs += dnrg;
				TR << ps.pdb_info()->chain(ni) << "-" << ps.pdb_info()->number(ni) << dnrg.dot(weights) << std::endl;
			}

			TR << std::endl;
		}

		TR << std::endl;
	}

	tot_w_nrg = tot_u_nrgs.dot(weights);
	TR << std::endl;
	TR << "##### TOTAL ENERGY: " <<  tot_w_nrg << std::endl;

	return 0;

} // try
catch ( utility::excn::EXCN_Base const & e ) {
	std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
}

} // main

