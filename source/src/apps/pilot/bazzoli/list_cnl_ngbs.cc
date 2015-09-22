// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @brief Lists the residues that are neighbor to a given constellation.
///
/// @param[in] -s <PDBFIL>, where <PDBFIL> is the path to the PDB file
///  containing the protein.
///
/// @param[in] -cnl_resfile <CNLFIL>, where <CNLFIL> is the path to a file
///  enumerating the residues that form the constellation. The file has the
///  following format:
///  I1 C1\n
///  ...
///  IN CN\n ,
///  where Ii and Ci are the residue index and the chain ID, respectively, of
///  the ith residue forming the constellation (i=1,...,N).
///
/// @details The ith output line contains the identifier of the ith residue in
///  the pose that is neighbor to the constellation (i=1,...,M, where M is the
///  number of neighbors).
///
/// @details In the current implementation, a residue is neighbor to the
///  constellation iff it is neighbor to at least one residue in the
///  constellation; in particular, residue B is neighbor to constellation
///  residue A iff protocols::neighbor::in_ngbat_sphere(A, B, pose) returns
///  true.
///
/// @details A chain identifier in <CNLFIL> that is equal to '_' indicates a
///  chain identifier of ' ' in <PDBFIL>.
///
/// @author Andrea Bazzoli (bazzoli@ku.edu)


#include <protocols/neighbor/Neighborhood.hh>
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/after_opts.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <fstream>
#include <string>

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.list_cnl_ngbs" );

using namespace basic::options::OptionKeys;
using core::Size;
using core::pose::Pose;
using utility::vector1;

OPT_KEY(String, cnl_resfile)


////////////////////////////////////////////////////////////////////////////////
//                                    MAIN                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Returns the residue number of a residue in a pose.
///
/// @parm[in] pdbnum residue number of the residue in its PDB file.
/// @parm[in] pdbchn chain identifier of the residue in the PDB file.
/// @parm[in] ps pose that the residue has been loaded into.
///
Size get_pose_resnum(int const pdbnum, char const pdbchn, Pose& ps) {

	for ( Size j = 1; j <= ps.total_residue(); ++j ) {
		if ( ( ps.pdb_info()->chain(j) == pdbchn ) && (ps.pdb_info()->number(j) == pdbnum) ) {
			return j;
		}
	}

	// residue not found
	TR << "ERROR!! Could not find residue" << pdbnum << " and chain " << pdbchn << std::endl;
	exit(1);
}


/// @brief: prints the identifiers of a set of residues
///
/// @param[in]: vec indexes of the residues in the pose
/// @param[in]: ps the pose
/// @param[in]: tr output tracer
///
/// @details: the ith output line contains the identifier of the residue
///   having the ith index in vec (i=1,...,vec.size()).
///
void print_res_ids(vector1<Size> const& vec, Pose const& ps, basic::Tracer& tr) {

	for ( Size i=1; i<=vec.size(); ++i ) {
		Size ri = vec[i];
		tr << ri << '(' << ps.pdb_info()->chain(ri) << ps.pdb_info()->number(ri) <<
			ps.pdb_info()->icode(ri) << ')' << std::endl;
	}
}


int main( int argc, char * argv [] )
{

	try {

		NEW_OPT( cnl_resfile, "set of residues forming the constellation", "cnl_resfile.txt" );

		devel::init(argc, argv);

		// create pose from pdb
		core::pose::Pose ps;
		std::string const input_pdb_name( basic::options::start_file() );
		core::import_pose::pose_from_pdb( ps, input_pdb_name );

		// load constellation: cnl[i] contains the pose index of the ith residue in
		// the constellation input file (i=1,...,N, where N is the number of residues
		// in the file)
		std::string const cnlf = basic::options::option[cnl_resfile];
		std::ifstream cnlfs(cnlf.c_str());
		vector1<Size> cnl;
		int ri;
		char rc;
		while ( cnlfs >> ri >> rc ) {
			if ( rc == '_' ) {
				rc = ' ';
			}
			cnl.push_back(get_pose_resnum(ri, rc, ps));
		}

		// find neighbors and list them on screen
		protocols::neighbor::Neighborhood n(cnl, ps,
			protocols::neighbor::in_ngbat_sphere);

		print_res_ids(n.get(), ps, TR);

	}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
}
