// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <core/conformation/Conformation.hh>
#include <core/scoring/rms_util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <numeric/angle.functions.hh>
#include <numeric/conversions.hh>

#include <basic/options/option.hh>

#include <ObjexxFCL/string.functions.hh>
#include <basic/options/option_macros.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <sstream>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <basic/Tracer.hh>
#include <devel/init.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/format.hh>

using namespace std;

OPT_KEY( IntegerVector, decoys_frag_len )

int main( int argc, char** argv ) {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( decoys_frag_len, "length of fragment to be test", utility::vector1<Size>() );

	try {
		//init
		devel::init( argc, argv );

		//decoys
		utility::vector1<pose::PoseOP> models = import_pose::poseOPs_from_files( option[OptionKeys::in::file::s](), core::import_pose::PDB_file);
		//native
		pose::Pose native;
		import_pose::pose_from_file( native, option[ in::file::native ]() , core::import_pose::PDB_file);
		util::switch_to_residue_type_set( native, chemical::CENTROID );

		Size N_frag( option[ decoys_frag_len ]().size() );
		for (Size n=1; n<=N_frag; n++) {
				//for each frag length
				Size len_pose(native.size());
				Size len_frag = option[ decoys_frag_len ]()[n];
				std::stringstream fn;
				utility::file::FileName pdbfn(option[ in::file::native ]());
				fn << pdbfn.base() << "." << len_frag << "mers.frag.rms";
				utility::io::ozstream outfn( fn.str() );

				for (Size fi=1; fi <= len_pose-len_frag+1; ++fi) {
					// for each position
					for (Size pi=1; pi <= models.size(); ++pi) {
						// for each decoys
						Real rmsd = scoring::CA_rmsd( native, *(models[pi]), fi, fi+len_frag-1 );
						outfn << len_frag << "\t" << fi << "\t" << pi << "\t" << rmsd << "\t0\t0" << endl;
					}
					outfn << endl;
				}
				outfn.close();
		}
	}
	catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

