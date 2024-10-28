#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/import_pose.hh>

#include <core/chemical/Atom.hh>
#include <core/chemical/AtomType.hh> 
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>

#include <protocols/electron_density/util.hh>
#include <protocols/electron_density/DensitySymmInfo.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>

#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/electron_density/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/random/random.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>

#include <basic/Tracer.hh>
#include <cmath>
#include <limits>
#include <numeric>

#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

//for debugging coordinates
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>

static basic::Tracer TR("amuenks.density_correlation_levels");

using namespace basic::options;
using namespace basic::options::OptionKeys;

int
main( int argc, char * argv[] ){
    try {
        
        devel::init(argc, argv);

        core::pose::Pose p;
        utility::vector1< utility::file::FileName > fname;
        fname = basic::options::option[in::file::s]();
        core::import_pose::pose_from_file(p, fname[1], core::import_pose::PDB_file);
        core::Size nres = p.size();

        core::scoring::electron_density::ElectronDensity density;
        density.readMRCandResize( basic::options::option[edensity::mapfile](),
                                  basic::options::option[edensity::mapreso](),
                                  basic::options::option[edensity::grid_spacing]() );
       
	protocols::electron_density::SetupForDensityScoringMoverOP dockindens( new protocols::electron_density::SetupForDensityScoringMover );
        dockindens->apply( p );

	core::scoring::electron_density::getDensityMap().set_nres( nres );
	core::scoring::electron_density::getDensityMap().setScoreWindowContext( true );
	core::scoring::electron_density::getDensityMap().setWindow( 1.0 );

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	scorefxn->set_weight( core::scoring::elec_dens_window, 1.0 );
	(*scorefxn)(p);

	core::conformation::Residue const &lig( p.residue( nres ) );
	
	TR << "Ligand is " << lig.name() << std::endl;
	numeric::xyzVector< core::Real > lig_com(0,0,0);
	core::Size nAtms = 0;
	
	//p.dump_pdb("correlation_app_dump.pdb");	
	core::Real pose_cc = core::scoring::electron_density::getDensityMap().matchPose( p );

	TR << "Pose density correlation: " << pose_cc << std::endl;

	std::string name = lig.name();

	for ( core::Size j=1; j<=lig.natoms(); ++j ) {
		lig_com += lig.xyz(j);
	}
	nAtms = lig.natoms();

	lig_com /= nAtms;

	core::Real pocket_cc = core::scoring::electron_density::getDensityMap().matchPoseInPocket( p, lig_com, 10 );
	TR << "Pocket density correlation " << name << " " << nres << ": " << pocket_cc << std::endl;
	core::Real lig_cc = core::scoring::electron_density::getDensityMap().matchRes( lig.seqpos(), lig, p, nullptr, false );
	TR << "Ligand density correlation " << name << " " << nres << ": " << lig_cc << std::endl;
	core::Real lig_cc_no_fg_matchres = core::scoring::electron_density::getDensityMap().matchRes( lig.seqpos(), lig, p, nullptr, false, false, true );
        TR << "Ligand density correlation no fg " << name << " " << nres << ": " << lig_cc_no_fg_matchres << std::endl;
	p.delete_residue_slow(nres);
	core::Real pocket_cc_no_lig = core::scoring::electron_density::getDensityMap().matchPoseInPocket( p, lig_com, 10 );
        TR << "Pocket density correlation no ligand " << name << " " << nres << ": " << pocket_cc_no_lig << std::endl;

	/*for ( core::Size i = 1; i <= p.size(); ++i ) {	
		if ( p.residue(i).is_ligand() ) {
			if ( !p.residue(i).is_metal() && !p.residue(i).is_virtual_residue() ) {
				core::conformation::Residue const &lig( p.residue( i ) );
				std::string name = lig.name();

				for ( core::Size j=1; j<=lig.natoms(); ++j ) {
					lig_com += lig.xyz(j);
				}
				nAtms = lig.natoms();
			       
				lig_com /= nAtms;

				core::Real pocket_cc = core::scoring::electron_density::getDensityMap().matchPoseInPocket( p, lig_com, 10 );
				TR << "Pocket density correlation " << name << " " << i << ": " << pocket_cc << std::endl;
				core::Real lig_cc = core::scoring::electron_density::getDensityMap().matchRes( lig.seqpos(), lig, p, nullptr, false );
				TR << "Ligand density correlation " << name << " " << i << ": " << lig_cc << std::endl;

				pose.delete_residue_slow(p.residue(i), p.residue(i));
			}
		}
	}*/

    } catch ( utility::excn::Exception const & e ){
        e.display();
        return -1;
    }
    return 0;
}
