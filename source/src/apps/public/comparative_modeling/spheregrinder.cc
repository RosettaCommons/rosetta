#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <core/scoring/constraints/util.hh>

#include <core/id/AtomID.hh>
#include <core/scoring/rms_util.hh>
#include <numeric/model_quality/rms.hh>
#include <ObjexxFCL/format.hh>
#include <utility/excn/Exceptions.hh>

#include <core/types.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace ObjexxFCL::format;

static THREAD_LOCAL basic::Tracer TR( "app.public.comparative_modeling.spheregrinder" );

utility::vector1< Vector >
get_sphere_xyz( pose::Pose const &pose,
	utility::vector1< id::AtomID > const &sphere_atomid
)
{
	utility::vector1< core::Size > resnos;

	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		Size resno = pose.pdb_info()->number(i);
		resnos.push_back( resno );
		//TR.Debug << i << "-> " << resno << std::endl;
	}

	utility::vector1< Vector > xyz( sphere_atomid.size() );
	for ( Size i = 1; i <= sphere_atomid.size(); ++i ) {
		id::AtomID atomid = sphere_atomid[i];
		// pdbresno
		//Size pdbno = atomid.rsd();
		//TR.Debug << i << " " << pdbno << " " << atomid.atomno() << " " << resnos.contains(pdbno) << std::endl;

		if ( !resnos.contains( atomid.rsd() ) ) continue;
		Size resno2 = resnos.index( atomid.rsd() );

		//TR.Debug << "ID: " << resno2 << " " << atomid.atomno() << std::endl;
		xyz[i] = pose.residue( resno2 ).xyz( atomid.atomno() );
		//if( atomid.atomno() == 1 ){
		//TR.Debug << "crd: " << resno2 << " " << atomid.atomno() << " " << xyz[i][0] << std::endl;
		//}
	}
	return xyz;
}

Real
SphereGrinder( pose::Pose const &pose_ref,
	pose::Pose const &pose,
	Real dcut2,
	std::string mode )
{
	// get sphere residue first
	utility::vector1< utility::vector1< id::AtomID > > sphere_atomid( pose_ref.total_residue() );
	utility::vector1< utility::vector1< id::AtomID > > sphere_atomid_ref( pose_ref.total_residue() );

	for ( Size ires = 1; ires <= pose_ref.total_residue(); ++ires ) {
		if ( !pose_ref.residue( ires ).has(" CA ") ) continue;
		Size iresno = pose_ref.pdb_info()->number(ires);

		if ( mode == "all" ) {
			for ( Size iatm = 1; iatm <= pose_ref.residue(ires).nheavyatoms(); ++iatm ) {
				id::AtomID id1 = id::AtomID(iatm,iresno);
				sphere_atomid[ires].push_back( id1 );
			}
		} else { // TODO
			sphere_atomid[ires].push_back( id::AtomID(pose_ref.residue(ires).atom_index( "CA" ), iresno ));
		}

		for ( Size jres = ires+1; jres <= pose_ref.total_residue(); ++jres ) {
			if ( !pose_ref.residue( jres ).has(" CA ") ) continue;
			Size jresno = pose_ref.pdb_info()->number(jres);

			if ( mode == "all" ) {
				for ( Size iatm = 1; iatm <= pose_ref.residue(ires).nheavyatoms(); ++iatm ) {
					for ( Size jatm = 1; jatm <= pose_ref.residue(jres).nheavyatoms(); ++jatm ) {

						id::AtomID id1 = id::AtomID(iatm,iresno);
						id::AtomID id2 = id::AtomID(jatm,jresno);
						if ( sphere_atomid[ires].contains( id2 ) ) continue;

						Vector const &xyz1 = pose_ref.residue( ires ).xyz(iatm);
						Vector const &xyz2 = pose_ref.residue( jres ).xyz(jatm);
						Real dis2 = xyz1.distance_squared( xyz2 );
						if ( dis2 < dcut2 ) {
							//sphere_atomid_ref[ires].push_back( id::AtomID(jatm,jres) );
							//sphere_atomid_ref[jres].push_back( id::AtomID(iatm,ires) );
							sphere_atomid[ires].push_back( id2 );
							sphere_atomid[jres].push_back( id1 );
						}
					}
				}

			} else if ( mode == "CA" ) {
				Vector const &xyz1 = pose_ref.residue( ires ).xyz("CA");
				Vector const &xyz2 = pose_ref.residue( jres ).xyz("CA");
				Real dis2 = xyz1.distance_squared( xyz2 );
				if ( dis2 < dcut2 ) {
					//sphere_atomid_ref[ires].push_back( id::AtomID(pose_ref.residue(jres).atom_index( "CA" ), jres ));
					//sphere_atomid_ref[jres].push_back( id::AtomID(pose_ref.residue(ires).atom_index( "CA" ), ires ));
					sphere_atomid[ires].push_back( id::AtomID(pose_ref.residue(jres).atom_index( "CA" ), jresno ));
					sphere_atomid[jres].push_back( id::AtomID(pose_ref.residue(ires).atom_index( "CA" ), iresno ));
				}
			}
		}
	}

	Size n2( 0 ), n4( 0 ), ncount( 0 );
	for ( Size ires = 1; ires <= pose_ref.total_residue(); ++ires ) {
		if ( sphere_atomid[ires].size() == 0 ) continue;
		ncount++;
		//TR.Debug << "look for" << ires << std::endl;
		//utility::vector1< Vector > xyz_ref = get_sphere_xyz( pose_ref, sphere_atomid_ref[ires] );
		utility::vector1< Vector > xyz_ref = get_sphere_xyz( pose_ref, sphere_atomid[ires] );
		utility::vector1< Vector > xyz = get_sphere_xyz( pose, sphere_atomid[ires] );
		Real sphere_rmsd = numeric::model_quality::calc_rms( xyz_ref, xyz );
		TR.Debug << "res/nsphere/rmsd/SGres" << " " << pose_ref.pdb_info()->number(ires) << " " << sphere_atomid[ires].size() << " " << sphere_rmsd << " "
			<< (sphere_rmsd<2.0?1:0)+(sphere_rmsd<4.0?1:0) << std::endl;
		if ( sphere_rmsd <= 2.0 ) n2++;
		if ( sphere_rmsd <= 4.0 ) n4++;
	}

	//std::cout << "n2/n4/ncount: " << n2 << " " << n4 << " "<< ncount << std::endl;
	return Real(n2+n4)*0.5/Real(ncount);
}

int main( int argc, char *argv [] )
{
	try {
		devel::init( argc, argv );

		chemical::ResidueTypeSetCOP rsd_set =
			chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

		scoring::ScoreFunctionCOP sfxn
			= scoring::ScoreFunctionFactory::create_score_function( option[ score::weights] );

		//Real dcut = option[ macro::dcut ]();
		Real dcut = 6.0;
		bool do_score = option[ score::weights ].user();

		protocols::simple_moves::ConstraintSetMoverOP loadCsts( new protocols::simple_moves::ConstraintSetMover );

		if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
			loadCsts->constraint_file( core::scoring::constraints::get_cst_fa_file_option() );
		}

		pose::Pose ref_pose;
		import_pose::pose_from_file( ref_pose, *rsd_set, option[ in::file::native ](), import_pose::PDB_file);
		std::cout << "SUMMARY: Tag   score   SG    cst  cart_bonded" << std::endl;

		if ( option[ in::file::silent ].user() ) {
			core::io::silent::SilentFileOptions sopt;
			core::io::silent::SilentFileData sfd( sopt );
			sfd.read_file( option[ in::file::silent ](1) );

			for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(); iter != sfd.end(); ++iter ) {
				pose::Pose pose;
				iter->fill_pose( pose, *rsd_set );

				std::string tag = iter->decoy_tag();
				TR.Debug << "Evaluating " << tag <<std::endl;
				Real SG = SphereGrinder( ref_pose, pose, dcut*dcut, "all" );

				Real score, cst, cart;
				if ( do_score ) {
					loadCsts->apply(pose);
					score = sfxn->score( pose );
					scoring::EnergyMap const emap = pose.energies().total_energies();
					cst = emap[ scoring::atom_pair_constraint ];
					cart = emap[ scoring::cart_bonded ];

				} else {
					score = iter->get_energy( "score" );
					cst = iter->get_energy( "atom_pair_constraint" );
					cart = iter->get_energy( "cart_bonded" );
				}

				std::cout << "SUMMARY: " << tag << " " << score << " " << SG <<  " " << cst << " " << " " << cart << std::endl;
			}

		} else {
			core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();
			while ( input.has_another_pose() ) {
				pose::Pose pose;
				input.fill_pose( pose, *rsd_set );
				if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) loadCsts->apply(pose);
				Real score = sfxn->score( pose );
				scoring::EnergyMap const emap = pose.energies().total_energies();
				Real cst = emap[ scoring::atom_pair_constraint ];

				Real SG = SphereGrinder( ref_pose, pose, dcut*dcut, "all" );

				std::cout << pose.pdb_info()->name() << " " << score << " " << SG << " " << cst << std::endl;
				//std::cout << tags[i] << " " << score << " " <<SG << std::endl;
			}
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
