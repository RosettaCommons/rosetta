#include <basic/options/option.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/types.hh>

#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>

#include <protocols/relax/FastRelax.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <core/scoring/rms_util.hh>
#include <numeric/model_quality/rms.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentFileData.hh>

#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>

#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <stdio.h>

using namespace core;
using namespace ObjexxFCL::format;

static THREAD_LOCAL basic::Tracer TR( "app.public.comparative_modeling.avrg_silent" );

utility::vector1< Real >
calculate_variations( utility::vector1< utility::vector1< Real > > const deviation ){
	utility::vector1< Real > rmsd_vars;
	rmsd_vars.resize(deviation.size());

	for ( Size ires = 1; ires <= deviation.size(); ++ires ) {
		utility::vector1< Real > const &dev_res = deviation[ires];

		Real rmsd( 0.0 );
		for ( Size i = 1; i <= dev_res.size(); ++i ) {
			rmsd += dev_res[i]*dev_res[i];
		}
		rmsd /= (Real)(dev_res.size());
		rmsd = std::sqrt( rmsd );

		rmsd_vars[ires] = rmsd;
	}

	return rmsd_vars;
}

void
smoothen_values( utility::vector1< Real > &CAvars,
	core::Size const winsize = 9 ){

	Size const pos_shift( (winsize - 1)/2 );
	utility::vector1< Real > const CAvars0( CAvars );

	utility::vector1< Real > w( winsize );
	for ( Size i_pos = 1; i_pos <= pos_shift+1; ++i_pos ) {
		w[winsize-i_pos] = i_pos*0.1;
		w[i_pos] = i_pos*0.1;
	}

	//Clean
	CAvars.resize( CAvars0.size(), 0.0 );

	for ( Size i_res = 1; i_res <= CAvars.size(); ++i_res ) {
		Real valsum( 0.0 );
		Real wsum( 0.0 );
		for ( Size i_w = 1; i_w <= winsize; ++i_w ) {
			Size const resno( i_res+i_w-pos_shift-1);
			if ( resno < 1 || resno > CAvars.size() ) continue;

			Real wval = w[i_w]*CAvars0[resno];
			valsum += wval;
			wsum += w[i_w];
		}
		CAvars[i_res] = valsum/wsum;
	}
}

void
add_deviations( pose::Pose ref_pose,
	utility::vector1< Vector > &crds,
	utility::vector1< id::AtomID > const & atomid,
	utility::vector1< utility::vector1< Real > > &deviation
)
{
	//core::Size i( 0 );

	for ( core::Size i = 1; i <= atomid.size(); ++i ) {
		id::AtomID const &id = atomid[i];
		Vector const &xyz1 = ref_pose.xyz( id );
		Vector const &xyz2 = crds[i];
		Size ires = id.rsd();
		Size iatm = id.atomno();
		if ( ref_pose.residue( ires ).atom_name( iatm ) == " CA " ) {
			Real const dis = xyz1.distance( xyz2 );
			deviation[ires].push_back( dis );
		}
	}
}

utility::vector1< Vector >
get_avrgcrd( utility::vector1< utility::vector1< Vector > > const &crds,
	Size const n,
	utility::vector1< bool > const &exclid
)
{
	Size ncrds( 0 );
	for ( Size icrd = 1; icrd <= crds.size(); ++icrd ) {
		if ( !exclid[icrd] ) ncrds++;
	}

	utility::vector1< Vector > avrgcrd( n, Vector( 0.0 ) );
	for ( core::Size i = 1; i <= n; ++i ) {
		for ( Size icrd = 1; icrd <= crds.size(); ++icrd ) {
			if ( !exclid[icrd] ) avrgcrd[i] += crds[icrd][i];
		}
		avrgcrd[i] /= ncrds;
	}
	return avrgcrd;
}

void report( utility::vector1< Real > const &CAvars )
{
	for ( core::Size ires = 1; ires <= CAvars.size(); ++ires ) {
		TR << "CAvar: " << I(4, int(ires) )
			<< " " << F(8,3,CAvars[ires])
			<< std::endl;
	}
}

utility::vector1< Vector >
get_aligned_crd( pose::Pose pose_ref,
	pose::Pose pose,
	core::Real &rmsd,
	core::Real &Sd,
	utility::vector1< id::AtomID > const &atomid )
{
	scoring::calpha_superimpose_pose( pose, pose_ref );

	utility::vector1< Vector > crds( atomid.size() );

	//core::Size i( 0 );
	rmsd = 0.0;
	Sd = 0.0;
	//for( Size ires = 1; ires <= pose_ref.total_residue(); ++ires ){
	//  if( !pose_ref.residue( ires ).has(" CA ") ) continue;
	for ( core::Size i = 1; i <= atomid.size(); ++i ) {
		Vector const &xyz1 = pose.xyz( atomid[i] );
		Vector const &xyz0 = pose_ref.xyz( atomid[i] );
		Real d2 = xyz1.distance_squared( xyz0 );
		rmsd += d2;
		Sd += 1.0/(1.0+d2/4.0);
		crds[i] = xyz1;
	}

	rmsd = std::sqrt( rmsd/crds.size() );
	Sd = 1.0 - Sd/crds.size();
	return crds;
}

utility::vector1< id::AtomID >
get_atomindex( pose::Pose const &pose_ref )
{
	utility::vector1< id::AtomID > atomindex;
	for ( Size ires = 1; ires <= pose_ref.total_residue(); ++ ires ) {
		if ( !pose_ref.residue( ires ).has(" CA ") ) continue;
		id::AtomID id1( pose_ref.residue(ires).atom_index( "N" ), ires );
		id::AtomID id2( pose_ref.residue(ires).atom_index( "CA" ), ires );
		id::AtomID id3( pose_ref.residue(ires).atom_index( "C" ), ires );
		id::AtomID id4( pose_ref.residue(ires).atom_index( "O" ), ires );

		atomindex.push_back( id1 );
		atomindex.push_back( id2 );
		atomindex.push_back( id3 );
		atomindex.push_back( id4 );
		if ( pose_ref.residue(ires).has( "H" ) ) {
			id::AtomID id5( pose_ref.residue(ires).atom_index( "H" ), ires );
			atomindex.push_back( id5 );
		}
	}
	return atomindex;
}

pose::Pose
get_avrgpose( utility::vector1< Vector > const &avrgcrd,
	pose::Pose const &pose_ref,
	utility::vector1< Real > const &CAvar,
	utility::vector1< id::AtomID > const & atomid,
	utility::vector1< Size > & constres,
	Real const varcut )
{
	// sidechains are brought from closest pose initially
	pose::Pose pose_tmp( pose_ref ), pose_avrg( pose_ref );

	protocols::moves::MoverOP tocen
		( new protocols::simple_moves::SwitchResidueTypeSetMover( chemical::CENTROID ) );
	protocols::moves::MoverOP tofa
		( new protocols::simple_moves::SwitchResidueTypeSetMover( chemical::FA_STANDARD ) );

	// hacky way to preserve trans/rot... make a temporary pose including avrgcrd info to superimpose output pose into
	for ( core::Size i = 1; i <= atomid.size(); ++i ) {
		pose_tmp.set_xyz( atomid[i], avrgcrd[i] );
	}

	scoring::calpha_superimpose_pose( pose_avrg, pose_tmp );

	// start from pose_ref, replace less-varying crd into avrgcrd (large varying will keep pose_ref crd)
	for ( core::Size i = 1; i <= atomid.size(); ++i ) {
		id::AtomID const &id = atomid[i];
		Size ires = id.rsd();
		Size iatm = id.atomno();
		if ( CAvar[ires] < varcut ) { // replace
			pose_avrg.set_xyz( id, avrgcrd[i] );
			if ( iatm == 1 ) constres.push_back( ires );
		}
	}

	// hacky way to recover reasonable side-chains..
	tocen->apply( pose_avrg );
	tofa->apply( pose_avrg );

	TR << "List of cst residues: ";
	for ( core::Size i = 1; i <= constres.size(); ++i ) TR << " " << constres[i];
	TR << std::endl;

	return pose_avrg;
}

bool myComparison( const std::pair< Real, Size > &a, const std::pair< Real, Size > &b )
{
	return a.first < b.first;
}

Size
get_close_and_lowE(
	utility::vector1< Real > const &scores,
	utility::vector1< bool > const &excluded
)
{
	Size closest( 0 );
	//Size Ntop( 5 );
	std::vector< std::pair< Real, Size > > id_sortby_energy;
	for ( Size i = 1; i <= scores.size(); ++i ) {
		id_sortby_energy.push_back( std::make_pair( scores[i], i ) );
	}

	std::sort( id_sortby_energy.begin(), id_sortby_energy.end(), myComparison );

	// first get index sorted by score
	for ( Size i = 0; i < id_sortby_energy.size(); ++i ) {
		TR << i << " " << id_sortby_energy[i].first << " " << id_sortby_energy[i].second << std::endl;
		if ( !excluded[ id_sortby_energy[i].second ] ) {
			closest = id_sortby_energy[i].second;
			break;
		}
	}

	return closest;
}


void
relax_with_restraints_on_constres( pose::Pose &pose_avrg,
	utility::vector1< Size > const & constres,
	scoring::ScoreFunctionCOP sfxn
)
{
	using namespace core::scoring::constraints;

	for ( core::Size i = 1; i <= constres.size(); ++i ) {
		core::Size resno( constres[i] );
		core::Size iatm( pose_avrg.residue( resno ).atom_index(" CA ") );
		core::id::AtomID atomID( iatm, resno );
		core::Vector const &xyz = pose_avrg.xyz( atomID );

		core::scoring::func::FuncOP fx( new BoundFunc( 0.0, 1.0, 0.5, "" ) );
		pose_avrg.add_constraint( ConstraintCOP( ConstraintOP(
			new CoordinateConstraint( atomID, atomID, xyz, fx )))
		);
	}

	// cart2
	std::vector< std::string > cmdlines;
	cmdlines.push_back( "repack" );
	cmdlines.push_back( "accept_to_best" );
	cmdlines.push_back( "switch:cartesian" );
	cmdlines.push_back( "repeat 2" );
	cmdlines.push_back( "ramp_repack_min 0.02  0.001    1.0  50");
	cmdlines.push_back( "ramp_repack_min 0.25  0.001    0.5  50");
	cmdlines.push_back( "ramp_repack_min 0.55  0.001    0.1 100");
	cmdlines.push_back( "ramp_repack_min 1.0   0.00001  0.1 200");
	cmdlines.push_back( "accept_to_best" );
	cmdlines.push_back( "endrepeat" );

	scoring::ScoreFunctionOP sfxn_loc = sfxn->clone();
	sfxn_loc->set_weight( core::scoring::coordinate_constraint, 10.0 );
	protocols::relax::FastRelax relax( sfxn_loc );
	relax.constrain_relax_to_start_coords( true );
	relax.set_script_from_lines( cmdlines );
	relax.min_type( "lbfgs_armijo_nonmonotone" );
	relax.apply( pose_avrg );
}

bool
reavrg( utility::vector1< Vector > &avrgcrd,
	utility::vector1< utility::vector1< Vector > > const &crds,
	Real const dcut,
	bool const use_Sd,
	utility::vector1< Size > &exclid
)
{
	Size nsel( 0 );
	for ( Size icrd = 1; icrd <= crds.size(); ++icrd ) {
		Real rmsd( 0.0 ), Sd( 0.0 );
		for ( core::Size i = 1; i <= avrgcrd.size(); ++i ) {
			Real d2 = avrgcrd[i].distance_squared( crds[icrd][i] );
			rmsd += d2;
			Sd += 1.0/(1.0+d2/4.0);
		}
		rmsd = std::sqrt( rmsd/avrgcrd.size() );
		Sd = 1.0 - Sd/avrgcrd.size();
		Real d( rmsd );
		if ( use_Sd ) d = Sd;

		if ( d > dcut ) exclid[icrd] = true;
		if ( !exclid[icrd] ) nsel++;
	}

	TR << "N selected for re-averaging: " << nsel
		<< " / " << crds.size() << std::endl;

	bool valid_target( true );
	if ( nsel > 10 ) { // should be at least 10 structures to average with
		avrgcrd = get_avrgcrd( crds, avrgcrd.size(), exclid );
	} else {
		valid_target = false;
	}

	return valid_target;
}

utility::vector1< Size >
get_exclid( utility::vector1< Real > const &dists,
	Real &dcut,
	bool dcut_dynamic )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1< Size > exclid;
	if ( option[ in::file::template_pdb ].user() ) {

		if ( dcut_dynamic ) {
			//dcut = 0.3;
			// start from dcut with increment of 0.1... hard coded for Sd only for now
			while ( true ) {
				Size nsel( 0 );
				std::cout << "Dynamic dcut, applying dcut " << dcut << std::endl;
				for ( Size i = 1; i <= dists.size(); ++i ) if ( dists[i] < dcut ) nsel++;
				if ( nsel > 20 ) break;
				dcut += 0.1;
			}
		} else {
			std::cout << "Using input dcut, applying dcut " << dcut << std::endl;
		}

		// otherwise use input dcut
		for ( Size i = 1; i <= dists.size(); ++i ) {
			if ( dists[i] > dcut ) {
				exclid.push_back( true );
			} else {
				exclid.push_back( false );
			}
		}
	} else {
		exclid.resize( dists.size(), false );
	}

	return exclid;
}

void
calc_rmsf_and_avrg()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	chemical::ResidueTypeSetCOP rsd_set
		= chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	scoring::ScoreFunctionCOP sfxn
		= scoring::ScoreFunctionFactory::create_score_function( option[ score::weights] );

	core::Real dcut = option[ cm::similarity_cut ]();

	bool use_Sd( true );
	bool dcut_dynamic( true );
	bool norelax( false );
	core::Real varcut( 2.0 );

	TR << "Options read (dcut/dynamic/varcut/Sscore) " << dcut << " " << dcut_dynamic << " " << use_Sd << " " << norelax << std::endl;

	pose::Pose pose_ref, pose;
	utility::vector1< Vector > avrgcrd;
	utility::vector1< utility::vector1< Vector > > crds;
	utility::vector1< utility::vector1< Real > > deviation;
	io::silent::SilentFileOptions sopt;
	io::silent::SilentFileData sfd( sopt );
	utility::vector1< id::AtomID > atomindex;
	utility::vector1< std::string > tags;
	utility::vector1< Real > scores;
	utility::vector1< Real > dists;

	if ( option[ in::file::template_pdb ].user() ) {
		import_pose::pose_from_file( pose_ref, *rsd_set, option[ in::file::template_pdb ](1) , core::import_pose::PDB_file);
		atomindex = get_atomindex( pose_ref );
	}

	sfd.read_file( *(option[ in::file::silent ]().begin()) );

	core::Size i( 0 );
	for ( io::silent::SilentFileData::iterator iter = sfd.begin();
			iter != sfd.end(); ++iter, ++i ) {
		iter->fill_pose( pose, *rsd_set );
		if ( i == 0 && !option[ in::file::template_pdb ].user() ) { // in case refpdb is not user-defined
			pose_ref = pose;
			atomindex = get_atomindex( pose_ref );
		}
		if ( pose.total_residue() != pose_ref.total_residue() ) continue;

		core::Real rmsd( 99.9 ), Sd( 0.0 ), d;
		utility::vector1< Vector > alcrd = get_aligned_crd( pose_ref, pose, rmsd, Sd, atomindex );
		TR << "istruct/rmsd/Sd: " << i << " " << rmsd << " " << Sd << std::endl;

		if ( use_Sd ) {
			d = Sd;
		} else {
			d = rmsd;
		}

		crds.push_back( alcrd );
		tags.push_back( iter->decoy_tag() );
		scores.push_back( iter->get_energy( "score" ) );
		dists.push_back( d );
	}

	utility::vector1< Size > exclid = get_exclid( dists, dcut, dcut_dynamic );

	if ( crds.size() == 0 ) {
		TR << "not enough input structure given!" << std::endl;
		return;
	}

	// get crd deviations
	avrgcrd.resize( atomindex.size(), Vector( 0.0 ) );
	deviation.resize( pose_ref.total_residue() );

	for ( core::Size i = 1; i <= crds.size(); ++i ) {
		if ( !exclid[i] ) add_deviations( pose_ref, crds[i], atomindex, deviation );
	}

	utility::vector1< Real > CAvar = calculate_variations( deviation );
	smoothen_values( CAvar );

	report( CAvar );

	avrgcrd = get_avrgcrd( crds, atomindex.size(), exclid );

	// do average again by excluding rmsd > dcut
	bool is_valid_target = reavrg( avrgcrd, crds, dcut, use_Sd, exclid );

	// terminate if it's not a valid target.
	if ( !is_valid_target ) {
		TR << "FAIL: Not enough structure to average!" << std::endl;
		return;
	}

	// get closest structure to avrgcrd
	pose::Pose closest;
	core::Size istruct( 1 );
	istruct = get_close_and_lowE( scores, exclid );

	sfd.get_structure( tags[istruct] ).fill_pose( closest, *rsd_set );

	TR << "Closest structure to avrg is assigned as " << tags[istruct] << "; variable structures will be brought from this structure." << std::endl;

	// steal coordinate from refcrd
	utility::vector1< Size > constres;
	pose::Pose pose_avrg = get_avrgpose( avrgcrd, closest, CAvar, atomindex, constres, varcut );

	// relax
	pose_avrg.dump_pdb( option[ out::prefix ]() + ".prerelax.pdb" );
	if ( !norelax ) relax_with_restraints_on_constres( pose_avrg, constres, sfxn );

	// add_bfactor
	for ( core::Size ires = 1; ires <= CAvar.size(); ++ires ) {
		for ( core::Size iatm = 1; iatm <= pose_avrg.residue(ires).natoms(); ++iatm ) {
			pose_avrg.pdb_info()->bfactor( ires, iatm, CAvar[ires]*39.47841760435743 ); // 8*pi^2
		}
	}

	pose_avrg.dump_pdb( option[ out::prefix ]() + ".relaxed.pdb" );

	return;
}

int main( int argc, char *argv [] )
{
	try {
		devel::init(argc, argv);
		calc_rmsf_and_avrg();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
