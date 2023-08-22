

// core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <devel/init.hh>
#include <basic/options/option_macros.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <basic/Tracer.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/conformation/Residue.hh>



// utility headers
#include <utility/vector1.hh>
#include <utility/file/FileName.fwd.hh>

// c++ headers
#include <string>

#include <protocols/jd2/JobDistributor.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/symmetry/SetupForSymmetryMover.hh>
#include <protocols/ligand_docking/GALigandDock/GridHash3D.hh>
#include <protocols/ligand_docking/GALigandDock/GALigandDock.hh>
#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>
#include <protocols/ligand_docking/GALigandDock/util.hh>



// namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace hbonds;
using namespace pose;
using namespace ObjexxFCL;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

OPT_1GRP_KEY(IntegerVector, settings, hb_res_list )


static basic::Tracer TR("apps.ligand_docking_hbonds");

void compute(core::pose::Pose const& pose, core::scoring::hbonds::HBondSet const& hbset,
	utility::vector1<core::Size> const& ligids,
	utility::vector1<core::Size> const& resids, core::Size & n_polars,
	core::Size & n_unsats, bool const& include_bb, std::string const& hb_metric
){
	using namespace protocols::ligand_docking::ga_ligand_dock;
	utility::vector1<hbDon> lig_hbDons, rec_hbDons;
	utility::vector1<hbAcc> lig_hbAccs, rec_hbAccs;


	for ( auto ligid : ligids ) {
		core::conformation::Residue const& ires = pose.residue(ligid);
		// ligand acceptors
		for ( auto anum=ires.accpt_pos().begin(), anume=ires.accpt_pos().end(); anum!=anume; ++anum ) {
			hbAcc acc;
			core::Size bnum = ires.atom_base( *anum );
			core::Size b0num = ires.abase2( *anum );
			acc.A =  ires.xyz( *anum );
			acc.B = ires.xyz( bnum );
			acc.B_0 = ires.xyz( b0num );
			acc.acctype = core::scoring::hbonds::get_hb_acc_chem_type( *anum, ires );
			lig_hbAccs.push_back(acc);
		}
		// ligand donors
		for ( auto hnum=ires.Hpos_polar().begin(), hnume=ires.Hpos_polar().end(); hnum!=hnume; ++hnum ) {
			hbDon don;
			core::Size dnum = ires.atom_base( *hnum );
			don.H =  ires.xyz( *hnum );
			don.D = ires.xyz( dnum );
			don.dontype = core::scoring::hbonds::get_hb_don_chem_type( dnum, ires );
			lig_hbDons.push_back(don);
		}
	}

	for ( auto resid : resids ) {
		core::conformation::Residue const& ires = pose.residue(resid);
		// receptor acceptors
		for ( auto anum=ires.accpt_pos().begin(), anume=ires.accpt_pos().end(); anum!=anume; ++anum ) {
			if ( !include_bb && !(*anum>ires.last_backbone_atom() && *anum<=ires.nheavyatoms()) ) continue;
			core::id::AtomID aid( *anum, resid );
			TR <<"atomid: " << aid <<  " hbonds: " << hbset.nhbonds(aid) << std::endl;
			if ( hbset.nhbonds(aid)>0 ) continue;
			hbAcc acc;
			core::Size bnum = ires.atom_base( *anum );
			core::Size b0num = ires.abase2( *anum );
			acc.A =  ires.xyz( *anum );
			acc.B = ires.xyz( bnum );
			acc.B_0 = ires.xyz( b0num );
			acc.acctype = core::scoring::hbonds::get_hb_acc_chem_type( *anum, ires );
			rec_hbAccs.push_back(acc);
		}
		// receptor donors
		for ( auto hnum=ires.Hpos_polar().begin(), hnume=ires.Hpos_polar().end(); hnum!=hnume; ++hnum ) {
			if ( !include_bb && !(*hnum>ires.first_sidechain_hydrogen()) ) continue;
			core::id::AtomID aid( *hnum, resid );
			TR <<"atomid: " << aid <<  " hbonds: " << hbset.nhbonds(aid) << std::endl;
			if ( hbset.nhbonds(aid)>0 ) continue;
			hbDon don;
			core::Size dnum = ires.atom_base( *hnum );
			don.H =  ires.xyz( *hnum );
			don.D = ires.xyz( dnum );
			don.dontype = core::scoring::hbonds::get_hb_don_chem_type( dnum, ires );
			rec_hbDons.push_back(don);
		}
	}
	n_polars = rec_hbAccs.size() + rec_hbDons.size();
	n_unsats = 0;
	core::Real maxHbdis = 4.2;
	core::Real maxHbdis2 = maxHbdis*maxHbdis;
	core::Real hb_energy_cutoff = -0.3;

	core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();
	core::scoring::hbonds::HBondOptions const & hbopt = sf->energy_method_options().hbond_options();
	core::scoring::hbonds::HBondDatabaseCOP hb_database = core::scoring::hbonds::HBondDatabase::get_database( hbopt.params_database_tag() );
	//compute the number of unsatisfied hbonds of receptor donors
	bool is_satisfied(false);
	for ( auto rec_hbDon : rec_hbDons ) {
		for ( auto lig_hbAcc : lig_hbAccs ) {
			is_satisfied = is_hb_satisfied(sf, hb_database, hbopt, lig_hbAcc, rec_hbDon, maxHbdis2, hb_energy_cutoff, hb_metric );
			if ( is_satisfied ) break;
		}
		if ( !is_satisfied ) n_unsats++;
	}
	//compute the number of unsatisfied hbonds of receptor acceptors
	for ( auto rec_hbAcc : rec_hbAccs ) {
		for ( auto lig_hbDon : lig_hbDons ) {
			is_satisfied = is_hb_satisfied(sf, hb_database, hbopt, rec_hbAcc, lig_hbDon, maxHbdis2, hb_energy_cutoff, hb_metric );
			if ( is_satisfied ) break;
		}
		if ( !is_satisfied ) n_unsats++;
	}
}


class LigandDockingHbondsReporter : public protocols::moves::Mover {
public:
	LigandDockingHbondsReporter()= default;

	void apply( core::pose::Pose & pose) override {

		utility::vector1<core::Size> res_list;
		if ( option[settings::hb_res_list].user() ) {
			TR << "The fixed residues are:" << std::endl;
			res_list = option[settings::hb_res_list]();
		}



		bool include_bb(false);
		std::string hb_metric("default");
		utility::vector1<core::Size> ligids;
		protocols::ligand_docking::ga_ligand_dock::get_ligand_resids( pose, ligids );

		core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();
		core::pose::Pose receptor( pose );
		core::Size startid, endid;
		startid = ( ligids[1] <= ligids.back())? ligids[1] : ligids.back() ;
		endid = ( ligids[1] > ligids.back())? ligids[1] : ligids.back() ;
		receptor.delete_residue_range_slow(startid, endid);

		receptor.update_residue_neighbors();
		sf->score(receptor);
		core::scoring::hbonds::HBondSet set1;
		fill_hbond_set( receptor, false, set1, false );
		core::Size n_polars(0);
		core::Size n_unsats(0);
		compute( pose, set1, ligids, res_list, n_polars, n_unsats, include_bb, hb_metric );
		core::pose::setPoseExtraScore( pose, "n_polars", n_polars);
		core::pose::setPoseExtraScore( pose, "n_unsats", n_unsats);

		core::Size n_hbonds_total(0), n_hbonds_max1(0);
		protocols::ligand_docking::ga_ligand_dock::compute_nhbonds( pose, ligids, res_list, n_hbonds_total, n_hbonds_max1, include_bb, hb_metric );
		core::pose::setPoseExtraScore( pose, "n_hbonds_total", n_hbonds_total);
		core::pose::setPoseExtraScore( pose, "n_hbonds_max1", n_hbonds_max1);

	}

	std::string get_name() const override {
		return "LigandDockingHbondsReporter";
	}
};

int
main( int argc, char * argv [] ) {
	using namespace protocols::moves;
	using namespace protocols;
	using namespace protocols::jd2;
	using namespace protocols::symmetry;

	try {
		NEW_OPT(settings::hb_res_list, "receptor residue list", 0);

		devel::init(argc, argv);

		SequenceMoverOP seq( new SequenceMover() );

		seq->add_mover( utility::pointer::make_shared< LigandDockingHbondsReporter >() );

		protocols::jd2::JobDistributor::get_instance()->go( seq );
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
