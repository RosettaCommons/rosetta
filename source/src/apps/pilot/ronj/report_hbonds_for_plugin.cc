// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/ronj/report_hbonds_for_plugin.cc
/// @brief  report all data about hydrogen bonds, formatted in a special way for a PyMOL plugin
/// @author Matthew O'Meara
/// @author Ron Jacak

#include <core/chemical/AtomType.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/option.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/types.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/Tracer.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <devel/init.hh>

// utility headers
#include <utility/io/izstream.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/file/FileName.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <string>
#include <sstream>
#include <boost/algorithm/string.hpp>  // for string split

// option keys includes

//Auto Headers
#include <protocols/moves/Mover.hh>
#include <utility/vector1.hh>


// namespaces
using namespace core;
using namespace scoring;
using namespace ObjexxFCL::format;

using utility::vector1;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.ronj.report_hbonds_for_plugin" );


// The following macros expand to the following.  They basically add these options to the Options system.
// #define OPT_1GRP_KEY( type, grp, key )
// namespace core { namespace options { namespace OptionKeys { namespace grp { basic::options::type##OptionKey const key( #grp":"#key ); }}}}
//OPT_1GRP_KEY( Boolean, HBondReporter, relax )
OPT_1GRP_KEY( String, HBondReporter, relevant_chains )


class HBondReporter : public protocols::moves::Mover {

public:


	HBondReporter() : allowNonProtein_( true ), relevant_chains_("*") {
		scfxn = core::scoring::ScoreFunctionOP( new ScoreFunction );
		scfxn->set_weight( core::scoring::hbond_lr_bb, 1.17 );
		scfxn->set_weight( core::scoring::hbond_sr_bb, 0.585 );
		scfxn->set_weight( core::scoring::hbond_bb_sc, 1.17 );
		scfxn->set_weight( core::scoring::hbond_sc, 1.1 );
	}


	virtual ~HBondReporter(){};


	virtual
	std::string get_name() const {
		return "HBondReporter";
	}


	/// @brief
	/// What is this for?
	///
	void load_job_data( Pose & pose ){

		std::vector< std::string > lines, tokens;
		std::stringstream chains;

		// the following is true if, say, the structure came from a silent file
		if ( !pose.pdb_info() ) { return; }

		if ( basic::options::option[ basic::options::OptionKeys::HBondReporter::relevant_chains ].user() ) {

			utility::io::izstream relevant_chains_file( basic::options::option[ basic::options::OptionKeys::HBondReporter::relevant_chains ]() );
			if ( !relevant_chains_file ) {
				TR << "No job data file to open. " << basic::options::option[ basic::options::OptionKeys::HBondReporter::relevant_chains ]() << std::endl;
				return;
			}

			utility::file::FileName pose_filename = pose.pdb_info()->name();
			pose_filename.to_local_name();

			std::string line;
			while ( getline( relevant_chains_file, line ) ) {
				//TR << "For line '" << line << "'" << std::endl;
				boost::algorithm::split( tokens, line, boost::algorithm::is_any_of("\t") );
				//TR << "\t tokens[0]="<< tokens[0] << " tokens[1]=" <<tokens[1] <<std::endl;
				//TR << "\t local pdb_filename " << pose_filename << std::endl;
				if ( tokens[0] == pose_filename ) {
					relevant_chains_ = tokens[1];
					TR << "Restricting pdb '" << pose.pdb_info()->name() << "' to chains " << relevant_chains_ << std::endl;
					//frag3_file = tokens[2];
					return;
				}
			}
		}

		TR << "Job data for PDB '"<< pose.pdb_info()->name() << "' not found!" << std::endl;
	}


	/// @brief
	/// Returns a string representing the name of this pose.  Uses the pdb_info object or the job distributor tag if
	/// no pdb_info object exists.
	///
	std::string pose_name( pose::Pose& pose ) {
		std::string name = "No_Name_Found";
		if ( pose.pdb_info() ) {
			name = pose.pdb_info()->name();
		} else if ( pose.data().has( pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
			name = static_cast< basic::datacache::CacheableString const & >( pose.data().get( pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ).str();
		} else {
			name = protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag();
		}
		return name;
	}


	/// @brief
	/// Runs a "simple multi relax" protocol on the whole pose, or if the relevant chains flag was used
	/// on just the relevant chains. But only if the relax flag was specified.
	/// I want this executable/mover to only print Hbond info.  If a user wants to relax their structure or do st else
	/// to it, they need to run that protocol first and then use this one to get Hbond info.
	///
	/*
	void decoyify( pose::Pose& pose ) {

	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	mm->clear();

	if ( relevant_chains_ == "*" || relevant_chains_.empty() ) {
	TR << "Information about relevant chains not found for " << pose_name( pose ) << "!" << std::endl;
	mm->set_bb( true );
	mm->set_chi( true );
	} else if ( !relevant_chains_.empty() ) {
	for ( Size jj=1; jj <= pose.size(); ++jj ) {
	if ( relevant_chains_.find( pose.pdb_info()->chain( jj ), 0 ) != std::string::npos ) {
	mm->set_bb( jj, true );
	mm->set_chi( jj, true );
	}
	}
	}

	protocols::relax::FastRelax relax( scfxn );
	TR << "Initial score: " << (*scfxn)(pose) << std::endl;
	relax.apply(pose);
	TR << "Decoy score: " << (*scfxn)(pose) << std::endl;
	}
	*/


	/// @brief
	/// Given a pose and a residue position, returns the bfactor for that residue.
	/// Why do we care about bfactors?
	///
	Real bfactor( pose::Pose& pose, Size resNum, Size atm, bool is_backbone ) {
		if ( !pose.pdb_info() ) { return 9999; }
		if ( is_backbone ) {
			return pose.pdb_info()->temperature( resNum, atm );
		} else {
			core::Real max_temp = -1000;
			for ( Size ii=1; ii <= pose.residue( resNum ).natoms(); ii++ ) {
				max_temp = std::max( max_temp, pose.pdb_info()->temperature( resNum, ii ) );
			}
			return max_temp;
		}
	}


	/// @brief
	///
	///
	///
	virtual void
	apply( pose::Pose& pose ) {

		TR << "getting hbonds for " << pose_name(pose) << std::endl;

		//load_job_data( pose );

		//if ( basic::options::option[ basic::options::OptionKeys::HBondReporter::relax ].value() ) { decoyify( pose ); }

		scoring::hbonds::HBondSet set1;
		(*scfxn)(pose);
		pose.update_residue_neighbors();
		set1.setup_for_residue_pair_energies( pose, false, false );

		// compute sasa pack
		core::Real const probe_radius(1.4);
		id::AtomID_Map< core::Real > atom_sasa;
		utility::vector1< core::Real > residue_sasa;
		scoring::calc_per_atom_sasa( pose, atom_sasa, residue_sasa, probe_radius );

		// compute dssp
		core::scoring::dssp::Dssp dssp( pose );
		dssp.insert_ss_into_pose( pose );

		core::Real hbond_energies = 0; // to get average energy for structure;
		Size nhbonds_used = 0;

		std::cout << "don_chain don_resname don_resnum don_pdbresnum don_atomname - ";
		std::cout << "acc_chain acc_resname acc_resnum acc_pdbresnum acc_atomname ";
		std::cout << "energy weight ";
		std::cout << "don_nb don_sasa don_bfactor donSS - acc_nb acc_sasa acc_bfactor accSS ";
		std::cout << std::endl;

		for ( Size i=1; i <= set1.nhbonds(); ++i ) {

			scoring::hbonds::HBond bond = set1.hbond( i );

			// need to access donor and acc residues as well as ints
			Size hatm = bond.don_hatm();
			//scoring::hbonds::HBEvalType type = bond.eval_type();  // unused ~Labonte
			Size don_resnum = bond.don_res(); Size acc_resnum = bond.acc_res();

			// get acc and donor residues from sequence numbers
			conformation::Residue accRes = pose.residue( acc_resnum ); std::string acc_resname = accRes.name3();
			conformation::Residue donRes = pose.residue( don_resnum ); std::string don_resname = donRes.name3();

			Size const datm( donRes.atom_base( hatm ) );
			//Vector const & hatm_xyz( donRes.atom( hatm ).xyz() );  // unused ~Labonte
			//Vector const & datm_xyz( donRes.atom( datm ).xyz() );  // unused ~Labonte

			Size aatm = bond.acc_atm();
			//Size const base( accRes.atom_base( aatm ) );
			//Size const base_of_base( accRes.atom_base( base ) );

			//Size const base2( accRes.abase2( aatm ) );
			//Size const dbase( donRes.atom_base( datm ) );
			//Vector const & aatm_xyz = accRes.atom( aatm ).xyz();  // unused ~Labonte
			//Vector const & base_xyz = accRes.atom( base ).xyz();  // unused ~Labonte
			//Vector const & base_of_base_xyz = accRes.atom( base_of_base ).xyz();  // unused ~Labonte

			//Vector const & base2_xyz = accRes.atom( base2 ).xyz();  // unused ~Labonte
			//Vector const & dbase_xyz = donRes.atom( dbase ).xyz();  // unused ~Labonte

			//Size donType = donRes.atom( datm ).type();
			//Size accType = accRes.atom( aatm ).type();

			bool isProtein = true;
			if ( !donRes.is_protein() || !accRes.is_protein() ) { isProtein = false; }
			if ( !( (set1.allow_hbond(i) && isProtein) || (set1.allow_hbond(i) && allowNonProtein_) ) ) {
				continue;
			}

			//const std::string donElement = donAtomType.element();
			//const std::string accElement = accAtomType.element();
			const std::string & don_atomname = donRes.atom_name( datm );
			const std::string & acc_atomname = accRes.atom_name( aatm );
			//const std::string & h_atomname = donRes.atom_name( hatm );
			// usual atomic numbers from periodic table
			//int donElemNum, accElemNum;  // unused ~Labonte

			char don_chain, acc_chain;
			int don_pdbresnum, acc_pdbresnum;
			char don_pdbicode, acc_pdbicode;

			if ( pose.pdb_info() ) {
				don_chain = pose.pdb_info()->chain( don_resnum ); acc_chain = pose.pdb_info()->chain( acc_resnum );
				don_pdbresnum = pose.pdb_info()->number( don_resnum ); acc_pdbresnum = pose.pdb_info()->number( acc_resnum );
				don_pdbicode = pose.pdb_info()->icode( don_resnum ); acc_pdbicode = pose.pdb_info()->icode( acc_resnum );
			} else {
				don_chain = 'A'; acc_chain = 'A';
				don_pdbresnum = don_resnum; acc_pdbresnum = acc_resnum;
				don_pdbicode = ' '; acc_pdbicode = ' ';
			}

			if ( relevant_chains_ != "*" ) {
				if ( ( relevant_chains_.find( don_chain, 0 ) == std::string::npos ) || (relevant_chains_.find( acc_chain, 0 ) == std::string::npos ) ) {
					continue;
				}
			}

			Size const donNbrs( set1.nbrs( don_resnum ) ); Size const accNbrs( set1.nbrs( acc_resnum ) );
			Real const donSasa( atom_sasa( don_resnum, datm ) ); Real const accSasa( atom_sasa( acc_resnum, aatm ) );

			Real donBfactor, accBfactor;
			if ( ! pose.pdb_info() ) {
				donBfactor = 9999; accBfactor = 9999;
			} else {
				donBfactor = pose.pdb_info()->temperature( don_resnum, datm ); accBfactor = pose.pdb_info()->temperature( acc_resnum, aatm );
			}

			char const & donSS = pose.secstruct( don_resnum ); char const & accSS = pose.secstruct( acc_resnum );

			//scoring::hbonds::Deriv deriv = bond.deriv();
			core::Real energy = bond.energy();
			//scoring::hbonds::hb_energy_deriv( type, datm_xyz, hatm_xyz, aatm_xyz, base_xyz, base2_xyz, energy, true, deriv );
			core::Real weight = bond.weight();

			hbond_energies += energy;

			// converts the atom element type (O, N, S, etc) to an integer - not sure for what so I'm commenting this out --ronj
			//if ( ! donElement.compare( "O" ) ) { donElemNum = 8; }
			//else if ( ! donElement.compare("N") ) { donElemNum = 7; }
			//else if ( ! donElement.compare("S") ) { donElemNum = 16; }
			//else { donElemNum = -1; }

			//if ( ! accElement.compare( "O" ) ) { accElemNum = 8; }
			//else if ( ! accElement.compare("N") ) { accElemNum = 7; }
			//else if ( ! donElement.compare("S") ) { accElemNum = 16; }
			//else { accElemNum = -1; }

			// the "type" of the hbond
			//std::cout << type << "\t";

			// chemical information about donor and acceptor
			std::cout << don_chain << " " << don_resname << " " << I(4,don_resnum) << " " << I(4,don_pdbresnum) << don_pdbicode << " " << don_atomname << " - ";
			std::cout << acc_chain << " " << acc_resname << " " << I(4,acc_resnum) << " " << I(4,acc_pdbresnum) << acc_pdbicode << " " << acc_atomname << "; ";

			// Energy and environmental weight
			std::cout << "hbE: " << F(7,4,energy) << ", weight: " << F(5,3,weight) << "; ";

			// misc parameters
			std::cout << I(2,donNbrs) << " " << F(5,2,donSasa) << " " << F(4,1,donBfactor) << " " << donSS << " - ";
			std::cout << I(2,accNbrs) << " " << F(5,2,accSasa) << " " << F(4,1,accBfactor) << " " << accSS << " ";

			//std::cout << h_atomname << " ";

			// don't care about the hbond energy function derivatives either --ronj
			// first and second derivative
			//std::cout << deriv.first[  0 ] << "\t" << deriv.first[  1 ] << "\t" << deriv.first[  2 ] << "\t";
			//std::cout << deriv.second[ 0 ] << "\t" << deriv.second[ 1 ] << "\t" << deriv.second[ 2 ] << "\t";

			// coordinates of atoms
			/*std::cout << datm_xyz.x()         << "\t" << datm_xyz.y()         << "\t" << datm_xyz.z()         << "\t";
			std::cout << hatm_xyz.x()         << "\t" << hatm_xyz.y()         << "\t" << hatm_xyz.z()         << "\t";
			std::cout << aatm_xyz.x()         << "\t" << aatm_xyz.y()         << "\t" << aatm_xyz.z()         << "\t";
			std::cout << base_xyz.x()         << "\t" << base_xyz.y()         << "\t" << base_xyz.z()         << "\t";
			std::cout << base_of_base_xyz.x() << "\t" << base_of_base_xyz.y() << "\t" << base_of_base_xyz.z() << "\t";
			std::cout << base2_xyz.x()        << "\t" << base2_xyz.y()        << "\t" << base2_xyz.z()        << "\t";
			std::cout << dbase_xyz.x()        << "\t" << dbase_xyz.y()        << "\t" << dbase_xyz.z()        << "\t";*/

			std::cout << std::endl;

			nhbonds_used++;
		}

		TR << "DONE getting hbonds. Number hydrogen bonds found: " << set1.nhbonds() << std::endl;

	}

private:
	bool allowNonProtein_;
	core::scoring::ScoreFunctionOP scfxn;
	std::string relevant_chains_;

};


int main( int argc, char* argv[] ) {

	try {


		//NEW_OPT( HBondReporter::relax, "Perform relaxation before dumping hbonds", false );
		NEW_OPT( HBondReporter::relevant_chains, "relevant_chains", "*" );

		devel::init(argc, argv);
		protocols::jd2::JobDistributor::get_instance()->go( protocols::moves::MoverOP( new HBondReporter ) );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

