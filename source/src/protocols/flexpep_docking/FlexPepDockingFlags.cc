// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (C) 199x-2008 Hebrew University, Jerusalem
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file   FlexPepDockingFlags.cc
///
/// @brief flags structure for FlexPepDocking protocols
/// @date January 1, 2009
/// @author Barak Raveh

#include <protocols/flexpep_docking/FlexPepDockingFlags.hh>

#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/flexPepDocking.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <utility/exit.hh>
#include <core/pose/util.hh>
#include <iostream>
#include <fstream>
#include <string>

#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "protocols.flexPepDocking.FlexPepDockingFlags" );

using namespace protocols::flexpep_docking;


//////////////////////////////////////////////////
/// ctr
///
/// @brief
/// initialize all flags from cmd-line options
//////////////////////////////////////////////////
protocols::flexpep_docking::FlexPepDockingFlags::FlexPepDockingFlags
()
{
	using namespace basic::options;

	// a-priori all flags are invalid
	valid_chain_bounds_ = false;
	valid_receptor_chain_ = false;
	valid_peptide_chain_ = false;
	valid_ref_start_struct_ = false;
	receptor_anchor_pos = -1; // -1 == invalid

	pep_fold_only = option[ OptionKeys::flexPepDocking::pep_fold_only ]();
  if (option[ OptionKeys::flexPepDocking::ref_startstruct].user() ) {
		ref_start_struct_ =
			option[ OptionKeys::flexPepDocking::ref_startstruct]();
		valid_ref_start_struct_ = true;
	}
	min_only = option[ OptionKeys::flexPepDocking::flexPepDockingMinimizeOnly ]();
	random_phi_psi_pert = false;
	if (option[ OptionKeys::flexPepDocking::random_phi_psi_preturbation ]())
		random_phi_psi_pert = true;// TODO: redundant, default value is enough
	random_phi_psi_pert_size =
		option[ OptionKeys::flexPepDocking::random_phi_psi_preturbation ]();
	extend = option[ OptionKeys::flexPepDocking::extend_peptide ]();
	randomRBstart = false; // TODO: redundant, default value is enough
	if (option[ OptionKeys::flexPepDocking::random_trans_start ]() ||
		option[ OptionKeys::flexPepDocking::random_rot_start ]())
		randomRBstart = true;
	lowres_abinitio = option[ OptionKeys::flexPepDocking::lowres_abinitio ]();
	lowres_preoptimize = option[ OptionKeys::flexPepDocking::lowres_preoptimize ]();
  pep_refine = option[ OptionKeys::flexPepDocking::pep_refine ]();
	rbMCM = option[ OptionKeys::flexPepDocking::rbMCM ](); // obsolete
	torsionsMCM = option[ OptionKeys::flexPepDocking::torsionsMCM ]();  // obsolete
	// the next section should be removed after we get rid of obsolee rmMCM and torsionsMCM completely
	if(pep_refine || lowres_preoptimize || lowres_abinitio)
		{ // overrides old rbMCM and torsionsMCM
			bool explicitFalseMCMs = 
				( option[ OptionKeys::flexPepDocking::rbMCM ].user() && !rbMCM) ||
				( option[ OptionKeys::flexPepDocking::torsionsMCM ].user() && !torsionsMCM);
			runtime_assert_msg(!explicitFalseMCMs,
												 "pep_refine / lowres_preoptimize / lowres_abinitio are not compatible with explicitly setting -rbMCM or -torsionsMCM to false");
			rbMCM = true;
			torsionsMCM = true;
		}

	rb_trans_size = option[ OptionKeys::flexPepDocking::random_trans_start ]();
	rb_rot_size = option[ OptionKeys::flexPepDocking::random_rot_start ]();
	peptide_loop_model = option[ OptionKeys::flexPepDocking::peptide_loop_model ]();
	smove_angle_range =
		option[ OptionKeys::flexPepDocking::smove_angle_range ]();
	design_peptide = option [ OptionKeys::flexPepDocking::design_peptide ]();
	backrub_opt = option[ OptionKeys::flexPepDocking::backrub_peptide ]();
	boost_fa_atr = option[ OptionKeys::flexPepDocking::boost_fa_atr ]();
	ramp_fa_rep = option [ OptionKeys::flexPepDocking::ramp_fa_rep ]();
	ramp_rama = option [ OptionKeys::flexPepDocking::ramp_rama ]();
	rep_ramp_cycles =
		option[ OptionKeys::flexPepDocking::rep_ramp_cycles ](); //10 default
	mcm_cycles =
		option[ OptionKeys::flexPepDocking::mcm_cycles ](); //8 default
	score_only = option[ OptionKeys::flexPepDocking::flexpep_score_only ]();
	use_cen_score = option[ OptionKeys::flexPepDocking::use_cen_score ]();
	min_receptor_bb = option[ OptionKeys::flexPepDocking::min_receptor_bb ];
	ppk_only = option[ OptionKeys::flexPepDocking::flexpep_prepack ]();
	no_prepack1 = option[ OptionKeys::flexPepDocking::flexpep_noprepack1 ]();
	no_prepack2 = option[ OptionKeys::flexPepDocking::flexpep_noprepack2 ]();
  score_filter = option[ OptionKeys::flexPepDocking::score_filter ]();
	hb_filter  = option[ OptionKeys::flexPepDocking::hb_filter ]();
	hotspot_filter  = option[ OptionKeys::flexPepDocking::hotspot_filter ]();
	frag3_weight  = option[ OptionKeys::flexPepDocking::frag3_weight ]();
	frag5_weight  = option[ OptionKeys::flexPepDocking::frag5_weight ]();
	frag9_weight  = option[ OptionKeys::flexPepDocking::frag9_weight ]();
	pSer2Asp_centroid = option[ OptionKeys::flexPepDocking::pSer2Asp_centroid ]();
	pSer2Glu_centroid = option[ OptionKeys::flexPepDocking::pSer2Glu_centroid ]();
	dumpPDB_abinitio = option[ OptionKeys::flexPepDocking::dumpPDB_abinitio ]();
	dumpPDB_lowres = option[ OptionKeys::flexPepDocking::dumpPDB_lowres ]();
	dumpPDB_hires = option[ OptionKeys::flexPepDocking::dumpPDB_hires ]();


  if ( option[ OptionKeys::flexPepDocking::receptor_chain ].user() )
		{
			this->set_receptor_chain
				(option[ OptionKeys::flexPepDocking::receptor_chain ]() );
			// TODO: validate string size!
		}
	if ( option[ OptionKeys::flexPepDocking::peptide_chain ].user() )
		{
			this->set_peptide_chain
				(option[ OptionKeys::flexPepDocking::peptide_chain ]().at(0) );
			// TODO: validate string size!
		}
	// params file
	if(option[ OptionKeys::flexPepDocking::params_file ].user() )
		{
			runtime_assert_msg(! (valid_peptide_chain_ || valid_receptor_chain_),
												 "Params file and cmd-line peptide or receptor chain are mutually exclusive");
			params_file = option[ OptionKeys::flexPepDocking::params_file ]();
			updateChainsAndAnchors_fromParamsFile(params_file);
		}
	runtime_assert_msg(! (pep_fold_only && (min_receptor_bb || valid_receptor_chain_) ),
										 "The flag -pep_fold_only is incompatible with receptor flags like -min_receptor_bb and -receptor chain");
	int mut_ex_opts = 0;
	if(min_only) mut_ex_opts++;
	if(ppk_only) mut_ex_opts++;
	if(torsionsMCM || rbMCM || lowres_preoptimize)
		mut_ex_opts++;
	runtime_assert_msg(mut_ex_opts <= 1,
		"More than one mutually exclusive flags selected");

}

bool FlexPepDockingFlags::is_ligand_present( core::pose::Pose const& pose ) const
 {
		 if(valid_chain_bounds_) {
				 if (pose.total_residue() > (Size)(peptide_nres_ + receptor_nres_))
						 return true;
				 else
						 return false;
		 }
		 std::cout << "ERROR: chain bounds invalid" << std::endl;
		 exit(-1);
 }

std::string FlexPepDockingFlags::receptor_chain() const
	{
		if(valid_chain_info() && !pep_fold_only) return receptor_chain_;
		std::cout << "ERROR: receptor chain invalid" << std::endl;
		exit(-1);
	}


char FlexPepDockingFlags::peptide_chain() const
{
	if(valid_chain_info()) return peptide_chain_;
	std::cout << "ERROR: peptide chain invalid" << std::endl;
	exit(-1);
}


int FlexPepDockingFlags::receptor_first_res() const
{
	if(valid_chain_info() && !pep_fold_only) return receptor_first_res_;
	std::cout << "ERROR: chain bounds invalid: receptor_first " << std::endl;
	exit(-1);
}


int FlexPepDockingFlags::receptor_last_res() const
{
	if(valid_chain_info() && !pep_fold_only)
		return receptor_first_res_ + receptor_nres_ - 1;
	std::cout << "ERROR: chain bounds invalid: receptor_last " << std::endl;
	exit(-1);
}


int FlexPepDockingFlags::receptor_nres() const
{
	if(valid_chain_info()) return receptor_nres_;
	if(pep_fold_only) return 0;
	std::cout << "ERROR: chain bounds invalid: receptor_nres " << std::endl;
	exit(-1);
}


int FlexPepDockingFlags::peptide_first_res() const
{
	if(valid_chain_info()) return peptide_first_res_;
	std::cout << "ERROR: chain bounds invalid: pep_first " << std::endl;
	exit(-1);
}


int FlexPepDockingFlags::peptide_last_res() const
{
	if(valid_chain_info())
		return peptide_first_res_ + peptide_nres_ - 1;
	std::cout << "ERROR: chain bounds invalid: pep_last " << std::endl;
	exit(-1);
}


int FlexPepDockingFlags::peptide_nres() const
{
	if(valid_chain_info()) return peptide_nres_;
	std::cout << "ERROR: chain bounds invalid: pep_nres " << std::endl;
	exit(-1);
}


std::string FlexPepDockingFlags::ref_start_struct() const
{
	if (valid_ref_start_struct()) return ref_start_struct_;
	std::cout << "ERROR: ref start structure not specified" << std::endl;
	exit(-1);
}


bool FlexPepDockingFlags::valid_ref_start_struct() const
{
	return valid_ref_start_struct_;
}


// if needed, set default chain ids for receptor and peptide
// and update parameters of chain length
//
// Default chains:
// if chains of receptor and peptide are invalid,
// use the first and second chain in the pose (respectively)
// if only one is invalid, take the first chain that differs from it
void
FlexPepDockingFlags::updateChains
( core::pose::Pose const& pose )
{
	// TODO: extend for multichains?
	// update chain ids if needed
	core::pose::PDBInfoCOP pdbinfo = pose.pdb_info();
	if(!pdbinfo){
		TR.Warning << "Missing PDBInfo in input pose - generating default PDBInfo from pose object" << std::endl;
		pdbinfo = core::pose::PDBInfoCOP( core::pose::PDBInfoOP( new core::pose::PDBInfo(pose) ) );
 	}
	core::Size resi = 1;
	// get receptor chain if needed
	if(!valid_receptor_chain_ && ! pep_fold_only)
		{
			receptor_chain_ = " ";
            receptor_chain_.at(0) = pdbinfo->chain(resi);
			// make sure receptor is different from peptide
            if(valid_peptide_chain_ && (receptor_chain_.at(0) == peptide_chain_)){
				do{	resi++;	}
				while(resi <= pose.total_residue() && pdbinfo->chain(resi) == peptide_chain_);
                receptor_chain_.at(0) = pdbinfo->chain(resi);
			}
			valid_receptor_chain_ = true;
		}
	// get peptide chain if needed
	if(!valid_peptide_chain_)
		{
			if(pep_fold_only)
				{
					this->set_peptide_chain(pdbinfo->chain(resi));
				} 
			else // docking mode
				{
					// skip receptor chain
                	while(resi <= pose.total_residue() &&
                				  receptor_chain_.find(pdbinfo->chain(resi)) != std::string::npos ) {
                		resi++;
                	}
					this->set_peptide_chain(pdbinfo->chain(resi));
				}
		}

	// find receptor boundaries
	resi = 1;
	if(! pep_fold_only)
		{
			while(resi <= pose.total_residue() && receptor_chain_.find(pdbinfo->chain(resi)) == std::string::npos) resi++;
			receptor_first_res_ = resi;
			do{ resi++; }
			while(resi <= pose.total_residue() && receptor_chain_.find(pdbinfo->chain(resi)) != std::string::npos);
			receptor_nres_ = resi - receptor_first_res_;
		}
	else
		{
			receptor_first_res_ = 0;
			receptor_nres_ = 0;
		}

	// find peptide boundaries
	resi = 1;
	while(resi <= pose.total_residue() && pdbinfo->chain(resi) != peptide_chain_) resi++;
	peptide_first_res_ = resi;
	do{ resi++; }
	while(resi <= pose.total_residue() && pdbinfo->chain(resi) == peptide_chain_);
	peptide_nres_ = resi - peptide_first_res_;
	if (peptide_nres_ > 30) {
		TR.Warning << "peptide chain is longer than 30 residues" << std::endl;
	}
	TR << "Receptor chain: " << receptor_chain_ << std::endl;
	TR << "Receptor first res: " << receptor_first_res_ << std::endl;
	TR << "Receptor nres: " << receptor_nres_ << std::endl;
	TR << "Peptide chain: " << peptide_chain_ << std::endl;
	TR << "Peptide first res: " << peptide_first_res_ << std::endl;
	TR << "Peptide nres: " << peptide_nres_ << std::endl;

	// declare valid results
	valid_chain_bounds_ = true;
}


// calc default anchors
void
FlexPepDockingFlags::setDefaultAnchors
( core::pose::Pose& pose ) // TODO: pose should be const, fix RB_geometry 4 this
{
	// TODO: extend for multichains?
	using namespace basic::options;
	using namespace protocols::geometry;
	peptide_anchors.clear();
	peptide_cuts.clear();
	// peptide anchor - from cmd line, if not set - peptide c.o.m
	if ( option[ OptionKeys::flexPepDocking::peptide_anchor ].user() )
    {
			this->peptide_anchors[1] = option[ OptionKeys::flexPepDocking::peptide_anchor ];
		}
	else {
	this->peptide_anchors[1] = core::pose::residue_center_of_mass
		( pose,
			peptide_first_res(),
			peptide_last_res() );
	}
	TR << "Peptide anchor: " << peptide_anchors[1] << std::endl;
	TR << "# peptide anchors: " << peptide_anchors.size() << std::endl;
	TR << "# peptide cuts: " << peptide_cuts.size() << std::endl;
  runtime_assert_msg(peptide_anchors[1] >= peptide_first_res() &&
										 peptide_anchors[1] <= peptide_last_res(),
										 "Peptide anchor out of range");

	// receptor anchor - nearest residue in protein
	core::Vector
		pep_anchor_ca( pose.residue( peptide_anchors[1] ).atom( "CA" ).xyz() );
	if(valid_receptor_chain_)
		{
			this->receptor_anchor_pos =	core::pose::return_nearest_residue
				( pose,
					receptor_first_res(),
					receptor_last_res() ,
					pep_anchor_ca);
			TR << "Receptor anchor: " << receptor_anchor_pos << std::endl;
			runtime_assert_msg( receptor_anchor_pos >= receptor_first_res() &&
													receptor_anchor_pos <= receptor_last_res() ,
													"Receptor anchor out of range");

		}
}


// support for old params-file format
// TODO: also take chain-ids from params file
// TODO: currently ignores cuts
// TODO: add support for pep_fold_only
void
FlexPepDockingFlags::updateChainsAndAnchors_fromParamsFile
( std::string const& params_file )
{
 	using namespace std;
 	bool local_debug = true;
	TR << "Reading params from file [" << params_file << "]" << endl;

	// invalidate any previous info
	peptide_anchors.clear();
	peptide_cuts.clear();
	valid_chain_bounds_ = false;
  valid_receptor_chain_ = false;
  valid_peptide_chain_ = false;
	receptor_nres_ = -1;
	peptide_first_res_ = -1;
	peptide_nres_ = -1;
	receptor_anchor_pos = -1;


 	// TODO: read params file somewhere else, move to XML?
	receptor_first_res_ = 1; // this is an assumption in params file
 	std::ifstream data( params_file.c_str() );
 	std::string line;
 	while ( getline( data,line) ) {
 		std::istringstream l( line );
 		int param;
 		std::string tag;
 		l >> tag >> param;
 		if ( tag == "scaffold_anchor_pos" ) { //TODO: change format of params to "receptor_anchor_pos"?
 			receptor_anchor_pos = param;
 			TR << "# receptor_anchor_pos " << param << endl;
 		}
 		else if ( tag == "nres_scaffold" ) {
			receptor_nres_ = param;
 			TR << "# nres_receptor " << param << endl;
 			peptide_first_res_ = receptor_nres_ + 1;	// TODO: is that always so?
 		}
 		else if ( tag == "nres_peptide" ) {
			peptide_nres_ = param;
			TR << "# nres_peptide " << param << endl;
 		}
 		else if ( tag == "peptide_anchor" ) {
 			int njump = param;
 			int pep_anchor;
 			l >> pep_anchor;
 			TR << "# peptide_anchor " << njump << " anchor " << pep_anchor <<endl;
 			peptide_anchors[njump]=pep_anchor;
 		}
 	}

 	TR << "# finished reading params" << endl;

	// offset peptide anchors / cuts residues, by peptide_first_res_
	std::map<int,int>::iterator iter;
	for(iter = peptide_anchors.begin(); iter != peptide_anchors.end(); iter++)
		iter->second += (peptide_first_res_ - 1);
	for(iter = peptide_cuts.begin(); iter != peptide_cuts.end(); iter++)
		iter->second += (peptide_first_res_ - 1);

	// check results validity
	if(receptor_first_res_ != -1 &&
		receptor_nres_ != -1 &&
		peptide_first_res_ != -1 &&
		peptide_nres_ != -1)
		{
			valid_chain_bounds_ = true;
		}

 	if(local_debug){
 		TR << "nres_receptor: " <<  receptor_first_res_ << endl
 			 << "nres_peptide: " << peptide_nres_ << endl
 			 << "pep_begin_res: " << peptide_first_res_ << endl
 			 << "receptor_anchor_pos: " << receptor_anchor_pos << endl
			 << "number of peptide anchors: " << peptide_anchors.size() << endl
			 << "number of peptide cuts: " << peptide_cuts.size() << endl
			 << endl;
 	}

	if(!valid_chain_bounds_ || !valid_anchors()){
		TR << "Missing or invalid information in parameters file "
			 << params_file << std::endl;
		exit(-1);
	}

}


