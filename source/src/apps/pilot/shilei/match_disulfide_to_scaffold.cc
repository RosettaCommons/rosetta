/**********
	This file was written by Will Sheffler to search the PDB for short loops with a disulfide bridge at their base.
	This file was modified by Vikram K. Mulligan on 23 August 2012 to check for breaks and to score the degree of loop exposure.
	Modify by Lei Shi to grafting a disulfide geometry to an exisiting scaffold.
**********/

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/conformation/util.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <stdio.h>
// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
// #include <devel/init.hh>
//Added by VKM:
#include <core/pose/PDB_Info.hh>
#include <numeric/xyzVector.hh>
#include <core/scoring/sasa.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/disulfides/FullatomDisulfideEnergyContainer.hh>
#include <core/scoring/Energies.hh>

static basic::Tracer TR( "apps.pilot_apps.shilei.match_disulfide_to_scaffold" );

OPT_1GRP_KEY(Integer,match_disulfide_to_scaffold,cdsf_max_res)
OPT_1GRP_KEY(String,match_disulfide_to_scaffold,input_disulfide)

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( match_disulfide_to_scaffold::cdsf_max_res ,"loop size" , 27 );
	NEW_OPT( match_disulfide_to_scaffold::input_disulfide,"input pdb containting disulfide to be matched","dummy");
}

int main(int argc, char *argv[]) {
	try {
	register_options();
	devel::init(argc,argv);
	using namespace std;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;

	//input disulfide geometry
        std::string input_disulfide=basic::options::option[ basic::options::OptionKeys::match_disulfide_to_scaffold::input_disulfide];
	std::ifstream infile;
   	infile.open(input_disulfide.c_str(),std::ios::in);
   	if ( ! infile.good() ) {
      		utility_exit_with_message( "Unable to open file: " + input_disulfide);
    	}

	//read file and check if there is disulfide
	core::pose::Pose pose1;
	core::import_pose::pose_from_pdb(pose1,input_disulfide);

	core::scoring::ScoreFunctionOP disulfide_scorefxn( core::scoring::ScoreFunctionFactory::create_score_function( "empty" ) );
  	disulfide_scorefxn->set_weight( dslf_ss_dst, 1.0);
  	disulfide_scorefxn->set_weight( dslf_cs_ang, 1.0);
  	disulfide_scorefxn->set_weight( dslf_ss_dih, 1.0);
  	disulfide_scorefxn->set_weight( dslf_ca_dih, 1.0);
	(*disulfide_scorefxn)(pose1);
	TR << "disulfide energy: " << pose1.energies().total_energies().dot( disulfide_scorefxn->weights() ) << std::endl;

	if ( pose1.energies().total_energies().dot( disulfide_scorefxn->weights() ) >= 0.0 ) {
		utility_exit_with_message(input_disulfide+" does not seem to have disulfide");
	}

	//need to find out which residues are in the disulfide
	core::scoring::disulfides::FullatomDisulfideEnergyContainerCOP dec = new core::scoring::disulfides::FullatomDisulfideEnergyContainer( pose1 );
	core::Size dslf_i=0,dslf_j=0;
	for (core::Size i=1; i<=pose1.total_residue(); i++) {
		if (pose1.residue(i).aa() != core::chemical::aa_cys) continue;
		if (dec->residue_forms_disulfide(i)) {
			dslf_i=i;
			dslf_j=dec->other_neighbor_id(i);
			TR << "found first disulfide between " << dslf_i << " and " << dslf_j << std::endl;
			break;
		}
	}

	if (dslf_i==0 || dslf_j==0) {
		utility_exit_with_message("could not find pairing disulfide");
	}

	/*
	//input matching scaffold
	utility::vector1<string> fnames = option[in::file::s]();

	for(int ifile=1; ifile <= (int)fnames.size(); ++ifile) { //Loop through all files provided
		string fn = utility::file_basename(fnames[ifile]);

		cout << "BEGIN " << fn << endl;
		core::pose::Pose pose;
		core::import_pose::pose_from_pdb(pose,fnames[ifile]); //Import the current PDB file
		for(int ir=1; ir <= (int)pose.n_residue(); ++ir) { //Loop through all residues
			if(pose.residue(ir).aa() != core::chemical::aa_cys) continue; //If the current residue isn't a cys, go on to the next.
			for(int jr=ir+3; jr < ir+option[cdsf_max_res](); ++jr) {
				if(jr > (int)pose.n_residue()) break;
				if(pose.residue(jr).aa() != core::chemical::aa_cys) continue;
				if(!core::conformation::is_disulfide_bond( pose.conformation(),ir,jr)) continue;
				bool terminus = false;
				bool gaps = false;
				bool baddistance = false;
				bool missingCN = false;
				numeric::xyzVector<double> Cpos;
				numeric::xyzVector<double> Npos;
				double CNdistance = 0;

				printf("Found disulfide-bridged loop from residue %i to residue %i.  Checking for gaps.\n", ir, jr);
				for(int i=ir; i < jr; ++i) {
					terminus |= pose.residue(i).is_terminus(); //Check up to aa just before the second member of the disulfide for a C-terminus
					//The above uses the bitwise OR operator |= .

					//Added by VKM: check that numbering in the pdb file is continuous:
					gaps|= (pose.pdb_info()->number(i+1) != (pose.pdb_info()->number(i)+1));
					//printf("%i\t%i\t%i\t%i\n", pose.pdb_info()->number(i+1), pose.pdb_info()->number(i), (int)(pose.pdb_info()->number(i+1) != (pose.pdb_info()->number(i)+1)), (int)gaps); //DELETE ME

					//Added here by VKM: check whether residue i is connected to residue i+1:
					gaps |= !pose.residue(i).is_bonded(i+1); //If consecutive residues aren't bonded, there's a gap.
					//if(pose.residue(i).is_bonded(i+1)) printf("Residue %i is bonded to residue %i.\n", i, i+1);
					//else printf("There's a gap between residues %i and %i.", i, i+1);
					//if(pose.residue(i).is_bonded(i+2)) printf("Residue %i is bonded to residue %i.\n", i, i+2);
					//else printf("There's a gap between residues %i and %i.", i, i+2);

					//check whether the distace from the nth C to the n+1st N is compatible with a peptide  bond:
					if(pose.residue(i).has("C") && pose.residue(i+1).has("N")) {
						Cpos=pose.residue(i).xyz("C");
						Npos=pose.residue(i+1).xyz("N");
						CNdistance=Cpos.distance(Npos);
						baddistance|=(CNdistance > 1.38 || CNdistance < 1.28);
						//printf("CN distance for peptide bonds between residues %i and %i is %f.\n", i, i+1, CNdistance);
					}
					else missingCN = true; //If either the C or the N is missing, this is a bad loop.

				}
				if(terminus) {
					printf ("C-terminus found.  This isn't a loop.  Discarding and moving on.\n");
					continue; //If there's a terminus in this loop, give up on it.
				}
				if(gaps) {
					printf ("Gap found.  This is a bad loop with missing data.  Discarding and moving on.\n");
					continue; //If there's a gap in this loop, give up on it.
				}
				if(missingCN) {
					printf ("Peptide bond atom(s) missing.  This is a bad loop with missing data.  Discarding and moving on.\n");
					continue; //If atoms are missing, give up on this loop.
				}
				if(baddistance) {
					printf("Bad peptide bond length detected.  This is a bad loop, probably with missing data.  Discarding and moving on.\n");
					continue; //If there are bad peptide bond lengths between adjacent residues, give up on this loop.
				}

				core::pose::Pose tmp; //If we haven't given up on this loop, create a new pose called "tmp"
				tmp.append_residue_by_jump(pose.residue(ir),1); //Add the initial residue to the pose.
				for(int i=ir+1; i <= jr; ++i) tmp.append_residue_by_bond(pose.residue(i)); //Add all subsequent residues in the loop to the pose.

				int exhb=0,inhb=0,inhbsc=0;
				get_hb_info(pose,ir,jr,exhb,inhb,inhbsc); //Get some stats about the residues in the loop.

				using namespace ObjexxFCL::format;
				//Indicate a hit, and write out some stats:
				cout << "HIT " << " " << I(4,ir) << " " << I(4,jr) << " " << I(4,exhb) << " " << I(4,inhb) << " " << I(4,inhbsc) << " " << outfile << endl;
				//tmp.dump_pdb(outfile); //Write out the "tmp" pose containing just the loop.
				loopnames.push_back(outfile); //Add the file name to the list of output file names.

			//Added by VKM: calculate solvent-accessible surface area (SASA) per atom for this pose, and store some stats:
				//Add elements to the SASA vectors:
				loop_sasa_alone.push_back(0);
				loop_sasa_instruct.push_back(0);
				loop_sasa_ratio.push_back(0);
				//Calculate the SASA for the whole structure (to get the loop SASA in the context of the structure):
				core::id::AtomID_Map<double> atom_sasa;
				utility::vector1<double> rsd_sasa;
				core::scoring::calc_per_atom_sasa(pose, atom_sasa, rsd_sasa, 2.0, true);
				//Get the loop SASA in the context of the structure:
				for (int vi = ir+1; vi < jr; ++vi) //Loop through all residues BETWEEN the cysteines forming the disulfide.
				{
					loop_sasa_instruct[loop_sasa_instruct.size()]+=rsd_sasa[vi];
					printf("Residue %i in structure:\t%f\n", vi, rsd_sasa[vi]); //DELETE ME
				}
				printf("Loop SASA in structure:\t%f\n", loop_sasa_instruct[loop_sasa_instruct.size()]); //DELETE ME
				//Calculate the SASA for the loop in isolation:
				core::id::AtomID_Map<double> atom_sasa_alone;
				utility::vector1<double> rsd_sasa_alone;
				core::scoring::calc_per_atom_sasa(tmp, atom_sasa_alone, rsd_sasa_alone, 2.0, true);
				for (int vi = 2; vi < (int)tmp.n_residue(); ++vi) //Loop through all residues BETWEEN the cysteines.
				{
					loop_sasa_alone[loop_sasa_alone.size()]+=rsd_sasa_alone[vi];
					printf("Residue %i in isolation:\t%f\n", vi, rsd_sasa_alone[vi]); //DELETE ME
				}
				printf("Loop SASA in isolation:\t%f\n", loop_sasa_alone[loop_sasa_alone.size()]); //DELETE ME
				//Calculate ratio:
				loop_sasa_ratio[loop_sasa_ratio.size()]=loop_sasa_instruct[loop_sasa_instruct.size()]/loop_sasa_alone[loop_sasa_alone.size()];
				printf("Loop SASA ratio:\t%f\n", loop_sasa_ratio[loop_sasa_ratio.size()]); //DELETE ME

				//Dump the PDB to a file ONLY if the loop is particularly exposed.
				if(loop_sasa_ratio[loop_sasa_ratio.size()] > 0.85) {
					printf("85% exposed loop!  Writing to file.\n");
					tmp.dump_pdb(outfile); }
			}
		}
		cout << "DONE " << fn << endl;
	}
	*/

        } catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
        }

	return 0;
}

