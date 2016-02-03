/// @brief   Pilot application to generate shape complementarity scores for biosensor design project
/// @author  Kyle Barlow

// Unit headers

// Project headers
#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <core/scoring/sc/ShapeComplementarityCalculator.hh>

static THREAD_LOCAL basic::Tracer TR( "ShapeScoring" );

//using namespace core::pack
//using namespace core::pack::task;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

typedef core::Size Size;
using utility::vector1;

void calc_shape_complementarity(core::pose::Pose &pose){

  vector1<Size> chain_lengths(10); // Maximum of 10 chains in pose supported
  Size number_of_chains = 0;
  Size ligand_chain_num = 0; // Will be set to the chain number of the ligand-containing chain
  Size protein_chain_a_num = 0; // Will be set to chain number of one of the protein chains
  Size protein_chain_b_num = 0; // Will be set to chain number of one of the protein chains

  for(utility::vector1<Size>::iterator it = chain_lengths.begin(); it != chain_lengths.end(); ++it) {
    *it = 0;
  }

  for (Size i=1; i<=pose.n_residue(); ++i) {
    //TR << pose.chain(i) ;
    chain_lengths[pose.chain(i)]+=1;
  }
  //TR << std::endl;

  for(utility::vector1<Size>::iterator it = chain_lengths.begin(); it != chain_lengths.end(); ++it) {
    if (*it>0)
      number_of_chains++;
  }

  runtime_assert(number_of_chains==3); // Specific to biosensor design

  for(Size i = 1; i <= chain_lengths.size(); i++) {
    if (chain_lengths[i]==1) {
      runtime_assert(ligand_chain_num==0); // Ensure that another chain of length 1 hasn't already been found
      ligand_chain_num=i;
    }
  }

  runtime_assert(ligand_chain_num!=0); // Ensure at least one chain of length 1 was found

  TR << "Ligand is on chain " << ligand_chain_num << std::endl;

  if (ligand_chain_num==1) {
    protein_chain_a_num = 2;
    protein_chain_b_num = 3;
  }
  else if (ligand_chain_num==2) {
    protein_chain_a_num = 1;
    protein_chain_b_num = 3;
  }
  else if (ligand_chain_num==3) {
    protein_chain_a_num = 1;
    protein_chain_b_num = 2;
  }

  core::scoring::sc::ShapeComplementarityCalculator scc_chain_a_ligand;
  scc_chain_a_ligand.Init();
  for (Size i=1; i<=pose.n_residue(); ++i) {
    if (pose.chain(i)==protein_chain_a_num)
      scc_chain_a_ligand.AddResidue(0,pose.residue(i));
    else if (pose.chain(i)==ligand_chain_num)
      scc_chain_a_ligand.AddResidue(1,pose.residue(i));
  }

  core::scoring::sc::ShapeComplementarityCalculator scc_chain_b_ligand;
  scc_chain_b_ligand.Init();
  for (Size i=1; i<=pose.n_residue(); ++i) {
    if (pose.chain(i)==protein_chain_b_num)
      scc_chain_b_ligand.AddResidue(0,pose.residue(i));
    else if (pose.chain(i)==ligand_chain_num)
      scc_chain_b_ligand.AddResidue(1,pose.residue(i));
  }

  if (scc_chain_a_ligand.Calc()) {
    TR << "CHAIN_A_LG_SC " << scc_chain_a_ligand.GetResults().sc << std::endl;
    TR << "CHAIN_A_LG_AREA " << scc_chain_a_ligand.GetResults().area << std::endl;
    TR << "CHAIN_A_LG_DIST " << scc_chain_a_ligand.GetResults().distance << std::endl;
  }

  if (scc_chain_b_ligand.Calc()) {
    TR << "CHAIN_B_LG_SC " << scc_chain_b_ligand.GetResults().sc << std::endl;
    TR << "CHAIN_B_LG_AREA " << scc_chain_b_ligand.GetResults().area << std::endl;
    TR << "CHAIN_B_LG_DIST " << scc_chain_b_ligand.GetResults().distance << std::endl;
  }

  return;
}

int main(int argc, char *argv[])
{
  // Initialize core.
  devel::init(argc, argv);

  // check if input PDB file is provided from command line option -s
  if ( !basic::options::option[ basic::options::OptionKeys::in::file::s ].user() ) {
    // exit if no PDB file is found
    utility_exit_with_message("Input PDB file not found");
    }

  // initialize pose
  core::pose::Pose p;
  // get name of pdb file from command line option -s
  std::string pdb_file = basic::options::option[ basic::options::OptionKeys::in::file::s ]()[1];
  // load pdb file into pose
  core::import_pose::pose_from_file( p, pdb_file , core::import_pose::PDB_file);
  // initialize score function
  ////core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();
  // score pose
  ////(*score_fxn)(p);
  // output score
  ////score_fxn->show(TR, p);

  calc_shape_complementarity(p);

  TR.flush();
}

