#include <devel/enzdes/EnzdesRemodelProtocol.hh>
#include <utility/file/FileName.hh>

#include <core/types.hh>

#include <protocols/enzdes/EnzdesFlexBBProtocol.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <protocols/motifs/Motif.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <protocols/motifs/MotifSearch.hh>
#include <protocols/motifs/motif_utils.hh>

#include <protocols/dna/util.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/ScoreFunction.hh>

#ifndef INCLUDED_devel_enzdes_motif_loop_build_hh
#define INCLUDED_devel_enzdes_motif_loop_build_hh

namespace devel {
namespace dna {

class MotifLoopBuild : public devel::enzdes::EnzdesRemodelMover
{
public:
	//For parser

	void apply(core::pose::Pose & pose);

	//Class functionality
	MotifLoopBuild();
	~MotifLoopBuild();

	void irc_build(core::pose::Pose & pose);
	void place_motifs( core::pose::Pose & pose,
		utility::vector1< core::Size > & flex_pos,
		protocols::motifs::MotifLibrary & motif_lib);
private:

	void build_inv_rots(core::pose::Pose  pose );

	protocols::motifs::MotifLibrary motif_lib;

	void mutate_dna(core::pose::Pose & pose);

	utility::vector1< core::Size > dna_design_pos_;
	core::scoring::ScoreFunctionOP  scorefxnOP;

	core::pose::Pose pose_;
	core::pose::Pose mut_pose_;

	std::map< std::string, core::pack::rotamer_set::RotamerSetOP >
		rot_map;
};


}}

#endif
