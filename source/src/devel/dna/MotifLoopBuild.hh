#include <devel/enzdes/EnzdesRemodelProtocol.hh>

#include <core/types.hh>

#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <protocols/motifs/MotifLibrary.hh>


#include <core/scoring/ScoreFunction.fwd.hh>

#ifndef INCLUDED_devel_enzdes_motif_loop_build_hh
#define INCLUDED_devel_enzdes_motif_loop_build_hh

namespace devel {
namespace dna {

class MotifLoopBuild : public devel::enzdes::EnzdesRemodelMover
{
public:
	//For parser

	void apply(core::pose::Pose & pose) override;

	//Class functionality
	MotifLoopBuild();
	~MotifLoopBuild() override;

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
