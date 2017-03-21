// TODO: document

#include "pathways_planners.h"

namespace pathways
{
// returns the angle differnce of to-from in the range [-180,180]
double Linear_planner_iterator::diff_angle(double from, double to) const
{
	double diff = to - from;
	if (diff < - 180.0) diff += 2*180.0;
	else if (diff > 180.0) diff = diff - 2*180.0;
	return diff;
}

// calculates a vector of the delta-values for each step in the iteration,
// s.t. the maximal step size of any given DOF is max_step_size
// NOTE: if step size does not allow even one step, than reduce step size to allow at least this step
void Linear_planner_iterator::calc_step_vector(double max_step_size)
{
	using namespace std;
	debug_assert(max_step_size > 0);
	// calculate difference vector between from and to
	double max_diff = 0;
	_v_step_delta.clear(); _v_step_delta.resize(_v_from.size());
	for(unsigned int i = 0; i < _v_step_delta.size();i++)
	{
		double cur_diff = diff_angle(_v_from[i], _v_to[i]);
		if ( std::abs(cur_diff) > max_diff)
			max_diff = std::abs(cur_diff);
		_v_step_delta[i] = cur_diff;
	}
	// make sure there are at least one steps (so reduce max_step_size to force that if neede)
	if(max_step_size > 0.5 * max_diff)
		max_step_size = 0.5 * max_diff; // 0.49 to avoid numerical problems
	// now normalize vector s.t. |max_diff| will reduce to max_step_size
	double normalize_factor = (max_diff > 0) ? (max_step_size / max_diff) : 0; // well max_diff==0 should hopefully never happen
	for(unsigned int i = 0; i < _v_step_delta.size();i++)
	{
		_v_step_delta[i] *= normalize_factor;
	}
	_total_steps_num = int(round( std::abs(max_diff / max_step_size)));

}
} // namespace pathways
