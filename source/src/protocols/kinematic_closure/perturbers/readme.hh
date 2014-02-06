// #error Don't include; for documentation only.

/// @namespace protocols::kinematic_closure::perturbers
///
/// @brief Algorithms for sampling new backbone torsions.
///
/// @details The kinematic closure algorithm allows the non-pivot torsions to 
/// be set without any restrictions.  The six pivot torsions are then used to 
/// ensure the backbone remains closed.  There is a great deal of flexibility 
/// in how the non-pivot torsions are chosen.  This namespace provides a number 
/// of different algorithms to sample these torsions.  These algorithms can be 
/// mixed in matched in any way you'd like using samplers::KicMover.
