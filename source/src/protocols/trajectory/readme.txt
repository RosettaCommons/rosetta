This namespace contains classes meant to work with long trajectories.  Although 
I am not an expert on the JD2, I don't believe it is suitable for recording 
long trajectories because output is only performed at the end of the job.  This 
means that the trajectory must be saved throughout the simulation, which may 
require a prohibitive amount of memory.

I don't really intend for the files in this namespace to be used in rosetta 
scripts.  Instead, I intend for them to be used internally within protocols 
written in C++.

Right now, all the classes in this namespace deal with databases, but I imagine 
that this could be expanded in the future.  For example, it might be nice to be 
able to read and write XTC trajectories.  If the need arises, I would advise 
making inheritance hierarchies that look something like this:

  | TrajectoryWriter (base writer class)
  |   DbTrajectoryWriter
  |   XtcTrajectoryWriter
  |   ...
  
  | TrajectoryReader (base reader class)
  |   DbTrajectoryReader
  |   XtcTrajectoryReader
  |   ...

The protocols::canonical_sampling namespace contains a class called  
PDBTrajectoryRecorder.  This class writes trajectories into a PDB format which 
is suitable for viewing in pymol.  I think this class would make more sense in 
this namespace.  If someone's interested in moving it, also change it's name to 
conform to the naming conventions using in this namespace: PdbTrajectoryWriter.  
A PdbTrajectoryReeader class seems unnecessary.
