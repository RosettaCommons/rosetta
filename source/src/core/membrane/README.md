Rosetta Membrane Protein Framework
====================================
Authors: 
Rebecca Faye Alford (rfalford12@gmail.com)
Julia Koehler Leman (julia.koehler1982@gmail.com)
Jeffrey J. Gray (jgray@jhu.edu)

Last Updated: 7/26/14

Database: 
------------------------------
  = MEM.params - Membrane position residue with default position normal=(0, 0, 1), center=(0, 0, 0)
  = EMB.params - Embedding position residue with default position normal=(0, 0, 1), center=(0, 0, 0)

MEM and EMB params have the residue type property MEMBRANE and are described in the fullatom and centroid residue type sets. 

Core Library: 
------------------------------
core/conformation/membrane: 
  = Membrane Exception Hierarchy - JD2 exceptions for when things go wrong. Shoud be used in the core layer, 
    and are particularly useful for designing testable IO code. 
  = Span: Describe a single transmembrane spanning region in the pose
  = SpanningTopology: Describe the total spanning topology of the protein (stored as a colleciton of Spans)
  = LipidAccInfo: Describe per-residue lipid exposed and buried surface areas (predicted by run_lips.pl 
    server)
  = MembraneParams: Enum for cleaning up code related to membrane position coordinate access
  = MembranePlanes: Store the position of the membrane planes as a series of virtual residues. Will display 
    as points - later used by the PyMol viewer for drawing CGO planes in visualization
  = MembraneInfo: Top Level container for membrane data (non-coordinate derived). 

core/membrane/geometry: 
  = EmbeddingDef: Define the embedding of a part of the pose - can be a single span, chain or whole pose. 
    Described by an embedding normal/center parameter
  = Embedding: Class for computing the embedding parameter. Can be from topology, sequence, or structure. 

Edits to Conformation: 
  = Methods for setting up a MembraneInfo object in conformation
  = Mehtods for accessing membrane positon parameters
  = Methods for dealing with coordinate derived membrane information

Protocols Library: 
-----------------------------
protocols/membrane: 
  = AddMembraneMover: Initialize the membrane protein framework. This mover adds a membrane position residue, 
    spanning topology, and optionally a lipophilicity object to the pose. Can be accessed via the command line, rosetta scripts or custom constructors (PyRosetta)
  = MembranePositionFromTopology: Move the membrane position to one based on the spanning topology and CA 
    coordinates in the pose. 
  = SetMembranePositionMover: Set the membrane normal/center position in the pose as custom positions
  = RandomMembranePositionMover: Randomnly rotate/translate the membrane position (first step to sampling)

PyMOL Mover Updates: 
-----------------------------
 = ShowMembranePlanesMover (in protocols/membrane/visualize) - Will add residues defining a top and bottom 
   plane to represent the membrane planes. This feature is designed specifically for visualization, but should not affect simulation results
 = PyMOLMover: Added flags to send membrane center, normal and plane points over to the PyMolPyRosetta server
 = PyMolPyRosettaServer: Recieves membrane message, will draw planes representing the membrane for real-time
   visualization

   This feature can be turned on using the flag -membrane_new:view_in_pymol when an application is supported
   by the framework (see AddmembraneMover).


