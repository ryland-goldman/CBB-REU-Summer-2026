## 1 Introduction

G4beamline is a particle tracking and simulation program based on the Geant4 toolkit [1] that is
specifically designed to easily simulate beamlines and related systems. It is flexible enough to simulate
complex beamlines like the MICE muon beam, the Neutrino Factory Study 2 SFOFO muon-cooling
channel, complex helical cooling channels, and many others. Because of its simple and straightforward
method of specifying the system to be simulated, it is also well suited for quickly answering questions
about particle interactions and tracking (e.g. “On average, how much energy does a 150 MeV proton
lose in a 1 mm Al window?”, “How large does the multiple-scattered beam grow 20 meters downstream
of the window?”). As a radioactive decay source, an isotropic source, and a Cosmic Ray “muon beam”
are included, the notion of “beamline” can be rather more general than usual.

The primary advantage of using G4beamline is that its description of the simulation is commensurate
with the complexity of the system being simulated, instead of being a significantly more complicated
C++ program. Most users need not face the challenges of learning C++ programming and the details of
the Geant4 toolkit – to use G4beamline there is no need to: a) know C++, b) learn the many aspects of
the Geant4 toolkit API, c) face the non-trivial challenges of installing the Geant4 toolkit and all its
required libraries, and d) learn how to solve any problems that arise while linking a complicated and
very large program. All of that is done during the production of the G4beamline distribution. Users with
special needs can download and install the source distribution, and learn how to build the program,
which will permit them to add their own C++ code and custom commands to G4beamline – this can be
much simpler than the direct use of Geant4 and its libraries.

The basic structure of a G4beamline simulation is to first define beamline elements (magnets, beam
pipes, windows, RF cavities, etc.), including their geometry, materials, fields, etc., and then to place
them into the world, usually along the beam direction. As bending magnets can be modeled, the “beam
direction” can change – see “Centerline Coordinates” below; it remains simple to place elements along
the nominal beam centerline. It should be noted that a G4beamline simulation is much closer to
specifying a real beamline that it is to the abstractions and approximations used in most acceleratorphysics codes. All descriptions and configurations are contained in a single ASCII input file, which also
provides values for various program parameters, specification of the initial beam, etc.

The tracking of particles through the simulated system is as accurate and realistic as the Geant4 toolkit
implements. The input file selects from any of the Geant4 physics lists, and can set values for the
various Geant4 tracking-accuracy parameters. This permits users to make trade-offs between CPU time
and simulation accuracy. Similarly, G4beamline permits the specification of magnetic map parameters,
permitting a trade-off between memory usage (and the CPU time to generate the map) and simulation
accuracy.

While G4beamline can make it rather simple to specify a simulation, it cannot substitute for knowledge
and experience about the problem domain or about particle-tracking simulations in general. Like all
computer programs, G4beamline is prone to “garbage in, garbage out”, especially when used by
unskilled users. It is strongly suggested that you use visualization to verify the geometry of your
simulation and that a handful of particles are tracked properly through it. Whenever possible you should
arrange to track through a simple geometry that you can compare to independent results, to make sure
that what you think is happening actually does occur in the simulation.
