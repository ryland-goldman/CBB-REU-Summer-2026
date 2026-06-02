# To-Do List: Linac Simulation

## Major Outstanding Items

### 1. CESR Control System Integration *(highest priority)*

Internal bookkeeping is coded, but no actual database calls are implemented. Loading parameters from the control system also requires calibrations — converting integer "computer units" to real machine units for each parameter. Some of these mappings will be non-linear, which complicates matters. Only linear mapping bookkeeping has been started.

### 2. Solenoid 1

Not yet in the simulation — more particles are lost in section 1 than necessary as a result. The Poisson input file is ready, and the GUI already has the parameter slots. Needed:
- Coil layout (check schematics 6043-03, 60-64)
- Number of windings and coil spacing
- Once known: generate fields, add GPT input line

Known: inner radius 5.25", outer 8.75", full box length 111.75"

### 3. Quad Sizes and Strengths

Lengths and positions in sections 1–8 are currently educated guesses. To fix:
- Physically measure quad lengths
- Measure current-to-field (A → T/m) conversion for each quad
- Update BMAD input file (straightforward once data is obtained)

**Current-to-field conversion is unknown for all quads** — this is critical for any quantitative BMAD results.

### 4. Verify Sections 2–8 in BMAD

GPT results are well-debugged; BMAD sections have not been carefully checked. Things to verify:
- Does changing a parameter actually change the simulation result?
- Is quad polarity correct?
- Does "relative phase = 0" do what the user would expect?
- Does the definition of reference energy in BMAD match expectations?

### 5. Electron Snout *(second highest priority)*

Not simulated at all. High beam transmission through section 8 but poor energy spread will translate into bad transmission in the snout — but quantitatively "how bad" requires simulation. Existing code for the positron snout may be adaptable.

### 6. Comparison to Reality

Compare simulation output to measured beam properties. Should wait until items 1–5 above are substantially complete. Key questions:
- Does the simulation match reality?
- Where does it break down?
- What can be tested to verify it?

---

## Secondary Items

- **Linac sections 2–8:** Change to "traveling wave" cavity type in BMAD/Tao, which affects their focusing. Requires Dave Sagan to add this to BMAD/Tao.
- **BMAD pipe boundaries:**
  - Find a way to extract pipe boundaries from BMAD without parsing nested input files; add to trajectory plots
  - Check whether boundaries actually work in BMAD (behavior is inconsistent)
  - Some pipe radii appear incorrect — needs verification (GPT values are from drawings and are reliable)
- **Port Matlab functionality to Java:**
  - Save/load current GUI state (avoids re-running long simulations)
  - Export graphs or graph data for external plotting programs

---

## Possible Future Extensions

- **Connect to Synchrotron / CESR:** BMAD layouts exist but are probably not directly compatible. Approach: save linac output, start a separate synchrotron simulation.
- **Wakefields / multi-bunch effects / HOMs:** GPT allows field maps to be superimposed, making HOM fields straightforward to add.
- **Break cylindrical symmetry in GPT:** Simulate couplers, kicker magnets, or HOMs with transverse kicks. Worth the added complexity?
