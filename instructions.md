Your task is to refactor the code used in ../Cornell. Structure it as follows:

- config
   - cathode.yaml
   - gun.yaml
   - injector.yaml
   - linac1.yaml
   - linac2.yaml
   - linac3.yaml
   - linac4-8.yaml
- sim
   - main.py
   - cathode.py
   - gun.py
   - injector.py
   - linac1-3.py
   - linac4-8.py
   - helpers
      - buildfields.py
      - loadparticles.py
      - tqdmwrapper.py
      - tools.py  # misc, feel free to add more files in this directory
   - plot
      - cathode.py
      - gun.py
      - injector.py
      - linac1-3.py
      - linac4-8.py
- logs
   - diags
      - cathode
      - gun
      - injector
      - linac1-3
      - linac4-8
   - plots
      - cathode
      - gun
      - injector
      - linac1-3
      - linac4-8
   - pipeline
      - log_[date].log
- docs
   - cathode.md
   - gun.md
   - injector.md
   - linac1-3.md
   - linac4-8.md
- fieldmaps
   - gdf
      - [import from ../Cornell/fieldmaps]
   - h5
      - [these will be built by sim/helpers/buildfields.py]
- README.md
- requirements.txt

**Do**
- Hardcode as many options as possible in the YAMLs so that they can be adjusted by me later. You may calculate values (such as phases) in the scratchpad, but then hardcode the values.
- Make sure the produced figures are equivalent to the old version.
- Make code easy to read. Helper functions that don't relate to the physics should go in sim/helpers. 

**Do not**
- List exact output numbers (e.g. 'beam accelerates to xx MeV with xx emittance') in the docs. 
- Overly complicate code. Comment important sections but be succinct and only do it when absolutely necessary. Most information should go in the README rather than the code.

You may make substantial design changes from the old repository, but ask first.
