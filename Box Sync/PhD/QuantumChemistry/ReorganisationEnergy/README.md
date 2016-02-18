# Calculate Reorganisation Energy 
     ''This follows the roughly nomenclature and working of JKP's thesis. 'we follow the method of Sakanou' - JKP Thesis 2.2 Reorganisation Energy. Instrusctions originally written by Florian Steiner and edited by Beth Rice''



1. Run 'polaron_reorg_optimisation_jobs_from_chk.sh' on a checkpoint for an arb molecule
--> generates two geometry opt jobs (neutral + charged)

2. Once complete, run 'polaron_reorg_energy_jobs_from_geom_chks.sh' on the two checkpoints returned.
--> generates 4 energy jobs "*_ion_E.chk" / "_neutral_E.chk"

3. Once everything's returned, run 'calc_reorg_energy.sh' on the root of the file names to calculate your energy in eV.
