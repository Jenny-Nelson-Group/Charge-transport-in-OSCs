for dir in C60 M B  B-isomers    T  T-isomers
do
 cd "${dir}"

 #steepest descent minim
 grompp -f ../steep.mdp -c *.pdb -p *.top
 mdrun

 #now the md...
 grompp -f ../md.mdp -c confout.gro -p *.top
 mdrun

 # this md.mdp is setup to do simulated annealing with a cycle time of 1ns
 #  .'. dump frames with dt=1000ps to get all the 'room temp' structures
 echo 0 | trjconv -f traj.xtc -dt 1000 -sep -nzero 3 -o ${dir}_SimulAnneal_.pdb
 # nb: note the rediculous echo 0 to select 'system'

 cd -
done
