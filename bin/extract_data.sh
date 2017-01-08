# Really primitive, but it works!
# Grep the "Dipole moment" line and provide +1 line of context

grep "Dipole moment" -A1 *.log | tee dipoles.dat
