# This script runs tests on the nn_mwpc program

# Save the current git hash of the program
# NOTE that this does not check that 
# (i) libnn_mwpc.so is compiled from the same commit 
# (ii) that there are uncommited files.
# This is up to the user to verify
rm test_logfiles/git_hash.out
git rev-parse HEAD >> test_logfiles/git_hash.out


# Compile nn_mwpc and move lib-file to this directory
(cd ../src/cpp_src; make clean)
(cd ../src/cpp_src; make -j so_2)
(cp ../src/cpp_src/libnn_mwpc.so .)


# Run the potential test
(make clean)
(make POT && ./POT)

# Test the nijmegen phase shifts
(make clean)
(make NPHASE && ./NPHASE)

# Test WPC LO
(make clean)
(make WPC && ./WPC)

# Test obaservables with nijmegen potential
(make clean)
(make OBS && ./OBS)

(python3 plot_obs.py DSG)
(python3 plot_obs.py PB)
(python3 plot_obs.py AYY)
(python3 plot_obs.py CKK)
(python3 plot_obs.py DT)
(python3 plot_obs.py PT)
(python3 plot_obs.py AZZ)
