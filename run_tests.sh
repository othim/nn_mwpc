# Runnig this script will recompile the code and make the ./obs
# excecutable. It will then run tests of the code.

(cd src/cpp_src; make clean)

(cd src/cpp_src; make obs)

# Run  phase shift test with nihjmegen potential. The flag test is telling the
# program that it is a test and just minimal output is printed for easy reading.
(cd src/cpp_src; ./obs phase test)

(cd src/cpp_src; ./obs test_1S0 no_test)

(cd src/cpp_src; ./obs DIAG test)

(cd src/cpp_src; ./obs SGT  no_test)
(cd src/cpp_src; ./obs SGTT no_test)
(cd src/cpp_src; ./obs SGTL no_test)

# Run tests for observables. Run Test of DSG to check also with Andreas code.
(cd src/cpp_src; ./obs DSG no_test)
(cd src/python_src; python3 plot_obs2.py DSG)

(cd src/cpp_src; ./obs PB no_test)
(cd src/python_src; python3 plot_obs.py PB)

(cd src/cpp_src; ./obs AYY no_test)
(cd src/python_src; python3 plot_obs.py AYY)

(cd src/cpp_src; ./obs CKK no_test)
(cd src/python_src; python3 plot_obs.py CKK)

