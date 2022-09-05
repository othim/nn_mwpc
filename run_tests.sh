# *****************************************************************************
# *****************************************************************************
# ****** Runnig this script will recompile the code and make the ./obs ********
# ****** excecutable. It will then run tests of the code. *********************
# *****************************************************************************
# *****************************************************************************



# Run tests with the LO WPC potential with Andreas settings and constants
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
(cd src/cpp_src; make clean)
# The compile flag -DANDREAS_CONST will define a variable ANDREAS_CONST
# cor the C++ preprocessor that will take the constants used to produced
# the test data and define change Mn <-> Mp in the rel.kinematics.
(cd src/cpp_src; make obs FLAGS=-DANDREAS_CONST)
(cd src/cpp_src; ./obs WPC_p_all no)
(cd src/cpp_src; ./obs WPC_PB no)

# Run test that checks my code with the Nijmegen potential agains Andreas code
# with the Nijmegen potential. The error corves should be spot on each other.
(cd src/cpp_src; ./obs DSG no_test)
(cd src/python_src; python3 plot_obs2.py DSG)

# Run tests with the LO WPC potential my settings 
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
(cd src/cpp_src; make clean)
(cd src/cpp_src; make obs)
(cd src/cpp_src; ./obs WPC_p_all no)
(cd src/cpp_src; ./obs WPC_PB no)
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# Run  phase shift test with nihjmegen potential. The flag test is telling the
# program that it is a test and just minimal output is printed for easy reading.
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
(cd src/cpp_src; make clean)

# The compile flag -DNIJM_CONST will define a variable NIJM_CONST
# cor the C++ preprocessor that will take the constants use in the 
# Nijmegen potential and define change Mn <-> Mp in the rel.kinematics.
(cd src/cpp_src; make obs FLAGS=-DNIJM_CONST)

(cd src/cpp_src; ./obs phase nijm)
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

(cd src/cpp_src; ./obs DT no_test)
(cd src/python_src; python3 plot_obs.py DT)

(cd src/cpp_src; ./obs PT no_test)
(cd src/python_src; python3 plot_obs.py PT)
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------




