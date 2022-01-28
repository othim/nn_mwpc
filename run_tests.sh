# Runnig this script will recompile the code and make the ./obs
# excecutable. It will then run tests of the code.

(cd src/cpp_src; make clean)

(cd src/cpp_src; make obs)

# Run  phase shift test with nihjmegen potential. The flag test is telling the
# program that it is a test and just minimal output is printed for easy reading.
(cd src/cpp_src; ./obs phase test)
