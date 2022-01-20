# Running this script with 'bash install.sh' will compile the code with the
# pybind interface activated (flag PY_BIND=-DPYBIND). It will then remove
# if there is an old nn_mwpc.so file and then copy the newly compiled .so
# file to this directory and install the python library.

(cd src/cpp_src; make clean)

(cd src/cpp_src; make so PY_BIND=-DPYBIND)

rm nn_mwpc.so

cp src/cpp_src/nn_mwpc.so .

pip install .

