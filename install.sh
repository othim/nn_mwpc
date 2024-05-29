# Running this script with 'bash install.sh' will compile the code with the
# pybind interface activated (flag PY_BIND=-DPYBIND). It will then remove
# if there is an old nn_mwpc.so file and then copy the newly compiled .so
# file to this directory and install the python library.

(cd src/cpp_src; make clean)

(cd src/cpp_src; make -j so PY_BIND=-DPYBIND)

rm libnn_mwpc.so

cp src/cpp_src/libnn_mwpc.so .
mv libnn_mwpc.so nn_mwpc.so

pip install .

git rev-parse HEAD > nn_mwpc_install_hash.txt

python3 python_module/run_tests_on_module.py
