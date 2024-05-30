(cd ../src/cpp_src; make clean)
(cd ../src/cpp_src; make -j so)

cp ../src/cpp_src/libnn_mwpc.so .
