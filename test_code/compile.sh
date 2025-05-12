(cd ../src/cpp_src; make clean)
(cd ../src/cpp_src; make -j so)

cp ../src/cpp_src/libnn_mwpc.so .

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/chalmers/users/toliver/phd/nn_mwpc/test_code

export LD_PRELOAD=$LD_PRELOAD:/chalmers/users/toliver/phd/nn_mwpc/test_code/nn_mwpc.so
