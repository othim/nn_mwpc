
(cd src/potential_codes/cdbonn; make cleanall)

(cd src/potential_codes/cdbonn; make)

(cd src/potential_codes/cdbonn; cp cdbonn_interface.o ../cpp_src)
(cd src/potential_codes/cdbonn; cp cdbonn.o ../cpp_src)

