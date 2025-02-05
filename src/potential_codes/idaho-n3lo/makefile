# the compiler to use
CC  := g++
CXX := g++-9
XLF := gfortran

SRCFORT_90_FILES := idaho-n3lo-interface.f90
SRCFORT_77_FILES := idaho-n3lo.f
OBJFORT_90_FILES := idaho-n3lo-interface.o
OBJFORT_77_FILES := idaho-n3lo.o

#On mac
LDFLAGS := -lgfortran -L/usr/local/Cellar/gcc/11.2.0/lib/gcc/11/
# On linux
#LDFLAGS := -lgfortran 

FORTFLAGS_90 := -O2 -fdefault-real-8 -fdefault-double-8 -cpp -ffree-form -ffree-line-length-1000 -fPIC -std=gnu
FORTFLAGS_77 := -O2 -fdefault-real-8 -fdefault-double-8 -fPIC -static  -std=gnu

run: $(OBJFORT_90_FILES) $(OBJFORT_77_FILES)
	#$(CC) $^ -o $@ $(LDFLAGS) $(LDLIBS) 

idaho-n3lo-interface.o: idaho-n3lo-interface.f90
	$(XLF) $(FORTFLAGS_90) -c $< -o $@
idaho-n3lo.o: idaho-n3lo.f
	$(XLF) $(FORTFLAGS_77) -c $< -o $@

clean:
	rm $(OBJFORT_90_FILES) $(OBJFORT_77_FILES)
