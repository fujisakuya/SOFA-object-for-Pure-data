current:
	echo make pd_darwin, pd_linux, pd_nt, pd_irix5, or pd_irix6

clean: ; rm -f -r conv_MacBinaries conv~.pd_darwin *.o

# ----------------------- Windows-----------------------
# note; you will certainly have to edit the definition of VC to agree with
# whatever you've got installed on your machine:

VC="D:\Program Files\Microsoft Visual Studio\Vc98"

pd_nt: obj1.dll obj2.dll obj3.dll obj4.dll obj5.dll dspobj~.dll

.SUFFIXES: .obj .dll

PDNTCFLAGS = /W3 /WX /DNT /DPD /nologo

PDNTINCLUDE = /I. /I\tcl\include /I..\..\src /I$(VC)\include

PDNTLDIR = $(VC)\lib
PDNTLIB = $(PDNTLDIR)\libc.lib \
	$(PDNTLDIR)\oldnames.lib \
	$(PDNTLDIR)\kernel32.lib \
	..\..\bin\pd.lib 

.c.dll:
	cl $(PDNTCFLAGS) $(PDNTINCLUDE) /c $*.c
	link /dll /export:$*_setup $*.obj $(PDNTLIB)

# override explicitly for tilde objects like this:
dspobj~.dll: dspobj~.c; 
	cl $(PDNTCFLAGS) $(PDNTINCLUDE) /c $*.c
	link /dll /export:dspobj_tilde_setup $*.obj $(PDNTLIB)

# ----------------------- LINUX i386 -----------------------

pd_linux: obj1.l_ia64 obj2.l_ia64 obj3.l_ia64 obj4.l_ia64 \
    obj5.l_ia64 dspobj~.l_ia64

pd_linux32: obj1.l_i386 obj2.l_i386 obj3.l_i386 obj4.l_i386 \
    obj5.l_i386 dspobj~.l_i386

.SUFFIXES: .l_i386 .l_ia64

LINUXCFLAGS = -DPD -O2 -funroll-loops -fomit-frame-pointer \
    -Wall -W -Wshadow -Wstrict-prototypes -Werror \
    -Wno-unused -Wno-parentheses -Wno-switch 
LINUXINCLUDE =  -I/Users/julian/Downloads/pd-0.47-1/src -I/usr/local/include

.c.l_i386:
	cc $(LINUXCFLAGS) $(LINUXINCLUDE) -o $*.o -c $*.c
	ld -shared -o $*.l_i386 $*.o -lc -lm
	strip --strip-unneeded $*.l_i386
	rm $*.o

.c.l_ia64:
	cc $(LINUXCFLAGS) $(LINUXINCLUDE) -fPIC -o $*.o -c $*.c
	ld -shared -o $*.l_ia64 $*.o -lc -lm
	strip --strip-unneeded $*.l_ia64
	rm $*.o

# ----------------------- Mac OSX -----------------------

pd_darwin: conv~.pd_darwin

.SUFFIXES: .pd_darwin


DARWININCLUDE =  -I/Applications/Pd-0.51-2.app/Contents/Resources/src -I/usr/local/include

DARWINCFLAGS = -DPD -O3 -Ofast -Wall -W -Wshadow -Wstrict-prototypes \
    -Wno-unused -Wno-parentheses -Wno-switch -Wno-unknown-pragmas -arch i386 -arch x86_64

DARWINLIBS = -L/usr/local/lib -lfftw3 -lc -lm -bundle -undefined suppress -arch i386 -arch x86_64 -flat_namespace

.c.pd_darwin:
	cc $(DARWINCFLAGS) $(DARWININCLUDE) -o $*.o -c $*.c
	cc $(DARWINLIBS)  -o $*.pd_darwin $*.o 
	rm -f $*.o

	mkdir -p conv_MacBinaries/libs
	cp conv~.pd_darwin conv_MacBinaries
	cp conv~-meta.pd conv_MacBinaries
	cp conv~-help.pd conv_MacBinaries
	./pastyroll-osx.sh libs conv_MacBinaries/conv~.pd_darwin
	zip -r conv_MacBinaries.zip conv_MacBinaries


#-------Documentation-----
doc:
	/Applications/Doxygen.app/Contents/Resources/doxygen  Doxyfile.cnf
	/usr/local/bin/astyle/astyle --style=java --suffix=none --errors-to-stdout *.c
        /usr/local/bin/astyle/astyle --style=java --suffix=none --errors-to-stdout *.h


