MYDIR = ./
FCOMPILER = gfortran -c -O3 
OBJS = burgersEq.o fftpackt.o

burgersEq:	$(OBJS)
		gfortran -O3 -o solveBurgersEq $(OBJS) 

burgersEq.o: 	$(MYDIR)burgersEq.f90
		$(FCOMPILER) $(MYDIR)burgersEq.f90

fftpackt.o:	$(MYDIR)fftpackt.f90
		$(FCOMPILER) $(MYDIR)fftpackt.f90
