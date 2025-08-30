MYDIR = ./
FCOMPILER = gfortran -c -O3 
OBJS = burgersWENOEq.o fftpackt.o

burgersWENOEq:	$(OBJS)
		gfortran -O3 -o solveBurgersWENOEq $(OBJS) 

burgersWENOEq.o: 	$(MYDIR)burgersWENOEq.f90
		$(FCOMPILER) $(MYDIR)burgersWENOEq.f90

fftpackt.o:			$(MYDIR)fftpackt.f90
		$(FCOMPILER) $(MYDIR)fftpackt.f90
