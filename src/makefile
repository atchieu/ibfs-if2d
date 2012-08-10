#Makefile for ibfs on neumann.caltech.edu
#Colonius

IBFSFILES = parameters.o grid.o myfft.o variables.o operators.o main.o
POSTFILES = parameters.o grid.o myfft.o variables.o operators.o post.o
COMPILE = g95
FLAGS1 = -O4 
FLAGS2 = -lfftw3 -o
FLAGS3 = -I/home/atchieu/Install/include -L/home/atchieu/Install/lib
# FLAGS3 = -I/usr/local/include -L/usr/local/lib

all:   	$(IBFSFILES) $(POSTFILES) 
	$(COMPILE)  $(IBFSFILES) $(FLAGS1) $(FLAGS2) ../bin/ibfs $(FLAGS3)
	$(COMPILE)  $(POSTFILES) $(FLAGS1) $(FLAGS2) ../bin/ibfspost $(FLAGS3)

ibfs: $(IBFSFILES) 
	$(COMPILE)  $(IBFSFILES) $(FLAGS1) $(FLAGS2) ../bin/ibfs $(FLAGS3)

post: $(POSTFILES) 
	$(COMPILE)  $(POSTFILES) $(FLAGS1) $(FLAGS2) ../bin/ibfspost $(FLAGS3)

clean:
	/bin/rm -f *.o *.mod core* *~

%.o : %.f90
	$(COMPILE) $(FLAGS1) $(FLAGS3) -c $<
