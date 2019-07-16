#############################################################################
# Name:       MD-jeep
#             the Branch & Prune algorithm for the DMDGP - makefile
# Author:     Antonio Mucherino, Leo Liberti, Carlile Lavor, Nelson Maculan
# Sources:    ansi C
# License:    GNU General Public License v.2
# History:    May 01 2010  v.0.1  first release
#############################################################################


OBJ =  main.o  bp.o  matrices.o  printfile.o  pruningtest.o  utils.o

mdjeep: $(OBJ)
	gcc -O3 -o mdjeep $(OBJ) -lm

.c.o:
	gcc -O3 -c $<

