##############################################################################################
# Name:       MD-jeep
#             the Branch & Prune algorithm for discretizable Distance Geometry - main program
# Author:     A. Mucherino, L. Liberti, D.S. Goncalves, C. Lavor, N. Maculan
# Sources:    ansi C
# License:    GNU General Public License v.2
# History:    May 01 2010  v.0.1  first release
#             May 10 2014  v.0.2  second release (no files added)
##############################################################################################


OBJ =  main.o  bp.o  matrices.o  printfile.o  pruningtest.o  utils.o

mdjeep: $(OBJ)
	gcc -O3 -o mdjeep $(OBJ) -lm

.c.o:
	gcc -O3 -c $<

