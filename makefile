################################################################################################################
# Name:       MD-jeep
#             the Branch & Prune algorithm for discretizable Distance Geometry - makefile
# Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
# Sources:    ansi C
# License:    GNU General Public License v.3
# History:    May 01 2010  v.0.1    first release
#             May 10 2014  v.0.2    second release (no files added)
#             Jul 28 2019  v.0.3.0  third release (added files vertex.c, distance.c, objfun.c, spg.c, splitime.c)
#             Mar 21 2020  v.0.3.1  no new files added
#################################################################################################################


OBJ= main.o bp.o vertex.o distance.o matrices.o pruningtest.o objfun.o spg.o utils.o printfile.o

mdjeep: $(OBJ) splitime.o
	gcc -O3 -o mdjeep $(OBJ) splitime.o -lm

splitime.o: splitime.c
	gcc -O2 -std=c99 -c splitime.c

.c.o:
	gcc -O3 -c $<

clean:
	\rm *.o

