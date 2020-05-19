/***************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - splitime
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    Jul 28 2019  v.0.3.0  introduced in this version
              Mar 21 2020  v.0.3.1  no changes
              May 19 2020  v.0.3.2  no changes
****************************************************************************************************/

#include "bp.h"

// separating the given time (the lap between the two given input arguments) 
// in hours, minutes, seconds, milliseconds, etc
char* splitime(struct timeval start,struct timeval end)
{
   int flag;
   int time_s,time_ms;
   char *timestr,*timetmp;

   time_s = end.tv_sec - start.tv_sec;
   time_ms = end.tv_usec - start.tv_usec;

   timestr = (char*)calloc(200,sizeof(char));
   timetmp = (char*)calloc(200,sizeof(char));

   flag = 0;  sprintf(timestr," ");
   while (time_ms > 0.0)
   {
      if (flag == 0)  sprintf(timetmp," %3d\u03BCs",time_ms%1000);
                else  sprintf(timetmp," %3dms",time_ms);
      time_ms = time_ms / 1000.0;
      strcat(timetmp,timestr);
      strcpy(timestr,timetmp);
      flag++;
   };

   flag = 0;
   while(time_s > 0.0)
   {
      if (flag == 0)  sprintf(timetmp," %2ds",time_s%60);
 else if (flag == 1)  sprintf(timetmp," %2dm",time_s%60);
                else  sprintf(timetmp," %2dh",time_s%60);
      strcat(timetmp,timestr);
      strcpy(timestr,timetmp);
      time_s = time_s / 60.0;
      flag++;
   };

   free(timetmp);

   return timestr;
};

