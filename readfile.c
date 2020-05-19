/*************************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - functions to read input files
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    May 19 2020  v.0.3.2 introduced in this version
*************************************************************************************************************/

#include "bp.h"

size_t errlen = 120;

// this function performs a preliminary analysis on an input text file
// -> the analysis consists in finding out how long every word (wordlen, wrt the given separator) 
//    and every line (linelen) are 
// -> the input needs to be a valid file pointer
// -> the returning value is the number of lines that are found in the file
size_t textFileAnalysis(FILE *input,char sep,size_t *wordlen,size_t *linelen)
{
   char c;
   size_t wlen,llen,nlines;
   unsigned long w,l;

   // moving the file pointer at the beginning of the file
   rewind(input);

   // performing the analysis
   w = 0;  wlen = 0;
   l = 0;  llen = 0;  nlines = 0;
   c = fgetc(input);

   // general loop
   while (c != EOF)
   {
      if (isSeparator(c,sep))
      {
         if (wlen < w)  wlen = w;
         w = 0;  l++;
      }
      else if (c == '\n')
      {
         nlines++;
         if (wlen < w)  wlen = w;
         if (llen < l)  llen = l;
         w = 0;  l = 0;
      }
      else if (c != '\r')
      {
         w++;  l++;
      };
      c = fgetc(input);
   };

   // last word / last line
   if (w != 0)  if (wlen < w)  wlen = w;
   if (l != 0)
   {
      nlines++;
      if (llen < l)  llen = l;
   };

   // ending
   (*wordlen) = wlen + 1;
   (*linelen) = llen + 1;
   return nlines;
};

// this function reads the MDfile and stores the information in the OPTION and INFORMATION structures
// -> the input needs to be a valid file pointer
// -> the returning value is a NULL char on success;
// -> it is otherwise a char pointer to a string describing the error
//   (the string is allocated for errlen characters)
char* readMDfile(FILE *input,OPTION *op,INFORMATION *info)
{
   int last;
   int count;
   size_t nlines,wordlen,linelen;
   char *c,*line;
   char *error;

   // initial check
   if (input == NULL)  return strdup("mdjeep: error: the pointer to the MDfile is invalid");

   // preliminary analysis of text file
   nlines = textFileAnalysis(input,' ',&wordlen,&linelen);
   if (nlines == 0)  return strdup("mdjeep: error: the MDfile seems to be empty");

   // moving the file pointer back to the beginning
   rewind(input);

   // memory allocation
   line = (char*)calloc(linelen+1,sizeof(char));
   error = (char*)calloc(errlen,sizeof(char));

   // initializing the mandatory variables (in info and op, some of them double-checked after reading the file)
   info->filename = NULL;
   info->format = 0UL;
   info->sep = ' ';  // default
   info->start = NULL;
   info->method = -1;
   info->refinement = -1;
   op->r = 5.0;  // default (for bp)
   op->eps = 0.001;  // default (for bp)
   op->maxtime = 3600;  // default (for bp)
   op->maxit = -1;
   op->eta = 0.99;  // default (for spg)
   op->gam = 1.e-4;  // default (for spg)
   op->epsobj = 1.e-7;   // default (for spg)
   op->epsg = 1.e-8;     // default (for spg)
   op->epsalpha = 1.e-12; // default (for spg)
   op->mumin = 1.e-12;    // default (for spg)
   op->mumax = 1.e+12;    // default (for spg)

   // reading the MDfile
   last = -1;  count = 0;
   do
   {
      c = fgets(line,linelen+1,input);
      if (c != NULL)
      {
         count++;
         removEndingChars(c);
         if (strlen(c) > 0)
         {
            if (c[0] != '#')  // lines starting with # are comments
            {
               if (!strncmp(c,"instance",8))
               {
                  last = 0;
                  c = nextColon(c+8);
                  if (c == NULL)
                  {
                     sprintf(error,"mdjeep: error while reading MDfile: key-word 'instance' needs to be followed by ':'");
                     free(line);  return error;
                  };
                  c = nextNonBlank(c+1);
                  if (c == NULL)
                  {
                     sprintf(error,"mdjeep: error while reading MDfile: missing instance name at line %d",count);
                     free(line);  return error;
                  };
                  info->name = strdup(c);  // instance name
               }
               else if (!strncmp(c,"method",6))
               {
                  last = 1;
                  c = nextColon(c+6);
                  if (c == NULL)
                  {
                     sprintf(error,"mdjeep: error while reading MDfile: key-word 'method' needs to be followed by ':'");
                     free(line);  return error;
                  };
                  c = nextNonBlank(c+1);
                  if (c == NULL)
                  {
                     sprintf(error,"mdjeep: error while reading MDfile: no method specified after method key-word at line %d",count);
                     free(line);  return error;
                  };
                  if (!strcmp(c,"bp"))
                  {
                     info->method = 0;  // bp
                  }
                  else if (!strcmp(c,"spg"))
                  {
                     info->method = 1;  // spg
                  }
                  else
                  {
                     sprintf(error,"mdjeep: error while reading MDfile: '%s' is an unknown method",c);
                     free(line);  return error;
                  };
               }
               else if (!strncmp(c,"refinement",10))
               {
                  last = 2;
                  c = nextColon(c+10);
                  if (c == NULL)
                  {
                     sprintf(error,"mdjeep: error while reading MDfile: key-word 'refinement' needs to be followed by ':'");
                     free(line);  return error;
                  };
                  c = nextNonBlank(c+1);
                  if (c == NULL)
                  {
                     sprintf(error,"mdjeep: error while reading MDfile: no refinement method specified after refinement key-word at line %d",count);
                     free(line);  return error;
                  };
                  if (!strcmp(c,"bp"))
                  {
                     info->refinement = 0;  // bp
                  }
                  else if (!strcmp(c,"spg"))
                  {
                     info->refinement = 1;  // spg
                  };
               }
               else if (!strncmp(c,"with",4))
               {
                  if (last == -1)
                  {
                     sprintf(error,"mdjeep: error while reading MDfile: key-word 'with' at line %d does not refer to any previous field",count);
                     free(line);  return error;
                  };
                  c = nextNonBlank(c+4);
                  if (c == NULL)
                  {
                     sprintf(error,"mdjeep: error while reading MDfile: key-word 'with' found at line %d but no attribute specified",count);
                     free(line);  return error;
                  };
                  if (last == 0)  // instance attributes
                  {
                     if (!strncmp(c,"file",4))  // file name
                     {
                        c = nextColon(c+4);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: 'with file' at line %d needs to be followed by ':'",count);
                           free(line);  return error;
                        };
                        c = nextNonBlank(c+1);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: unexpected end of line after 'with file:' at line %d",count);
                           free(line);  return error;
                        };
                        info->filename = strdup(c);
                     }
                     else if (!strncmp(c,"format",6))  // file format
                     {
                        c = nextColon(c+6);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: 'with format' at line %d needs to be followed by ':'",count);
                           free(line);  return error;
                        };
                        c++;
                        info->format = readFormat(c);
                        if (info->format < 6UL)
                        {
                           if (info->format == 0UL)
                              sprintf(error,"mdjeep: error while reading MDfile: format specified at line %d seems to be empty",count);
                           else if (info->format == 1UL)
                              sprintf(error,"mdjeep: error while reading MDfile: specified format at line %d is too long",count);
                           else if (info->format == 2UL)
                              sprintf(error,"mdjeep: error while reading MDfile: multiple use of format elements at line %d",count);
                           else if (info->format == 3UL)
                              sprintf(error,"mdjeep: error while reading MDfile: format element number missing (can be either 1 or 2) at line %d",count);
                           else if (info->format == 4UL)
                              sprintf(error,"mdjeep: error while reading MDfile: unknown format element at line %d",count);
                           else
                              sprintf(error,"mdjeep: error while reading MDfile: elements Id1, Id2, lb and/or ub are missing at line %d",count);
                           free(line);  return error;
                        };
                     }
                     else if (!strncmp(c,"separator",9))  // separator (default is blank char and tab)
                     {
                        c = nextColon(c+9);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: 'with separator' at line %d needs to be followed by ':'",count);
                           free(line);  return error;
                        };
                        c = nextNonBlank(c+1);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: unexpected end of line after 'with separator:' at line %d",count);
                           free(line);  return error;
                        }
                        if (strlen(c) < 3 || c[0] != '\'' || c[2] != '\'')
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: separator needs to be enclosed between two apostrophes (' ')");
                           free(line);  return error;
                        };
                        info->sep = c[1];
                     }
                     else
                     {
                        sprintf(error,"mdjeep: unknown attribute for 'instance' at line %d",count);
                        free(line);  return error;
                     };
                  }
                  else if (last == 1 || last == 2)  // main or refinement method
                  {
                     if (!strncmp(c,"resolution",10))  // resolution (bp)
                     {
                        if (last == 1 && info->method != 0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: resolution is not an attribute of selected method");
                           free(line);  return error;
                        };
                        if (last == 2 && info->refinement != 0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: resolution is not an attribute of refinement method");
                           free(line);  return error;
                        };
                        c = nextColon(c+10);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: 'with resolution' at line %d needs to be followed by ':'",count);
                           free(line);  return error;
                        };
                        c = nextNonBlank(c+1);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: unexpected end of line after 'with resolution:' at line %d",count);
                           free(line);  return error;
                        };
                        if (!isReal(c))
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified resolution at line %d is not a real number",count);
                           free(line);  return error;
                        };
                        op->r = atof(c);
                        if (op->r <= 0.0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified resolution at line %d is non-positive",count);
                           free(line);  return error;
                        };
                     }
                     else if (!strncmp(c,"tolerance",9))  // tolerance (bp)
                     {
                        if (last == 1 && info->method != 0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: tolerance is not an attribute of selected method");
                           free(line);  return error;
                        };
                        if (last == 2 && info->refinement != 0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: tolerance is not an attribute of refinement method");
                           free(line);  return error;
                        };
                        c = nextColon(c+9);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: 'with tolerance' at line %d needs to be followed by ':'",count);
                           free(line);  return error;
                        };
                        c = nextNonBlank(c+1);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: unexpected end of line after 'with tolerance:' at line %d",count);
                           free(line);  return error;
                        };
                        if (!isReal(c))
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified tolerance at line %d is not a real number",count);
                           free(line);  return error;
                        };
                        op->eps = atof(c);
                        if (op->eps <= 0.0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified tolerance at line %d is non-positive",count);
                           free(line);  return error;
                        };
                     }
                     else if (!strncmp(c,"maxtime",7))  // maxtime (bp)
                     {
                        if (last == 1 && info->method != 0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: maxtime is not an attribute of selected method");
                           free(line);  return error;
                        };
                        if (last == 2 && info->refinement != 0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: maxtime is not an attribute of refinement method");
                           free(line);  return error;
                        };
                        c = nextColon(c+7);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: 'with maxtime' at line %d needs to be followed by ':'",count);
                           free(line);  return error;
                        };
                        c = nextNonBlank(c+1);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: unexpected end of line after 'with maxtime:' at line %d",count);
                           free(line);  return error;
                        };
                        if (!isInteger(c))
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified maxtime value at line %d is not given in seconds",count);
                           free(line);  return error;
                        };
                        op->maxtime = atoi(c);
                        if (op->maxtime <= 0.0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified maxtime value at line %d is non-positive",count);
                           free(line);  return error;
                        };
                     }
                     else if (!strncmp(c,"startpoint",10))  // startpoint (spg)
                     {
                        if (info->method != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: startpoint at line %d is not an attribute of the method",count);
                           free(line);  return error;
                        };
                        if (last == 2)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: startpoint at line %d cannot be set up when spg is refinement method",count);
                           free(line);  return error;
                        };
                        c = nextColon(c+10);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: 'with startpoint' at line %d needs to be followed by ':'",count);
                           free(line);  return error;
                        };
                        c = nextNonBlank(c+1);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: unexpected end of line after 'with startpoint:' at line %d",count);
                           free(line);  return error;
                        };
                        info->start = strdup(c);
                     }
                     else if (!strncmp(c,"maxit",5))  // maxit (spg)
                     {
                        if (last == 1 && info->method != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: maxit is not an attribute of selected method");
                           free(line);  return error;
                        };
                        if (last == 2 && info->refinement != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: maxit is not an attribute of refinement method");
                           free(line);  return error;
                        };
                        c = nextColon(c+5);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: 'with maxit' at line %d needs to be followed by ':'",count);
                           free(line);  return error;
                        };
                        c = nextNonBlank(c+1);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: unexpected end of line after 'with maxit:' at line %d",count);
                           free(line);  return error;
                        };
                        if (!isInteger(c))
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified maxit value at line %d is not an integer number",count);
                           free(line);  return error;
                        };
                        op->maxit = atoi(c);
                        if (op->maxit <= 0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified maxit value at line %d is non-positive",count);
                           free(line);  return error;
                        };
                     }
                     else if (!strncmp(c,"eta",3))  // eta (spg)
                     {
                        if (last == 1 && info->method != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: eta is not an attribute of selected method");
                           free(line);  return error;
                        };
                        if (last == 2 && info->refinement != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: eta is not an attribute of refinement method");
                           free(line);  return error;
                        };
                        c = nextColon(c+3);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: 'with eta' at line %d needs to be followed by ':'",count);
                           free(line);  return error;
                        };
                        c = nextNonBlank(c+1);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: unexpected end of line after 'with eta:' at line %d",count);
                           free(line);  return error;
                        };
                        if (!isReal(c))
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified eta value at line %d is not a real number",count);
                           free(line);  return error;
                        };
                        op->eta = atof(c);
                        if (op->eta < 0.80 || op->eta > 1.0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified eta value at line %d is out of bounds [0.8,1.0)",count);
                           free(line);  return error;
                        };
                     }
                     else if (!strncmp(c,"gamma",5))  // gamma (spg)
                     {
                        if (last == 1 && info->method != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: gamma is not an attribute of selected method");
                           free(line);  return error;
                        };
                        if (last == 2 && info->refinement != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: gamma is not an attribute of refinement method");
                           free(line);  return error;
                        };
                        c = nextColon(c+5);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: 'with gamma' at line %d needs to be followed by ':'",count);
                           free(line);  return error;
                        };
                        c = nextNonBlank(c+1);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: unexpected end of line after 'with gamma:' at line %d",count);
                           free(line);  return error;
                        };
                        if (!isReal(c))
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified gamma value at line %d is not a real number",count);
                           free(line);  return error;
                        };
                        op->gam = atof(c);
                        if (op->gam < 0.0 || op->gam >= 1.0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified gamma value at line %d is not valid",count);
                           free(line);  return error;
                        };
                     }
                     else if (!strncmp(c,"epsobj",6))  // epsobj (spg)
                     {
                        if (last == 1 && info->method != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: epsobj is not an attribute of selected method");
                           free(line);  return error;
                        };
                        if (last == 2 && info->refinement != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: epsobj is not an attribute of refinement method");
                           free(line);  return error;
                        };
                        c = nextColon(c+6);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: 'with epsobj' at line %d needs to be followed by ':'",count);
                           free(line);  return error;
                        };
                        c = nextNonBlank(c+1);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: unexpected end of line after 'with epsobj:' at line %d",count);
                           free(line);  return error;
                        };
                        if (!isReal(c))
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified epsobj value at line %d is not a real number",count);
                           free(line);  return error;
                        };
                        op->epsobj = atof(c);
                        if (op->epsobj < 0.0 || op->epsobj >= 1.0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified epsobj value at line %d is not valid",count);
                           free(line);  return error;
                        };
                     }
                     else if (!strncmp(c,"epsg",4))  // epsg (spg)
                     {
                        if (last == 1 && info->method != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: epsg is not an attribute of selected method");
                           free(line);  return error;
                        };
                        if (last == 2 && info->refinement != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: epsg is not an attribute of refinement method");
                           free(line);  return error;
                        };
                        c = nextColon(c+4);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: 'with epsg' at line %d needs to be followed by ':'",count);
                           free(line);  return error;
                        };
                        c = nextNonBlank(c+1);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: unexpected end of line after 'with epsg:' at line %d",count);
                           free(line);  return error;
                        };
                        if (!isReal(c))
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified epsg value at line %d is not a real number",count);
                           free(line);  return error;
                        };
                        op->epsg = atof(c);
                        if (op->epsg < 0.0 || op->epsg >= 1.0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified epsg value at line %d is not valid",count);
                           free(line);  return error;
                        };
                     }
                     else if (!strncmp(c,"epsalpha",8))  // epsalpha (spg)
                     {
                        if (last == 1 && info->method != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: epsalpha is not an attribute of selected method");
                           free(line);  return error;
                        };
                        if (last == 2 && info->refinement != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: epsalpha is not an attribute of refinement method");
                           free(line);  return error;
                        };
                        c = nextColon(c+8);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: 'with epsalpha' at line %d needs to be followed by ':'",count);
                           free(line);  return error;
                        };
                        c = nextNonBlank(c+1);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: unexpected end of line after 'with epsalpha:' at line %d",count);
                           free(line);  return error;
                        };
                        if (!isReal(c))
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified epsalpha value at line %d is not a real number",count);
                           free(line);  return error;
                        };
                        op->epsalpha = atof(c);
                        if (op->epsalpha < 0.0 || op->epsalpha >= 1.0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified epsalpha value at line %d is not valid",count);
                           free(line);  return error;
                        };
                     }
                     else if (!strncmp(c,"mumin",5))  // mumin (spg)
                     {
                        if (last == 1 && info->method != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: mumin is not an attribute of selected method");
                           free(line);  return error;
                        };
                        if (last == 2 && info->refinement != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: mumin is not an attribute of refinement method");
                           free(line);  return error;
                        };
                        c = nextColon(c+5);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: 'with mumin' at line %d needs to be followed by ':'",count);
                           free(line);  return error;
                        };
                        c = nextNonBlank(c+1);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: unexpected end of line after 'with mumin:' at line %d",count);
                           free(line);  return error;
                        };
                        if (!isReal(c))
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified mumin value at line %d is not a real number",count);
                           free(line);  return error;
                        };
                        op->mumin = atof(c);
                        if (op->mumin < 0.0 || op->mumin > 1.0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified mumin value at line %d is not valid",count);
                           free(line);  return error;
                        };
                     }
                     else if (!strncmp(c,"mumax",5))  // mumax (spg)
                     {
                        if (last == 1 && info->method != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: mumax is not an attribute of selected method");
                           free(line);  return error;
                        };
                        if (last == 2 && info->refinement != 1)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: mumax is not an attribute of refinement method");
                           free(line);  return error;
                        };
                        c = nextColon(c+5);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: 'with mumax' at line %d needs to be following by ':'",count);
                           free(line);  return error;
                        };
                        c = nextNonBlank(c+1);
                        if (c == NULL)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: unexpected end of line after 'with mumax:' at line %d",count);
                           free(line);  return error;
                        };
                        if (!isReal(c))
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified mumax value at line %d is not a real number",count);
                           free(line);  return error;
                        };
                        op->mumax = atof(c);
                        if (op->mumax < 1.0)
                        {
                           sprintf(error,"mdjeep: error while reading MDfile: specified mumax value at line %d is not valid",count);
                           free(line);  return error;
                        };
                     }
                     else
                     {
                        if (last == 1)
                           sprintf(error,"mdjeep: unknown attribute for 'method' at line %d",count);
                        else
                           sprintf(error,"mdjeep: unknown attribute for 'refinement' at line %d",count);
                        free(line);  return error;
                     };
                  }
                  else
                  {
                     sprintf(error,"mdjeep: error while reading MDfile: internal error");
                     free(line);  return error;
                  };
               }
               else
               {
                  sprintf(error,"mdjeep: error while reading MDfile: syntax error at line %d",count);
                  free(line);  return error;
               };
            };
         };
      };
   }
   while (c != NULL);

   // have all mandatory information been loaded? are they coherent?
   if (info->filename == NULL)
   {
      sprintf(error,"mdjeep: error while reading MDfile: instance file name not specified in the MDfile");
      free(line);  return error;
   };
   if (info->format == 0UL)
   {
      sprintf(error,"mdjeep: error while reading MDfile: file format not specified");
      free(line);  return error;
   };
   if (info->method == -1)
   {
      sprintf(error,"mdjeep: error while reading MDfile: main method not specified (can be 'bp' or 'spg')");
      free(line);  return error;
   };
   if (info->method == 0 && info->refinement == 0)
   {
      sprintf(error,"mdjeep: error while reading MDfile: bp cannot be invoked as a refinement method for itself");
      free(line);  return error;
   };
   if (info->method == 1 && info->refinement == 0)
   {
      sprintf(error,"mdjeep: error while reading MDfile: spg cannot use bp as a refinement method");
      free(line);  return error;
   };
   if (info->method == 1 && info->refinement == 1)
   {
      sprintf(error,"mdjeep: error while reading MDfile: spg cannot be invoked as a refinement method for itself");
      free(line);  return error;
   };
   if (info->method == 1 && info->start == NULL)
   {
      sprintf(error,"mdjeep: error while reading MDfile: startpoint attribute not set up, impossible to run spg without starting point");
      free(line);  return error;
   };
   if (info->method == 1 && op->maxit == -1)
   {
      sprintf(error,"mdjeep: error while reading MDfile: maxit attribute needs to be specified when spg is the main method");
      free(line);  return error;
   };
   if ((info->method == 1 || info->refinement == 1) && (op->mumin >= op->mumax))
   {
      sprintf(error,"mdjeep: error while reading MDfile: mumin is greater than or equal to mumax");
      free(line);  return error;
   };

   // ending
   free(line);  free(error);
   return NULL;
};

// verifying whether a text file contains a valid list of distances (instance file)
// (it is verified that all entry types are constant line by line)
// -> the input needs to be a valid file pointer, sep is the separator in the text file
// -> typelist is the binary encoding of the format (if returning value is true, it is 0UL otherwise)
// -> memory is a pre-allocated char array which can entirely hold all words in the distance file 
//    (see previous function, msize is the size of this memory, which is the longest line length detected in the file)
// -> the returning value is a boolean indicating whether the distance file is valid or not
bool isDistanceFileValid(FILE *input,char sep,unsigned long *typelist,size_t msize,char *memory)
{
   unsigned long t1,t2;

   // moving the file pointer at the beginning of the file
   rewind(input);

   // resetting typelist
   (*typelist) = 0UL;

   // reading the first (nonempty) line of the file
   t1 = 0UL;
   while (t1 == 0UL && memory != NULL)
   {
      memory = fgets(memory,msize,input);
      if (memory == NULL)  return false;
      removEndingChars(memory);
      t1 = detectTypes(memory,sep);
   };
   if (t1 == 0UL)  return false;

   // comparing the first line to the others
   while (memory != NULL)
   {
      memory = fgets(memory,msize,input);
      if (memory != NULL)
      {
         removEndingChars(memory);
         t2 = detectTypes(memory,sep);
         if (t2 != 0UL)  if (t1 != t2)  return false;
      };
   };

   // the file is valid
   (*typelist) = t1;  // list of types encoded in binary format
   return true;
};

// reading the file format specified through an array of characters
// -> the array of char c needs to be a valid pointer (with allocated memory)
// -> the binary output is the format specified with the following code for the format elements: 
//        011? = Id, 100? = groupId, 101? = Name, 110? = groupName, ? indicates vertex 0 or 1
//        1110 = lb, 1111 = ub, 0000 = ignore
//    an error code is returned in the following cases:
//         0UL = the char array is empty
//        0001 = the format is too long (exceeds the bit size 8*sizeof(unsigned long))
//        0010 = the same format element appears more than once in the format
//        0011 = element vertex number is missing or not corresponding to neither 1 or 2
//        0100 = unknown format element in the char array
//        0101 = one of the following elements is missing: id1, id2, lb and/or ub
unsigned long readFormat(const char *c)
{
   int i,k,maxk;
   bool id1,groupid1,name1,groupname1;
   bool id2,groupid2,name2,groupname2;
   bool lb,ub;
   size_t len;
   unsigned long cf;
   unsigned long format = 0UL;

   // nothing can be done if the char array is empty
   if (isLastChar(c[0]))  return format;

   // computing the length of the char array
   len = strlen(c) + 1;  // we also count '\0'
   if (len == 1)  return format;

   // initializing control variables
   id1 = false;  groupid1 = false;  name1 = false;  groupname1 = false;
   id2 = false;  groupid2 = false;  name2 = false;  groupname2 = false;
   lb = false;  ub = false;

   // computing the maximum number of bits
   maxk = 8*sizeof(unsigned long);

   // reading the char array and coding the information in binary format
   i = 0;  k = 0;
   while (!isLastChar(c[i]))
   {
      if (isSeparator(c[i],' ') || isNewLineDelimiter(c[i]))
      {
         i++;
      }
      else if (!strncasecmp(c+i,"id",2))
      {
         cf = 6UL;  // 011?
         if (c[i+2] == '1')
         {
            if (id1)  return 2UL;
            id1 = true;
         }
         else if (c[i+2] == '2')
         {
            cf = cf + 1UL;
            if (id2)  return 2UL;
            id2 = true;
         }
         else return 3UL;
         k = k + 4;
         if (k > maxk)  return 1UL;
         format = (format << 4) | cf;
         i = i + 3;
      }
      else if (!strncasecmp(c+i,"groupid",7))
      {
         cf = 8UL;  // 100?
         if (c[i+7] == '1')
         {
            if (groupid1)  return 2UL;
            groupid1 = true;
         }
         else if (c[i+7] == '2')
         {
            cf = cf + 1UL;
            if (groupid2)  return 2UL;
            groupid2 = true;
         }
         else return 3UL;
         k = k + 4;
         if (k > maxk)  return 1UL;
         format = (format << 4) | cf;
         i = i + 8;
      }
      else if (!strncasecmp(c+i,"name",4))
      {
         cf = 10UL;  // 101?
         if (c[i+4] == '1')
         {
            if (name1)  return 2UL;
            name1 = true;
         }
         else if (c[i+4] == '2')
         {
            cf = cf + 1UL;
            if (name2)  return 2UL;
            name2 = true;
         }
         else return 3UL;
         k = k + 4;
         if (k > maxk)  return 1UL;
         format = (format << 4) | cf;
         i = i + 5;
      }
      else if (!strncasecmp(c+i,"groupname",9))
      {
         cf = 12UL;  // 110?
         if (c[i+9] == '1')
         {
            if (groupname1)  return 2UL;
            groupname1 = true;
         }
         else if (c[i+9] == '2')
         {
            cf = cf + 1UL;
            if (groupname2)  return 2UL;
            groupname2 = true;
         }
         else return 3UL;
         k = k + 4;
         if (k > maxk)  return 1UL;
         format = (format << 4) | cf;
         i = i + 10;
      }
      else if (!strncasecmp(c+i,"lb",2))
      {
         cf = 14UL;  // 1110
         if (lb)  return 2UL;
         lb = true;
         k = k + 4;
         if (k > maxk)  return 1UL;
         format = (format << 4) | cf;
         i = i + 2;
      }
      else if (!strncasecmp(c+i,"ub",2))
      {
         cf = 15UL;  // 1111
         if (ub)  return 2UL;
         ub = true;
         k = k + 4;
         if (k > maxk)  return 1UL;
         format = (format << 4) | cf;
         i = i + 2;
      }
      else if (!strncasecmp(c+i,"ignore",6))
      {
         cf = 0UL;  // 0000
         k = k + 4;
         if (k > maxk)  return 1UL;
         format = (format << 4) | cf;
         i = i + 6;
      }
      else
      {
         return 4UL;
      };
   };

   // any missing label among the necessary ones?
   if (!id1 || !id2 || !lb || !ub)  return 5UL;

   // the format is correct
   return format;
};

// verifying the index range for the vertices in the distance list
// -> the input needs to be a valid file pointer, sep is the separator
// -> format is the expected input line format in binary
// -> memory is a char array of length msize, pre-allocated and able to contain an entire line of the input file
// -> n0 is the smallest identifier found in the file (output, pointer)
// -> the returning value is the number of identified vertices (it is 0 if an error occurs)
int numberOfVerticesInFile(FILE *input,char sep,unsigned long format,int *n0,size_t msize,char *memory)
{
   int i,id;
   int nmin,nmax;
   int nformat,nf;
   char tmp,*pointer;
   unsigned long f,cf;
   size_t l;

   // moving the file pointer at the beginning of the file
   rewind(input);

   // counting the bits used for the format
   f = format;
   nformat = 0;
   while (f != 0UL)
   {
      nformat = nformat + 4;
      f = f >> 4;
   };
   if (nformat < 4)  return 0;

   // reading the identifiers in the file
   nmin = -1;  nmax = 0;
   do
   {
      memory = fgets(memory,msize,input);
      if (memory != NULL)
      {
         l = removEndingChars(memory);
         if (l > 0)
         {
            i = 0;  nf = nformat - 4;
            pointer = memory + i;
            while (i <= l + 1)
            {
               if (isLastChar(memory[i]) || isSeparator(memory[i],sep))
               {
                  if (!isLastChar(pointer[0]))
                  {
                     cf = (format >> nf) & 15UL;
                     tmp = memory[i];  memory[i] = '\0';
                     if (!isLastChar(pointer[0]))
                     {
                        if (cf == 6UL || cf == 7UL)
                        {
                           if (!isInteger(pointer))  return 0;
                           id = atoi(pointer);
                           if (nmax < id)  nmax = id;
                           if (nmin == -1)  nmin = nmax;
                           if (nmin > id)  nmin = id;
                        };
                        nf = nf - 4;
                     };
                     memory[i] = tmp;  i++;
                     while (!isLastChar(memory[i]) && isSeparator(memory[i],sep))  i++;
                     pointer = memory + i;
                  }
                  else break;
               }
               else i++;
            };
            if (nf >= 0 && nf != nformat - 4)  return 0;
         };
      };
   }
   while (memory != NULL);

   // ending
   if (nmin == -1)  return 0;
   (*n0) = nmin;
   return nmax - nmin + 1;
};

// reading instance file (distance list) in specified input format
// -> the input needs to be a valid file pointer, sep is the separator
// -> n is the total number of vertices, n0 is the smallest vertex rank
// -> v is the array of type VERTEX (output, memory needs to be pre-allocated)
// -> format is the specified line format in binary
// -> memory is a char array of length msize, pre-allocated and able to contain an entire line of the input file
// -> returning value is: -1 on succeed
//                        -2 if a format error occurred
//                        -3 the presence of a distance from a vertex to itself is detected
//                        -4 some vertex ranks are missing in the file
//                        -5 some lower bounds are larger than the upper bounds
//                   0 <= id if this vertex was found for the second time in the file but with different attributes 
//                           (groupid,name,groupname)
int readDistanceFile(FILE *input,char sep,int n,int n0,VERTEX *v,unsigned long format,size_t msize,char *memory)
{
   int i,j,k;
   int nf,nformat;
   int id1,id2,gid1,gid2;
   char *name1,*name2,*gname1,*gname2;
   char tmp,*pointer;
   size_t l;
   double lb,ub;
   unsigned long f,cf;

   // moving the file pointer at the beginning of the file
   rewind(input);

   // counting the bits used for the format
   f = format;
   nformat = 0;
   while (f != 0UL)
   {
      nformat = nformat + 4;
      f = f >> 4;
   };
   if (nformat < 4)  return -2;

   // initializing all vertices with Id = -1 (which means: not defined)
   for (i = 0; i < n; i++)  v[i].Id = -1;

   /* reading the file */
   do
   {
      // reading the next line of the file
      memory = fgets(memory,msize,input);
      if (memory != NULL)
      {
         // variable initialization
         id1 = -1;  id2 = -1;
         gid1 = 0;  gid2 = 0;
         name1 = strdup("(no name)");
         name2 = strdup("(no name)");
         gname1 = strdup("(no group name)");
         gname2 = strdup("(no group name)");
         lb = -1.0;  ub = -1.0;

         // cleaning the string ending chars
         l = removEndingChars(memory);
         if (l > 0)
         {
            i = 0;  nf = nformat - 4;
            pointer = memory + i;
            while (i <= l + 1)
            {
               if (isLastChar(memory[i]) || isSeparator(memory[i],sep))
               {
                  if (!isLastChar(pointer[0]))
                  {
                     cf = (format >> nf) & 15UL;
                     tmp = memory[i];  memory[i] = '\0';
                     if (!isLastChar(pointer[0]))
                     {
                        if (cf == 6UL)
                        {
                           if (!isInteger(pointer))  return -2;
                           id1 = atoi(pointer);
                        }
                        else if (cf == 7UL)
                        {
                           if (!isInteger(pointer))  return -2;
                           id2 = atoi(pointer);
                        }
                        else if (cf == 8UL)
                        {
                           if (!isInteger(pointer))  return -2;
                           gid1 = atoi(pointer);
                        }
                        else if (cf == 9UL)
                        {
                           if (!isInteger(pointer))  return -2;
                           gid2 = atoi(pointer);
                        }
                        else if (cf == 10UL)
                        {
                           free(name1);
                           name1 = strdup(pointer);
                        }
                        else if (cf == 11UL)
                        {
                           free(name2);
                           name2 = strdup(pointer);
                        }
                        else if (cf == 12UL)
                        {
                           free(gname1);
                           gname1 = strdup(pointer);
                        }
                        else if (cf == 13UL)
                        {
                           free(gname2);
                           gname2 = strdup(pointer);
                        }
                        else if (cf == 14UL)
                        {
                           if (!isReal(pointer))  return -2;
                           lb = atof(pointer);
                        }
                        else if (cf == 15UL)
                        {
                           if (!isReal(pointer))  return -2;
                           ub = atof(pointer);
                        }
                        else if (cf != 0UL)
                        {
                           return false;
                        };
                        nf = nf - 4;
                     };
                     memory[i] = tmp;  i++;
                     while (!isLastChar(memory[i]) && isSeparator(memory[i],sep))  i++;
                     pointer = memory + i;
                  }
                  else break;
               }
               else i++;
            };

            /* inserting the data in the VERTEX array */
            if (id1 != -1 && id2 != -1 && lb != -1.0 && ub != -1.0)
            {
               // inserting vertex with Id1
               i = id1 - n0;
               if (v[i].Id == -1)
                  initVertex(&v[i],id1,gid1,name1,gname1);
               else
                  if (v[i].groupId != gid1 || strcmp(v[i].Name,name1) || strcmp(v[i].Group,gname1))  return id1;

               // inserting vertex with Id2
               j = id2 - n0;
               if (v[j].Id == -1)
                  initVertex(&v[j],id2,gid2,name2,gname2);
               else
                  if (v[j].groupId != gid2 || strcmp(v[j].Name,name2) || strcmp(v[j].Group,gname2))  return id2;

               // inserting the distance values
               if (i > j)
               {
                  k = i;  i = j;  j = k;
               }
               else if (i == j)  return -3;
               if (lb > ub)  return -5;
               if (v[j].ref == NULL)
                  v[j].ref = initReference(i,lb,ub);
               else if (getReference(v,i,j) == NULL)
                  addDistance(v[j].ref,i,lb,ub);

               // deallocating memory allocated by strdup
               free(name1);  free(name2);
               free(gname1);  free(gname2);
            };
         };
      };
   }
   while (memory != NULL);

   // verifying that all vertices were defined
   for (i = 0; i < n; i++)  if (v[i].Id == -1)  return -4;

   // ending
   return -1;
};

// reading a 3D conformation from a text file (fixed format, list of triplets of coordinates)
// -> the 3D conformation is the starting point necessary for SPG when invoked as main method
// -> it is supposed that the memory has been allocated for X to keep at least n triplets of coordinates
// -> the returning value is the number of loaded triplets of coordinates
int readStartingPoint(FILE *input,int n,double **X)
{
   int i;

   // moving the file pointer at the beginning of the file
   rewind(input);

   // reading the coordinates from the file
   i = 0;
   while (i < n && EOF != fscanf(input,"%lf%lf%lf",&X[0][i],&X[1][i],&X[2][i]))  i++;

   // ending
   return i;
};

