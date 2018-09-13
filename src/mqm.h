/**********************************************************************
 *
 * mqm.h
 *
 * Copyright (c) 1996-2011 by
 * Ritsert C Jansen, Danny Arends, Pjotr Prins and Karl W Broman
 *
 * initial MQM C code written between 1996-2002 by Ritsert C. Jansen
 * improved for the R-language by Danny Arends, Pjotr Prins and Karl W. Broman
 *
 * Modified by Pjotr Prins and Danny Arends
 * last modified Feb 2011
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU General Public License,
 *     version 3, as published by the Free Software Foundation.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but without any warranty; without even the implied warranty of
 *     merchantability or fitness for a particular purpose.  See the GNU
 *     General Public License, version 3, for more details.
 *
 *     A copy of the GNU General Public License, version 3, is available
 *     at http://www.r-project.org/Licenses/GPL-3
 *
 * C functions for the R/qtl package
 *
 **********************************************************************/


#ifndef __MQM_H
  #define __MQM_H

  #include <R.h>
  #include <R_ext/Utils.h>
  #include "standalone.h"
  #include "util.h"
  #include "mqmdatatypes.h"
  #include "mqmprob.h"
  #include "mqmmixture.h"
  #include "mqmregression.h"
  #include "mqmaugment.h"
  #include "mqmeliminate.h"
  #include "mqmmapqtl.h"
  #include "mqmscan.h"

  #ifdef STANDALONE
    // Running mqm stand alone (without R)
    extern FILE* redirect_info;  // Redirect output for testing
    extern int debuglevel;  // Redirect output for testing

    #define message(type, format, ...) { \
      fprintf(redirect_info,"%s: ",type); \
      fprintf(redirect_info, format, ## __VA_ARGS__); \
      fprintf(redirect_info,"\n"); }

    #define fatal(s, ...) { message("FATAL",s, ## __VA_ARGS__); exit(127); }

    #define debug_trace(format, ...) { \
      if(debuglevel > 0){ \
        fprintf(redirect_info,"TRACE "); \
        fprintf(redirect_info,"%s %d:",__FILE__,__LINE__); \
        fprintf(redirect_info,format, ## __VA_ARGS__); \
      }\
    }
  #else
    #ifdef ENABLE_C99_MACROS
      #define message(type, format, ...) { Rprintf(format, ## __VA_ARGS__);Rprintf("\n");}
      #define fatal(s, ...) { message("FATAL",s, ## __VA_ARGS__); Rf_error(s); }
      #define debug_trace(format, ...) { }
    #else
      void message(const char*, ...);
      void fatal(const char*, ...);
      void debug_trace(const char*, ...);
    #endif // ENABLE_C99_MACROS
  #endif // !STANDALONE

  #ifdef ENABLE_C99_MACROS
    #define info(format, ...) { message("INFO",format, ## __VA_ARGS__); }
    #define verbose(format, ...) if (verbose) { info(format, ## __VA_ARGS__); }
  #else
    void info(const char*, ...);
    void verbose(const char*, ...);
  #endif // !ENABLE_C99_MACROS

#endif // MQM_H
