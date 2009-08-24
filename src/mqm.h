/**********************************************************************
 *
 * mqm.h
 *
 * Master include file for common mqm.headers
 *
 * copyright (c) 2009 Ritsert Jansen, Danny Arends, Pjotr Prins and Karl Broman
 *
 * last modified July, 2009
 * first written July, 2009
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
 **********************************************************************/


#ifndef __MQM_H
  #define __MQM_H

  #include <R.h>
  // #include <R_ext/PrtUtil.h>
  // #include <R_ext/RS.h> 
  // #include <R_ext/Utils.h>
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
  extern FILE* redirect_info;  // Redirect output for testing
  #define message(type, format, ...) { \
    fprintf(redirect_info,"%s: ",type); \
    fprintf(redirect_info, format, ## __VA_ARGS__); \
    fprintf(redirect_info,"\n"); } 
  // #define warning(s) { message("WARNING",s); }
  #define fatal(s) { message("FATAL",s); exit(127); }
#else
  #define message(type,s) { Rprintf("%s: %s\n",type,s); }
  // #define warning(s) { Rf_warning(s); }
  #define fatal(s) { message("FATAL",s); Rf_error(s); }
#endif

#ifdef NDEBUG
  #define info(format, ...) 
#else
  #define info(format, ...) { message("INFO",format, ## __VA_ARGS__); }
  #define verbose(format, ...) if (verbose) { info(format, ## __VA_ARGS__); }
#endif

#endif // MQM_H
