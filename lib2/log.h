/**
 * Logging functions
 *
 * @version $Revision::       $$Date::             $
 * @copyright Copyright (c) 2004 to the present, Mauricio Villegas <mauvilsa@upv.es>
 */

#ifndef __MV_LOG_H__
#define __MV_LOG_H__

#include <stdio.h>

void _proxy_log(FILE *log_file, int log_level, int level, const char *fmt, ...)
    __attribute__((format (printf, 4, 5)));
#define proxy_log(log_file, log_level, level, fmt, ...) _proxy_log(log_file, log_level, level, fmt"\n", ##__VA_ARGS__)

#endif
