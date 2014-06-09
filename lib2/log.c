/**
 * Logging functions
 *
 * @version $Revision::       $$Date::             $
 * @copyright Copyright (c) 2004 to the present, Mauricio Villegas <mauvilsa@upv.es>
 */

#include "log.h"

#include <stdarg.h>

void _proxy_log(FILE *log_file, int log_level, int level, const char *fmt, ...) {
  va_list arg;

  if( ( ! log_level ) || level > log_level )
    return;

  va_start(arg, fmt);
  vfprintf(log_file, fmt, arg);
  va_end(arg);
}
