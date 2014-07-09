/**
 * Compatibility definitions for unimplemented shared memory functions
 *
 * @version $Revision::       $$Date::             $
 * @copyright Copyright (c) 2004 to the present, Mauricio Villegas <mauvilsa@upv.es>
 */

#ifndef __MV_MEMSHNO_H__
#define __MV_MEMSHNO_H__

#include "mem.h"

#define shmalloc_I1v(D,V,c,sh)       malloc_I1v(D,V,c)

#define shmalloc_I1m(R,C,M,c,sh)     malloc_I1m(R,C,M,c)
#define shmalloc_F1m(R,C,M,c,sh)     malloc_F1m(R,C,M,c)

#define shmattach_I1v(D,V,sh)        1
#define shmattach_F1m(R,C,M,sh)      1

#define ashmget(s,D,c,V,sh)          mem((s)*(D),c,V)
#define mshmget(s,R,C,c,M,sh)        mmem(R,C,s,c,M)
#define vrshmget(s,nz,R,C,c,M,r,sh)  vrmem(s,nz,R,C,c,M,r)

#define shmfree(p,sh,id,od)          nullfree(p)

#endif
