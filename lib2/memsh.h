/**
 * Functions for allocating shared memory
 *
 * @version $Revision::       $$Date::             $
 * @copyright Copyright (c) 2004 to the present, Mauricio Villegas <mauvilsa@upv.es>
 */

#ifndef __MV_MEMSH_H__
#define __MV_MEMSH_H__

#include "mem.h"

int shmalloc_I1v( int D, I1** _vec, char clr, int* _sid );

int shmalloc_I1m( int R, int C, I1*** _mat, char clr, int* _sid );
int shmalloc_F1m( int R, int C, F1*** _mat, char clr, int* _sid );

int shmattach_I1v( int* _D, I1** _vec, int sid );

int shmattach_F1m( int* _R, int* _C, F1*** _mat, int sid );


int ashmget( int size, int D, char clr, char** _vec, int* _shmid );
int ashmat( int size, int* _D, char** _vec, int shmid );
int mshmget( int size, int R, int C, char clr, char*** _mat, int* _shmid );
int mshmat( int size, int* _R, int* _C, char*** _mat, int shmid );
int vrshmget( int size, int tnnz, int* R, int C, char clr, char*** _mat, int** _R, int* _shmid );
int vrshmat( int size, int** _R, int* _C, char*** _mat, int shmid );
int shmfree( void *ptr, void *shptr, int shmid, int onlydetach );

#endif
