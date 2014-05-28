/**
 * Functions for allocating memory
 *
 * @version $Revision::       $$Date::             $
 * @copyright Copyright (c) 2004 to the present, Mauricio Villegas <mauvilsa@upv.es>
 */

#ifndef __MV_MEM_H__
#define __MV_MEM_H__

#include <sys/ipc.h>
#include <sys/shm.h>

#define bfree(M,b) free((char**)M-b)

#define malloc_intv(size,vec,clr) mem((size)*sizeof(int),clr,(char**)(vec))
#define malloc_doublev(size,vec,clr) mem((size)*sizeof(double),clr,(char**)(vec))

//void nfree(void** p);
int mem(int size,int clr,char** _p);
int mmem(int R,int C,int size,int clr,char*** _mat);
int vrmem(int* R,int C,int size,int clr,char*** _mat);
int bmem(int R,int C,int size,int clr,int brd,char*** _mat);
int tmem(int D,int size,int clr,char*** _mat);
int mclone(char** mat, int R, int C, int size, char*** _clon);

#define shmalloc_floatm(R,C,M,clr,sid) mshmget((R),(C),sizeof(float),(clr),(char***)(M),(sid))
#define shmattach_floatm(R,C,M,sid) mshmat((R),(C),sizeof(float),(char***)(M),(sid))

int mshmget( int R, int C, int size, int clr, char*** _mat, int* _shmid );
int mshmat( int* _R, int* _C, int size, char*** _mat, int shmid );
void mshmdt( void* p );

#endif
