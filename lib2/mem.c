/**
 * Functions for allocating memory
 *
 * @version $Revision::       $$Date::             $
 * @copyright Copyright (c) 2004 to the present, Mauricio Villegas <mauvilsa@upv.es>
 */

#include "mem.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//#define bfree(M,b) free((char**)M-b)
void bfree( void* mat, int brd )
  { free((char**)mat-brd); }

int clone_graym( gray** mat, int R, int C, gray*** _clon )
  { return mclone((char**)mat,R,C,sizeof(gray),(char***)_clon); }

int malloc_I1v( int size, I1** _vec, char clr )
  { return mem(size*sizeof(I1),clr,(char**)_vec); }
int malloc_F1v( int size, F1** _vec, char clr )
  { return mem(size*sizeof(F1),clr,(char**)_vec); }
int malloc_D1v( int size, D1** _vec, char clr )
  { return mem(size*sizeof(D1),clr,(char**)_vec); }
int malloc_iI1v( int size, iI1** _vec, char clr )
  { return mem(size*sizeof(iI1),clr,(char**)_vec); }
int malloc_iF1v( int size, iF1** _vec, char clr )
  { return mem(size*sizeof(iF1),clr,(char**)_vec); }

int malloc_graym( int imW, int imH, gray*** _im, char clr )
  { return mmem(imW,imH,sizeof(gray),clr,(char***)_im); }
int malloc_pixelm( int imW, int imH, pixel*** _im, char clr )
  { return mmem(imW,imH,sizeof(pixel),clr,(char***)_im); }
int malloc_I1m( int R, int C, I1*** _mat, char clr )
  { return mmem(R,C,sizeof(I1),clr,(char***)_mat); }
int malloc_F1m( int R, int C, F1*** _mat, char clr )
  { return mmem(R,C,sizeof(F1),clr,(char***)_mat); }


void nullfree(void* ptr) {
  void** p = (void**)ptr;
  if( p[0] != NULL ) {
    free(p[0]);
    p[0] = NULL;
  }
  else
    fprintf(stderr,"nullfree: warning: tried to free a null pointer\n");
}

int mclone(char** mat, int R, int C, int size, char*** _clon) {
  int err = mmem(R,C,size,0,_clon);
  if(err!=EXIT_SUCCESS)
    return err;

  memcpy(_clon[0][0],mat[0],R*C*size);

  return EXIT_SUCCESS;
}

int mem(int size,char clr,char** _p) {
  char *p;

  if((p=(char*)malloc(size))==NULL)
    return EXIT_FAILURE;
  if(clr)
    memset(p,0,size);
  *_p=p;

  return EXIT_SUCCESS;
}

int mmem(int R,int C,int size,char clr,char*** _mat) {
  int c;
  char **mat,*vec;

  mat=(char**)malloc(C*sizeof(char*)+R*C*size);
  if(mat==NULL)
    return EXIT_FAILURE;
  vec=(char*)(mat+C);
  for(c=0;c<C;c++,vec+=R*size)
    mat[c]=vec;
  if(clr)
    memset(mat[0],0,R*C*size);

  *_mat=mat;
  return EXIT_SUCCESS;
}

int bmem(int R,int C,int size,char clr,int brd,char*** _mat) {
  int c;
  char **mat,*vec;

  mat=(char**)malloc((C+2*brd)*sizeof(char*)+(R+2*brd)*(C+2*brd)*size);
  if(mat==NULL)
    return EXIT_FAILURE;
  vec=(char*)(mat+C+2*brd)+brd*size;
  for(c=0;c<C+2*brd;c++,vec+=(R+2*brd)*size)
    mat[c]=vec;
  if(clr)
    memset((char*)(mat+C+2*brd),0,(R+2*brd)*(C+2*brd)*size);

  *_mat=mat+brd;
  return EXIT_SUCCESS;
}

int tmem(int D,int size,char clr,char*** _mat) {
  int d;
  char **mat,*vec;

  mat=(char**)malloc(D*sizeof(char*)+size*D*(D+1)/2);
  if(mat==NULL)
    return EXIT_FAILURE;
  vec=(char*)(mat+D);
  for(d=0;d<D;d++,vec+=d*size)
    mat[d]=vec;
  if(clr)
    memset(mat[0],0,size*D*(D+1)/2);

  *_mat=mat;
  return EXIT_SUCCESS;
}

void vrmem_index( int size, int* R, int C, char** mat ) {
  int c;
  char *ptr = mat[0]+R[0]*size;
  for( c=1; c<C; ptr+=R[c]*size, c++ )
    mat[c] = ptr;
}

int vrmem( int size, int tnnz, int* R, int C, char clr, char*** _mat, int** _R ) {
  char **mat, *ptr;

  ptr = (char*)malloc(C*(sizeof(int)+sizeof(char*))+tnnz*size);
  if( ptr == NULL )
    return EXIT_FAILURE;
  if( R != NULL )
    memcpy(ptr,R,C*sizeof(int));
  else
    memset(ptr,0,C*sizeof(int));
  R = (int*)ptr;
  mat = (char**)(ptr+C*sizeof(int));
  mat[0] = (char*)(mat+C);
  vrmem_index(size,R,C,mat);
  if( clr )
    memset(mat[0],0,tnnz*size);

  *_R = R;
  *_mat = mat;
  return EXIT_SUCCESS;
}
