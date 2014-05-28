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

/*void nfree(void** p) {
  if(p[0]!=NULL) {
    free(p[0]);
    p[0]=NULL;
  } else
    fprintf(stderr,"nfree: warning: tryed to free a null pointer\n");
}*/

int mem(int size,int clr,char** _p)
{
  char *p;

  if((p=(char*)malloc(size))==NULL)
    return EXIT_FAILURE;
  if(clr)
    memset(p,0,size);
  *_p=p;

  return EXIT_SUCCESS;
}

int mmem(int R,int C,int size,int clr,char*** _mat)
{
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

int vrmem(int* R,int C,int size,int clr,char*** _mat)
{
  int c;
  char **mat,*vec;

  int TR=0;
  for(c=0;c<C;c++)
    TR+=R[c];

  mat=(char**)malloc(C*sizeof(char*)+TR*C*size);
  if(mat==NULL)
    return EXIT_FAILURE;
  vec=(char*)(mat+C);
  for(c=0;c<C;vec+=R[c]*size,c++)
    mat[c]=vec;
  if(clr)
    memset(mat[0],0,TR*C*size);

  *_mat=mat;
  return EXIT_SUCCESS;
}

int bmem(int R,int C,int size,int clr,int brd,char*** _mat)
{
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

int tmem(int D,int size,int clr,char*** _mat)
{
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

int mclone(char** mat, int R, int C, int size, char*** _clon) {
  int err = mmem(R,C,size,0,_clon);
  if(err!=EXIT_SUCCESS)
    return err;

  memcpy(_clon[0][0],mat[0],R*C*size);

  return EXIT_SUCCESS;
}

int mshmget( int R, int C, int size, int clr, char*** _mat, int* _shmid ) {
  int shmid;
  char **mat, *vec;

  shmid = shmget( IPC_PRIVATE, 4*sizeof(int)+R*C*size, 0660|IPC_CREAT|IPC_EXCL );
  if( shmid == -1 )
    return EXIT_FAILURE;

  vec = (char*)shmat( shmid, NULL, 0 );
  ((int*)vec)[0] = R;
  ((int*)vec)[1] = C;
  ((int*)vec)[2] = size;
  ((int*)vec)[3] = 'm';
  shmdt( vec );

  if( mshmat( NULL, NULL, size, &mat, shmid ) )
    return EXIT_FAILURE;

  if( clr )
    memset( mat[0], 0, R*C*size );

  *_mat = mat;
  *_shmid = shmid;
  return EXIT_SUCCESS;
}

int mshmat( int* _R, int* _C, int size, char*** _mat, int shmid ) {
  int c, R, C;
  char **mat, *vec;

  vec = (char*)shmat( shmid, NULL, 0 );
  if( ((int*)vec)[3] != 'm' ||
      ((int*)vec)[2] != size )
      return EXIT_FAILURE;
  R = ((int*)vec)[0];
  C = ((int*)vec)[1];

  mat = (char**)malloc(C*sizeof(char*));
  if( mat == NULL ) {
    shmdt( vec );
    return EXIT_FAILURE;
  }

  vec += 4*sizeof(int);
  for( c=0; c<C; c++, vec += R*size )
    mat[c] = vec;

  if( _R != NULL )
    *_R = R;
  if( _C != NULL )
    *_C = C;
  *_mat = mat;
  return EXIT_SUCCESS;
}

void mshmdt( void* p ) {
  int *vec = (int*)p;
  switch( *(vec-1) ) {
    case 'm':
      shmdt( vec-4 );
      break;
    default:
      fprintf(stderr,"mshmdt: warning: unknown type to detach '%c'\n",(char)(*(vec-1)));
  }
}
