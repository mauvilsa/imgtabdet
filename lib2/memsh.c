/**
 * Functions for allocating shared memory
 *
 * @version $Revision::       $$Date::             $
 * @copyright Copyright (c) 2004 to the present, Mauricio Villegas <mauvilsa@upv.es>
 */

#include "memsh.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/ipc.h>
#include <sys/shm.h>

int shmalloc_I1v( int D, I1** _vec, char clr, int* _sid )
  { return ashmget(sizeof(int),D,clr,(char**)_vec,_sid); }

int shmalloc_I1m( int R, int C, I1*** _mat, char clr, int* _sid )
  { return mshmget(sizeof(I1),R,C,clr,(char***)_mat,_sid); }
int shmalloc_F1m( int R, int C, F1*** _mat, char clr, int* _sid )
  { return mshmget(sizeof(F1),R,C,clr,(char***)_mat,_sid); }

int shmattach_I1v( int* _D, I1** _vec, int sid )
  { return ashmat(sizeof(I1),_D,(char**)_vec,sid); }

int shmattach_F1m( int* _R, int* _C, F1*** _mat, int sid )
  { return mshmat(sizeof(F1),_R,_C,(char***)_mat,sid); }


int ashmget( int size, int D, char clr, char** _vec, int* _shmid ) {
  if( _shmid == NULL )
    return mem(size*D,clr,_vec);

  int shmid = shmget( IPC_PRIVATE, 3*sizeof(int)+size*D, 0660|IPC_CREAT|IPC_EXCL );
  if( shmid == -1 )
    return EXIT_FAILURE;

  int *ptr = (int*)shmat( shmid, NULL, 0 );
  ptr[0] = D;
  ptr[1] = size;
  ptr[2] = 'a';
  shmdt( ptr );

  char *vec;
  if( ashmat( size, NULL, &vec, shmid ) )
    return EXIT_FAILURE;

  if( clr )
    memset( vec, 0, size*D );

  *_vec = vec;
  *_shmid = shmid;
  return EXIT_SUCCESS;
}

int ashmat( int size, int* _D, char** _vec, int shmid ) {
  int *ptr = (int*)shmat( shmid, NULL, 0 );
  if( ptr[2] != 'a' ||
      ptr[1] != size ) {
      shmdt( ptr );
      return EXIT_FAILURE;
  }
  if( _D != NULL )
    *_D = ptr[0];
  *_vec = (char*)(ptr+3);
  return EXIT_SUCCESS;
}

int mshmget( int size, int R, int C, char clr, char*** _mat, int* _shmid ) {
  if( _shmid == NULL )
    return mmem(R,C,size,clr,_mat);

  int shmid = shmget( IPC_PRIVATE, 4*sizeof(int)+R*C*size, 0660|IPC_CREAT|IPC_EXCL );
  if( shmid == -1 )
    return EXIT_FAILURE;

  int *ptr = (int*)shmat( shmid, NULL, 0 );
  ptr[0] = R;
  ptr[1] = C;
  ptr[2] = size;
  ptr[3] = 'm';
  shmdt( ptr );

  char **mat;
  if( mshmat( size, NULL, NULL, &mat, shmid ) )
    return EXIT_FAILURE;

  if( clr )
    memset( mat[0], 0, R*C*size );

  *_mat = mat;
  *_shmid = shmid;
  return EXIT_SUCCESS;
}

int mshmat( int size, int* _R, int* _C, char*** _mat, int shmid ) {
  int c, R, C;
  char **mat, *ptr;

  ptr = (char*)shmat( shmid, NULL, 0 );
  if( ((int*)ptr)[3] != 'm' ||
      ((int*)ptr)[2] != size ) {
      shmdt( ptr );
      return EXIT_FAILURE;
  }
  R = ((int*)ptr)[0];
  C = ((int*)ptr)[1];

  mat = (char**)malloc(C*sizeof(char*));
  if( mat == NULL ) {
    shmdt( ptr );
    return EXIT_FAILURE;
  }

  ptr += 4*sizeof(int);
  for( c=0; c<C; c++, ptr += R*size )
    mat[c] = ptr;

  if( _R != NULL )
    *_R = R;
  if( _C != NULL )
    *_C = C;
  *_mat = mat;
  return EXIT_SUCCESS;
}

int vrshmget( int size, int tnnz, int* R, int C, char clr, char*** _mat, int** _R, int* _shmid ) {
  int shmid;

  if( _shmid == NULL )
    return vrmem(size,tnnz,R,C,clr,_mat,_R);

  shmid = shmget( IPC_PRIVATE, (4+C)*sizeof(int)+tnnz*size, 0660|IPC_CREAT|IPC_EXCL );
  if( shmid == -1 )
    return EXIT_FAILURE;

  int *ptr = (int*)shmat( shmid, NULL, 0 );
  ptr[0] = tnnz;
  ptr[1] = C;
  ptr[2] = size;
  ptr[3] = 's';
  if( R != NULL )
    memcpy(ptr+4,R,C*sizeof(int));
  else
    memset(ptr+4,0,C*sizeof(int));
  shmdt( ptr );

  if( vrshmat( size, _R, NULL, _mat, shmid ) )
    return EXIT_FAILURE;

  if(clr)
    memset(_mat[0][0],0,tnnz*size);

  *_shmid = shmid;
  return EXIT_SUCCESS;
}

int vrshmat( int size, int** _R, int* _C, char*** _mat, int shmid ) {
  int C, *R;
  char **mat, *ptr;

  ptr = (char*)shmat( shmid, NULL, 0 );
  if( ((int*)ptr)[3] != 's' ||
      ((int*)ptr)[2] != size ) {
      shmdt( ptr );
      return EXIT_FAILURE;
  }
  C = ((int*)ptr)[1];
  //tnnz = ((int*)ptr)[0];

  mat = (char**)malloc(C*sizeof(char*));
  if( mat == NULL ) {
    shmdt( ptr );
    return EXIT_FAILURE;
  }
  ptr += 4*sizeof(int);
  R = (int*)ptr;
  mat[0] = ptr+C*sizeof(int);
  vrmem_index(size,R,C,mat);

  if( _C != NULL )
    *_C = C;
  *_R = R;
  *_mat = mat;
  return EXIT_SUCCESS;
}

int shmfree( void *ptr, void *shptr, int shmid, int onlydetach ) {
  if( shmid != 0 ) {
    int *shp = (int*)(((void**)shptr)[0]);
    switch( *(shp-1) ) {
      case 'm':
      case 's':
        //fprintf(stderr,"detach %c\n",(char)(*(shp-1)));
        shmdt( shp-4 );
        break;
      case 'a':
        ptr = NULL;
        //fprintf(stderr,"detach %c, ptr=%p\n",(char)(*(shp-1)),ptr);
        shmdt( shp-3 );
        break;
      default:
        fprintf(stderr,"mshmdt: warning: unknown type to detach '%c'\n",(char)(*(shp-1)));
        return EXIT_FAILURE;
    }
    if( ! onlydetach )
      shmctl( shmid, IPC_RMID, 0 );
  }
  if( ptr != NULL )
    nullfree(ptr);
  return EXIT_SUCCESS;
}
