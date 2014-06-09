/**
 * Functions for reading and writing matrices
 *
 * @version $Revision::       $$Date::             $
 * @copyright Copyright (c) 2004 to the present, Mauricio Villegas <mauvilsa@upv.es>
 */

/* todo: change fixed length buffer with getline */

#include "iomat.h"

#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define BUFFSZ 256

int fscanMATHead( FILE* file, int* _R, int* _C, char* _name, char* _type, int* _elements, char* _format );
int fscanMAT( char type, int elements, FILE* file, char*** _matrix, int* _R, int* _C, char* _name, int* _shmid );
int fprintMAT( char type, int elements, FILE* file, char** matrix, int R, int C, char* name, char format );
int freadMAT( char type, int elements, char* fname, char*** _matrix, int* _R, int* _C, char* _name, int* _shmid );
int fwriteMAT( char type, int elements, char* fname, char** matrix, int R, int C, char* name, char format );
int fscanVEC( char type, int elements, FILE* file, char* vec, int D, FILE* mess );

int fscanMAT_I1(FILE* F,I1*** M,int* R,int* C,char* N,int* s)
  { return fscanMAT('I',1,F,(char***)M,R,C,N,s); }
int fscanMAT_F1(FILE* F,F1*** M,int* R,int* C,char* N,int* s)
  { return fscanMAT('F',1,F,(char***)M,R,C,N,s); }

int freadMAT_I1(char* F,I1*** M,int* R,int* C,char* N,int* s)
  { return freadMAT('I',1,F,(char***)M,R,C,N,s); }
int freadMAT_F1(char* F,F1*** M,int* R,int* C,char* N,int* s)
  { return freadMAT('F',1,F,(char***)M,R,C,N,s); }

int fprintMAT_I1(FILE* F,I1** M,int R,int C,char* N,char f)
  { return fprintMAT('I',1,F,(char**)M,R,C,N,f); }
int fprintMAT_F1(FILE* F,F1** M,int R,int C,char* N,char f)
  { return fprintMAT('F',1,F,(char**)M,R,C,N,f); }

int fwriteMAT_I1(char* F,I1** M,int R,int C,char* N,char f)
  { return fwriteMAT('I',1,F,(char**)M,R,C,N,f); }
int fwriteMAT_F1(char* F,F1** M,int R,int C,char* N,char f)
  { return fwriteMAT('F',1,F,(char**)M,R,C,N,f); }

int fscanVEC_F1(FILE* F,F1* V,int D,FILE* l)
  { return fscanVEC('F',1,F,(char*)V,D,l); }


static int type_size( char type ) {
  switch(type) {
  case 'G': return sizeof(char);
  case 'I': return sizeof(int);
  case 'F': return sizeof(float);
  case 'D': return sizeof(double);
  }
  return 0;
}

static int fscanf_tuple( char type, int elements, int size, FILE* file, char* ptr ) {
  if( elements == 1 )
    switch(type) {
    case 'G': return fscanf( file, scanf_G, (gray*)ptr );
    case 'I': return fscanf( file, scanf_I, (int*)ptr );
    case 'F': return fscanf( file, scanf_F, (float*)ptr );
    case 'D': return fscanf( file, scanf_D, (double*)ptr );
    }

  int e, z = 0;
  e = fscanf( file, " (" );
  for( e=elements-1; e>1; e--, ptr+=size )
    switch(type) {
    case 'G': z += fscanf( file, scanf_G ",", (gray*)ptr );   break;
    case 'I': z += fscanf( file, scanf_I ",", (int*)ptr );    break;
    case 'F': z += fscanf( file, scanf_F ",", (float*)ptr );  break;
    case 'D': z += fscanf( file, scanf_D ",", (double*)ptr ); break;
    }
  switch(type) {
  case 'G': z += fscanf( file, scanf_G "," scanf_G ")", (gray*)ptr,   (gray*)(ptr+size) );   break;
  case 'I': z += fscanf( file, scanf_I "," scanf_I ")", (int*)ptr,    (int*)(ptr+size) );    break;
  case 'F': z += fscanf( file, scanf_F "," scanf_F ")", (float*)ptr,  (float*)(ptr+size) );  break;
  case 'D': z += fscanf( file, scanf_D "," scanf_D ")", (double*)ptr, (double*)(ptr+size) ); break;
  }
  return z;
}

static void fprintf_tuple( char type, int elements, int size, FILE* file, char* ptr ) {
  if( elements == 1 )
    switch(type) {
    case 'G': fprintf( file, " " printf_G, *((gray*)ptr) );   break;
    case 'I': fprintf( file, " " printf_I, *((int*)ptr) );    break;
    case 'F': fprintf( file, " " printf_F, *((float*)ptr) );  break;
    case 'D': fprintf( file, " " printf_D, *((double*)ptr) ); break;
    }
  else {
    int e;
    fprintf( file, " (" );
    for( e=elements-1; e>1; e--, ptr+=size )
      switch(type) {
      case 'G': fprintf( file, printf_G ",", *((gray*)ptr) );   break;
      case 'I': fprintf( file, printf_I ",", *((int*)ptr) );    break;
      case 'F': fprintf( file, printf_F ",", *((float*)ptr) );  break;
      case 'D': fprintf( file, printf_D ",", *((double*)ptr) ); break;
      }
    switch(type) {
    case 'G': fprintf( file, printf_G "," printf_G ")", *((gray*)ptr),   *((gray*)(ptr+size)) );   break;
    case 'I': fprintf( file, printf_I "," printf_I ")", *((int*)ptr),    *((int*)(ptr+size)) );    break;
    case 'F': fprintf( file, printf_F "," printf_F ")", *((float*)ptr),  *((float*)(ptr+size)) );  break;
    case 'D': fprintf( file, printf_D "," printf_D ")", *((double*)ptr), *((double*)(ptr+size)) ); break;
    }
  }
}

int fscanMATHead( FILE* file, int* _R, int* _C, char* _name, char* _type, int* _elements, char* _format ) {
  char format = 0;
  char type = '-';
  int elements = 1;
  int R=0, C=0, z=0;

  int ch = fgetc(file);
  ungetc(ch,file);

  /* Binary: "M{G|I|F|D}{int#ROWS}{int#COLS}{int#ELEMS}{GZIP_DATA}" */
  if( ch == 'M' ) {
    //long offset = ftell(file);
    ch = fgetc(file);
    type = fgetc(file);
    if( ! ( type=='G' || type=='I' || type=='F' || type=='D' ) )
      return EXIT_FAILURE;
    z += fread(&R,sizeof(int),1,file);
    z += fread(&C,sizeof(int),1,file);
    z += fread(&elements,sizeof(int),1,file);
    ch = fgetc(file);
    //if( _format != NULL ) // leave file position after the header
      ungetc(ch,file);
    //else
    //  fseek(file,offset,SEEK_SET);
    if( z!=3 || ch!=0x1F )
      return EXIT_FAILURE;
    format = format_MGZ;
  }

  /* ASCII with header: "% {#ROWS} {#COLS} [{#TYPE} [{#ELEMS}]]" */
  else if( ch == '%' ) {
    char ttype = type;
    int telements = elements;
    char buffer[BUFFSZ];
    if( fgets(buffer,BUFFSZ,file) != buffer )
      return EXIT_FAILURE;
    if( ! ( sscanf(buffer,"%% %d %d %c %d%*[\n]",&R,&C,&ttype,&telements) == 4 ||
            sscanf(buffer,"%% %d %d %c%*[\n]",&R,&C,&ttype) == 3 ||
            sscanf(buffer,"%% %d %d%*[\n]",&R,&C) == 2 ) )
    //if( fscanf(file,"%% %d %d\n",&R,&C) != 2 )
      return EXIT_FAILURE;
    type = ttype;
    elements = telements;
    format = format_ASCII;
  }

  /* ASCII without header */
  else if( ch==' ' || ch=='.' || ch=='-' || (ch>='0'&&ch<='9') ) {
    long offset = ftell(file);
    int z=0;
    while( ! feof(file) ) {
      int r=0;
      char ch='$';
      while( ch != '\n' ) {
        ch = fgetc(file);
        while( ch==' ' || ch=='\t' || ch=='\r' )
          ch = fgetc(file);
        if( feof(file) )
          break;
        if( ch != '\n' ) {
          r++;
          while( ! ( ch==' ' || ch=='\t' || ch=='\n' || ch=='\r' ) )
            ch = fgetc(file);
        }
      }
      z += r;
      if( C==0 )
        R = r;
      C++;
    }
    C--;
    fseek(file,offset,SEEK_SET);
    if( C==0 || R==0 || z!=C*R )
      return EXIT_FAILURE;

    format = format_ASCII;
  }

  /* Octave style matrix */
  else {
    char buffer[BUFFSZ];
    long offset = ftell(file);
    do {
      //if( _format != NULL ) // leave file position after the header
        offset = ftell(file);
      if( fgets(buffer,BUFFSZ,file) == NULL )
        return EXIT_FAILURE;
      if( ! strncmp(buffer,"# type: ",8) ) {
        if( ! strncmp(buffer+8,"matrix",6) )
          elements = 1;
        else if( ! strncmp(buffer+8,"complex matrix",14) )
          elements = 2;
        else if( sscanf(buffer+8,"%d-tuple matrix",&elements) != 1 )
          return EXIT_FAILURE;
        format = format_MAT;
      }
      else if( ! strncmp(buffer,"# rows: ",8) )
        R = atoi(buffer+8);
      else if( ! strncmp(buffer,"# columns: ",11) )
        C = atoi(buffer+11);
      else if( ! strncmp(buffer,"# name: ",8) ) {
        if( _name != NULL ) {
          strncpy(_name,buffer+8,BUFFSZ-8);
          if( strrchr(_name,'\n') != NULL )
            *strrchr(_name,'\n') = '\0';
          if( strrchr(_name,'\r') != NULL )
            *strrchr(_name,'\r') = '\0';
        }
      }
      else if( strncmp(buffer,"# Created by ",13) )
        buffer[0] = '\0';
    } while( buffer[0] == '#' );
    fseek(file,offset,SEEK_SET);
  }

  if( _R != NULL )
    _R[0] = R;
  if( _C != NULL )
    _C[0] = C;
  if( _format != NULL )
    _format[0] = format;
  if( _elements != NULL )
    _elements[0] = elements;
  if( _type != NULL )
    _type[0] = type;

  return EXIT_SUCCESS;
}

int fscanMAT( char type, int elements, FILE* file, char*** _matrix, int* _R, int* _C, char* _name, int* _shmid ) {
  int felements, R, C, RC, r, c, Csize, Esize, z=0;
  int size = type_size(type);
  char *p, *mat;
  char fformat, ftype;

  /* read header */
  if( fscanMATHead(file,&R,&C,_name,&ftype,&felements,&fformat) ) {
    fprintf(stderr,"fscanMAT: error: problems reading matrix header\n");
    return EXIT_FAILURE;
  }

  /* check data is of correct type */
  if( size==0 ||
      felements!=elements ||
      ( ftype!='-' && type!=ftype ) ) {
    fprintf(stderr,"fscanMAT: error: unexpected matrix type M%c%d\n",ftype,felements);
    return EXIT_FAILURE;
  }

  /* reserve memory */
  if( mshmget(size,R,C,0,_matrix,_shmid) ) {
    fprintf(stderr,"fscanMAT: error: unable to reserve memory\n");
    return EXIT_FAILURE;
  }

  if( _R != NULL )
    *_R = R;
  if( _C != NULL )
    *_C = C;

  RC = R*C;
  mat = _matrix[0][0];
  gzFile zfile;

  /* load matrix depending on each format */
  switch(fformat) {

  case format_MGZ:
    z = fflush(file);
    //fprintf(stderr,"fflush returned %d\n",z);
  /*if((file=fopen("xid.mgz+","rb"))==NULL) {
    fprintf(stderr,"freadMAT: error: unable to open file\n");
    return EXIT_FAILURE;
  }
  fseeko(file,14,SEEK_SET);*/
    if( (zfile=gzdopen(dup(fileno(file)),"rb")) == NULL ) {
      fprintf(stderr,"fscanMAT: error: unable to open stream for gzip decompression\n");
      return EXIT_FAILURE;
    }
    z = gzread(zfile,mat,size*elements*RC);
    //fprintf(stderr,"gzread returned %d\n",z);
    z /= size;
    gzclose(zfile);
    break;

  case format_ASCII:
    Esize = size*elements;
    for( c=C-1, p=mat; c>=0; c-- )
      for( r=R-1; r>=0; r--, p+=Esize )
        z += fscanf_tuple(type,elements,size,file,p);
    break;

  case format_MAT:
    Esize = size*elements;
    Csize = R*Esize;
    for( r=R-1; r>=0; r--, mat+=Esize )
      for( c=C-1, p=mat; c>=0; c--, p+=Csize )
        z += fscanf_tuple(type,elements,size,file,p);
    break;

  default:
    fprintf(stderr,"fscanMAT: error: unknown matrix format\n");
    return EXIT_FAILURE;
  }

  if( z != elements*RC ) {
    fprintf(stderr,"R=%d C=%d e=%d RC=%d z=%d\n",R,C,elements,RC,z);
    fprintf(stderr,"fscanMAT: error: unexpected matrix format\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int fprintMAT( char type, int elements, FILE* file, char** matrix, int R, int C, char* name, char format ) {
  int r, c, Csize, Esize;
  int size = type_size(type);
  char *p, *mat=matrix[0];
  gzFile zfile;

  if( C>1 && ( (matrix[1]-mat) != size*elements*R ) ) {
    fprintf(stderr,"fprintMAT: error: matrix must be contiguous in memory\n");
    return EXIT_FAILURE;
  }

  switch(format) {

  case format_MGZ:
    fprintf(file,"M%c",type);
    fwrite(&R,sizeof(int),1,file);
    fwrite(&C,sizeof(int),1,file);
    fwrite(&elements,sizeof(int),1,file);
    fflush(file);
    if( (zfile=gzdopen(dup(fileno(file)),"wb")) == NULL ) {
      fprintf(stderr,"fprintMAT: error: unable to open stream for gzip compression\n");
      return EXIT_FAILURE;
    }
    c = size*elements*R*C;
    r = gzwrite(zfile,mat,c);
    gzclose(zfile);
    if( r != c ) {
      fprintf(stderr,"fprintMAT: error: problem writing to gzip compression stream\n");
      return EXIT_FAILURE;
    }
    break;

  case format_ASCII:
    if( elements==1 )
      fprintf(file,"%% %d %d %c\n",R,C,type);
    else
      fprintf(file,"%% %d %d %c %d\n",R,C,type,elements);
  case format_ASCIINOHEAD:
    Esize = size*elements;
    for(c=C-1,p=mat;c>=0;c--) {
      for(r=R-1;r>=0;r--,p+=Esize)
        fprintf_tuple(type,elements,size,file,p);
      fprintf(file,"\n");
    }
    break;

  case format_MAT:
    fprintf(file,"# name: %s\n",name!=NULL?name:"data");
    if( elements <= 2 )
      fprintf(file,"# type: %smatrix\n",elements==1?"":"complex");
    else
      fprintf(file,"# type: %d-tuple matrix\n",elements);
    fprintf(file,"# rows: %d\n# columns: %d\n",R,C);
    Esize = size*elements;
    Csize = R*Esize;
    for(r=R-1;r>=0;r--,mat+=Esize) {
      for(c=C-1,p=mat;c>=0;c--,p+=Csize)
        fprintf_tuple(type,elements,size,file,p);
      fprintf(file,"\n");
    }
    break;

  default:
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int freadMAT( char type, int elements, char* fname, char*** _matrix, int* _R, int* _C, char* _name, int* _shmid ) {
  FILE *file;
  if((file=fopen(fname,"rb"))==NULL) {
    fprintf(stderr,"freadMAT: error: unable to open file %s\n",fname);
    return EXIT_FAILURE;
  }
  int err = fscanMAT(type,elements,file,_matrix,_R,_C,_name,_shmid);
  fclose(file);
  return err;
}

int fwriteMAT( char type, int elements, char* fname, char** matrix, int R, int C, char* name, char format ) {
  FILE *file;
  if((file=fopen(fname,"wb"))==NULL) {
    fprintf(stderr,"fwriteMAT: error: unable to open file %s\n",fname);
    return EXIT_FAILURE;
  }
  int err = fprintMAT(type,elements,file,matrix,R,C,name,format);
  fclose(file);
  return err;
}

int fscanVEC( char type, int elements, FILE* file, char* vec, int D, FILE* mess ) {
  int err = 0;
  int d = 0, dd = 1;
  int size = type_size(type);
  int ch = fgetc(file);

  if( feof(file) ) {
    if( mess != NULL )
      fprintf(mess,"fscanVEC: warning: reached end of file\n");
    return EXIT_FAILURE;
  }

  ungetc(ch,file);
  if( ch=='M' ) {
    char **mat;
    err = fscanMAT(type,elements,file,&mat,&d,&dd,NULL,NULL);
    if( (!err) && d==D && dd==1 )
      memcpy(vec,mat[0],D*size*elements);
    if( !err )
      free(mat);
  }

  else {
    d=0;
    ch='-';
    char cval[BUFFSZ];
    char dummy[16];
    while(1) {
      ch=fgetc(file);
      if(ch=='\n' || feof(file))
        break;
      else if(ch==' ' || ch=='\t')
        continue;
      else {
        ungetc(ch,file);
        int p=0;
        do {
          ch=fgetc(file);
          if(p<BUFFSZ)
            cval[p++]=ch;
        } while(!(ch==' ' || ch=='\t' || ch=='\n' || ch=='\r'));
        cval[p-1]='\0';
        ungetc(ch,file);
        char *valp=vec+d*size;
        if(d>=D)
          valp = dummy;
        int rf=0;
        switch(type) {
        case 'G': rf = sscanf(cval,scanf_G,(gray*)valp);   break;
        case 'I': rf = sscanf(cval,scanf_I,(int*)valp);    break;
        case 'F': rf = sscanf(cval,scanf_F,(float*)valp);  break;
        case 'D': rf = sscanf(cval,scanf_D,(double*)valp); break;
        }
        if(rf==1)
          d++;
        else
          err--;
      }
    }
  }

  if( err || d!=D || dd!=1 ) {
    if( mess != NULL ) {
      if( err )
        fprintf(mess,"fscanVEC: error: unable to read vector\n");
      if( d!=D || dd!=1 )
        fprintf(mess,"fscanVEC: error: vector is not of the right dimensionality\n");
    }
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
