/**
 * Functions for reading and writing matrices
 *
 * @version $Revision::       $$Date::             $
 * @copyright Copyright (c) 2004 to the present, Mauricio Villegas <mauvilsa@upv.es>
 */

#ifndef __MV_IOMAT_H__
#define __MV_IOMAT_H__

#include "memsh.h"

#include <stdio.h>

#define format_ASCII       'a'
#define format_ASCIINOHEAD 'A'
#define format_MAT         'm'
#define format_LIBSVM      'S'
#define format_MGZ         'Z'

#define scanf_G   "%hhd"
#define scanf_I   "%d"
#define scanf_F   "%f"
#define scanf_D   "%lf"

#define printf_G  "%d"
#define printf_I  "%d"
#define printf_F  "%.9g"
#define printf_D  "%.15g"

int fscanMAT_I1(FILE* F,I1*** M,int* R,int* C,char* N,int* s);
int fscanMAT_F1(FILE* F,F1*** M,int* R,int* C,char* N,int* s);

int freadMAT_I1(char* F,I1*** M,int* R,int* C,char* N,int* s);
int freadMAT_F1(char* F,F1*** M,int* R,int* C,char* N,int* s);

int fprintMAT_I1(FILE* F,I1** M,int R,int C,char* N,char f);
int fprintMAT_F1(FILE* F,F1** M,int R,int C,char* N,char f);

int fwriteMAT_I1(char* F,I1** M,int R,int C,char* N,char f);
int fwriteMAT_F1(char* F,F1** M,int R,int C,char* N,char f);

int fscanVEC_F1(FILE* F,F1* V,int D,FILE* l);

#endif
