/**
 * Functions for reading and writing images
 *
 * @version $Revision::       $$Date::             $
 * @copyright Copyright (c) 2004 to the present, Mauricio Villegas <mauvilsa@upv.es>
 */

#ifndef __MV_IMGPROC_CV_H__
#define __MV_IMGPROC_CV_H__

#include "imgio_magick.h"

#ifdef CV20
#include <opencv/cv.h>
#else
#include <opencv2/imgproc/imgproc_c.h>
#endif

int getpixels_magick_cv8u( Img* img, IplImage* cvimg );
int setpixels_magick_cv8u( Img* img, IplImage* cvimg );

int ccomp_contours_opencv( Img *img, CvMemStorage** _storage, CvSeq** _contours, int method, double summeps );
int join_ccomp_opencv( Img *img );

#endif
