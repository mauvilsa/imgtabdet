/**
 * Functions for reading and writing images
 *
 * @version $Revision::       $$Date::             $
 * @copyright Copyright (c) 2004 to the present, Mauricio Villegas <mauvilsa@upv.es>
 */

#ifndef __MV_IMGPROC_CV_H__
#define __MV_IMGPROC_CV_H__

#include "imgio_magick.h"

#include <opencv2/core/version.hpp>
#if CV_MAJOR_VERSION == 2
#include <opencv2/imgproc/imgproc.hpp>
#else
#include <opencv2/imgproc.hpp>
#endif
#include <opencv2/imgproc/imgproc_c.h>

int getpixels_magick_cv8u1_mat( Img* img, cv::Mat& cvimg );
int getpixels_magick_cv8u( Img* img, IplImage* cvimg );
int setpixels_magick_cv8u( Img* img, IplImage* cvimg );

int ccomp_contours_opencv( Img *img, CvMemStorage** _storage, CvSeq** _contours, int method, double summeps );
int join_ccomp_opencv( Img *img );

#endif
