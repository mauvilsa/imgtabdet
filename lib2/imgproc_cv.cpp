/**
 * Functions for reading and writing images using MagickCore
 *
 * @version $Revision$$Date::             $
 * @copyright Copyright (c) 2014 to the present, Mauricio Villegas <mauvilsa@upv.es>
 */

#include "imgproc_cv.h"

#ifdef CV20
#include <opencv/cxflann.h>
#else
#include <opencv2/flann/flann.hpp>
#endif

int getpixels_magick_cv8u( Img* img, IplImage* cvimg ) {
  ExceptionInfo *exception = AcquireExceptionInfo();
  const PixelPacket *pixs =
    GetVirtualPixels( img->image, 0, 0, img->width, img->height, exception );

  if( exception->severity != UndefinedException )
    CatchException( exception );
  DestroyExceptionInfo( exception );

  gray *ptr = (gray*)cvimg->imageData;

  int x,y,n;
  for( y=0,n=0; y<cvimg->imageSize; y+=cvimg->widthStep )
    for( x=0; x<img->width; x++,n++ )
#if MAGICKCORE_QUANTUM_DEPTH == 16
      ptr[y+x] = GetPixelGray(pixs+n) >> 8;
#elif MAGICKCORE_QUANTUM_DEPTH == 8
      ptr[y+x] = GetPixelGray(pixs+n);
#endif

  return SUCCESS;
}

int setpixels_magick_cv8u( Img* img, IplImage* cvimg ) {
  gray *gimg = NULL;
  if( malloc_grayv( img->width*img->height, &gimg, FALSE ) )
    return FAILURE;

  int x,y,n;
  for( y=0,n=0; y<cvimg->imageSize; y+=cvimg->widthStep )
    for( x=0; x<img->width; x++,n++ )
      gimg[n] = cvimg->imageData[y+x];

  int err = setpixels_magick_graym( img, gimg );

  free(gimg);

  return err;
}

int ccomp_contours_opencv( Img *img, CvMemStorage** _storage, CvSeq** _contours, int method, double summeps ) {
  IplImage* cvimg = cvCreateImage( cvSize(img->width,img->height), IPL_DEPTH_8U, 1 );
  getpixels_magick_cv8u( img, cvimg );

  if( *_storage == NULL )
    *_storage = cvCreateMemStorage(0);

  cvFindContours( cvimg, *_storage, _contours, sizeof(CvContour), CV_RETR_EXTERNAL, method, cvPoint(0,0) );

  if( summeps > 0.0 ) {
    CvMemStorage* storage = cvCreateMemStorage(0);
    *_contours = cvApproxPoly( *_contours, sizeof(CvContour), storage, CV_POLY_APPROX_DP, summeps, TRUE );
    cvReleaseMemStorage( _storage );
    *_storage = storage;
  }

  cvReleaseImage( &cvimg );

  return SUCCESS;
}

int join_ccomp_opencv( Img *img ) {
  IplImage* cvimg = cvCreateImage( cvSize(img->width,img->height), IPL_DEPTH_8U, 1 );
  getpixels_magick_cv8u( img, cvimg );

  //cvSaveImage( "cvout1.png", cvimg, 0 );

  CvSeq* contours = NULL;
  CvMemStorage* storage = cvCreateMemStorage(0);

  cvFindContours( cvimg, storage, &contours, sizeof(CvContour),
                  CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE, cvPoint(0,0) );

  getpixels_magick_cv8u( img, cvimg );

  //const char* attrs[] = { "recursive", "1", 0 };
  //cvSave("contours.xml", contours, NULL, NULL, cvAttrList(attrs, 0));
  //eval convert 1_01.png $(cat contours.xml | tr '\n' ' ' | xmlstarlet sel -t -m //data -v . -n | sed 's|^ *| |; s| *$||; s|   *| |g;' | awk '{printf(" -draw \"polyline");for(n=1;n<=NF;n+=2)printf(" %s,%s",$n,$(n+1));printf("\"");}') tmp.png

  int x,y,n;
  int C = 0;
  int N = 0;

  CvSeq* contn = contours;
  while( contn != NULL ) {
    C ++;
    N += contn->total;
    contn = contn->h_next;
  }

  //logger( 0, "C: %d", C );
  //logger( 0, "N: %d", N );

  CvSeq **conts = (CvSeq**)malloc(C*sizeof(CvSeq*));
  for( x=1,conts[0]=contours; x<C; x++ )
    conts[x] = conts[x-1]->h_next;

  int* ccid = NULL;
  malloc_I1v( C+N, &ccid, FALSE );
  int *nnccid = ccid + C;

  for( x=0; x<C; x++ )
    ccid[x] = x;

  for( x=0; x<C; x++ )
    if( ccid[x] == x ) { // test if not already joined
      int nx=0, ny=0;

      /// copy x-th component coordinates for test ///
      int xx;
      for( xx=0; xx<C; xx++ )
        nx += x == ccid[xx] ? conts[xx]->total : 0 ;

      if( nx == N )
        break;

      //cv::Mat ftest = cv::Mat::zeros( nx, 2, CV_32S );
      cv::Mat ftest = cv::Mat::zeros( nx, 2, CV_32F );
      int nn;
      for( xx=0,nn=0; xx<C; xx++ )
        if( x == ccid[xx] )
          for( n=0; n<conts[xx]->total; n++,nn++ ) {
            int* p = (int*)cvGetSeqElem(conts[xx], n);
            //ftest.at<int>(nn,0) = p[0];
            //ftest.at<int>(nn,1) = p[1];
            ftest.at<float>(nn,0) = p[0];
            ftest.at<float>(nn,1) = p[1];
          }

      /// copy all other component coordinates for reference ///
      //cv::Mat fref = cv::Mat::zeros( N-nx, 2, CV_32S );
      cv::Mat fref = cv::Mat::zeros( N-nx, 2, CV_32F );
      for( y=0,nn=0; y<C; y++ )
        if( x != ccid[y] )
          for( n=0; n<conts[y]->total; n++,nn++,ny++ ) {
            int* p = (int*)cvGetSeqElem(conts[y], n);
            //fref.at<int>(nn,0) = p[0];
            //fref.at<int>(nn,1) = p[1];
            fref.at<float>(nn,0) = p[0];
            fref.at<float>(nn,1) = p[1];
            nnccid[ny] = y;
          }

      /// find nearest neighbors ///
      cv::Mat nnidx( N-nx, 1, CV_32S );
      cv::Mat nndst( N-nx, 1, CV_32F );
      cv::flann::Index flann_index( fref,
        cv::flann::LinearIndexParams() );//,
        //cv::flann::KDTreeIndexParams(4),
//        cvflann::FLANN_DIST_EUCLIDEAN );

      flann_index.knnSearch( ftest, nnidx, nndst, 1, cv::flann::SearchParams() );

      float dst = nndst.at<float>(0,0);
      int i = 0;
      for( n=1; n<nx; n++ )
        if( dst > nndst.at<float>(n,0) ) {
          dst = nndst.at<float>(n,0);
          i = n;
        }

      /// relabel joined components ///
      int j = nnidx.at<int>(i,0);
      int c = nnccid[j];
      while( c != ccid[c] )
        //{ fprintf(stderr," %d=>%d",c,ccid[c]);
        c = ccid[c];
        //} fprintf(stderr,"\n");
      for( xx=0; xx<C; xx++ )
        if( x == ccid[xx] )
          //{ fprintf(stderr," %d(%d)->%d",xx,ccid[xx],c);
          ccid[xx] = c;
          //} fprintf(stderr,"\n");

      //logger( 0, "x=%d dst=%g p1=%g,%g p2=%g,%g", x, dst, ftest[2*i], ftest[2*i+1], fref[2*j], fref[2*j+1] );

      /// join components with a line ///
      //cvLine(cvimg, cvPoint(ftest.at<int>(i,0),ftest.at<int>(i,1)), cvPoint(fref.at<int>(j,0),fref.at<int>(j,1)), cvScalar(255,255,255,0), 1, 8, 0 );
      cvLine(cvimg, cvPoint(ftest.at<float>(i,0),ftest.at<float>(i,1)), cvPoint(fref.at<float>(j,0),fref.at<float>(j,1)), cvScalar(255,255,255,0), 1, 8, 0 );

      ftest.release();
      fref.release();
      nnidx.release();
      nndst.release();
      flann_index.release();
    }

  //cvSaveImage( "cvout2.png", cvimg, 0 );

  if( setpixels_magick_cv8u( img, cvimg ) )
    return FAILURE;

  free(conts);
  free(ccid);

  cvReleaseMemStorage( &storage );
  cvReleaseImage( &cvimg );

  return SUCCESS;
}
