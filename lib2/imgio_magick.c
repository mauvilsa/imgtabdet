/**
 * Functions for reading and writing images using MagickCore
 *
 * @version $Revision$$Date::             $
 * @copyright Copyright (c) 2014 to the present, Mauricio Villegas <mauvilsa@upv.es>
 */

#include "imgio_magick.h"

#include <flann/flann.h>

//#define CM_PER_IN 2.54
#define IN_PER_CM 0.39370078740157480314

void free_Img( Img *img ) {
  DestroyImageInfo( img->info );
  DestroyImage( img->image );
  free(img);
}

Img* create_img( ImageInfo* image_info, Image* image ) {
  ExceptionInfo *exception = AcquireExceptionInfo();
  unsigned char is_gray = IsGrayImage( image, exception );
  unsigned char is_opaque = IsOpaqueImage( image, exception );

  if( exception->severity != UndefinedException ) {
    CatchException( exception );
    DestroyExceptionInfo( exception );
    return NULL;
  }
  DestroyExceptionInfo( exception );

  //IdentifyImage(image,stderr,MagickTrue);

  Img* img = (Img*)calloc( 1, sizeof(Img) );
  if( img == NULL )
    return NULL;

  /* @todo convert to true color sRGB if not gray and not true color sRGB
     @todo change orientation to 0 if different */

  //if( image->orientation != UndefinedOrientation )

  if( image->units == PixelsPerInchResolution ) {
    image->units = PixelsPerCentimeterResolution;
    image_info->units = PixelsPerCentimeterResolution;
    image->x_resolution *= IN_PER_CM;
    image->y_resolution *= IN_PER_CM;
  }

  img->width = image->columns;
  img->height = image->rows;
  img->res_x = image->units != UndefinedResolution ? image->x_resolution : 0;
  img->res_y = image->units != UndefinedResolution ? image->y_resolution : 0;
  img->is_gray = is_gray;
  img->is_opaque = is_opaque;
  img->info = image_info;
  img->image = image;

  return img;
}

int readimg_magick( char* fname, Img** _img, FILE* logfile ) {
  ExceptionInfo *exception = AcquireExceptionInfo();
  ImageInfo *image_info = CloneImageInfo( NULL );
  CopyMagickString( image_info->filename, fname, MaxTextExtent );
  Image *image = ReadImage( image_info, exception );

  if( exception->severity != UndefinedException ) {
    CatchException( exception );
    DestroyImageInfo( image_info );
    DestroyExceptionInfo( exception );
    return FAILURE;
  }
  DestroyExceptionInfo( exception );

  if( image == NULL ) {
    DestroyImageInfo( image_info );
    fprintf( logfile, "%s: error: unable to read image: %s\n", "readimg_magick", fname );
    return FAILURE;
  }

  *_img = create_img( image_info, image );

  if( *_img == NULL ) {
    DestroyImageInfo( image_info );
    DestroyImage( image );
    fprintf( logfile, "%s: error: problems creating Img structure\n", "readimg_magick" );
    return FAILURE;
  }

  return SUCCESS;
}

int scanimg_magick( FILE* file, Img** _img, FILE* logfile ) {
  ExceptionInfo *exception = AcquireExceptionInfo();
  ImageInfo *image_info = CloneImageInfo( (ImageInfo*)NULL );
  image_info->file = file;
  Image *image = ReadImage( image_info, exception );

  if( exception->severity != UndefinedException ) {
    CatchException( exception );
    DestroyImageInfo( image_info );
    DestroyExceptionInfo( exception );
    return FAILURE;
  }
  DestroyExceptionInfo( exception );

  if( image == (Image*)NULL ) {
    DestroyImageInfo( image_info );
    fprintf( logfile, "%s: error: unable to read image\n", "scanimg_magick" );
    return FAILURE;
  }

  *_img = create_img( image_info, image );

  if( *_img == NULL ) {
    DestroyImageInfo( image_info );
    DestroyImage( image );
    fprintf( logfile, "%s: error: problems creating Img structure\n", "scanimg_magick" );
    return FAILURE;
  }

  return SUCCESS;
}

int writeimg_magick( char* fname, Img* img, FILE* logfile ) {
  CopyMagickString( img->image->filename, fname, MaxTextExtent );
  if( WriteImage( img->info, img->image ) == MagickFalse ) {
    fprintf( logfile, "%s: error: unable to write image: %s\n", "writeimg_magick", fname );
    return FAILURE;
  }

  return SUCCESS;
}

int printimg_magick( FILE* file, char* format, Img* img, FILE* logfile ) {
  char out[16];
  sprintf(out,"%s:out",format);

  img->info->file = file;
  CopyMagickString( img->image->filename, out, MaxTextExtent );

  if( WriteImage( img->info, img->image ) == MagickFalse ) {
    fprintf( logfile, "%s: error: unable to write image\n", "printimg_magick" );
    return FAILURE;
  }

  return SUCCESS;
}

int getpixels_magick_graym( Img* img, gray** gimg ) {
  ExceptionInfo *exception = AcquireExceptionInfo();
  const PixelPacket *pixs =
    GetVirtualPixels( img->image, 0, 0, img->width, img->height, exception );

  if( exception->severity != UndefinedException ) {
    CatchException( exception );
    DestroyExceptionInfo( exception );
    return FAILURE;
  }
  DestroyExceptionInfo( exception );

  int n;
  for( n=img->width*img->height-1; n>=0; n-- )
    gimg[0][n] = GetPixelGray(pixs+n) >> 8; // will only work with Q16

  return SUCCESS;
}

int setpixels_magick_graym( Img* img, gray** gimg ) {
  ExceptionInfo *exception = AcquireExceptionInfo();

  /// @todo this discards all previous metadata, how to preserve it? ///
  if( img->image != NULL )
    DestroyImage( img->image );
  img->image = ConstituteImage( img->width, img->height, "I", CharPixel, gimg[0], exception );

/*
  /// @todo only works by creating new image, why ??? ///
  if( img->image != NULL )
    DestroyImage( img->image );
  MagickPixelPacket black;
  QueryMagickColor("black", &black, exception);
  if( exception->severity != UndefinedException ) {
    CatchException( exception );
    DestroyExceptionInfo( exception );
    return FAILURE;
  }
  img->image = NewMagickImage( img->info, img->width, img->height, &black );

  PixelPacket *pixs =
    QueueAuthenticPixels( img->image, 0, 0, img->width, img->height, exception );
  //  GetAuthenticPixels( img->image, 0, 0, img->width, img->height, exception );

  if( exception->severity != UndefinedException ) {
    CatchException( exception );
    DestroyExceptionInfo( exception );
    return FAILURE;
  }

  int n;
  for( n=img->width*img->height-1; n>=0; n-- )
    SetPixelGray( pixs+n, ((int)gimg[0][n]) << 8 ); // will only work with Q16

  if( SyncAuthenticPixels( img->image, exception ) != MagickTrue )
    return FAILURE;
*/

  if( exception->severity != UndefinedException ) {
    CatchException( exception );
    DestroyExceptionInfo( exception );
    return FAILURE;
  }

  DestroyExceptionInfo( exception );

  return SUCCESS;
}

int getpixels_magick_cv8u( Img* img, IplImage* cvimg ) {
  ExceptionInfo *exception = AcquireExceptionInfo();
  const PixelPacket *pixs =
    GetVirtualPixels( img->image, 0, 0, img->width, img->height, exception );

  if( exception->severity != UndefinedException ) {
    CatchException( exception );
    DestroyExceptionInfo( exception );
    return FAILURE;
  }
  DestroyExceptionInfo( exception );

  gray *ptr = (gray*)cvimg->imageData;

  int x,y,n;
  for( y=0,n=0; y<cvimg->imageSize; y+=cvimg->widthStep )
    for( x=0; x<img->width; x++,n++ )
      ptr[y+x] = GetPixelGray(pixs+n) >> 8; // will only work with Q16

  return SUCCESS;
}

int setpixels_magick_cv8u( Img* img, IplImage* cvimg ) {
  ExceptionInfo *exception = AcquireExceptionInfo();

  /// @todo this discards all previous metadata, how to preserve it? ///
  /*if( img->image != NULL )
    DestroyImage( img->image );
  img->image = ConstituteImage( img->width, img->height, "I", CharPixel, (gray*)cvimg->imageData, exception );*/


  /// @todo only works by creating new image, why ??? ///
  if( img->image != NULL )
    DestroyImage( img->image );
  MagickPixelPacket black;
  QueryMagickColor("black", &black, exception);
  if( exception->severity != UndefinedException ) {
    CatchException( exception );
    DestroyExceptionInfo( exception );
    return FAILURE;
  }
  img->image = NewMagickImage( img->info, img->width, img->height, &black );

  PixelPacket *pixs =
    QueueAuthenticPixels( img->image, 0, 0, img->width, img->height, exception );
  //  GetAuthenticPixels( img->image, 0, 0, img->width, img->height, exception );

  if( exception->severity != UndefinedException ) {
    CatchException( exception );
    DestroyExceptionInfo( exception );
    return FAILURE;
  }

  gray *ptr = (gray*)cvimg->imageData;

  int x,y,n;
  for( y=0,n=0; y<cvimg->imageSize; y+=cvimg->widthStep )
    for( x=0; x<img->width; x++,n++ )
      SetPixelGray( pixs+n, ((int)ptr[y+x]) << 8 ); // will only work with Q16

  if( SyncAuthenticPixels( img->image, exception ) != MagickTrue )
    return FAILURE;

  if( exception->severity != UndefinedException ) {
    CatchException( exception );
    DestroyExceptionInfo( exception );
    return FAILURE;
  }

  DestroyExceptionInfo( exception );

  return SUCCESS;
}

int togray_magick( Img* img ) {
  if( img->is_gray )
    return SUCCESS;

  TransformImageColorspace( img->image, GRAYColorspace );

  ExceptionInfo *exception = AcquireExceptionInfo();
  img->is_gray = IsGrayImage( img->image, exception );

  if( exception->severity != UndefinedException ) {
    CatchException( exception );
    DestroyExceptionInfo( exception );
    return FAILURE;
  }
  DestroyExceptionInfo( exception );

  return img->is_gray ? SUCCESS : FAILURE;
}

int add_border_magick( Img* img, int size, PixelPacket color ) {
  ExceptionInfo *exception = AcquireExceptionInfo();
  const RectangleInfo bord = { size, size, 0, 0 };

  img->image->border_color = color;

  Image *bimg = BorderImage( img->image, &bord, exception );
  if( exception->severity != UndefinedException ) {
    CatchException( exception );
    DestroyExceptionInfo( exception );
    return FAILURE;
  }
  DestroyExceptionInfo( exception );

  DestroyImage( img->image );
  img->image = bimg;
  img->width = bimg->columns;
  img->height = bimg->rows;

  return SUCCESS;
}
