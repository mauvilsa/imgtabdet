/**
 * Tool that tries to find line grids (i.e. a table) in an image
 *
 * @version $Revision$$Date::             $
 * @copyright Copyright (c) 2016 to the present, Mauricio Villegas <mauvilsa@upv.es>
 */

/*** Includes *****************************************************************/
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <limits>
#include <algorithm>
#include <set>

#include "log.h"
#include "iomat.h"
#include "imgio_magick.h"

#ifdef ENABLE_OPENCV_LSD
#include "imgproc_cv.h"
#include <opencv2/line_descriptor.hpp>
using namespace cv::line_descriptor;
#endif

#include <opencv2/core/version.hpp>
#if CV_MAJOR_VERSION == 2
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp>
#else
#include <opencv2/imgproc.hpp>
#include <opencv2/core.hpp>
#endif

using namespace cv;
using namespace std;

/*** Definitions **************************************************************/
static char tool[] = "imgtabdet";
static char revnum[] = "$Revision$";
static char revdate[] = "$Date$";

#define SEGDIR_DOWN  0
#define SEGDIR_UP    1
#define SEGDIR_LEFT  2
#define SEGDIR_RIGHT 3
char const *dirname[] = { "down", "up", "left", "right" };

#define DIR_VERT 0
#define DIR_HORZ 1

#define TABL_TOP    0
#define TABL_BOTTOM 1
#define TABL_LEFT   2
#define TABL_RIGHT  3
char const *tabname[] = { "top", "bottom", "left", "right" };

#define SORT_PAIR_SAME_NOSUM 256

char *ifn = NULL;
char const *gb_units = "pixels";

#define OUT_ASCII 0
#define OUT_XMLPAGE 1

int gb_format = OUT_ASCII;

bool gb_opencv_lsd = false;

double gb_pair_max_width = 5.0; // @todo make it depend on resolution
double gb_join_max_sep = 100.0; // @todo make it depend on resolution
double gb_join_max_parallel = 4.0; // @todo make it depend on resolution
double gb_join_min_dot = 0.999657324976; // 1.5 degrees
double gb_join_min_short_dot = 0.999048221582; // 2.5 degrees
//double gb_join_min_sep_dot = 0.965925826289; // 15.0 degrees
double gb_join_min_sep_dot = 0.939692620786; // 20.0 degrees

double gb_margrow_top[2] = { 0.0, 0.0 };
double gb_margrow_bottom[2] = { 0.0, 0.0 };
double gb_margcol_left[2] = { 0.0, 0.0 };
double gb_margcol_right[2] = { 0.0, 0.0 };

double gb_twopage = 0.0;
int tabdef[2][2];

//char gb_procimgs = FALSE;
double gb_density = FALSE;

FILE *logfile = NULL;
int verbosity = 1;

/*** Functions ****************************************************************/
void print_usage( FILE *file ) {
  fprintf( file, "Description: Tries to find a table structure in an image using the line segment detector (LSD) algorithm\n" );
  fprintf( file, "Usage: %s [options] <in_img>\n", tool );
  fprintf( file, "\n" );
#ifdef ENABLE_OPENCV_LSD
  fprintf( file, "  -O              Use the opencv LSD implementation (def.=%s)\n", strbool(gb_opencv_lsd) );
#endif
  fprintf( file, "  -F FORMAT       Output format among: 'ascii', 'xmlpage' (def.=ascii)\n" );
  fprintf( file, "  -T tabdef       Table definition, {rows}x{cols}[,{rows}x{cols}] (def.=auto)\n" );
  fprintf( file, "  -D xsplit       Set double page with split at x=xsplit (def.=%s)\n", strbool(gb_twopage) );
  fprintf( file, "  -S sep          Maximum length of join segments (def.=%g %s)\n", gb_join_max_sep, gb_units );
  fprintf( file, "  -P pdist        Maximum parallel separation to consider joining (def.=%g %s)\n", gb_join_max_parallel, gb_units );
  fprintf( file, "  -A angle        Maximum angle between segments to consider joining (def.=%g degrees)\n", acos(gb_join_min_dot)*180.0/M_PI );
  fprintf( file, "  -a sangle       Maximum angle between segments to consider joining when join is shorter than the smallest segment (def.=%g degrees)\n", acos(gb_join_min_short_dot)*180.0/M_PI );
  fprintf( file, "  -J jangle       Maximum angle between segments and separation segment to consider joining (def.=%g degrees)\n", acos(gb_join_min_sep_dot)*180.0/M_PI );
  fprintf( file, "  -M mwidth       Maximum width for lines (def.=%g %s)\n", gb_pair_max_width, gb_units );

  fprintf( file, "  -m ({pg}[tblr]{fact})+  Add margin rows/columns as a percentage of height/width (def.=false)\n" );

  //fprintf( file, "  -u (mm|pixels)  Units for ALL distance parameters (def.=%s)\n", gb_units );
  fprintf( file, "  -d density      Specify the image density in pixels per cm (def.=%s)\n", strbool(gb_density) );
  //fprintf( file, "  -p              Save images %s_*.png of processing steps (def.=%s)\n", tool, strbool(gb_procimgs) );
  fprintf( file, "  -l lfile        Logging to 'lfile' (def.=stderr)\n" );
  fprintf( file, "  -V (-|+|level)  Set verbosity level (def.=%d)\n", verbosity );
  fprintf( file, "  -h              Print this usage information and exit\n" );
  fprintf( file, "  -v              Print version and exit\n" );
  fprintf( file, "\n");
  print_svn_rev( file );
}

/**
 * Function that creates a temporal file using the mktemp command
 */
void mktemp( const char* tempbase, char *tempname ) {
  char cmd[FILENAME_MAX];
  sprintf( cmd, "mktemp %s_XXXXXXXX", tempbase );
  FILE *p = popen( cmd, "r" );
  if( p != NULL ) {
    sprintf( cmd, "%%%ds\n", FILENAME_MAX-1 );
    if( fscanf( p, cmd, tempname ) != 1 )
      tempname[0] = '\0';
    pclose(p);
  }
}

/**
 * Function that returns the index of the minimum value in a vector
 */
template<typename T>
int argMin( const vector<T> vect ) {
  if( vect.size() == 0 )
    return -1;
  return min_element( vect.begin(), vect.end() ) - vect.begin() ;
}
template int argMin( const vector<double> vect );
template int argMin( const vector<float> vect );

/**
 * Function that returns the indexes of the sorted input but only for values within the limit if provided and finite, i.e. not NaN, Inf or -Inf. Optionally it can give a vector with the sorted values.
 */
template<typename T>
vector<int> sortIdxLimit( const vector<T> &vect, int flags, double limit = numeric_limits<double>::quiet_NaN(), vector<T> *_sorted = NULL, int keep_at_least = 0 ) {
  if( vect.size() == 0 ) {
    if( _sorted != NULL )
      *_sorted = vector<T>(0);
    return vector<int>(0);
  }
  vector<int> idx( vect.size() );
  vector<int> fidx( vect.size() );
  /*printf( "srt: limit=%g flags=%d\n", limit, flags );
  printf( " orig:" );
  for( size_t n=0; n<vect.size(); n++ )
    printf( " %g", vect[n] );
  printf( "\n" );*/
  sortIdx( vect, idx, flags );
  bool ascending = flags & CV_SORT_DESCENDING ? false : true ;
  if( isnan(limit) )
    limit = ascending ?
      numeric_limits<double>::infinity() :
      -numeric_limits<double>::infinity() ;
  int num = 0;
  for( int n=0; n<int(vect.size()); n++ )
    if( isfinite(double(vect[idx[n]])) && 
        ( ( ascending && double(vect[idx[n]]) < limit ) ||
          ( ! ascending && double(vect[idx[n]]) > limit ) ||
          n < keep_at_least ) )
      fidx[num++] = idx[n];
  if( num < int(vect.size()) )
    fidx.resize(num);
  /*printf( " srt_idx:" );
  for( size_t n=0; n<fidx.size(); n++ )
    printf( " %d", fidx[n] );
  printf( "\n" );*/
  if( _sorted != NULL ) {
    vector<T> sorted( num );
    for( int n=0; n<num; n++ )
      sorted[n] = vect[fidx[n]];
    /*printf( " srt_val:" );
    for( size_t n=0; n<sorted.size(); n++ )
      printf( " %g", sorted[n] );
    printf( "\n" );*/
    *_sorted = sorted;
  }
  return fidx;
}
template vector<int> sortIdxLimit( const vector<double> &vect, int flags, double limit = numeric_limits<double>::quiet_NaN(), vector<double> *_sorted = NULL, int keep_at_least = 0 );

/**
 * Function that returns the pairs of indexes of the sorted sum of two input vectors but only for finite values within the limit if provided. Optionally it can give a vector with the sorted sum values.
 */
template<typename T>
vector<Point2i> pairSumSortIdxLimit( const vector<T> &vect1, const vector<T> &vect2, int flags, double limit = numeric_limits<double>::quiet_NaN(), vector<T> *_sorted = NULL ) {
  vector<T> vsum;
  vector<Point2i> pidx;
  bool same_nosum = ( &vect1 == &vect2 && flags&SORT_PAIR_SAME_NOSUM ) ? true : false;
  //printf( "pairSumSortIdxLimit: p1=%p p2=%p same_nosum=%s\n", &vect1, &vect2, strbool(same_nosum) );
  for( size_t n=0; n<vect1.size(); n++ )
    for( size_t m=0; m<vect2.size(); m++ ) {
      double sum = ( m==n && same_nosum ) ? vect1[n] : vect1[n] + vect2[m];
      if( isfinite(double(sum)) ) {
        vsum.push_back( sum );
        pidx.push_back( Point2i(n,m) );
      }
    }
  vector<int> idx = sortIdxLimit( vsum, flags&0xff, limit, _sorted );
  vector<Point2i> fidx( idx.size() );
  for( size_t n=0; n<idx.size(); n++ )
    fidx[n] = pidx[idx[n]];
  //for( size_t n=0; n<fidx.size(); n++ )
  //  printf( " %zu(%d,%d): %g %g -> %g\n", n, fidx[n].x, fidx[n].y, vect1[fidx[n].x], vect2[fidx[n].y], (*_sorted)[n] );
  return fidx;
}
template vector<Point2i> pairSumSortIdxLimit( const vector<double> &vect1, const vector<double> &vect2, int flags, double limit = numeric_limits<double>::quiet_NaN(), vector<double> *_sorted = NULL );


/**
 * Function that returns the index of the first occurrence value in vector or -1 if value not found
 */
template<typename T>
int findValue( vector<T> vect, T value ) {
  vector<int>::iterator iter = find( vect.begin(), vect.end(), value );
  return iter == vect.end() ? -1 : iter - vect.begin() ;
}
template int findValue( vector<int>, int value );

/**
 * Function that checks if a point is within a line segment
 */
template<typename T>
bool pointInSegment( const Point_<T> &segm_start, const Point_<T> &segm_end, const Point_<T> &point ) {
  Point_<T> segm = segm_end - segm_start;
  Point_<T> start_point = segm_start - point;
  Point_<T> end_point = segm_end - point;
  return 1.0001*segm.ddot(segm) >= start_point.ddot(start_point) + end_point.ddot(end_point);
}
template bool pointInSegment<float>( const Point2f &segm_start, const Point2f &segm_end, const Point2f &point );
template bool pointInSegment<double>( const Point2d &segm_start, const Point2d &segm_end, const Point2d &point );

/**
 * Function that finds the intersection point between two lines defined by pairs of points or returns false if no intersection
 */
template<typename T>
bool intersection( const Point_<T> line1_point1, const Point_<T> line1_point2, const Point_<T> line2_point1, const Point_<T> line2_point2, Point_<T> &ipoint ) {
  Point_<T> x = line2_point1 - line1_point1;
  Point_<T> direct1 = line1_point2 - line1_point1;
  Point_<T> direct2 = line2_point2 - line2_point1;

  double cross = direct1.x*direct2.y - direct1.y*direct2.x;
  if( fabs(cross) < /*EPS*/1e-8 )
    return false;

  double t1 = (x.x * direct2.y - x.y * direct2.x)/cross;
  ipoint = line1_point1 + direct1 * t1;

  return true;
}
template bool intersection( const Point2f line1_point1, const Point2f line1_point2, const Point2f line2_point1, const Point2f line2_point2, Point2f &ipoint );
template bool intersection( const Point2d line1_point1, const Point2d line1_point2, const Point2d line2_point1, const Point2d line2_point2, Point2d &ipoint );

/**
 * Function that computes the vector from a point to a line defined by two points
 */
template<typename T>
void pointToLine( const Point_<T> &line_point1, const Point_<T> &line_point2, const Point_<T> &point, Point_<T> &point_to_line ) {
  Point_<T> line_norm = line_point2 - line_point1;
  double line_length = norm(line_norm);
  line_norm = line_norm * (1.0/line_length);
  Point_<T> point1_to_point = line_point1 - point;
  double point1_to_point_proj = point1_to_point.ddot(line_norm);
  point_to_line = point1_to_point - point1_to_point_proj*line_norm;
}
template void pointToLine<double>( const Point2d &line_point1, const Point2d &line_point2, const Point2d &point, Point_<double> &point_to_line );
template void pointToLine<float>( const Point2f &line_point1, const Point2f &line_point2, const Point2f &point, Point_<float> &point_to_line );

/**
 * Function that computes a point on a line that extends beyond a segment a factor of its length
 */
template<typename T>
void extendSegment( const Point_<T> &segment1, const Point_<T> &segment2, double factor, Point_<T> &point ) {
  Point_<T> segment = segment2 - segment1;
  double length = norm(segment);
  segment.x /= length;
  segment.y /= length;
  point = ( factor < 0 ? segment1 : segment2 ) + factor * length * segment;
}
template void extendSegment( const Point2f &segment1, const Point2f &segment2, double factor, Point2f &point );

/**
 * Function that fits a set of points to a line segment
 */
template<typename T>
void fitSegment( const vector<Point_<T> > &points, vector<Point_<T> > &segment, int distType, double param, double reps, double aeps ) {
  Vec4f line;
  fitLine(points,line,distType,param,reps,aeps);

  Point_<T> line_direct( line[0], line[1] );
  Point_<T> line_point1( line[2], line[3] );
  Point_<T> line_point2 = line_point1 + line_direct;
  Point_<T> aux;
  double proj_min = 0.0;
  double proj_max = 0.0;
  for( size_t n=0; n<points.size(); n++ ) {
    pointToLine( line_point1, line_point2, points[n], aux );
    aux = aux + points[n];
    double proj = line_direct.ddot( aux - line_direct );
    if( n==0 || proj < proj_min ) {
      proj_min = proj;
      segment[0] = aux;
    }
    if( n==0 || proj > proj_max ) {
      proj_max = proj;
      segment[1] = aux;
    }
  }

  line_point1 = segment[0] - points[0];
  line_point2 = segment[1] - points[0];
  if( line_point1.ddot(line_point1) > line_point2.ddot(line_point2) ) {
    aux = segment[0];
    segment[0] = segment[1];
    segment[1] = aux;
  }
}
template void fitSegment<float>( const vector<Point2f> &points, vector<Point2f> &segment, int distType, double param, double reps, double aeps );

/**
 * Function that samples points along a poly-segment
 */
template<typename T>
vector<Point_<T> > samplePolySegment( const vector<Point_<T> > &polysegm, const double samp_sep )
  {
  if( polysegm.size() == 0 )
    return polysegm;

  double totlgth = 0.0;
  double lgth[polysegm.size()-1];
  for( int n=0; n<int(polysegm.size())-1; n++ ) {
    Point_<T> segm = polysegm[n+1] - polysegm[n];
    lgth[n] = norm(segm);
    totlgth += lgth[n];
  }
  double sampsep = totlgth/double(int(0.5+totlgth/samp_sep));

  vector<Point_<T> > sampled;
  sampled.push_back( polysegm[0] );
  double samppos = sampsep;
  for( int n=0; n<int(polysegm.size())-1; n++ ) {
    double dx = ( polysegm[n+1].x - polysegm[n].x ) / lgth[n];
    double dy = ( polysegm[n+1].y - polysegm[n].y ) / lgth[n];
    while( true ) {
      if( n == int(polysegm.size())-2 &&
          lgth[n] - samppos < 0.5*sampsep )
        break;
      if( samppos > lgth[n] ) {
        samppos -= lgth[n];
        break;
      }
      Point_<T> sample( polysegm[n].x+dx*samppos, polysegm[n].y+dy*samppos );
      sampled.push_back( sample );
      samppos += sampsep;
    }
  }
  sampled.push_back( polysegm[polysegm.size()-1] );
  return sampled;
}
template vector<Point2f> samplePolySegment( const vector<Point2f> &polysegm, const double samp_sep );

/**
 * Function that checks if two segments overlap
 */
template<typename T>
bool segmentsOverlap( const Point_<T> &segm1_start, const Point_<T> &segm1_end, const Point_<T> &segm2_start, const Point_<T> &segm2_end ) {
  return ( pointInSegment(segm1_start,segm1_end,segm2_start) ||
           pointInSegment(segm1_start,segm1_end,segm2_end) ||
           pointInSegment(segm2_start,segm2_end,segm1_start) ||
           pointInSegment(segm2_start,segm2_end,segm1_end) );
}
template bool segmentsOverlap( const Point_<float> &segm1_start, const Point_<float> &segm1_end, const Point_<float> &segm2_start, const Point_<float> &segm2_end );
template bool segmentsOverlap( const Point_<double> &segm1_start, const Point_<double> &segm1_end, const Point_<double> &segm2_start, const Point_<double> &segm2_end );

/**
 * Function for selecting segments contained in a given range
 */
template<typename T>
void filterSegmentsRange( const vector<vector<Point_<T> > > segments, const int xmin, const int ymin, const int width, int const height, vector<vector<Point_<T> > > *_filtered, vector<int> *_indexes ) {

  vector<int> indexes;
  vector<vector<Point_<T> > > filtered;

  int xmax = xmin + width - 1;
  int ymax = ymin + height - 1;

  for( size_t n=0; n<segments.size(); n++ ) {
    Point_<T> p1 = segments[n][0];
    Point_<T> p2 = segments[n][1];
    if( p1.x >= xmin && p1.x <= xmax &&
        p2.x >= xmin && p2.x <= xmax &&
        p1.y >= ymin && p1.y <= ymax &&
        p2.y >= ymin && p2.y <= ymax ) {
      if( _filtered != NULL )
        filtered.push_back( segments[n] );
      if( _indexes != NULL )
        indexes.push_back( n );
    }
  }

  if( _filtered != NULL )
    *_filtered = filtered;
  if( _indexes != NULL )
    *_indexes = indexes;
}
template void filterSegmentsRange( const vector<vector<Point2f> > segments, const int xmin, const int ymin, const int width, int const height, vector<vector<Point2f> > *_filtered, vector<int> *_indexes );

/**
 * Function that distributes segments among the 4 directions: down, up, left, right
 */
template<typename T>
void segmentsPerDirection( const vector<vector<Point_<T> > > segments, vector<vector<vector<Point_<T> > > > *_perdirection, vector<vector<int> > *_indexes ) {

  vector<vector<vector<Point_<T> > > > perdirection( 4 );
  vector<vector<int> > indexes( 4 );
  Point_<T> right( 1.0, 0.0 );

  for( size_t n=0; n<segments.size(); n++ ) {
    Point_<T> p1 = segments[n][0];
    Point_<T> p2 = segments[n][1];
    Point_<T> segmdir = p2 - p1;
    double cos_theta = segmdir.ddot(right)/norm(segmdir);

    int direction;
    if( cos_theta > M_SQRT1_2 )
      direction = SEGDIR_RIGHT;
    else if( cos_theta < -M_SQRT1_2 )
      direction = SEGDIR_LEFT;
    else if( segmdir.y > 0 )
      direction = SEGDIR_DOWN;
    else
      direction = SEGDIR_UP;

    logger( 2, "info: n=%d: %s %d: %.1f,%.1f %.1f,%.1f", int(n), dirname[direction], int(perdirection[direction].size())+1, p1.x, p1.y, p2.x, p2.y );

    if( _perdirection != NULL )
      perdirection[direction].push_back( segments[n] );
    if( _indexes != NULL )
      indexes[direction].push_back( n );
  }

  if( _perdirection != NULL )
    *_perdirection = perdirection;
  if( _indexes != NULL )
    *_indexes = indexes;
}
template void segmentsPerDirection( const vector<vector<Point2f> > segments, vector<vector<vector<Point2f> > > *_perdirection, vector<vector<int> > *_indexes );

/**
 * Function that checks if two segments are aligned
 */
template<typename T>
bool segmentsAligned( const vector<Point_<T> > segm1, const vector<Point_<T> > segm2 ) {
  Point2d pi1( segm1[0].x, segm1[0].y );
  Point2d pi2( segm1[1].x, segm1[1].y );
  Point2d pj1( segm2[0].x, segm2[0].y );
  Point2d pj2( segm2[1].x, segm2[1].y );
  Point2d si = pi2 - pi1;
  Point2d sj = pj2 - pj1;
  Point2d sij = pj1 - pi2;
  double mag_i = norm(si);
  double mag_j = norm(sj);

  /// Distance between end of i and start of j ///
  double sep_ij = norm(sij);

  /// Dot product between segments ///
  double dot = si.ddot(sj)/(mag_i*mag_j);
  double dotsep = 0.5*( si.ddot(sij)/(mag_i*sep_ij) + sj.ddot(sij)/(mag_j*sep_ij) );

  /// Average of: ///
  /// - Distance between end point of segment i and line j ///
  /// - Distance between start point of segment j and line i ///
  /*Point2d i_to_j;
  Point2d j_to_i;
  pointToLine( pj1, pj2, pi2, i_to_j );
  pointToLine( pi1, pi2, pj1, j_to_i );
  double dst_ij_ji = 0.5*( norm(i_to_j) + norm(j_to_i) );*/

  double minlgth = mag_i < mag_j ? mag_i : mag_j ;

  /// ??? ///
  double gb_join_rect_fact = 0.3;
  double min_sep_dot = gb_join_min_sep_dot;
  if( sep_ij <= gb_join_max_parallel ) {
    // @todo create a function that does this line {inter,extra}polation ?
    double slope = gb_join_min_sep_dot / ( gb_join_max_parallel*(1.0-gb_join_rect_fact) );
    double ycept = -slope * gb_join_max_parallel * gb_join_rect_fact;
    min_sep_dot = slope * sep_ij + ycept;
  }

  /// Consider aligned depending on restrictions ///
  if( ( ( dot > gb_join_min_dot ) ||
        ( dot > gb_join_min_short_dot && sep_ij < minlgth ) ) &&
      dotsep > min_sep_dot /*&&
      sep_ij < gb_join_max_sep &&
      dst_ij_ji < gb_join_max_parallel*/ )
    return true;

  return false;
}
template bool segmentsAligned( const vector<Point2f> segm1, const vector<Point2f> segm2 );

/**
 * Function for joining almost aligned segments
 */
template<typename T>
void joinAlignedSegments( const vector<vector<Point_<T> > > segments, vector<vector<int> > *_polysegm_indexes, vector<double> *_polysegm_lengths ) {

  vector<vector<int> > polysegm_indexes;
  vector<double> polysegm_lengths;

  int num_segms = segments.size();
  double segm_norm[segments.size()];
  for( int i=0; i<num_segms; i++ ) {
    Point2d pi1( segments[i][0].x, segments[i][0].y );
    Point2d pi2( segments[i][1].x, segments[i][1].y );
    segm_norm[i] = norm( pi2 - pi1 );
  }

  vector<int> joined_to_index( num_segms, -1 );
  vector<double> joined_to_score_dist( num_segms, -1.0 );
  vector<double> joined_to_score_dot( num_segms, -1.0 );

  for( int i=0; i<num_segms; i++ ) {
    Point2d pi1( segments[i][0].x, segments[i][0].y );
    Point2d pi2( segments[i][1].x, segments[i][1].y );
    Point2d si = pi2 - pi1;
    double mag_i = segm_norm[i];

    for( int j=0; j<num_segms; j++ )
      if( i != j ) {
        Point2d pj1( segments[j][0].x, segments[j][0].y );
        Point2d pj2( segments[j][1].x, segments[j][1].y );
        Point2d sj = pj2 - pj1;
        Point2d sij = pj1 - pi2;
        double mag_j = segm_norm[j];

        /// Distance between end of i and start of j ///
        double sep_ij = norm(sij);

        /// Dot product between segments ///
        double dot = si.ddot(sj)/(mag_i*mag_j);
        double dotsep = 0.5*( si.ddot(sij)/(mag_i*sep_ij) + sj.ddot(sij)/(mag_j*sep_ij) );

        /// Average of: ///
        /// - Distance between end point of segment i and line j ///
        /// - Distance between start point of segment j and line i ///
        Point2d i_to_j;
        Point2d j_to_i;
        pointToLine( pj1, pj2, pi2, i_to_j );
        pointToLine( pi1, pi2, pj1, j_to_i );
        double dst_ij_ji = 0.5*( norm(i_to_j) + norm(j_to_i) );

        //logger( 4, "i=%d j=%d sep=%g dot=%g dotsep=%g dst_ij_ji=%g", i, j, sep_ij, dot, dotsep, dst_ij_ji );

        double minlgth = mag_i < mag_j ? mag_i : mag_j ;

        /// ??? ///
        double gb_join_rect_fact = 0.3;
        double min_sep_dot = gb_join_min_sep_dot;
        if( sep_ij <= gb_join_max_parallel ) {
          // @todo create a function that does this line {inter,extra}polation ?
          double slope = gb_join_min_sep_dot / ( gb_join_max_parallel*(1.0-gb_join_rect_fact) );
          double ycept = -slope * gb_join_max_parallel * gb_join_rect_fact;
          min_sep_dot = slope * sep_ij + ycept;
          //logger( 0, "i=%d j=%d min_sep_dot=%g dotsep=%g sep=%g", i, j, min_sep_dot, dotsep, sep_ij );
        }

        /// Candidate join depending on a few restrictions on the computed relations ///
        if( ( ( dot > gb_join_min_dot ) ||
              ( dot > gb_join_min_short_dot && sep_ij < minlgth ) ) &&
            dotsep > min_sep_dot &&
            sep_ij < gb_join_max_sep &&
            dst_ij_ji < gb_join_max_parallel ) {

          double score_dot = ( (mag_i+mag_j)*(1.0-dot) + (sep_ij)*(1.0-dotsep) )
                             / ( mag_i+mag_j+sep_ij );
          double score_dist = sep_ij + dst_ij_ji;

          // @todo better criterion combining dist and dot ?
          if( joined_to_score_dist[i] == -1.0 ||
              score_dist < joined_to_score_dist[i] ) {
            joined_to_index[i] = j;
            joined_to_score_dist[i] = score_dist;
            joined_to_score_dot[i] = score_dot;
          }
        }

      } // for( int j=0; j<num_segms; j++ ) if( i != j ) {
  } // for( int i=0; i<num_segms; i++ ) {

  /// Keep only one join per start point of segments ///
  for( int i=0; i<num_segms; i++ )
    if( joined_to_index[i] != -1 ) {
      int ibest = i;
      int j = joined_to_index[ibest];
      double score_dist = joined_to_score_dist[ibest];
      //double score_dot = joined_to_score_dot[ibest];
      for( int ii=i+1; ii<num_segms; ii++ )
        if( j == joined_to_index[ii] ) {
          if( score_dist > joined_to_score_dist[ii] ) {
            joined_to_index[ibest] = -1;
            //joined_to_score_dist[ibest] = -1.0;
            //joined_to_score_dot[ibest] = -1.0;
            ibest = ii;
            score_dist = joined_to_score_dist[ibest];
            //score_dot = joined_to_score_dot[ibest];
          }
          else {
            joined_to_index[ii] = -1;
            //joined_to_score_dist[ii] = -1.0;
            //joined_to_score_dot[ii] = -1.0;
          }
        }
    }

  /// Get indexes and lengths of joined poly-segments ///
  for( int i=0; i<num_segms; i++ ) {
    bool startsegm = true;
    for( int j=0; j<num_segms; j++ )
      if( i == joined_to_index[j] ) {
        startsegm = false;
        break;
      }
    if( ! startsegm )
      continue;

    vector<int> indexes;
    double length = 0.0;
    Point2d prev;

    //Point2d pi1( segments[i][0].x, segments[i][0].y );

    int j = i;
    do {
      indexes.push_back( j );
      if( i != j ) {
        Point2d sep = prev - Point2d( segments[j][0].x, segments[j][0].y );
        length += norm(sep);
      }
      length += segm_norm[j];
      prev = Point2d( segments[j][1].x, segments[j][1].y );
      j = joined_to_index[j];
    } while( j != -1 );

    polysegm_indexes.push_back( indexes );
    polysegm_lengths.push_back( length );
  }

  if( _polysegm_indexes != NULL )
    *_polysegm_indexes = polysegm_indexes;
  if( _polysegm_lengths != NULL )
    *_polysegm_lengths = polysegm_lengths;
}
template void joinAlignedSegments( const vector<vector<Point2f> > segments, vector<vector<int> > *_polysegm_indexes, vector<double> *_polysegm_lengths );

/**
 * Function that fits segments to poly-segments
 */
template<typename T>
void fitSegmentsToPolySegments( const vector<vector<Point_<T> > > segments, const vector<vector<int> > polysegm_indexes, vector<vector<Point_<T> > > *_fitted ) {

  vector<vector<Point_<T> > > fitted;

  for( size_t i=0; i<polysegm_indexes.size(); i++ )
    if( polysegm_indexes[i].size() == 1 )
      fitted.push_back( segments[polysegm_indexes[i][0]] );
    else {
      vector<Point_<T> > points( 2*polysegm_indexes[i].size() );
      for( size_t j=0,jj=0; j<polysegm_indexes[i].size(); j++,jj+=2 ) {
        int idx = polysegm_indexes[i][j];
        points[jj] = segments[idx][0];
        points[jj+1] = segments[idx][1];
      }
      vector<Point_<T> > fit( 2 );
      fitSegment( points, fit, CV_DIST_L2, 0, 0.01, 0.01 );
      fitted.push_back( fit );
    }

  *_fitted = fitted;
}
template void fitSegmentsToPolySegments( const vector<vector<Point2f> > segments, const vector<vector<int> > polysegm_indexes, vector<vector<Point2f> > *_fitted );

/**
 * Function that pairs-up close segments of opposing directions
 */
template<typename T>
void pairUpOpposingSegments( const vector<vector<vector<Point_<T> > > > perdirection, const vector<vector<double> > polysegm_lengths, const vector<vector<vector<int> > > polysegm_indexes, const vector<vector<vector<Point_<T> > > > polysegm_fit, vector<vector<vector<Point_<T> > > > *_pairup_fit ) {

  vector<vector<vector<Point2f> > > pairup_fit;

  int num_pairchecks = 0;
  //for( int k=SEGDIR_DOWN; k<=SEGDIR_RIGHT; k+=2 ) {
  for( size_t k=0; k<perdirection.size(); k+=2 ) {
    vector<vector<Point2f> > pairupk_fit;

    vector<double> polysegmki_lengths = polysegm_lengths[k];
    vector<double> polysegmkj_lengths = polysegm_lengths[k+1];
    vector<vector<Point2f> > polysegmki_fit = polysegm_fit[k];
    vector<vector<Point2f> > polysegmkj_fit = polysegm_fit[k+1];

    vector<int> srti( polysegmki_lengths.size() );
    vector<int> srtj( polysegmkj_lengths.size() );
    sortIdx( polysegmki_lengths, srti, CV_SORT_DESCENDING );
    sortIdx( polysegmkj_lengths, srtj, CV_SORT_DESCENDING );

    vector<bool> jpaired( polysegmkj_fit.size(), false );

    // @todo implement multi pairings, i.e. each pairing may have multiple poly-segments in both directions

    for( size_t ii=0; ii<polysegmki_fit.size(); ii++ ) {
      int i = srti[ii];
      Point2f pi1 = polysegmki_fit[i][0];
      Point2f pi2 = polysegmki_fit[i][1];

      vector<int> pairup;
      vector<double> pairup_dst;
      vector<vector<Point2f> > pairup_range;

      for( size_t jj=0; jj<polysegmkj_fit.size(); jj++ ) {
        int j = srtj[jj];
        if( jpaired[j] )
          continue;

        num_pairchecks ++;

        Point2f pj1 = polysegmkj_fit[j][0];
        Point2f pj2 = polysegmkj_fit[j][1];

        /// Get points on line i closest to extremes of segment j ///
        Point2f j1_to_i;
        Point2f j2_to_i;
        pointToLine( pi1, pi2, pj1, j1_to_i );
        pointToLine( pi1, pi2, pj2, j2_to_i );
        Point2f j1_on_i = pj1 + j1_to_i;
        Point2f j2_on_i = pj2 + j2_to_i;

        /// Check that segment j is at the left (if vertical) or above (if horizontal) of segment i ///
        Matx33d tri1( 1.0, pj1.x, pj1.y,
                          1.0, pj2.x, pj2.y,
                          1.0, j1_on_i.x, j1_on_i.y );
        Matx33d tri2( 1.0, pj1.x, pj1.y,
                          1.0, pj2.x, pj2.y,
                          1.0, j2_on_i.x, j2_on_i.y );
        double det1 = determinant( tri1 );
        double det2 = determinant( tri2 );
        if( det1 < 0.0 || det2 < 0.0 ) // sign indicates if clockwise or counterclockwise
          continue;

        /// Check that segments overlap ///
        if( ! segmentsOverlap(pi1,pi2,j1_on_i,j2_on_i) )
          continue;
        /*bool j1_on_i_check = pointInSegment(pi1,pi2,j1_on_i);
        bool j2_on_i_check = pointInSegment(pi1,pi2,j2_on_i);
        bool i1_on_j_check = pointInSegment(j1_on_i,j2_on_i,pi1);
        bool i2_on_j_check = pointInSegment(j1_on_i,j2_on_i,pi2);
        if( ! ( j1_on_i_check || j2_on_i_check || i1_on_j_check || i2_on_j_check ) )
          continue;

        /// If a point is outside segment i, get point on j at the extreme of i ///
        if( ! j1_on_i_check ) {
          j1_on_i = (j1_on_i-pi1).ddot(j1_on_i-pi1) < (j1_on_i-pi2).ddot(j1_on_i-pi2) ? pi1 : pi2 ;
          pointToLine( pj1, pj2, j1_on_i, j1_to_i );
          pj1 = j1_on_i + j1_to_i;
          j1_to_i = -j1_to_i;
        }
        if( ! j2_on_i_check ) {
          j2_on_i = (j2_on_i-pi1).ddot(j2_on_i-pi1) < (j2_on_i-pi2).ddot(j2_on_i-pi2) ? pi1 : pi2 ;
          pointToLine( pj1, pj2, j2_on_i, j2_to_i );
          pj2 = j2_on_i + j2_to_i;
          j2_to_i = -j2_to_i;
        }*/

        /// Check that distance between segments is small ///
        double dst_ij_ji = 0.5*( norm(j1_to_i) + norm(j2_to_i) );
        if( dst_ij_ji >= gb_pair_max_width )
          continue;

        /// Check already added pairings and if overlap keep only the closest one ///
        bool pairing_overlap = false;
        for( size_t pu=0; pu<pairup.size(); pu++ )
          if( pairup[pu] != -1 &&
              segmentsOverlap(j1_on_i,j2_on_i,pairup_range[pu][0],pairup_range[pu][1]) ) {
            //logger( 0, "pair-up overlap i=%d j=%d pu=%d dst_j=%g dst_pu=%g %s%s", i, j, pairup[pu], dst_ij_ji, pairup_dst[pu], dirname[k], dirname[k+1] );
            if( dst_ij_ji > pairup_dst[pu] ) {
              pairing_overlap = true;
              break;
            }
            else {
              jpaired[pairup[pu]] = false;
              pairup[pu] = -1;
            }
          }
        if( pairing_overlap )
          continue;

        //logger( 0, "%s,%s i=%d j=%d dst=%g j1_on_i=%s j2_on_i=%s det1=%g det2=%g", dirname[k], dirname[k+1], i, j, dst_ij_ji, strbool(j1_on_i_check), strbool(j2_on_i_check), det1, det2 );
        //logger( 0, "%d %s pairs-up with %d %s", i, dirname[k], j, dirname[k+1] );

        pairup.push_back( j );
        pairup_dst.push_back( dst_ij_ji );
        vector<Point2f> psegm;
        psegm.push_back( j1_on_i );
        psegm.push_back( j2_on_i );
        pairup_range.push_back( psegm );
        jpaired[j] = true;
      } // for( size_t jj=0; jj<polysegmkj_fit.size(); jj++ ) {

      /// Fit a segment to paired poly-segments ///
      if( pairup.size() > 0 ) {
        int n = 0;
        double half = 0.0;
        for( size_t j=0; j<pairup.size(); j++ ) {
          if( pairup[j] == -1 )
            continue;
          half += pairup_dst[j];
          n ++;
        }
        half /= n+n;

        vector<Point2f> points;
        vector<Point2f> midpoly;
        for( size_t si=0; si<polysegm_indexes[k][i].size(); si++ ) {
          vector<Point2f> segm = perdirection[k][polysegm_indexes[k][i][si]];
          //float *p = dirsegs[k]->ptr<float>(polysegm_indexes[k][i][si]);
          Point2f p1 = segm[0];
          Point2f p2 = segm[1];
          //Point2f p1( p[0], p[1] );
          //Point2f p2( p[2], p[3] );
          Point2f vdir = p2 - p1;
          double nrm = norm(vdir);
          double tmp = vdir.x/nrm;
          vdir.x = -vdir.y/nrm;
          vdir.y = tmp;
          midpoly.push_back( p1 + half*vdir );
          midpoly.push_back( p2 + half*vdir );
        }
        midpoly = samplePolySegment( midpoly, 4*gb_pair_max_width );
        points.insert( points.end(), midpoly.begin(), midpoly.end() );

        for( size_t j=0; j<pairup.size(); j++ ) {
          if( pairup[j] == -1 )
            continue;
          vector<Point2f> midpoly;
          for( size_t sj=0; sj<polysegm_indexes[k+1][pairup[j]].size(); sj++ ) {
            vector<Point2f> segm = perdirection[k+1][polysegm_indexes[k+1][pairup[j]][sj]];
            //float *p = dirsegs[k+1]->ptr<float>(polysegm_indexes[k+1][pairup[j]][sj]);
            Point2f p1 = segm[0];
            Point2f p2 = segm[1];
            //Point2f p1( p[0], p[1] );
            //Point2f p2( p[2], p[3] );
            Point2f vdir = p2 - p1;
            double nrm = norm(vdir);
            double tmp = vdir.x/nrm;
            vdir.x = -vdir.y/nrm;
            vdir.y = tmp;
            midpoly.push_back( p1 + half*vdir );
            midpoly.push_back( p2 + half*vdir );
          }
          midpoly = samplePolySegment( midpoly, 4*gb_pair_max_width );
          points.insert( points.end(), midpoly.begin(), midpoly.end() );
        }
        vector<Point2f> segment( 2 );
        fitSegment( points, segment, CV_DIST_L2, 0, 0.01, 0.01 );

        /// Have all pairing fits be direction down or right ///
        if( ( k == SEGDIR_DOWN && segment[0].y > segment[1].y ) ||
            ( k == SEGDIR_LEFT && segment[0].x > segment[1].x ) ) {
          Point2f tmp = segment[0];
          segment[0] = segment[1];
          segment[1] = tmp;
        }

        pairupk_fit.push_back( segment );
      }

      // @todo compute a score for pairings based on straightness, opposite directions overlap, pairing distance consistency, maybe between 0.5 and 1 and give unpaired poly-segments scores between 0 and 0.5

    } // for( size_t ii=0; ii<polysegmki_fit.size(); ii++ ) {

    pairup_fit.push_back( pairupk_fit );
  //} // for( int k=SEGDIR_DOWN; k<=SEGDIR_RIGHT; k+=2 ) {
  } // for( int k=0; k<perdirection.size(); k+=2 ) {

  *_pairup_fit = pairup_fit;
}
template void pairUpOpposingSegments( const vector<vector<vector<Point2f> > > perdirection, const vector<vector<double> > polysegm_lengths, const vector<vector<vector<int> > > polysegm_indexes, const vector<vector<vector<Point2f> > > polysegm_fit, vector<vector<vector<Point2f> > > *_pairup_fit );


/**
 * Compute for every segment a score for being either: table begin line, table in between line or table end line.
 */
void tableLineScores(
  const vector<vector<Point2f> > horzsegms,
  const vector<vector<Point2f> > vertsegms,
  vector<vector<double> > &_scores,
  vector<vector<vector<int> > > &_perpidxs,
  vector<vector<vector<int> > > &_perp1idxs,
  vector<vector<vector<int> > > &_perp2idxs,
  vector<vector<vector<double> > > &_perp1dsts,
  vector<vector<vector<double> > > &_perp2dsts
  ) {

  vector<vector<double> > sumlengths( 4 );
  vector<vector<double> > scores( 4 );
  vector<vector<vector<int> > > perpidxs( 4 );
  vector<vector<vector<int> > > perp1idxs( 4 );
  vector<vector<vector<int> > > perp2idxs( 4 );
  vector<vector<vector<double> > > perp1dsts( 4 );
  vector<vector<vector<double> > > perp2dsts( 4 );

  for( int k=0; k<4; k++ ) {
    int sz = ( k == TABL_LEFT || k == TABL_RIGHT ) ? vertsegms.size() : horzsegms.size() ;
    sumlengths[k] = vector<double>( sz, 0.0 );
    scores[k]     = vector<double>( sz, 0.0 );
    perpidxs[k]   = vector<vector<int> >( sz );
    perp1idxs[k]  = vector<vector<int> >( sz );
    perp2idxs[k]  = vector<vector<int> >( sz );
    perp1dsts[k]  = vector<vector<double> >( sz );
    perp2dsts[k]  = vector<vector<double> >( sz );
  }

  //vector<vector<vector<int> > > isects( 2 );
  //isects[DIR_HORZ]  = vector<vector<int> >( horzsegms.size() );
  //isects[DIR_VERT] = vector<vector<int> >( vertsegms.size() );

  for( size_t i=0; i<vertsegms.size(); i++ ) {
    Point2f vert1 = vertsegms[i][0];
    Point2f vert2 = vertsegms[i][1];

    for( size_t j=0; j<horzsegms.size(); j++ ) {
      Point2f horz1 = horzsegms[j][0];
      Point2f horz2 = horzsegms[j][1];
      double lgth_horz = norm( horz1 - horz2 );
      double lgth_vert = norm( vert1 - vert2 );

      /// Intersection and distances between segment and perpendicular ///
      Point2f isect;
      intersection( horz1, horz2, vert1, vert2, isect );
      double dst_isect_horz1 = norm( isect - horz1 );
      double dst_isect_horz2 = norm( isect - horz2 );
      double dst_isect_vert1 = norm( isect - vert1 );
      double dst_isect_vert2 = norm( isect - vert2 );
      bool isect_in_horz = pointInSegment( horz1, horz2, isect );
      bool isect_in_vert = pointInSegment( vert1, vert2, isect );

      /// Dist. from isect to closest segment endpoint, zero if isect in segment ///
      double dst_isect_horz = isect_in_horz ? 0.0 : min( dst_isect_horz1, dst_isect_horz2 );
      double dst_isect_vert = isect_in_vert ? 0.0 : min( dst_isect_vert1, dst_isect_vert2 );

      /// Skip if intersection far from segment along the segment ///
      /// or if perpendicular far from segment ///
      bool horz_candidate = true;
      if( ( dst_isect_horz > 0.2*gb_join_max_sep ) ||
          ( ! isect_in_vert && min( dst_isect_vert1, dst_isect_vert2 ) > 0.5*gb_join_max_sep ) )
        horz_candidate = false;
      bool vert_candidate = true;
      if( ( dst_isect_vert > 0.2*gb_join_max_sep ) ||
          ( ! isect_in_horz && min( dst_isect_horz1, dst_isect_horz2 ) > 0.5*gb_join_max_sep ) )
        vert_candidate = false;

      /// Table top stats ///
      if( horz_candidate && dst_isect_vert1 < dst_isect_vert2 ) {
        scores[TABL_TOP][j] += lgth_vert*( dst_isect_vert1 + dst_isect_horz );
        sumlengths[TABL_TOP][j] += lgth_vert;
        perpidxs[TABL_TOP][j].push_back( i );
        perp1dsts[TABL_TOP][j].push_back( dst_isect_vert1 + dst_isect_horz1 );
        perp2dsts[TABL_TOP][j].push_back( dst_isect_vert1 + dst_isect_horz2 );
      }
      /// Table bottom stats ///
      if( horz_candidate && dst_isect_vert1 > dst_isect_vert2 ) {
        scores[TABL_BOTTOM][j] += lgth_vert*( dst_isect_vert2 + dst_isect_horz );
        sumlengths[TABL_BOTTOM][j] += lgth_vert;
        perpidxs[TABL_BOTTOM][j].push_back( i );
        perp1dsts[TABL_BOTTOM][j].push_back( dst_isect_vert2 + dst_isect_horz1 );
        perp2dsts[TABL_BOTTOM][j].push_back( dst_isect_vert2 + dst_isect_horz2 );
      }
      /// Table left stats ///
      if( vert_candidate && dst_isect_horz1 < dst_isect_horz2 ) {
        scores[TABL_LEFT][i] += lgth_horz*( dst_isect_horz1 + dst_isect_vert );
        sumlengths[TABL_LEFT][i] += lgth_horz;
        perpidxs[TABL_LEFT][i].push_back( j );
        perp1dsts[TABL_LEFT][i].push_back( dst_isect_horz1 + dst_isect_vert1 );
        perp2dsts[TABL_LEFT][i].push_back( dst_isect_horz1 + dst_isect_vert2 );
      }
      /// Table right stats ///
      if( vert_candidate && dst_isect_horz1 > dst_isect_horz2 ) {
        scores[TABL_RIGHT][i] += lgth_horz*( dst_isect_horz2 + dst_isect_vert );
        sumlengths[TABL_RIGHT][i] += lgth_horz;
        perpidxs[TABL_RIGHT][i].push_back( j );
        perp1dsts[TABL_RIGHT][i].push_back( dst_isect_horz2 + dst_isect_vert1 );
        perp2dsts[TABL_RIGHT][i].push_back( dst_isect_horz2 + dst_isect_vert2 );
      }
      /// Update table middle stats ///
      //if( isect_in_horz && isect_in_vert ) {
      //  isects[DIR_HORZ][j].push_back( i );
      //  isects[DIR_VERT][i].push_back( j );
      //}
    } // for( size_t j=0; j<horzsegms.size(); j++ ) {
  } // for( size_t i=0; i<vertsegms.size(); i++ ) {

  /// Scores defined as perpendicuar-length-weighted intersection distance divided by number of close perpendiculars ///
  for( int k=0; k<4; k++ )
    for( size_t i=0; i<scores[k].size(); i++ ) {
      scores[k][i] = perpidxs[k][i].size() == 0 ?
        numeric_limits<double>::infinity() :
        scores[k][i] / ( sumlengths[k][i] * perpidxs[k][i].size() );
      vector<double> dsts1;
      vector<double> dsts2;
      //printf( "#### %s %zu ####\n", tabname[k], i );
      perp1idxs[k][i] = sortIdxLimit( perp1dsts[k][i], CV_SORT_ASCENDING, 0.3*gb_join_max_sep, &dsts1, 1 );
      perp2idxs[k][i] = sortIdxLimit( perp2dsts[k][i], CV_SORT_ASCENDING, 0.3*gb_join_max_sep, &dsts2, 1 );
      perp1dsts[k][i] = dsts1;
      perp2dsts[k][i] = dsts2;
      //printf( "%s %zu: sco=%g #idx1=%zu->%zu #idx2=%zu->%zu\n", tabname[k], i, scores[k][i], perp1dsts[k][i].size(), dsts1.size(), perp2dsts[k][i].size(), dsts2.size() );
      if( dsts1.size() == 0 || 
          ( dsts2.size() == 0 &&
            ( k == TABL_TOP || k == TABL_BOTTOM ) ) ) {
        scores[k][i] = numeric_limits<double>::infinity();
        continue;
      }
      for( int j=perp1idxs[k][i].size()-1; j>=0; j-- )
        perp1idxs[k][i][j] = perpidxs[k][i][ perp1idxs[k][i][j] ];
      for( int j=perp2idxs[k][i].size()-1; j>=0; j-- )
        perp2idxs[k][i][j] = perpidxs[k][i][ perp2idxs[k][i][j] ];
    }

  _scores = scores;
  _perpidxs = perpidxs;
  _perp1idxs = perp1idxs;
  _perp2idxs = perp2idxs;
  _perp1dsts = perp1dsts;
  _perp2dsts = perp2dsts;
}


/*** Program ******************************************************************/
int main( int argc, char *argv[] ) {
  logfile = stderr;
  int err = 0;
  FILE *ifd = NULL;
  char lsdtmp[FILENAME_MAX];
  char pgmtmp[FILENAME_MAX];

  /// Parse input arguments ///
  if( ! strncmp( "--", argc>1 ? argv[1] : "", 2 ) ) {
    print_svn_rev( logfile );
    return SUCCESS;
  }

  char const *optstr = "OF:T:D:S:P:A:a:J:M:m:u:d:pl:V:hv";
  while( getopt(argc,argv,optstr) != -1 );
  int nopts = optind;
  optind = 1;

  char *p;
  int n;
  while( ( n = getopt(argc,argv,optstr) ) != -1 )
    switch( n ) {
    case 'O':
#ifdef ENABLE_OPENCV_LSD
      gb_opencv_lsd = true;
#else
      die( "error: not compiled with opencv LSD support" );
#endif
      break;
    case 'F':
      if( ! strcmp(optarg,"ascii") )
        gb_format = OUT_ASCII;
      else if( ! strcmp(optarg,"xmlpage") )
        gb_format = OUT_XMLPAGE;
      else {
        print_usage( logfile );
        return FAILURE;
      }
      break;
    case 'T':
      p = optarg;
      for( int m=0; m<2; m++ ) {
        if( strchr(p,'x') == NULL )
          die( "error: expected table definition to be of form {rows}x{cols}[,{rows}x{cols}]" );
        tabdef[m][0] = atoi(p);
        p = strchr(p,'x') + 1;
        tabdef[m][1] = atoi(p);
        if( ( p=strchr(p,',') ) == NULL )
          m = 2;
        p ++;
      }
      break;
    case 'D':
      gb_twopage = atof(optarg);
      break;
    case 'S':
      gb_join_max_sep = atof(optarg);
      break;
    case 'P':
      gb_join_max_parallel = atof(optarg);
      break;
    case 'A':
      gb_join_min_dot = cos(atof(optarg)*M_PI/180.0);
      break;
    case 'a':
      gb_join_min_short_dot = cos(atof(optarg)*M_PI/180.0);
      break;
    case 'J':
      gb_join_min_sep_dot = cos(atof(optarg)*M_PI/180.0);
      break;
    case 'M':
      gb_pair_max_width = atof(optarg);
      break;
    case 'm':
      p = optarg;
      do {
        if( p != optarg )
          p ++;
        int pg;
        char side;
        double fact;
        if( 3 != sscanf( p, "%d%c%lf", &pg, &side, &fact ) )
          die( "error: unexpected margin column/row specification: %s", p );
        if( pg < 1 || pg > 2 )
          die( "error: margin page has to either 1 or 2: %s", p );
        if( fact < 0.0 || fact > 1.0 )
          die( "error: expected margin to be a value in the range [0,1]: %s", p );
        switch( side ) {
          case 't':
            gb_margrow_top[pg-1] = fact;
            break;
          case 'b':
            gb_margrow_bottom[pg-1] = fact;
            break;
          case 'l':
            gb_margcol_left[pg-1] = fact;
            break;
          case 'r':
            gb_margcol_right[pg-1] = fact;
            break;
          default:
            die( "error: unexpected margin column/row specification: %s", p );
        }
        p = strchr(p,',');
      } while( p != NULL );
      break;
    case 'u':
      if( ! strcmp(optarg,"mm") || ! strcmp(optarg,"pixels") )
        gb_units = optarg;
      else
        die( "error: units can only be mm or pixels" );
      break;
    case 'd':
      gb_density = atof(optarg);
      break;
    //case 'p':
    //  gb_procimgs = !gb_procimgs;
    //  break;
    case 'l':
      if( ! strcmp( optarg, "-" ) ||
          ! strcmp( optarg, "stdout" ) ||
          ! strcmp( optarg, "/dev/stdout" ) )
        logfile = stdout;
      else if( (logfile = fopen(optarg,"ab")) == NULL )
        die( "error: unable to %s file %s", "open log", optarg );
      break;
    case 'V':
      if( optarg[0] == '+' )
        verbosity ++;
      else if( optarg[0] == '-' )
        verbosity --;
      else
        verbosity = atoi( optarg );
      verbosity = verbosity < 0 ? 0 : verbosity;
      break;
    default:
      logger( 0, "error: incorrect input argument (-%c)", n );
      err = FAILURE;
    case 'h':
      print_usage( logfile );
      return err;
    case 'v':
      print_svn_rev( logfile );
      return err;
    }

  if( argc - nopts > 1 ) {
    logger( 0, "error: expected at most one non-option argument" );
    print_usage( logfile );
    return FAILURE;
  }

  if( argc - nopts > 0 )
    ifn = argv[optind++];

  if( ifn == NULL ) {
    logger( 0, "error: input image required" );
    print_usage( logfile );
    return FAILURE;
  }

  logger( 2, "input image: %s", ifn );

  /// read image ///
  MagickCoreGenesis( (char*)NULL, MagickFalse );

  Img* img = NULL;
  if( readimg_magick( ifn, &img, logfile ) )
    die( "error: unable to read image" );

  if( gb_density ) {
    set_density_magick( img, gb_density );
    logger( 2, "read image of size %dx%d pixels and provided density %.0f ppc", img->width, img->height, gb_density );
  }
  else if( img->res_x == 0 )
    logger( 2, "read image of size %dx%d pixels and undefined density", img->width, img->height );
  else
    logger( 2, "read image of size %dx%d pixels and density %.0fx%.0f ppc", img->width, img->height, img->res_x, img->res_y );

  if( gb_twopage )
    logger( 0, "double page, separation at %g", gb_twopage );
  if( gb_twopage && gb_twopage < 1.0 )
    gb_twopage *= img->width;

  /// convert physical units to pixels ///
  /*if( ( gb_dilate || gb_summeps || gb_rlsa[0] || gb_rlsa[1] || gb_rlsa[2] || gb_rlsa[3] ) &&
      ! strcmp(gb_units,"mm") ) {
    if( img->res_x == 0 )
      die( "error: image does not specify density which is required for mm units" );
    if( img->res_x != img->res_y )
      die( "error: expected image density to be the same for vertical and horizontal" );
    double fact = img->res_x / 10.0;
    gb_dilate *= fact;
    for( n=0; n<4; n++ )
      gb_rlsa[n] *= fact;
    gb_summeps *= fact;
  }*/

#ifdef ENABLE_OPENCV_LSD
  /// Use opencv's line_descriptor contrib module to get line segments ///
  vector<KeyLine> lines;
  if( gb_opencv_lsd ) {
    logger( 1, "LSD segments from opencv's line_descriptor contrib module" );
    Mat imageMat = Mat( img->height, img->width, CV_8UC1 );
    getpixels_magick_cv8u1_mat( img, imageMat );
    Mat mask = Mat::ones( imageMat.size(), CV_8UC1 );
    Ptr<LSDDetector> bd = LSDDetector::createLSDDetector();
    bd->detect( imageMat, lines, 2, 1, mask );
  }
#else
  if( gb_opencv_lsd ) {
    die( "error: not compiled with opencv LSD support" );
  }
#endif

  /// Use external lsd command to get line segments ///
  else {
    logger( 1, "LSD segments from external command" );

    //tmpnam( lsdtmp );
    mktemp( tool, lsdtmp );
    snprintf( pgmtmp, FILENAME_MAX-4, "%s.pgm", lsdtmp );

    char cmd[6+2*FILENAME_MAX];
    sprintf( cmd, "lsd %s %s", pgmtmp, lsdtmp );

    logger( 4, "executing: %s", cmd );
    writeimg_magick( pgmtmp, img, logfile );
    int retcode = system( cmd );

    if( retcode != 0 )
      die( "problems executing lsd external command" );

    if( (ifd=fopen(lsdtmp,"r")) == NULL )
      die( "problems reading lsd results file: %s", lsdtmp );
  }

  /// Print Page XML header ///
  bool xmlpage = gb_format == OUT_XMLPAGE ? true : false ;

  if ( xmlpage ) {
    char buf[80];
    time_t now = time(0);
    struct tm tstruct;
    tstruct = *localtime( &now );
    strftime( buf, sizeof(buf), "%Y-%m-%dT%X", &tstruct );

    printf( "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n" );
    printf( "<PcGts xmlns=\"http://schema.primaresearch.org/PAGE/gts/pagecontent/2013-07-15\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://schema.primaresearch.org/PAGE/gts/pagecontent/2013-07-15 http://schema.primaresearch.org/PAGE/gts/pagecontent/2013-07-15/pagecontent.xsd\">\n" );
    printf( "<Metadata>\n" );
    printf( "<Creator>%s</Creator>\n", tool );
    printf( "<Created>%s</Created>\n", buf );
    printf( "<LastChange>%s</LastChange>\n", buf );
    printf( "</Metadata>\n" );
    printf( "<Page imageFilename=\"%s\" imageWidth=\"%d\" imageHeight=\"%d\">\n", strrchr(ifn,'/') == NULL ? ifn : strrchr(ifn,'/')+1, img->width, img->height );
  }


  /// Load lsd segments list into vectors ///
  F1 invect[7];
  vector<vector<Point2f> > lsd_segments;
  n = 0;
  while( true ) {
    Point2f p1;
    Point2f p2;

#ifdef ENABLE_OPENCV_LSD
    /// Either from OpenCV ///
    if( gb_opencv_lsd ) {
      if( n == int(lines.size()) )
        break;
      p1.x = lines[n].startPointX;
      p1.y = lines[n].startPointY;
      p2.x = lines[n].endPointX;
      p2.y = lines[n].endPointY;
    }
#endif

    n++;
    /// Or from external lsd command ///
    if( ! gb_opencv_lsd ) {
      if( feof(ifd) )
        break;
      /// Read segment parameters ///
      if( fscanVEC_F1( ifd, invect, 7, logfile ) ) {
        if( ! feof(ifd) )
          die( "problems parsing input line %d of lsd results file %s", n, lsdtmp );
        continue;
      }
      p1.x = invect[0];
      p1.y = invect[1];
      p2.x = invect[2];
      p2.y = invect[3];
    }

    /// Compute slope and discard if not almost vertical/horizontal ///
    double abs_dx = fabs( p1.x - p2.x );
    double abs_dy = fabs( p1.y - p2.y );
    double slope = abs_dx > abs_dy ? abs_dy/abs_dx : abs_dx/abs_dy ;
    if( slope > 0.1 ) {
      logger( 2, "info: n=%d: skipping, slope=%g", n, slope );
      continue;
    }

    /// Append to vector ///
    vector<Point2f> lsd_segment;
    lsd_segment.push_back( p1 );
    lsd_segment.push_back( p2 );
    lsd_segments.push_back( lsd_segment );
  }

  /// Remove external LSD temporal files ///
  if( ifd != NULL ) {
    fclose( ifd );
    unlink( lsdtmp );
    unlink( pgmtmp );
  }

  /// Process each page in case there are multiple ///
  for( int pg=0; pg<(gb_twopage?2:1); pg++ ) {

    /// Filter segments for page ///
    vector<vector<Point2f> > filtered_segments;
    //vector<int> filtered_indexes;
    if( ! gb_twopage )
      filtered_segments = lsd_segments;
    else if( pg == 0 )
      filterSegmentsRange( lsd_segments, 0, 0, int(gb_twopage), img->height, &filtered_segments, NULL /*&filtered_indexes*/ );
    else
      filterSegmentsRange( lsd_segments, int(gb_twopage), 0, int(gb_twopage), img->height, &filtered_segments, NULL /*&filtered_indexes*/ );

    /// Sort segments by direction ///
    vector<vector<vector<Point2f> > > directed_segments;
    //vector<vector<int> > directed_indexes;
    segmentsPerDirection( filtered_segments, &directed_segments, NULL /*&directed_indexes*/ );

    /// Join almost aligned segments for each direction ///
    vector<vector<double> > polysegm_lengths;
    vector<vector<vector<int> > > polysegm_indexes;
    vector<vector<vector<Point2f> > > polysegm_fit;

    for( size_t k=0; k<directed_segments.size(); k++ ) {
      vector<double> polysegmk_lengths;
      vector<vector<int> > polysegmk_indexes;
      joinAlignedSegments( directed_segments[k], &polysegmk_indexes, &polysegmk_lengths );

      /// Fit a segment to each poly-segment ///
      vector<vector<Point2f> > polysegmk_fit;
      fitSegmentsToPolySegments( directed_segments[k], polysegmk_indexes, &polysegmk_fit );

      polysegm_lengths.push_back( polysegmk_lengths );
      polysegm_indexes.push_back( polysegmk_indexes );
      polysegm_fit.push_back( polysegmk_fit );
    }

    /// Pair-up poly-segments of opposing directions ///
    vector<vector<vector<Point2f> > > pairup_fits;
    pairUpOpposingSegments( directed_segments, polysegm_lengths, polysegm_indexes, polysegm_fit, &pairup_fits );

    /// Compute pair-up scores for table begin, end or in between ///
    vector<vector<double> > scores;
    vector<vector<vector<int> > > perpidxs;
    vector<vector<vector<int> > > perp1idxs;
    vector<vector<vector<int> > > perp2idxs;
    vector<vector<vector<double> > > perp1dsts;
    vector<vector<vector<double> > > perp2dsts;
    tableLineScores( pairup_fits[DIR_HORZ], pairup_fits[DIR_VERT], scores, perpidxs, perp1idxs, perp2idxs, perp1dsts, perp2dsts );


    /*for( int k=0; k<4; k++ )
    for( size_t j=0; j<scores[k].size(); j++ )
    if( scores[k][j] < 1e9 ) {
      printf( "%s %zu: score=%g\n", tabname[k], j, scores[k][j] );
      printf( "  perpidxs (%zu):", perpidxs[k][j].size() );
      for( size_t i=0; i<perpidxs[k][j].size(); i++ )
        printf( " %d", perpidxs[k][j][i] );
      printf( "\n" );
      printf( "  perp1 (%zu,%zu):", perp1idxs[k][j].size(), perp1dsts[k][j].size() );
      for( size_t i=0; i<perp1idxs[k][j].size(); i++ )
        printf( " %d->%.1f", perp1idxs[k][j][i], perp1dsts[k][j][i] );
      printf( "\n" );
      printf( "  perp2 (%zu,%zu):", perp2idxs[k][j].size(), perp2dsts[k][j].size() );
      for( size_t i=0; i<perp2idxs[k][j].size(); i++ )
        printf( " %d->%.1f", perp2idxs[k][j][i], perp2dsts[k][j][i] );
      printf( "\n" );
    }*/


    // @todo the following should eventually be moved to a function

    double score_limit = gb_join_max_sep;
    vector<double> scores_top;
    vector<double> scores_bottom;
    vector<double> scores_left;
    vector<double> scores_right;
    vector<double> scores_leftright;

    vector<Point2i> order_top = pairSumSortIdxLimit( scores[TABL_TOP], scores[TABL_TOP], CV_SORT_ASCENDING | SORT_PAIR_SAME_NOSUM, score_limit, &scores_top );
    vector<Point2i> order_bottom = pairSumSortIdxLimit( scores[TABL_BOTTOM], scores[TABL_BOTTOM], CV_SORT_ASCENDING | SORT_PAIR_SAME_NOSUM, score_limit, &scores_bottom );
    vector<Point2i> order_left = pairSumSortIdxLimit( scores[TABL_LEFT], scores[TABL_LEFT], CV_SORT_ASCENDING | SORT_PAIR_SAME_NOSUM, score_limit, &scores_left );
    vector<Point2i> order_right = pairSumSortIdxLimit( scores[TABL_RIGHT], scores[TABL_RIGHT], CV_SORT_ASCENDING | SORT_PAIR_SAME_NOSUM, score_limit, &scores_right );
    vector<Point2i> order_leftright = pairSumSortIdxLimit( scores_left, scores_right, CV_SORT_ASCENDING, score_limit, &scores_leftright );

    /*printf( "order_top.size()=%zu\n", order_top.size() );
    printf( "order_bottom.size()=%zu\n", order_bottom.size() );
    printf( "order_left.size()=%zu\n", order_left.size() );
    printf( "order_right.size()=%zu\n", order_right.size() );
    printf( "order_leftright.size()=%zu\n", order_leftright.size() );*/

    /// Precompute bottom candidates perpendicular scores ///
    vector<vector<double> > perpbottomscores( order_bottom.size() );
    vector<vector<Point2i> > perpbottomorder( order_bottom.size() );
    for( size_t b=0; b<order_bottom.size(); b++ ) {
      int bl = order_bottom[b].x;
      int br = order_bottom[b].y;
      vector<Point2f> segm_bl = pairup_fits[DIR_HORZ][bl];
      vector<Point2f> segm_br = pairup_fits[DIR_HORZ][br];
      if( bl != br && ! segmentsAligned( segm_bl, segm_br ) )
        continue;

      perpbottomorder[b] = pairSumSortIdxLimit( perp1dsts[TABL_BOTTOM][bl], perp2dsts[TABL_BOTTOM][br], CV_SORT_ASCENDING, gb_join_max_sep, &(perpbottomscores[b]) );
    }

    /// Loop through top candidates ///
    vector<vector<int> > border_candidates;
    vector<double> border_scores;

    for( size_t t=0; t<order_top.size(); t++ ) {
      int tl = order_top[t].x;
      int tr = order_top[t].y;
      vector<Point2f> segm_tl = pairup_fits[DIR_HORZ][tl];
      vector<Point2f> segm_tr = pairup_fits[DIR_HORZ][tr];
      if( tl != tr && ! segmentsAligned( segm_tl, segm_tr ) )
        continue;

      vector<double> perptopscores;
      vector<Point2i> perptoporder = pairSumSortIdxLimit( perp1dsts[TABL_TOP][tl], perp2dsts[TABL_TOP][tr], CV_SORT_ASCENDING, gb_join_max_sep, &perptopscores );

      /*printf( "tl=%d tr=%d sco=%g:\n", tl, tr, scores_top[t] );
      printf( "  perptoporder (%zu):\n", perptoporder.size() );
      for( size_t i=0; i<perptoporder.size(); i++ )
        printf( "    %zu (%d,%d) -> %.1f\n", i, perp1idxs[TABL_TOP][tl][ perptoporder[i].x ], perp2idxs[TABL_TOP][tr][ perptoporder[i].y ], perptopscores[i] );*/

      /// Loop through left and right candidates ///
      for( size_t s=0; s<order_leftright.size(); s++ ) {
        int lt = order_left[order_leftright[s].x].x;
        int lb = order_left[order_leftright[s].x].y;
        int rt = order_right[order_leftright[s].y].x;
        int rb = order_right[order_leftright[s].y].y;
        vector<Point2f> segm_lt = pairup_fits[DIR_VERT][lt];
        vector<Point2f> segm_lb = pairup_fits[DIR_VERT][lb];
        vector<Point2f> segm_rt = pairup_fits[DIR_VERT][rt];
        vector<Point2f> segm_rb = pairup_fits[DIR_VERT][rb];
        //printf( "  tl=%d tr=%d lt=%d lb=%d rt=%d rb=%d:\n", tl, tr, lt, lb, rt, rb );
        if( lt != lb && ! segmentsAligned( segm_lt, segm_lb ) )
          continue;
        if( rt != rb && ! segmentsAligned( segm_rt, segm_rb ) )
          continue;

        /// Loop through top corner candidates ///
        for( size_t i=0; i<perptoporder.size(); i++ ) {
          int topleft = perp1idxs[TABL_TOP][tl][ perptoporder[i].x ];
          int topright = perp2idxs[TABL_TOP][tr][ perptoporder[i].y ];
          if( topleft != lt || topright != rt )
            continue;

          double sco = scores_top[t] + scores_leftright[s] + perptopscores[i];
          //printf( "rnk: t(%zu)=%.1f s(%zu)=%.1f i(%zu)=%1f sco=%.1f top=(%d,%d) left=(%d,%d) right=(%d,%d)\n", t, scores_top[t], s, scores_leftright[s], i, perptopscores[i], sco, tl, tr, lt, lb, rt, rb );

          vector<int> border( 8, -1 );
          border[0] = tl;
          border[1] = tr;
          border[2] = lt;
          border[3] = lb;
          border[4] = rt;
          border[5] = rb;

          /// Loop through bottom candidates ///
          for( size_t b=0; b<order_bottom.size(); b++ ) {
            int bl = order_bottom[b].x;
            int br = order_bottom[b].y;
            /// Loop through bottom corner candidates ///
            for( size_t c=0; c<perpbottomorder[b].size(); c++ ) {
              if( perpbottomorder[b].size() == 0 )
                continue;
              int bottomleft = perp1idxs[TABL_BOTTOM][bl][ perpbottomorder[b][c].x ];
              int bottomright = perp2idxs[TABL_BOTTOM][br][ perpbottomorder[b][c].y ];

              if( bottomleft == lb && bottomright == rb ) {
                border[6] = bl;
                border[7] = br;
                sco += scores_bottom[b] + perpbottomscores[b][c];
                b = order_bottom.size();
                break;
              }
            } // for( size_t c=0; c<perpbottomorder[b].size(); i++ ) {
          } // for( size_t b=0; b<order_bottom.size(); b++ ) {

          /// If no bottom, check that sides are more or less the same length ///
          bool skip_candidate = false;
          if( border[6] < 0 ) {
            double lgth_left = norm( segm_lb[1] - segm_lt[0] );
            double lgth_right = norm( segm_rb[1] - segm_rt[0] );
            double side_ratio = lgth_left < lgth_right ? lgth_right/lgth_left : lgth_left/lgth_right ;
            if( side_ratio > 1.2 )
              skip_candidate = true;
          }
          if( skip_candidate )
            continue;

          /// Add table border candidate ///
          border_candidates.push_back( border );
          border_scores.push_back( sco );
        } // for( size_t i=0; i<perptoporder.size(); i++ ) {
      } // for( size_t s=0; s<order_leftright.size(); s++ ) {
    } // for( size_t t=0; t<order_top.size(); t++ ) {

    if( border_candidates.size() == 0 )
      logger( 0, "warning: unable to detect a table" );

    for( size_t i=0; i<border_candidates.size(); i++ )
      logger( 0, "bord_cand(%zu): top=(%d,%d) left=(%d,%d) right=(%d,%d) bottom=(%d,%d) sco=%g", i, border_candidates[i][0], border_candidates[i][1], border_candidates[i][2], border_candidates[i][3], border_candidates[i][4], border_candidates[i][5], border_candidates[i][6], border_candidates[i][7], border_scores[i] );


    /// Output results ///
    if( ! xmlpage ) {
      for( size_t k=0; k<directed_segments.size(); k++ ) {
        /// Print the filtered LSD segments for each direction ///
        for( size_t i=0; i<directed_segments[k].size(); i++ ) {
          Point2f p1 = directed_segments[k][i][0];
          Point2f p2 = directed_segments[k][i][1];
          printf( "%s %g,%g %g,%g\n", dirname[k], p1.x, p1.y, p2.x, p2.y );
        }
        for( size_t i=0; i<polysegm_indexes[k].size(); i++ ) {
          /// Print the join segments ///
          for( size_t j=1; j<polysegm_indexes[k][i].size(); j++ ) {
            int prev = polysegm_indexes[k][i][j-1];
            int curr = polysegm_indexes[k][i][j];
            Point2f prev2 = directed_segments[k][prev][1];
            Point2f curr1 = directed_segments[k][curr][0];
            printf( "%s_join %g,%g %g,%g %zu %zu\n", dirname[k], prev2.x, prev2.y, curr1.x, curr1.y, i, j );
          }
          /// Print the joined poly-segment ///
          printf( "%s_poly", dirname[k] );
          for( size_t j=0; j<polysegm_indexes[k][i].size(); j++ ) {
            int jj = polysegm_indexes[k][i][j];
            Point2f p1 = directed_segments[k][jj][0];
            Point2f p2 = directed_segments[k][jj][1];
            printf( " %g,%g %g,%g", p1.x, p1.y, p2.x, p2.y );
          }
          printf( "\n" );
          /// Print the fitted poly-segment ///
          printf( "%s_polyfit %g,%g %g,%g\n", dirname[k], polysegm_fit[k][i][0].x, polysegm_fit[k][i][0].y, polysegm_fit[k][i][1].x, polysegm_fit[k][i][1].y );
        }
      } // for( size_t k=0; k<directed_segments.size(); k++ ) {

      /// Print the paired-up poly-segments ///
      for( size_t k=0; k<directed_segments.size(); k+=2 ) {
        int kk = k/2;
        for( size_t i=0; i<pairup_fits[kk].size(); i++ ) {
          vector<Point2f> segment = pairup_fits[kk][i];
          printf( "pairup_%s%s %g,%g %g,%g\n", dirname[k], dirname[k+1], segment[0].x, segment[0].y, segment[1].x, segment[1].y );
        }
      }

    } // if( ! xmlpage ) {

    /// Loop through table candidates ///
    int num_table = 0;
    for( int i=0; i<min(1,int(border_candidates.size())); i++ ) {
      Point2f t1 = pairup_fits[DIR_HORZ][border_candidates[i][0]][0];
      Point2f t2 = pairup_fits[DIR_HORZ][border_candidates[i][1]][1];
      Point2f l1 = pairup_fits[DIR_VERT][border_candidates[i][2]][0];
      Point2f l2 = pairup_fits[DIR_VERT][border_candidates[i][3]][1];
      Point2f r1 = pairup_fits[DIR_VERT][border_candidates[i][4]][0];
      Point2f r2 = pairup_fits[DIR_VERT][border_candidates[i][5]][1];
      Point2f b1 = border_candidates[i][6]>=0 ? pairup_fits[DIR_HORZ][border_candidates[i][6]][0] : l2;
      Point2f b2 = border_candidates[i][7]>=0 ? pairup_fits[DIR_HORZ][border_candidates[i][7]][1] : r2;

      if( border_candidates[i][6] < 0 || border_candidates[i][7] < 0 )
        logger( 0, "warning: discarded possible bottom line" );

    // @todo This needs to be improved, if unjoined use not only extreme points but a new fit, when no bottom use a fit (possibly restricted to the expected angle) of the closest vertical pairup bottom endpoints plus an offset computed from the top endpoints
      Point2f topleft;
      Point2f topright;
      Point2f bottomleft;
      Point2f bottomright;
      intersection( t1, t2, l1, l2, topleft );
      intersection( t1, t2, r1, r2, topright );
      intersection( b1, b2, l1, l2, bottomleft );
      intersection( b1, b2, r1, r2, bottomright );

      /// Print table border ///
      if( ! xmlpage )
        printf( "border %g,%g %g,%g %g,%g %g,%g\n",
          topleft.x, topleft.y,
          topright.x, topright.y,
          bottomright.x, bottomright.y,
          bottomleft.x, bottomleft.y );
      else {
        char tabid[20];
        if( gb_twopage )
          sprintf( tabid, "pg%d_tb%d", pg+1, ++num_table );
        else
          sprintf( tabid, "tb%d", ++num_table );
        printf( "<TextRegion id=\"%s\" type=\"floating\">\n", tabid );
        printf( "<Coords points=\"%g,%g %g,%g %g,%g %g,%g\"/>\n",
          topleft.x, topleft.y,
          topright.x, topright.y,
          bottomright.x, bottomright.y,
          bottomleft.x, bottomleft.y );
        printf( "</TextRegion>\n" );
      }

      /// Determine column and row separators ///
      int t1i = border_candidates[i][0];
      int t2i = border_candidates[i][1];
      int l1i = border_candidates[i][2];
      int l2i = border_candidates[i][3];
      int r1i = border_candidates[i][4];
      int r2i = border_candidates[i][5];
      int b1i = border_candidates[i][6];
      int b2i = border_candidates[i][7];

      set<int> colsep;
      set<int> rowsep;

      for( size_t n=0; n<perpidxs[TABL_TOP][t1i].size(); n++ )
        colsep.insert( perpidxs[TABL_TOP][t1i][n] );
      if( t1i != t2i )
        for( size_t n=0; n<perpidxs[TABL_TOP][t2i].size(); n++ )
          colsep.insert( perpidxs[TABL_TOP][t2i][n] );
      if( b1i >= 0 ) {
        for( size_t n=0; n<perpidxs[TABL_BOTTOM][b1i].size(); n++ )
          colsep.insert( perpidxs[TABL_BOTTOM][b1i][n] );
        if( b1i != b2i )
          for( size_t n=0; n<perpidxs[TABL_BOTTOM][b2i].size(); n++ )
            colsep.insert( perpidxs[TABL_BOTTOM][b2i][n] );
      }

      for( size_t n=0; n<perpidxs[TABL_LEFT][l1i].size(); n++ )
        rowsep.insert( perpidxs[TABL_LEFT][l1i][n] );
      if( t1i != t2i )
        for( size_t n=0; n<perpidxs[TABL_LEFT][l2i].size(); n++ )
          rowsep.insert( perpidxs[TABL_LEFT][l2i][n] );
      for( size_t n=0; n<perpidxs[TABL_RIGHT][r1i].size(); n++ )
        rowsep.insert( perpidxs[TABL_RIGHT][r1i][n] );
      if( b1i != b2i )
        for( size_t n=0; n<perpidxs[TABL_RIGHT][r2i].size(); n++ )
          rowsep.insert( perpidxs[TABL_RIGHT][r2i][n] );

      colsep.erase( l1i );
      colsep.erase( l2i );
      colsep.insert( r1i );
      colsep.insert( r2i );
      rowsep.erase( t1i );
      rowsep.erase( t2i );
      rowsep.insert( b1i /*>= 0 ? b1i : -l2i*/ );
      rowsep.insert( b2i /*>= 0 ? b2i : -r2i*/ );

      double lgth_horz = 0.5*( norm(topright-topleft) + norm(bottomright-bottomleft) );
      double lgth_vert = 0.5*( norm(bottomleft-topleft) + norm(bottomright-topright) );

      /// Remove separators whose length does not agree with the table ///
      for( set<int>::iterator it=colsep.begin(); it != colsep.end(); ) {
        Point2f sep1 = pairup_fits[DIR_VERT][*it][0];
        Point2f sep2 = pairup_fits[DIR_VERT][*it][1];
        //logger( 0, "check colsep %d => %g", *it, lgth_vert/norm(sep2-sep1) );
        set<int>::iterator current_it = it++;
        if( lgth_vert/norm(sep2-sep1) > 1.1 )
          colsep.erase( current_it );
      }

      for( set<int>::iterator it=rowsep.begin(); it!=rowsep.end(); ) {
        Point2f sep1 = *it >= 0 ?
          pairup_fits[DIR_HORZ][*it][0] :
          pairup_fits[DIR_VERT][l2i][1] ;
        Point2f sep2 = *it >= 0 ?
          pairup_fits[DIR_HORZ][*it][1] :
          pairup_fits[DIR_VERT][r2i][1] ;
        //logger( 0, "check rowsep %d => %g", *it, lgth_horz/norm(sep2-sep1) );
        set<int>::iterator current_it = it++;
        if( lgth_horz/norm(sep2-sep1) > 1.1 )
          rowsep.erase( current_it );
      }

      /*fprintf( stderr, "colsep(%zu):", colsep.size() );
      for( set<int>::iterator it=colsep.begin(); it!=colsep.end(); it++ )
        fprintf( stderr, " %d", *it );
      fprintf( stderr, "\n" );
      fprintf( stderr, "rowsep(%zu):", rowsep.size() );
      for( set<int>::iterator it=rowsep.begin(); it!=rowsep.end(); it++ )
        fprintf( stderr, " %d", *it );
      fprintf( stderr, "\n" );*/

      /// Sort from left to right and top to bottom the column/row separators ///
      vector<Point2f> linecol1;
      vector<Point2f> linecol2;
      vector<Point2f> linerow1;
      vector<Point2f> linerow2;
      vector<double> linecoldst;
      vector<double> linerowdst;

      linecol1.push_back( topleft );
      linecol2.push_back( bottomleft );
      linecoldst.push_back( 0.0 );
      for( set<int>::iterator it=colsep.begin(); it!=colsep.end(); it++ ) {
        Point2f isect1;
        Point2f isect2;
        Point2f sep1 = pairup_fits[DIR_VERT][*it][0];
        Point2f sep2 = pairup_fits[DIR_VERT][*it][1];
        intersection( topleft, topright, sep1, sep2, isect1 );
        intersection( bottomleft, bottomright, sep1, sep2, isect2 );
        linecol1.push_back( isect1 );
        linecol2.push_back( isect2 );
        linecoldst.push_back( norm(sep1-topleft) + norm(sep2-bottomleft) );
      }

      linerow1.push_back( topleft );
      linerow2.push_back( topright );
      linerowdst.push_back( 0.0 );
      for( set<int>::iterator it=rowsep.begin(); it!=rowsep.end(); it++ ) {
        Point2f isect1;
        Point2f isect2;
        Point2f sep1 = *it >= 0 ?
          pairup_fits[DIR_HORZ][*it][0] :
          pairup_fits[DIR_VERT][l2i][1] ;
        Point2f sep2 = *it >= 0 ?
          pairup_fits[DIR_HORZ][*it][1] :
          pairup_fits[DIR_VERT][r2i][1] ;
        intersection( topleft, bottomleft, sep1, sep2, isect1 );
        intersection( topright, bottomright, sep1, sep2, isect2 );
        linerow1.push_back( isect1 );
        linerow2.push_back( isect2 );
        linerowdst.push_back( norm(sep1-topleft) + norm(sep2-topright) );
      }

      vector<int> linecolorder = sortIdxLimit( linecoldst, CV_SORT_ASCENDING );
      vector<int> lineroworder = sortIdxLimit( linerowdst, CV_SORT_ASCENDING );

      /// Create matrix of the table cells ///
      Mat cellpts( lineroworder.size(), linecolorder.size(), CV_32FC2 );

      /// Fill first and last rows and columns
      for( int col=0; col<cellpts.cols; col++ ) {
        cellpts.at<Point2f>( 0, col ) = linecol1[linecolorder[col]];
        cellpts.at<Point2f>( cellpts.rows-1, col ) = linecol2[linecolorder[col]];
      }
      for( int row=0; row<cellpts.rows; row++ ) {
        cellpts.at<Point2f>( row, 0 ) = linerow1[lineroworder[row]];
        cellpts.at<Point2f>( row, cellpts.cols-1 ) = linerow2[lineroworder[row]];
      }

      for( int col=1; col<cellpts.cols; col++ ) {
        Point2f col1 = cellpts.at<Point2f>( 0, col );
        Point2f col2 = cellpts.at<Point2f>( lineroworder.size()-1, col );
        for( int row=1; row<cellpts.rows; row++ ) {
          Point2f row1 = cellpts.at<Point2f>( row, 0 );
          Point2f row2 = cellpts.at<Point2f>( row, linecolorder.size()-1 );
          Point2f isect;
          intersection( col1, col2, row1, row2, isect );
          cellpts.at<Point2f>( row, col ) = isect;
          //fprintf( stderr, "=== cel=%zu,%zu lcol=%.0f,%.0f-%.0f,%.0f lrow=%.0f,%.0f-%.0f,%.0f isect=%.0f,%.0f\n", row, col, col1.x, col1.y, col2.x, col2.y, row1.x, row1.y, row2.x, row2.y, isect.x, isect.y );
        }
      }


      /// Add margin columns/rows ///
      if( gb_margcol_left[pg] || gb_margcol_right[pg] ) {
        int addcols = gb_margcol_left[pg] && gb_margcol_right[pg] ? 2 : 1 ;
        Mat new_cellpts( cellpts.rows, cellpts.cols+addcols, CV_32FC2 );
        copyMakeBorder( cellpts, new_cellpts, 0, 0, gb_margcol_left[pg]?1:0, gb_margcol_right[pg]?1:0, BORDER_CONSTANT, Scalar(-1.0,-1.0) );

        for( int row=0; row<cellpts.rows; row++ ) {
          Point2f margin;
          Point2f left = cellpts.at<Point2f>( row, 0 );
          Point2f right = cellpts.at<Point2f>( row, cellpts.cols-1 );
          if( gb_margcol_left[pg] ) {
            extendSegment( right, left, gb_margcol_left[pg], margin );
            new_cellpts.at<Point2f>( row, 0 ) = margin;
          }
          if( gb_margcol_right[pg] ) {
            extendSegment( left, right, gb_margcol_right[pg], margin );
            new_cellpts.at<Point2f>( row, new_cellpts.cols-1 ) = margin;
          }
        }

        cellpts = new_cellpts;
      }

      if( gb_margrow_top[pg] || gb_margrow_bottom[pg] ) {
        int addrows = gb_margrow_top[pg] && gb_margrow_bottom[pg] ? 2 : 1 ;
        Mat new_cellpts( cellpts.rows+addrows, cellpts.cols, CV_32FC2 );
        copyMakeBorder( cellpts, new_cellpts, gb_margrow_top[pg]?1:0, gb_margrow_bottom[pg]?1:0, 0, 0, BORDER_CONSTANT, Scalar(-1.0,-1.0) );

        for( int col=0; col<cellpts.cols; col++ ) {
          Point2f margin;
          Point2f top = cellpts.at<Point2f>( 0, col );
          Point2f bottom = cellpts.at<Point2f>( cellpts.rows-1, col );
          if( gb_margrow_top[pg] ) {
            extendSegment( bottom, top, gb_margrow_top[pg], margin );
            new_cellpts.at<Point2f>( 0, col ) = margin;
          }
          if( gb_margrow_bottom[pg] ) {
            extendSegment( top, bottom, gb_margrow_bottom[pg], margin );
            new_cellpts.at<Point2f>( new_cellpts.rows-1, col ) = margin;
          }
        }

        cellpts = new_cellpts;
      }


      /// Print table cells as regions ///
      if( xmlpage )
      for( int col=1; col<cellpts.cols; col++ )
        for( int row=1; row<cellpts.rows; row++ ) {
          char cellid[32];
          char tmpbuf[32];
          int colid = gb_margcol_left[pg] ? col-1 : col ;
          int rowid = gb_margrow_top[pg] ? row-1 : row ;
          sprintf( cellid, "tb%d_%d_%d", num_table, rowid, colid );
          if( gb_twopage ) {
            strncpy( tmpbuf, cellid, 32 );
            sprintf( cellid, "pg%d_%s", pg+1, tmpbuf );
          }
          if( ( row == 1 && gb_margrow_top[pg] ) ||
              ( row == cellpts.rows-1 && gb_margrow_bottom[pg] ) ||
              ( col == 1 && gb_margcol_left[pg] ) ||
              ( col == cellpts.cols-1 && gb_margcol_right[pg] ) ) {
            strncpy( tmpbuf, cellid, 32 );
            sprintf( cellid, "%s_marg", tmpbuf );
          }
          Point2f topleft = cellpts.at<Point2f>( row-1, col-1 );
          Point2f topright = cellpts.at<Point2f>( row-1, col );
          Point2f bottomright = cellpts.at<Point2f>( row, col );
          Point2f bottomleft = cellpts.at<Point2f>( row, col-1 );
          printf( "<TextRegion id=\"%s\" type=\"floating\">\n", cellid );
          printf( "<Coords points=\"%g,%g %g,%g %g,%g %g,%g\"/>\n",
            topleft.x, topleft.y,
            topright.x, topright.y,
            bottomright.x, bottomright.y,
            bottomleft.x, bottomleft.y );
          printf( "</TextRegion>\n" );
        }
    }
  } // for( int pg=0; pg<(gb_twopage?2:1); pg++ ) {


/*
for f in $(find ~/DataBases/HTR/qidenus/1stDATASET/3x10_GT_2016-01-20 -name '*.jpg'); do
  ff=$(echo "$f" | sed 's|.*_\([^_]*_[0-9][0-9][0-9][0-9]\)\.jpg$|\1|');
  convert "$f" -resize 20% testdata/$ff.png;
  echo $ff;
done

id=0003;
id=0048;
id=0043;
id=0018;
id=mean;

for id in $(find testdata/ -name '*_[0-9][0-9][0-9][0-9].png' | sed 's|.\+/||; s|\.png||;'); do
  echo $id;

bname=testdata/$id
out=testdata/out

./imgtabdet -V1 -D 0.5 -J 20 testdata/$id.png 2> testdata/err > $out;
#./imgtabdet -F xmlpage -V1 -D 0.5 -J 20 testdata/$id.png 2> testdata/err > $out;
#htrsh_pagexml_round < $out > $out.xml

func_draw () {
draw=( -fill none
  -stroke red
  $( awk '{ if($1=="left"||$1=="down") printf(" -draw line_%s_%s",$2,$3); }' $out )
  #$( awk -v OFS='_' '{ if($1=="left_polyfit"||$1=="down_polyfit") { $1=" -draw polyline"; printf("%s",$0); } }' $out )
  -stroke blue
  $( awk '{ if($1=="right"||$1=="up") printf(" -draw line_%s_%s",$2,$3); }' $out )
  #$( awk -v OFS='_' '{ if($1=="right_polyfit"||$1=="up_polyfit") { $1=" -draw polyline"; printf("%s",$0); } }' $out )
  -stroke green
  $( awk '{ if($1=="right_join"||$1=="up_join"||$1=="left_join"||$1=="down_join") printf(" -draw line_%s_%s",$2,$3); }' $out )
  -stroke orange
  $( awk -v OFS='_' '{ if($1=="pairup_downup"||$1=="pairup_leftright") { $1=" -draw polyline"; printf("%s",$0); } }' $out )
  #$( awk -v OFS='_' '{ if($1=="tab_top"||$1=="tab_left"||$1=="tab_right"||$1=="tab_bottom") { $1=" -draw polyline"; printf("%s",$0); } }' $out )
  -stroke cyan
  $( awk -v OFS='_' '{ if($1=="border") { $1=" -draw polygon"; printf("%s",$0); } }' $out )
  );
draw=("${draw[@]//_/ }");
convert $bname.png -colorspace rgb "${draw[@]}" ${bname}_tabdet.png;
}
func_draw;

mv testdata/err testdata/err_$id;
mv testdata/out testdata/out_$id;

done
*/

  /// Close Page XML ///
  if ( xmlpage )
    printf( "</Page>\n</PcGts>\n" );

  /// Release resources ///
  if( logfile != stdout && logfile != stderr )
    fclose( logfile );

  free_Img( img );

  MagickCoreTerminus();

  return SUCCESS;
}
