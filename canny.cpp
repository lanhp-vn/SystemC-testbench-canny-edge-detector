/* source: http://marathon.csee.usf.edu/edge/edge_detection.html */
/* URL: ftp://figment.csee.usf.edu/pub/Edge_Comparison/source_code/canny.src */

/* ECPS 203 Assignment 4 solution */

/* off-by-one bugs fixed by Rainer Doemer, 10/05/23 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "systemc.h"


#define VERBOSE 0

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0
#define BOOSTBLURFACTOR 90.0
#ifndef M_PI
#define M_PI 3.14159265356789202346
#endif
#define SIGMA 0.6
#define TLOW  0.3
#define THIGH 0.8

#define COLS 2704
#define ROWS 1520
#define SIZE COLS*ROWS
#define VIDEONAME "Engineering"
#define IMG_IN    "video/" VIDEONAME "%03d.pgm"
#define IMG_OUT   VIDEONAME "%03d_edges.pgm"
#define IMG_NUM   30 /* number of images processed (1 or more) */
#define AVAIL_IMG 30 /* number of different image frames (1 or more) */

#define WINSIZE 21
#define SET_STACK_SIZE set_stack_size(128*1024*1024)

typedef struct Image_s
{
	unsigned char img[SIZE];

	Image_s(void)
	{
	   for (int i=0; i<SIZE; i++)
	   { 
	      img[i] = 0;
	   }
	}

	Image_s& operator=(const Image_s& copy)
	{
	   for (int i=0; i<SIZE; i++)
	   { 
	      img[i] = copy.img[i];
	   }            
	   return *this;
	}

	operator unsigned char*()
	{
	   return img;
	}

	unsigned char& operator[](const int index)
	{
	   return img[index];
	}
} IMAGE;


SC_MODULE(Stimulus)
{
   sc_fifo_out<IMAGE> ImgOut;
   char infilename[70];
   IMAGE Image;
   int i=0, n=0;

   void main(void)
   {
      for(i=0; i<IMG_NUM; i++)
      {
         n = i % AVAIL_IMG;
         sprintf(infilename, IMG_IN, n+1);
         /****************************************************************************
         * Read in the image. This read function allocates memory for the image.
         ****************************************************************************/
         if(VERBOSE) printf("Reading the image %s.\n", infilename);
         if(read_pgm_image(infilename, Image, ROWS, COLS) == 0){
            fprintf(stderr, "Error reading the input image, %s.\n", infilename);
            exit(1);
         }
         ImgOut.write(Image);
      }
   }

   SC_CTOR(Stimulus)
   {
      SC_THREAD(main);
      SET_STACK_SIZE;
   }

   int read_pgm_image(const char *infilename, unsigned char *image, int rows, int cols)
   {
      FILE *fp;
      char buf[71];
      int r, c;

      /***************************************************************************
      * Open the input image file for reading if a filename was given. If no
      * filename was provided, set fp to read from standard input.
      ***************************************************************************/
      if(infilename == NULL) fp = stdin;
      else{
         if((fp = fopen(infilename, "r")) == NULL){
            fprintf(stderr, "Error reading the file %s in read_pgm_image().\n",
               infilename);
            return(0);
         }
      }

      /***************************************************************************
      * Verify that the image is in PGM format, read in the number of columns
      * and rows in the image and scan past all of the header information.
      ***************************************************************************/
      fgets(buf, 70, fp);
      if(strncmp(buf,"P5",2) != 0){
         fprintf(stderr, "The file %s is not in PGM format in ", infilename);
         fprintf(stderr, "read_pgm_image().\n");
         if(fp != stdin) fclose(fp);
         return(0);
      }
      do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */
      sscanf(buf, "%d %d", &c, &r);
      if(c != cols || r != rows){
         fprintf(stderr, "The file %s is not a %d by %d image in ", infilename, cols, rows);
         fprintf(stderr, "read_pgm_image().\n");
         if(fp != stdin) fclose(fp);
         return(0);
      }
      do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */

      /***************************************************************************
      * Read the image from the file.
      ***************************************************************************/
      if((unsigned)rows != fread(image, cols, rows, fp)){
         fprintf(stderr, "Error reading the  data in read_pgm_().\n");
         if(fp != stdin) fclose(fp);
         return(0);
      }

      if(fp != stdin) fclose(fp);
      return(1);
   }
};

SC_MODULE(Monitor)
{
   sc_fifo_in<IMAGE> ImgIn;
   IMAGE Edge;
   char outfilename[128];    /* Name of the output "edge" image */
   int i=0, n=0;

   void main(void)
   {
      for(i=0; i<IMG_NUM; i++)
      {
         ImgIn.read(Edge);
         /****************************************************************************
         * Write out the edge image to a file.
         ****************************************************************************/
         n = i % AVAIL_IMG;
         sprintf(outfilename, IMG_OUT, n+1);
         if(VERBOSE) printf("Writing the edge iname in the file %s.\n", outfilename);
         if(write_pgm_image(outfilename, Edge, ROWS, COLS, "", 255) == 0){
            fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
            exit(1);
         }
      }
   }

   SC_CTOR(Monitor)
   {
      SC_THREAD(main);
      SET_STACK_SIZE;
   }

   int write_pgm_image(const char *outfilename, unsigned char *image, int rows,
      int cols, const char *comment, int maxval)
   {
      FILE *fp;

      /***************************************************************************
      * Open the output image file for writing if a filename was given. If no
      * filename was provided, set fp to write to standard output.
      ***************************************************************************/
      if(outfilename == NULL) fp = stdout;
      else{
         if((fp = fopen(outfilename, "w")) == NULL){
            fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
               outfilename);
            return(0);
         }
      }

      /***************************************************************************
      * Write the header information to the PGM file.
      ***************************************************************************/
      fprintf(fp, "P5\n%d %d\n", cols, rows);
      if(comment != NULL)
         if(strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
      fprintf(fp, "%d\n", maxval);

      /***************************************************************************
      * Write the image data to the file.
      ***************************************************************************/
      if((unsigned)rows != fwrite(image, cols, rows, fp)){
         fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
         if(fp != stdout) fclose(fp);
         return(0);
      }

      if(fp != stdout) fclose(fp);
      return(1);
   }
   
};

SC_MODULE(DataIn)
{
   sc_fifo_in<IMAGE> ImgIn;
   sc_fifo_out<IMAGE> ImgOut;
   IMAGE Image;
   int i = 0;

   void main(void)
   {
      for(i=0; i<IMG_NUM; i++)
      {
         ImgIn.read(Image);
         ImgOut.write(Image);
      }
   }

   SC_CTOR(DataIn)
   {
      SC_THREAD(main);
      SET_STACK_SIZE;
   }
};

SC_MODULE(DUT)
{
   IMAGE Image, Edge;
   sc_fifo_in<IMAGE> ImgIn;
   sc_fifo_out<IMAGE> ImgOut;
   int i=0;

   void main(void)
   {
      for(i=0; i<IMG_NUM; i++)
      {
         ImgIn.read(Image);
         /****************************************************************************
         * Perform the edge detection. All of the work takes place here.
         ****************************************************************************/
         if(VERBOSE) printf("Starting Canny edge detection.\n");
         canny(Image, ROWS, COLS, SIGMA, TLOW, THIGH, Edge);
         ImgOut.write(Image);
      }
   }

   SC_CTOR(DUT)
   {
      SC_THREAD(main);
      SET_STACK_SIZE;
   }

   void canny(unsigned char *image, int rows, int cols, float sigma, float tlow, float thigh, unsigned char *edge)
   {
      unsigned char nms[SIZE]    /* Points that are local maximal magnitude. */
            = {0};
      short int smoothedim[SIZE] /* The image after gaussian smoothing.      */
            = {0},
               delta_x[SIZE]    /* The first derivative image, x-direction. */
            = {0},
               delta_y[SIZE]    /* The first derivative image, y-direction. */
            = {0},
               magnitude[SIZE]  /* The magnitude of the gadient image.      */
            = {0};

      /****************************************************************************
      * Perform gaussian smoothing on the image using the input standard
      * deviation.
      ****************************************************************************/
      if(VERBOSE) printf("Smoothing the image using a gaussian kernel.\n");
      gaussian_smooth(image, rows, cols, sigma, smoothedim);

      /****************************************************************************
      * Compute the first derivative in the x and y directions.
      ****************************************************************************/
      if(VERBOSE) printf("Computing the X and Y first derivatives.\n");
      derivative_x_y(smoothedim, rows, cols, delta_x, delta_y);

      /****************************************************************************
      * Compute the magnitude of the gradient.
      ****************************************************************************/
      if(VERBOSE) printf("Computing the magnitude of the gradient.\n");
      magnitude_x_y(delta_x, delta_y, rows, cols, magnitude);

      /****************************************************************************
      * Perform non-maximal suppression.
      ****************************************************************************/
      if(VERBOSE) printf("Doing the non-maximal suppression.\n");
      non_max_supp(magnitude, delta_x, delta_y, rows, cols, nms);

      /****************************************************************************
      * Use hysteresis to mark the edge pixels.
      ****************************************************************************/
      if(VERBOSE) printf("Doing hysteresis thresholding.\n");
      apply_hysteresis(magnitude, nms, rows, cols, tlow, thigh, edge);
   }

   void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols, short int *magnitude)
   {
      int r, c, pos, sq1, sq2;

      for(r=0,pos=0;r<rows;r++){
         for(c=0;c<cols;c++,pos++){
            sq1 = (int)delta_x[pos] * (int)delta_x[pos];
            sq2 = (int)delta_y[pos] * (int)delta_y[pos];
            magnitude[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
         }
      }

   }

   void derivative_x_y(short int *smoothedim, int rows, int cols, short int *delta_x, short int *delta_y)
   {
      int r, c, pos;

      /****************************************************************************
      * Compute the x-derivative. Adjust the derivative at the borders to avoid
      * losing pixels.
      ****************************************************************************/
      if(VERBOSE) printf("   Computing the X-direction derivative.\n");
      for(r=0;r<rows;r++){
         pos = r * cols;
         delta_x[pos] = smoothedim[pos+1] - smoothedim[pos];
         pos++;
         for(c=1;c<(cols-1);c++,pos++){
            delta_x[pos] = smoothedim[pos+1] - smoothedim[pos-1];
         }
         delta_x[pos] = smoothedim[pos] - smoothedim[pos-1];
      }

      /****************************************************************************
      * Compute the y-derivative. Adjust the derivative at the borders to avoid
      * losing pixels.
      ****************************************************************************/
      if(VERBOSE) printf("   Computing the Y-direction derivative.\n");
      for(c=0;c<cols;c++){
         pos = c;
         delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos];
         pos += cols;
         for(r=1;r<(rows-1);r++,pos+=cols){
            delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos-cols];
         }
         delta_y[pos] = smoothedim[pos] - smoothedim[pos-cols];
      }
   }

   void gaussian_smooth(unsigned char *image, int rows, int cols, float sigma, short int *smoothedim)
   {
      int r, c, rr, cc,     /* Counter variables. */
         windowsize,        /* Dimension of the gaussian kernel. */
         center;            /* Half of the windowsize. */
      float tempim[SIZE]    /* Buffer for separable filter gaussian smoothing. */
         = {0.0},
            kernel[WINSIZE] /* A one dimensional gaussian kernel. */
         = {0.0},
            dot,            /* Dot product summing variable. */
            sum;            /* Sum of the kernel weights variable. */

      /****************************************************************************
      * Create a 1-dimensional gaussian smoothing kernel.
      ****************************************************************************/
      if(VERBOSE) printf("   Computing the gaussian smoothing kernel.\n");
      make_gaussian_kernel(sigma, kernel, &windowsize);
      center = windowsize / 2;

      /****************************************************************************
      * Blur in the x - direction.
      ****************************************************************************/
      if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
      for(r=0;r<rows;r++){
         for(c=0;c<cols;c++){
            dot = 0.0;
            sum = 0.0;
            for(cc=(-center);cc<=center;cc++){
               if(((c+cc) >= 0) && ((c+cc) < cols)){
                  dot += (float)image[r*cols+(c+cc)] * kernel[center+cc];
                  sum += kernel[center+cc];
               }
            }
            tempim[r*cols+c] = dot/sum;
         }
      }

      /****************************************************************************
      * Blur in the y - direction.
      ****************************************************************************/
      if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
      for(c=0;c<cols;c++){
         for(r=0;r<rows;r++){
            sum = 0.0;
            dot = 0.0;
            for(rr=(-center);rr<=center;rr++){
               if(((r+rr) >= 0) && ((r+rr) < rows)){
                  dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
                  sum += kernel[center+rr];
               }
            }
            smoothedim[r*cols+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
         }
      }
   }

   void make_gaussian_kernel(float sigma, float *kernel, int *windowsize)
   {
      int i, center;
      float x, fx, sum=0.0;

      *windowsize = 1 + 2 * ceil(2.5 * sigma);
      center = (*windowsize) / 2;

      if(VERBOSE) printf("      The kernel has %d elements.\n", *windowsize);

      for(i=0;i<(*windowsize);i++){
         x = (float)(i - center);
         fx = pow(2.71828, -0.5*x*x/(sigma*sigma)) / (sigma * sqrt(6.2831853));
         kernel[i] = fx;
         sum += fx;
      }

      for(i=0;i<(*windowsize);i++) kernel[i] /= sum;

      if(VERBOSE){
         printf("The filter coefficients are:\n");
         for(i=0;i<(*windowsize);i++)
            printf("kernel[%d] = %f\n", i, kernel[i]);
      }
   }

   void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval, int cols)
   {
      short *tempmagptr;
      unsigned char *tempmapptr;
      int i;
      int x[8] = {1,1,0,-1,-1,-1,0,1},
         y[8] = {0,1,1,1,0,-1,-1,-1};

      for(i=0;i<8;i++){
         tempmapptr = edgemapptr - y[i]*cols + x[i];
         tempmagptr = edgemagptr - y[i]*cols + x[i];

         if((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval)){
            *tempmapptr = (unsigned char) EDGE;
            follow_edges(tempmapptr,tempmagptr, lowval, cols);
         }
      }
   }

   void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols, float tlow, float thigh, unsigned char *edge)
   {
      int r, c, pos, numedges, highcount, lowthreshold, highthreshold, hist[32768];
      short int maximum_mag=0;

      /****************************************************************************
      * Initialize the edge map to possible edges everywhere the non-maximal
      * suppression suggested there could be an edge except for the border. At
      * the border we say there can not be an edge because it makes the
      * follow_edges algorithm more efficient to not worry about tracking an
      * edge off the side of the image.
      ****************************************************************************/
      for(r=0,pos=0;r<rows;r++){
         for(c=0;c<cols;c++,pos++){
      if(nms[pos] == POSSIBLE_EDGE) edge[pos] = POSSIBLE_EDGE;
      else edge[pos] = NOEDGE;
         }
      }

      for(r=0,pos=0;r<rows;r++,pos+=cols){
         edge[pos] = NOEDGE;
         edge[pos+cols-1] = NOEDGE;
      }
      pos = (rows-1) * cols;
      for(c=0;c<cols;c++,pos++){
         edge[c] = NOEDGE;
         edge[pos] = NOEDGE;
      }

      /****************************************************************************
      * Compute the histogram of the magnitude image. Then use the histogram to
      * compute hysteresis thresholds.
      ****************************************************************************/
      for(r=0;r<32768;r++) hist[r] = 0;
      for(r=0,pos=0;r<rows;r++){
         for(c=0;c<cols;c++,pos++){
      if(edge[pos] == POSSIBLE_EDGE) hist[mag[pos]]++;
         }
      }

      /****************************************************************************
      * Compute the number of pixels that passed the nonmaximal suppression.
      ****************************************************************************/
      for(r=1,numedges=0;r<32768;r++){
         if(hist[r] != 0) maximum_mag = r;
         numedges += hist[r];
      }

      highcount = (int)(numedges * thigh + 0.5);

      /****************************************************************************
      * Compute the high threshold value as the (100 * thigh) percentage point
      * in the magnitude of the gradient histogram of all the pixels that passes
      * non-maximal suppression. Then calculate the low threshold as a fraction
      * of the computed high threshold value. John Canny said in his paper
      * "A Computational Approach to Edge Detection" that "The ratio of the
      * high to low threshold in the implementation is in the range two or three
      * to one." That means that in terms of this implementation, we should
      * choose tlow ~= 0.5 or 0.33333.
      ****************************************************************************/
      r = 1;
      numedges = hist[1];
      while((r<(maximum_mag-1)) && (numedges < highcount)){
         r++;
         numedges += hist[r];
      }
      highthreshold = r;
      lowthreshold = (int)(highthreshold * tlow + 0.5);

      if(VERBOSE){
         printf("The input low and high fractions of %f and %f computed to\n",
      tlow, thigh);
         printf("magnitude of the gradient threshold values of: %d %d\n",
      lowthreshold, highthreshold);
      }

      /****************************************************************************
      * This loop looks for pixels above the highthreshold to locate edges and
      * then calls follow_edges to continue the edge.
      ****************************************************************************/
      for(r=0,pos=0;r<rows;r++){
         for(c=0;c<cols;c++,pos++){
      if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)){
               edge[pos] = EDGE;
               follow_edges((edge+pos), (mag+pos), lowthreshold, cols);
      }
         }
      }

      /****************************************************************************
      * Set all the remaining possible edges to non-edges.
      ****************************************************************************/
      for(r=0,pos=0;r<rows;r++){
         for(c=0;c<cols;c++,pos++) if(edge[pos] != EDGE) edge[pos] = NOEDGE;
      }
   }

   void non_max_supp(short *mag, short *gradx, short *grady, int nrows, int ncols, unsigned char *result)
   {
      int rowcount, colcount,count;
      short *magrowptr,*magptr;
      short *gxrowptr,*gxptr;
      short *gyrowptr,*gyptr,z1,z2;
      short m00; short gx=0; short gy=0;
      float mag1,mag2; float xperp=0;float yperp=0;
      unsigned char *resultrowptr, *resultptr;

      /****************************************************************************
      * Zero the edges of the result image.
      ****************************************************************************/
      for(count=0,resultrowptr=result,resultptr=result+ncols*(nrows-1);
         count<ncols; resultptr++,resultrowptr++,count++){
         *resultrowptr = *resultptr = (unsigned char) 0;
      }

      for(count=0,resultptr=result,resultrowptr=result+ncols-1;
         count<nrows; count++,resultptr+=ncols,resultrowptr+=ncols){
         *resultptr = *resultrowptr = (unsigned char) 0;
      }

      /****************************************************************************
      * Suppress non-maximum points.
      ****************************************************************************/
      for(rowcount=1,magrowptr=mag+ncols+1,gxrowptr=gradx+ncols+1,
         gyrowptr=grady+ncols+1,resultrowptr=result+ncols+1;
         rowcount<nrows-1; /* bug fix 10/05/23, RD */
         rowcount++,magrowptr+=ncols,gyrowptr+=ncols,gxrowptr+=ncols,
         resultrowptr+=ncols){
         for(colcount=1,magptr=magrowptr,gxptr=gxrowptr,gyptr=gyrowptr,
            resultptr=resultrowptr;colcount<ncols-1; /* bug fix 10/05/23, RD */
            colcount++,magptr++,gxptr++,gyptr++,resultptr++){
            m00 = *magptr;
            if(m00 == 0){
               *resultptr = (unsigned char) NOEDGE;
            }
            else{
               xperp = -(gx = *gxptr)/((float)m00);
               yperp = (gy = *gyptr)/((float)m00);
            }

            if(gx >= 0){
               if(gy >= 0){
                     if (gx >= gy)
                     {
                           /* 111 */
                           /* Left point */
                           z1 = *(magptr - 1);
                           z2 = *(magptr - ncols - 1);

                           mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;

                           /* Right point */
                           z1 = *(magptr + 1);
                           z2 = *(magptr + ncols + 1);

                           mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                     }
                     else
                     {
                           /* 110 */
                           /* Left point */
                           z1 = *(magptr - ncols);
                           z2 = *(magptr - ncols - 1);

                           mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

                           /* Right point */
                           z1 = *(magptr + ncols);
                           z2 = *(magptr + ncols + 1);

                           mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp;
                     }
                  }
                  else
                  {
                     if (gx >= -gy)
                     {
                           /* 101 */
                           /* Left point */
                           z1 = *(magptr - 1);
                           z2 = *(magptr + ncols - 1);

                           mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;

                           /* Right point */
                           z1 = *(magptr + 1);
                           z2 = *(magptr - ncols + 1);

                           mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
                     }
                     else
                     {
                           /* 100 */
                           /* Left point */
                           z1 = *(magptr + ncols);
                           z2 = *(magptr + ncols - 1);

                           mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

                           /* Right point */
                           z1 = *(magptr - ncols);
                           z2 = *(magptr - ncols + 1);

                           mag2 = (z1 - z2)*xperp  + (m00 - z1)*yperp;
                     }
                  }
               }
               else
               {
                  if ((gy = *gyptr) >= 0)
                  {
                     if (-gx >= gy)
                     {
                           /* 011 */
                           /* Left point */
                           z1 = *(magptr + 1);
                           z2 = *(magptr - ncols + 1);

                           mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

                           /* Right point */
                           z1 = *(magptr - 1);
                           z2 = *(magptr + ncols - 1);

                           mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
                     }
                     else
                     {
                           /* 010 */
                           /* Left point */
                           z1 = *(magptr - ncols);
                           z2 = *(magptr - ncols + 1);

                           mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

                           /* Right point */
                           z1 = *(magptr + ncols);
                           z2 = *(magptr + ncols - 1);

                           mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
                     }
                  }
                  else
                  {
                     if (-gx > -gy)
                     {
                           /* 001 */
                           /* Left point */
                           z1 = *(magptr + 1);
                           z2 = *(magptr + ncols + 1);

                           mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

                           /* Right point */
                           z1 = *(magptr - 1);
                           z2 = *(magptr - ncols - 1);

                           mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
                     }
                     else
                     {
                           /* 000 */
                           /* Left point */
                           z1 = *(magptr + ncols);
                           z2 = *(magptr + ncols + 1);

                           mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

                           /* Right point */
                           z1 = *(magptr - ncols);
                           z2 = *(magptr - ncols - 1);

                           mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
                     }
                  }
               }

               /* Now determine if the current point is a maximum point */

               if ((mag1 > 0.0) || (mag2 > 0.0))
               {
                  *resultptr = (unsigned char) NOEDGE;
               }
               else
               {
                  if (mag2 == 0.0)
                     *resultptr = (unsigned char) NOEDGE;
                  else
                     *resultptr = (unsigned char) POSSIBLE_EDGE;
               }
         }
      }
   }
};

SC_MODULE(DataOut)
{
   sc_fifo_in<IMAGE> ImgIn;
   sc_fifo_out<IMAGE> ImgOut;

   IMAGE Image;
   int i = 0;

   void main(void)
   {
      for(i=0; i<IMG_NUM; i++)
      {
         ImgIn.read(Image);
         ImgOut.write(Image);
      }
   }

   SC_CTOR(DataOut)
   {
      SC_THREAD(main);
      SET_STACK_SIZE;
   }
};

SC_MODULE(Platform)
{
   DataIn din;
   DUT canny;
   DataOut dout;

   sc_fifo_in<IMAGE> ImgIn;
   sc_fifo_out<IMAGE> ImgOut;
   sc_fifo<IMAGE> q1, q2;

   void before_end_of_elaboration(void)
   {
      din.ImgIn.bind(ImgIn);
      din.ImgOut.bind(q1);
      canny.ImgIn.bind(q1);
      canny.ImgOut.bind(q2);
      dout.ImgIn.bind(q2);
      dout.ImgOut.bind(ImgOut);
   }

   SC_CTOR(Platform) : din("din"), dout("dout"), canny("canny"),  q1("q1", 1), q2("q2", 1)
   {}
};

SC_MODULE(Top)
{
   Stimulus stimulus;
   Monitor monitor;
   Platform platform;

   sc_fifo<IMAGE> q1, q2;

   void before_end_of_elaboration(void)
   {
      stimulus.ImgOut.bind(q1);
      platform.ImgIn.bind(q1);
      platform.ImgOut.bind(q2);
      monitor.ImgIn.bind(q2);
   }

   SC_CTOR(Top):stimulus("stimulus"), platform("platform"), monitor("monitor"), q1("q1", 1), q2("q2", 1)
   {}
};
Top top("top");

int sc_main(int argc, char*[])
{
   sc_start();
   return 0;
}