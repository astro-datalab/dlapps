/**
 *  DLPNG -- Utilities to convert to PNG images.
 */

#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "fitsio.h"
#include "../dlApps.h"


#define	PNG_OK		0
#define	PNG_ERR		1

/*  PNG bitmap/pixel structures.
 */
typedef struct { 			// A colored pixel
    uint8_t red;
    uint8_t green;
    uint8_t blue;
} pixel_t;

typedef struct { 			// A 2-D image
    pixel_t *pixels;
    size_t   width;
    size_t   height;
} bitmap_t;
    


/*  Public methods.
 */
int dl_FITS2PNG (char *fits_fname, char *png_fname);


/*  Private methods.
 */
static int png_saveToFile (bitmap_t *bitmap, const char *path);

/*  External scaling procedures.
 */
extern void dl_applyZScale (unsigned char **pix, int nx, int ny, int bitpix, 
	float z1, float z2);
extern void dl_imgZScale (unsigned char *im, int nx, int ny, int bitpix,
    float *z1, float *z2, float contrast, int opt_size, int len_stdline);




/**
 *  DL_FITS2PNG -- Convert the FITS image located in 'from' to a scaled
 *  PNG preview image in the 'to' file.
 */
int
dl_FITS2PNG (char *fits_fname, char *png_fname)
{
    bitmap_t img;
    int   i, bitpix, bytepix, naxis, npix, anynul, datatype;
    int   totpix, status = 0;
    long  naxes[7], first = 1;
    fitsfile *fptr;
    pixel_t  *p;
    uint8_t  *a;
    double   *array, nulval=0.0;


    /*  Open the input FITS image and get the dimensions.
     */
    fits_open_file (&fptr, fits_fname, READONLY, &status);
    if (status) {
        fits_report_error (stderr, status);
	return (ERR);
    }
    fits_get_img_param (fptr, 2, &bitpix, &naxis, naxes, &status);

    /*  Compute pixel and array sizes.
     */
    switch (bitpix) {
    case BYTE_IMG:      datatype = TBYTE;       break;
    case SHORT_IMG:     datatype = TSHORT;      break;
    case LONG_IMG:      datatype = TLONG;       break;
    case FLOAT_IMG:     datatype = TFLOAT;      break;
    case DOUBLE_IMG:    datatype = TDOUBLE;     break;
    }
    npix = totpix = naxes[0] * naxes[1];
    bytepix = abs(bitpix) / 8;

    /*  Allocate enough space for the entire array.  Note this may be larger
     *  than what we actually need for the raster.
     */
    array = (double *) calloc (npix, sizeof(double));

    /*  Read the image array from the file.
     */
    while (totpix > 0 && !status) {
       fits_read_img (fptr, datatype, first, npix, &nulval, array, &anynul,
	   &status);
       totpix = totpix - npix;
       first  = first  + npix;
    }


    /* Compute the optimal zscale values.  The image raster is modified
     * in-place and loaded into the PNG structure.
     */
    int  nx = naxes[0], ny = naxes[1], nsample=900;
    float  z1, z2, contrast = 0.25, np = nx * ny;

    nsample = max (1000, (int) (npix * 0.005));

    dl_imgZScale ((unsigned char *)array, nx, ny, bitpix, &z1, &z2, 
	contrast, nsample, 
        max (2, (int) ((float)ny / sqrt(np / (float)nsample))) );
    //printf ("z1 = %g  z2 = %g\n", z1, z2);
    a = (unsigned char *) array;
    dl_applyZScale (&a, nx, ny, bitpix, z1, z2);

    /*  Now create the output PNG file.  First setup the data structure with
     *  dimensions and the scaled 8-bit values.  Note we're assuming a
     *  gray-scale colormap.
     */
    img.width = naxes[0]; 			// Create an image
    img.height = naxes[1];
    img.pixels = calloc (sizeof (pixel_t), img.width * img.height);

    a = (unsigned char *) array;		// load png bitmap
    for (p=img.pixels, i=0; i < npix; i++, a++, p++)
	p->red = p->green = p->blue = *a;
    free (array);

    /* Write the image to a file 'img.png'.
     */
    png_saveToFile (&img, png_fname);

    return 0;
}


/**
 *  PNG_SAVETPFILE --  Write "bitmap" to a PNG file specified by "path".
 *
 *  @returns 		0 on success, non-zero on error. 
 */
static int
png_saveToFile (bitmap_t *bitmap, const char *path)
{
    FILE * fp;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    png_byte **row_pointers = NULL, *row = NULL;
    int   x, y, iy = 0;
    int pixel_size = 3;			// number of colors per pixel
    int depth = 8;			// depth of each pixel
    

    if (! (fp = fopen (path, "wb")) ) {
	fprintf (stderr, "Failed to open file '%s'\n", path);
        return PNG_ERR;
    }

    /*  Initialize the output PNG file.
     */
    if ((png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING,
	NULL, NULL, NULL)) == NULL) {
	    fprintf (stderr, "Error: png_create_write_struct returns NULL\n");
            fclose (fp);
            return PNG_ERR;
    }
    
    if ((info_ptr = png_create_info_struct (png_ptr)) == NULL) {
	fprintf (stderr, "Error: png_create_info_struct returns NULL\n");
        png_destroy_write_struct (&png_ptr, &info_ptr);
        fclose (fp);
        return PNG_ERR;
    }
    
    if (setjmp (png_jmpbuf (png_ptr))) { 	// set up error handling
	fprintf (stderr, "Error: png_jmpbuf failed to set\n");
        png_destroy_write_struct (&png_ptr, &info_ptr);
        fclose (fp);
	return PNG_ERR;
    }
    
    /*  Set image attributes.
     */
    png_set_IHDR (png_ptr, info_ptr, bitmap->width, bitmap->height, depth,
                  PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                  PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    
    /* Initialize rows of the PNG.  Note we y-flip the image at this point.
     */
    row_pointers = png_malloc (png_ptr, bitmap->height * sizeof (png_byte *));
    for (y = bitmap->height-1; y >= 0; y--) {
        row = png_malloc (png_ptr, sizeof(uint8_t) * bitmap->width*pixel_size);
        row_pointers[iy++] = row;
        for (x = 0; x < bitmap->width; ++x) {
            pixel_t * pixel = (bitmap->pixels + bitmap->width * y + x);
            *row++ = pixel->red;
            *row++ = pixel->green;
            *row++ = pixel->blue;
        }
    }
    
    /* Write the image data to the output file.
     */
    png_init_io (png_ptr, fp);
    png_set_rows (png_ptr, info_ptr, row_pointers);
    png_write_png (png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

    /* The routine has successfully written the file, clean up and return.
     */
    for (y = 0; y < bitmap->height; y++)
        png_free (png_ptr, row_pointers[y]);
    png_free (png_ptr, row_pointers);
    
    png_destroy_write_struct (&png_ptr, &info_ptr);
    fclose (fp);

    return (PNG_OK);
}


/**
 *  Task main().
 */
#ifdef PNG_TASK_TEST
int main (int argc, char *argv[])
{
    dl_FITS2PNG (argv[1], argv[2]);
}
#endif
