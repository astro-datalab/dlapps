#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "fitsio.h"

int main(int argc, char *argv[])
{
    /* FITS file pointers defined in fitsio.h */
    fitsfile *infptr, *outfptr, *phufptr;

    int status = 0, ii = 1, iteration = 0, single = 0, hdupos;
    int hdutype, bitpix, bytepix, naxis = 0, nkeys, datatype = 0, anynul;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    long first, totpix = 0, npix;
    double *array, bscale = 1.0, bzero = 0.0, nulval = 0.;
    char card[81], fpath[160], *ip;


    memset (fpath, 0, 160);
    strcpy (fpath, argv[1]);

    /* Open the input file and create output file */
    fits_open_file (&infptr, fpath, READONLY, &status);

    if ((ip = strchr (fpath, '[')))		// get PHU pointer
	*ip = '\0';
    fits_open_file (&phufptr, fpath, READONLY, &status);

    fits_create_file (&outfptr, argv[2], &status);

    if (status != 0) {    
        fits_report_error(stderr, status);
        return(status);
    }

    fits_get_hdu_num(infptr, &hdupos);  /* Get the current HDU position */

    /* Copy only a single HDU if a specific extension was given */ 
    if (hdupos != 1 || strchr(argv[1], '[')) single = 1;

    for (; !status; hdupos++)  /* Main loop through each extension */
    {

      fits_get_hdu_type(infptr, &hdutype, &status);

      if (hdutype == IMAGE_HDU) {
          /* get image dimensions and total number of pixels in image */
          for (ii = 0; ii < 9; ii++)
              naxes[ii] = 1;

          fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status);

          totpix = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4]
             * naxes[5] * naxes[6] * naxes[7] * naxes[8];
      }

      if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0) { 
          /* just copy tables and null images */
          fits_copy_hdu(infptr, outfptr, 0, &status);

      } else {
          /* Explicitly create new image, to support compression */
          fits_create_img(outfptr, bitpix, naxis, naxes, &status);

          /* copy all the PHU user keywords (not the structural keywords) */
          fits_get_hdrspace (phufptr, &nkeys, NULL, &status); 
          for (ii = 1; ii <= nkeys; ii++) {
              fits_read_record (phufptr, ii, card, &status);
              if (fits_get_keyclass (card) > TYP_CMPRS_KEY)
                  fits_write_record (outfptr, card, &status);
          }

          /* copy all the user keywords (not the structural keywords) */
          fits_get_hdrspace (infptr, &nkeys, NULL, &status); 
          for (ii = 1; ii <= nkeys; ii++) {
              fits_read_record (infptr, ii, card, &status);
              if (fits_get_keyclass (card) > TYP_CMPRS_KEY)
                  fits_write_record (outfptr, card, &status);
          }


          switch(bitpix) {
              case BYTE_IMG:
                  datatype = TBYTE;
                  break;
              case SHORT_IMG:
                  datatype = TSHORT;
                  break;
              case LONG_IMG:
                  datatype = TLONG;
                  break;
              case FLOAT_IMG:
                  datatype = TFLOAT;
                  break;
              case DOUBLE_IMG:
                  datatype = TDOUBLE;
                  break;
          }

          bytepix = abs(bitpix) / 8;

          npix = totpix;
          iteration = 0;

          /* try to allocate memory for the entire image */
          /* use double type to force memory alignment */
          array = (double *) calloc(npix, bytepix);

          /* if allocation failed, divide size by 2 and try again */
          while (!array && iteration < 10)  {
              iteration++;
              npix = npix / 2;
              array = (double *) calloc(npix, bytepix);
          }

          if (!array)  {
              printf("Memory allocation error\n");
              return(0);
          }

          /* turn off any scaling so that we copy the raw pixel values */
          fits_set_bscale(infptr,  bscale, bzero, &status);
          fits_set_bscale(outfptr, bscale, bzero, &status);

          first = 1;
          while (totpix > 0 && !status)
          {
             /* read all or part of image then write it back to output file */
             fits_read_img(infptr, datatype, first, npix, 
                     &nulval, array, &anynul, &status);

             fits_write_img(outfptr, datatype, first, npix, array, &status);
             totpix = totpix - npix;
             first  = first  + npix;
          }
          free(array);
      }

      if (single) break;  /* quit if only copying a single HDU */
      fits_movrel_hdu(infptr, 1, NULL, &status);  /* try to move to next HDU */
    }

    if (status == END_OF_FILE)  status = 0; /* Reset after normal error */

    fits_close_file(outfptr,  &status);
    fits_close_file(infptr, &status);
    fits_close_file(phufptr, &status);

    /* if error occurred, print out error message */
    if (status)
       fits_report_error(stderr, status);
    return(status);
}
