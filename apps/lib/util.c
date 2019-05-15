/**
 * UTIL.C -- Utility routines for the package.
 *
 *	   bswap2  (a, b, nbytes)
 *	   bswap4  (a, aoff, b, boff, nbytes)
 *	   bswap8  (a, aoff, b, boff, nbytes)
 *     is_swapped  ()
 */


#include <stdio.h>
#include <unistd.h>
#include <string.h>
#ifdef ULTRIX
#include <sys/types.h>
#endif
#include <unistd.h>

#ifndef AIXV3
#ifndef OSF1
typedef unsigned char uchar;
#endif
#endif



/* BSWAP2 - Move bytes from array "a" to array "b", swapping successive
 * pairs of bytes.  The two arrays may be the same but may not be offset
 * and overlapping.
 */
void
bswap2 (
  char    *a,         	                        // input array
  char    *b,         	                        // output array
  int     nbytes         	                // number of bytes to swap
)
{
        register char *ip=a, *op=b, *otop;
        register unsigned temp;

        /* Swap successive pairs of bytes.
         */
        for (otop = op + (nbytes & ~1);  op < otop;  ) {
            temp  = *ip++;
            *op++ = *ip++;
            *op++ = temp;
        }

        /* If there is an odd byte left, move it to the output array.
         */
        if (nbytes & 1)
            *op = *ip;
}


/* BSWAP4 - Move bytes from array "a" to array "b", swapping the four bytes
 * in each successive 4 byte group, i.e., 12345678 becomes 43218765.
 * The input and output arrays may be the same but may not partially overlap.
 */
void
bswap4 (
  char	*a,			                // input array
  int	aoff,			                // first byte in input array
  char	*b,			                // output array
  int	boff,			                // first byte in output array
  int	nbytes			                // number of bytes to swap
)
{
	register char	*ip, *op, *tp;
	register int	n;
	static	char temp[4];

	tp = temp;
	ip = (char *)a + aoff - 1;
	op = (char *)b + boff - 1;

	/* Swap successive four byte groups.
	 */
	for (n = nbytes >> 2;  --n >= 0;  ) {
	    *tp++ = *ip++;
	    *tp++ = *ip++;
	    *tp++ = *ip++;
	    *tp++ = *ip++;
	    *op++ = *--tp;
	    *op++ = *--tp;
	    *op++ = *--tp;
	    *op++ = *--tp;
	}

	/* If there are any odd bytes left, move them to the output array.
	 * Do not bother to swap as it is unclear how to swap a partial
	 * group, and really incorrect if the data is not modulus 4.
	 */
	for (n = nbytes & 03;  --n >= 0;  )
	    *op++ = *ip++;
}


/* BSWAP8 - Move bytes from array "a" to array "b", swapping the eight bytes
 * in each successive 8 byte group, i.e., 12345678 becomes 87654321.
 * The input and output arrays may be the same but may not partially overlap.
 */
void
bswap8 (
  char  *a,			                // input array
  int	aoff,			                // first byte in input array
  char  *b,			                // output array
  int	boff,			                // first byte in output array
  int	nbytes		                        // number of bytes to swap
)
{
	register char	*ip, *op, *tp;
	register int	n;
	static	char temp[8];

	tp = temp;
	ip = (char *)a + aoff - 1;
	op = (char *)b + boff - 1;

	/* Swap successive eight byte groups.
	 */
	for (n = nbytes >> 3;  --n >= 0;  ) {
	    *tp++ = *ip++;
	    *tp++ = *ip++;
	    *tp++ = *ip++;
	    *tp++ = *ip++;
	    *tp++ = *ip++;
	    *tp++ = *ip++;
	    *tp++ = *ip++;
	    *tp++ = *ip++;
	    *op++ = *--tp;
	    *op++ = *--tp;
	    *op++ = *--tp;
	    *op++ = *--tp;
	    *op++ = *--tp;
	    *op++ = *--tp;
	    *op++ = *--tp;
	    *op++ = *--tp;
	}

	/* If there are any odd bytes left, move them to the output array.
	 * Do not bother to swap as it is unclear how to swap a partial
	 * group, and really incorrect if the data is not modulus 8.
	 */
	for (n = nbytes & 03;  --n >= 0;  )
	    *op++ = *ip++;
}


/* IS_SWAPPED -- See if this is a byte-swapped (relative to Sun) machine.
 */
int
is_swapped ()
{
        union {
            char ch[4];
            int  i;
        } u;

        u.i = 1;
        return ((int) u.ch[0]);
}




/*  SSTRIP -- Strip leading and trailing spaces in a string.
 */
char *
sstrip (char *s)
{
    char *ip = s;

    if (!s || !*s)
        return s;

    /* Remove trailing spaces.  */
    for (ip=(s + strlen(s) - 1); *ip == ' ' && ip > s; ip--) *ip = '\0';

    /* Remove leading spaces.   */
    for (ip=s; *ip && *ip == ' '; ip++) ;

    return (ip);
}
