#include "Heegaard.h"
#include "Heegaard_Dec.h"
 
/*  sort elements "first" through "last"-1 */

static void qksort1(first, last)
{
	static 		int i;		/*  "static" to save stack space  */
	register 	int j;
 
	while (last - first > 1) 
		{
		i = first;
		j = last;
		for (;;) 
			{
			while (++i < last && qkst_compare(i, first) < 0) 	;
			while (--j > first && qkst_compare(j, first) > 0)	;
			if (i >= j)	break;
			qkst_swap(i, j);
		}
		qkst_swap(first, j);
		if (j - first < last - (j + 1)) 
			{
			qksort1(first, j);
			first = j + 1;			/*  qsort1(j + 1, last);  */
			}
		else 
			{
			qksort1(j + 1, last);
			last = j;				/*  qsort1(first, j);  */
			}
	}
}


/*  sort "nelems" elements, using user's "qkst_compare" and "qkst_swap" routines */

void qksort(unsigned int nelems)
{
	qksort1(0, nelems);
}
