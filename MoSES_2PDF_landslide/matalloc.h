/********************************************************
  This header file contains macros for allocating memory
  of matrices and vectors. The routines are:

  * NEW_ARRAY(array, type, nnn)
  * NEW_MATRIX(mat, type, nrow, ncol)
  * NEW_3DMATRIX(mat, type, nrow, ncol, ncomp)
  * FREE_ARRAY(array )
  * FREE_MATRIX(matrix)
  * FREE_3DMATRIX(matrix)
     
  *****************************************************/
#include <stdlib.h>
#include <stdio.h>


#define NEW_ARRAY(array, type, nnn)                        \
{                                                          \
  if (!( array = (type *) calloc(nnn, sizeof(type)) ))     \
    {                                                      \
      fprintf(stderr,"Allocation failure - aborting \n");  \
      exit(1);                                             \
    }                                                      \
}

     
#define NEW_MATRIX(mat, type, nrow, ncol)                  \
{                                                          \
  int pPp, qQq;                                            \
							   \
  if (!( mat = (type **) calloc(nrow, sizeof(type *)) ))   \
    {                                                      \
      fprintf(stderr,"Allocation failure - aborting \n");  \
      exit(1);                                             \
    }                                                      \
							   \
  if (!(mat[0] = (type *)                                  \
	calloc((nrow)*(ncol), sizeof(type)) ))             \
    {                                                      \
      fprintf(stderr,"Allocation failure - aborting \n");  \
      exit(1);                                             \
    }                                                      \
                                                           \
							   \
  for (qQq=1, pPp=ncol; qQq<nrow; qQq++)                   \
    {                                                      \
      mat[qQq] = mat[0]+pPp;                               \
      pPp += (ncol);                                       \
    }                                                      \
}


#define NEW_3DMATRIX(mat, type, nrow, ncol, ncomp)         \
{                                                          \
  int pPp, qQq, rRr;                                       \
							   \
  if (!( mat = (type ***) calloc(nrow, sizeof(type **)) )) \
    {                                                      \
      fprintf(stderr,"Allocation failure - aborting \n");  \
      exit(1);                                             \
    }                                                      \
							   \
  if (!(mat[0] = (type **)                                 \
	calloc((nrow)*(ncol), sizeof(type *)) ))           \
    {                                                      \
      fprintf(stderr,"Allocation failure - aborting \n");  \
      exit(1);                                             \
    }                                                      \
                                                           \
							   \
  for (qQq=1, pPp=ncol; qQq<nrow; qQq++)                   \
    {                                                      \
      mat[qQq] = mat[0]+pPp;                               \
      pPp += (ncol);                                       \
    }                                                      \
                                                           \
  if (!(mat[0][0] = (type *)                               \
	calloc((nrow)*(ncol)*(ncomp), sizeof(type)) ))     \
    {                                                      \
      fprintf(stderr,"Allocation failure - aborting \n");  \
      exit(1);                                             \
    }                                                      \
							   \
  for (qQq=0, rRr=0; qQq<nrow; qQq++)                      \
    for (pPp=0; pPp<ncol; pPp++)                           \
      {                                                    \
	mat[qQq][pPp] = mat[0][0]+rRr;                     \
	rRr += (ncomp);                                    \
      }                                                    \
}

#define FREE_ARRAY(array )                                 \
{                                                          \
  free(array);                                             \
}


#define FREE_MATRIX(matrix)                                \
{                                                          \
  free(*matrix);                                           \
  free(matrix);                                            \
}

#define FREE_3DMATRIX(matrix)                              \
{                                                          \
  free(**matrix);                                          \
  free(*matrix);                                           \
  free(matrix);                                            \
}
							     
