#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "matrix.h"

int solve_matrix (int n, double *a, double *b, double eps)
{ 
  for (int i = 0; i < n; i++)
    {
      double main_elt = a[i * 3 + 1];
      if (fabs (main_elt) < eps)
        return -1;
        
      //a[i * 3 + 1] = 1; 
      b[i] /= main_elt;
      if (i == n - 1)
        break;
      a[i * 3 + 2] /= main_elt;
      
      double mult = a[(i + 1) * 3 + 0];
      //a[(i + 1) * 3 + 0] = 0;
      a[(i + 1) * 3 + 1] -= a[i * 3 + 2] * mult;
        
      b[i + 1] -= b[i] * mult;
    }
    
  for (int i = n - 1; i > 0; i--)
    {
      double mult = a[(i - 1) * 3 + 2];
      //a[(i - 1) * 3 + 2] = 0;
      b[i - 1] -= b[i] * mult;
    }
  return 0;
}

