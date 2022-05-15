
#include <QPainter>
#include <stdio.h>

#include "window.h"
#include "matrix.h"
#include <iostream>

#define DEFAULT_A -10
#define DEFAULT_B 10
#define DEFAULT_N 10
#define EPS 1e-14

static
double f_0 (double x)
{
  (void)x;
  return 1;
}

static
double d_0 (double x)
{
  (void)x;
  return 0;
}

static
double f_1 (double x)
{
  return x;
}

static
double f_2 (double x)
{
  return x * x;
}

static
double d_2 (double x)
{
  return 2 * x;
}

static
double f_3 (double x)
{
  return x * x * x;
}

static
double d_3 (double x)
{
  return 3 * x * x;
}

static
double f_4 (double x)
{
  return x * x * x * x;
}


static
double d_4 (double x)
{
  return 4 * x * x * x;
}

static
double f_5 (double x)
{
  return exp (x);
}

static
double f_6 (double x)
{
  return 1.0 / (25 * x * x + 1);
}

static
double d_6 (double x)
{
  return (50 * x) * (-1.0 / ((25 * x * x + 1) * (25 * x * x + 1))) ;
}

void set_points (int n, double a, double b, double *x, double *y, double (*f) (double), double &f_max)
{
  f_max = 0;
  double step = (b - a) / (n - 1);
  for (int i = 0; i < n; i++)
    {
      x[i] = a + i * step;
      y[i] = f (x[i]);
      if (fabs (y[i]) > f_max)
        f_max = fabs (y[i]);
    }
}


void calculate_sum (int n, int k, double *x, double *y, double *fx)
{
  memset (fx, 0, n * (k + 1) * sizeof (double));
  for (int i = 0; i < n; i++)
    {
      fx [i * (k + 1)] = y [i];
    }
  for (int j = 1; j < k + 1; j++)
    {
      for (int i = 0; i < n - j; i++)
        {
          if (fabs (x[j + i] - x[i]) < EPS)
            return;
          fx[i * (k + 1) + j] = (fx[(i + 1) * (k + 1) + j - 1] - fx[i * (k + 1) + j - 1]) / (x[j + i] - x[i]);
        }
    }
}

int bin_search (int n, double *x, double point)
{
  int l = 0, r = n - 1;
  while (r - l > 1)
    {
      int m = (r + l) / 2;
      if (x[m] > point)
        r = m;
      else
        l = m;
    }
  return l;
}


double newton (int n, int k, double *x, double *y, double *fx, double point)
{
  (void)y;
  
  int i = bin_search (n, x, point);
  if (i > n - k - 1)
    i = n - k - 1;
   
  double res = fx[i * (k + 1) + k];
  
  //f(xi) + (x - x_i)(f(xi, xi+1) + (x - x_i+1)(f(xi, xi+1, xi+2) + (x - x_i+2)f(xi, xi+1, xi+2, xi+3)))
  
  for (int var_cnt = k + 1; var_cnt >= 2; var_cnt--)
    {
      res = res * (point - x[i + var_cnt - 2]) + fx[i * (k + 1) + var_cnt - 2];
    }
    
  return res;
}


void set_and_solve_system (int n, int k, double *a, double *b, double *x, double *fx, double (*d) (double))
{
  memset (a, 0, n * 3 * sizeof (double));
  for (int i = 1; i < n - 1; i++)
    {
      a[i * 3 + 0] = x[i + 1] - x[i];
      a[i * 3 + 1] = 2 * (x[i + 1] - x[i - 1]);
      a[i * 3 + 2] = x[i] - x[i - 1];
      b[i] = 3 * fx [(i - 1) * (k + 1) + 1] * (x[i + 1] - x[i]) 
               + 3 * fx [i * (k + 1) + 1] * (x[i] - x[i - 1]);
    } 
  a[0 * 3 + 1] = 1;
  b[0] = d (x[0]);
  a[(n - 1) * 3 + 1] = 1;
  b[n - 1] = d (x[n - 1]);
  
  solve_matrix (n, a, b, EPS);
}


double cube_spline (int n, int k, double *x, double *y, double *fx, double *d, double point)
{
  int i = bin_search (n, x, point);
  
  if (fabs (x[i + 1] - x[i]) < EPS)
    return 0;
  //double fx_in_1 = (i + 1 < n ? (y[i + 1] - y[i]) / (x[i + 1] - x[i]) : 0);
  double c_1i = y[i];
  double c_2i = d[i];
  double c_3i = (3 * fx[i * (k + 1) + 1] - 2 * d[i]  - d [i + 1]) / (x[i + 1] - x[i]);
  double c_4i = (d[i] + d [i + 1] - 2 * fx[i * (k + 1) + 1]) / ((x[i + 1] - x[i]) * (x[i + 1] - x[i]));
  
  double result = c_1i
                + c_2i * (point - x[i])
                + c_3i * (point - x[i]) * (point - x[i])
                + c_4i * (point - x[i]) * (point - x[i]) * (point - x[i]);
  
  return result;
}


Window::Window (QWidget *parent, double in_a, double in_b, int in_n, int in_func_id)
  : QWidget (parent)
{
  a = in_a;
  b = in_b;
  n = in_n;
  func_id = in_func_id;
  func_id--;
  k = 4;
  line_up = 0;
  scale = 1;
  p = 0;
  f_max = 0;
  
  redeclare_memory ();
  
  set_func ();
}

void Window::redeclare_memory () 
{
  if (x != nullptr)
    delete [] x;
  if (x != nullptr)
    delete [] y;
  if (x != nullptr)
    delete [] fx;
  if (x != nullptr)
    delete [] A;
  if (x != nullptr)
    delete [] B;
  x = new double[n];
  y = new double[n];
  fx = new double[n * (k + 1)];
  A = new double[n * 3];
  B = new double[n];
}

Window::~Window ()
{
  delete [] x;
  delete [] y;
  delete [] fx;
  delete [] A;
  delete [] B;
}

QSize Window::minimumSizeHint () const
{
  return QSize (100, 100);
}

QSize Window::sizeHint () const
{
  return QSize (1000, 1000);
}

/// change count of drawings
void Window::change_line_up ()
{
  line_up = (line_up + 1) % 4;
  
  update ();
}

void Window::increase_n ()
{
  n *= 2;
  redeclare_memory ();
  set_points (n, a, b, x, y, f, f_max);
  calculate_sum (n, k, x, y, fx);
  set_and_solve_system (n, k, A, B, x, fx, d);
  update ();
}

void Window::decrease_n ()
{
  if (n == 1)
    return;
  n /= 2;
  redeclare_memory ();
  set_points (n, a, b, x, y, f, f_max);
  calculate_sum (n, k, x, y, fx);
  set_and_solve_system (n, k, A, B, x, fx, d);
  update ();
}


void Window::increase_scale ()
{
  scale *= 2;
  update ();
}

void Window::decrease_scale ()
{
  scale /= 2;
  update ();
}

void Window::increase_p ()
{
  if (n == 1)
    return;
  p++;
  y[n / 2] += 0.1 * f_max;
  calculate_sum (n, k, x, y, fx);
  set_and_solve_system (n, k, A, B, x, fx, d);
  update ();
}

void Window::decrease_p ()
{
  if (n == 1)
    return;
  p--;
  y[n / 2] -= 0.1 * f_max;
  calculate_sum (n, k, x, y, fx);
  set_and_solve_system (n, k, A, B, x, fx, d);
  update ();
}


/// change current function for drawing
void Window::set_func ()
{
  func_id = (func_id + 1) % 7;
  p = 0;
  switch (func_id)
    {
      case 0:
        f_name = "f (x) = 1";
        f = f_0;
        d = d_0;
        break;
      case 1:
        f_name = "f (x) = x";
        f = f_1;
        d = f_0;
        break;
      case 2:
        f_name = "f (x) = x * x";
        f = f_2;
        d = d_2;
        break;
      case 3:
        f_name = "f (x) = x * x * x";
        f = f_3;
        d = d_3;
        break;
      case 4:
        f_name = "f (x) = x * x * x * x";
        f = f_4;
        d = d_4;
        break;
      case 5:
        f_name = "f (x) = exp (x)";
        f = f_5;
        d = f_5;
        break;
      case 6:
        f_name = "f (x) = 1 / (25 * x * x  + 1)";
        f = f_6;
        d = d_6;
        break;
    }  
  set_points (n, a, b, x, y, f, f_max);
  calculate_sum (n, k, x, y, fx);
  set_and_solve_system (n, k, A, B, x, fx, d);
 
  update ();
}

/// render graph
void Window::paintEvent (QPaintEvent * /* event */)
{
  a *= scale;
  b *= scale;  
  QPainter painter (this);
  double x1, x2, y1, y2;
  double max_y, min_y;
  double delta_y, delta_x = (b - a) / width();
  int wid = width ();
  QPen pen_black(Qt::black, 0, Qt::SolidLine); 
  QPen pen_red(Qt::red, 0, Qt::SolidLine); 
  QPen pen_green(Qt::green, 0, Qt::SolidLine); 
  QPen pen_blue(Qt::blue, 0, Qt::SolidLine); 

  painter.setPen (pen_black);

  // calculate min and max for current function
  max_y = min_y = 0;
  for (x1 = a; x1 - b < 1.e-6; x1 += delta_x)
    {
      y1 = f (x1);
      if (y1 < min_y)
        min_y = y1;
      if (y1 > max_y)
        max_y = y1;
    }

  delta_y = 0.01 * (max_y - min_y);
  min_y -= delta_y;
  max_y += delta_y;

  // save current Coordinate System
  painter.save ();

  // make Coordinate Transformations
  painter.translate (0.5 * width (), 0.5 * height ());
  double eps = 1e-14;
  if (fabs (max_y - min_y) < eps)
    {
      max_y = 1 + eps;
      min_y = 1 - eps;
    }
  painter.scale (width () / (b - a), -height () / (max_y - min_y));
  painter.translate (-0.5 * (a + b), -0.5 * (min_y + max_y));

  if (line_up != 3)
    {
      // draw approximated line for graph
      
      x1 = a;
      y1 = f (x1);
      int ind = 1;
      for (int q = 1; q <= wid; q++) 
        {
          x2 = a + delta_x * q;
          y2 = f (x2) + p * (ind == n / 2) * 0.1 * f_max;
          
          painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));

          x1 = x2, y1 = y2;
          ind++;
        }
      x2 = b;
      y2 = f (x2);
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
    }
    
  if (n <= 50 && (line_up == 0 || line_up == 2))
    {  
      //draw newton line
      painter.setPen (pen_green);
      
      x1 = a;
      y1 = newton (n, k, x, y, fx, x1);
      for (int q = 1; q <= wid; q++) 
        {
          x2 = a + delta_x * q;
          y2 = newton (n, k, x, y, fx, x2);
          
          painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));

          x1 = x2, y1 = y2;
        }
      x2 = b;
      y2 = newton (n, k, x, y, fx, x2);
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
    }
    
  if (line_up == 1 || line_up == 2)
    {
      //draw cube spline
      painter.setPen (pen_blue);
      
      x1 = a;
      y1 = cube_spline (n, k, x, y, fx, B, x1);
      for (int q = 1; q <= wid; q++) 
        {
          x2 = a + delta_x * q;
          y2 = cube_spline (n, k, x, y, fx, B, x2);
          
          painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));

          x1 = x2, y1 = y2;
        }
      x2 = b;
      y2 = cube_spline (n, k, x, y, fx, B, x2);
      painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
    }
    
  double residual_newton = 0, residual_spline = 0;
  if (line_up == 3)
    {
      eps = 1e-17;
      // calculate min and max for residual
      double max_rd = 0, min_rd = 0;
      int ind = 1;
      for (int q = 0; q <= wid; q++)
        {
          x1 = a + q * delta_x;
          y1 = fabs (f (x1) + p * (ind == n / 2) * 0.1 * f_max - cube_spline (n, k, x, y, fx, B, x1));
          ind++;
          if (y1 < min_rd)
            min_rd = y1;
          if (y1 > max_rd)
            max_rd = y1;
        }
      if (n <= 50)  
        {
          ind = 1;
          for (int q = 0; q <= wid; q++)
          {
            x1 = a + q * delta_x;
            y1 = fabs (f (x1) + p * (ind == n / 2) * 0.1 * f_max - newton (n, k, x, y, fx, x1));
            ind++;
            if (y1 < min_rd)
              min_rd = y1;
            if (y1 > max_rd)
              max_rd = y1;
          }
        }
      

      if (n <= 50)
        {
          painter.setPen (pen_green);
      
          x1 = a;
          y1 = fabs (f (x1) - newton (n, k, x, y, fx, x1));
          residual_newton = (residual_newton < y1 ? y1 : residual_newton);
          ind = 1;
          for (int q = 1; q <= wid; q++) 
            {
              x2 = a + delta_x * q;
              y2 = fabs (f (x2) + p * (ind == n / 2) * 0.1 * f_max - newton (n, k, x, y, fx, x2));
              residual_newton = (residual_newton < y2 ? y2 : residual_newton);
              if (fabs (max_rd - min_rd) < eps)
                y2 = 0;
              else
                y2 = (y2 - min_rd) / (max_rd - min_rd) / 1.2;
              ind++;
              painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));

              x1 = x2, y1 = y2;
            }
          x2 = b;
          y2 = fabs (f (x2) - newton (n, k, x, y, fx, x2));
          residual_newton = (residual_newton < y2 ? y2 : residual_newton);
          painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
        }
        
        painter.setPen (pen_blue);
      
        x1 = a;
        y1 = fabs (f (x1) - cube_spline (n, k, x, y, fx, B, x1));
        residual_spline = (residual_spline < y1 ? y1 : residual_spline);
        ind = 1;
        for (int q = 1; q <= wid; q++) 
          {
            x2 = a + delta_x * q;
            y2 = fabs (f (x2) + p * (ind == n / 2) * 0.1 * f_max - cube_spline (n, k, x, y, fx, B, x2));
            residual_spline = (residual_spline < y2 ? y2 : residual_spline);
            if (fabs (max_rd - min_rd) < eps)
              y2 = 0;
            else
              y2 = (y2 - min_rd) / (max_rd - min_rd) / 1.2;
            ind++;
            painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));

            x1 = x2, y1 = y2;
          }
        x2 = b;
        y2 = fabs (f (x2) - cube_spline (n, k, x, y, fx, B, x2));
        residual_spline = (residual_spline < y2 ? y2 : residual_spline);
        painter.drawLine (QPointF (x1, y1), QPointF (x2, y2));
        
        
        // draw axis
        painter.setPen (pen_red);
        painter.drawLine (a, 0, b, 0);
        painter.drawLine (0, 1, 0, -1);
    }
  
  if (line_up != 3)
  {
    // draw axis
    painter.setPen (pen_red);
    painter.drawLine (a, 0, b, 0);
    painter.drawLine (0, max_y, 0, min_y);
  }
  // restore previously saved Coordinate System
  painter.restore ();

  // render function name
  painter.setPen ("black");
  painter.drawText (0, 20, f_name);
  std::string s = "n = " + std::to_string (n);
  painter.drawText (0, 80, s.c_str ());
  s = "scale = " + std::to_string (scale);
  painter.drawText (0, 140, s.c_str ());
  double f_max = std::max (fabs (max_y), fabs (min_y));
  s = "max|F| = " + std::to_string (f_max);
  painter.drawText (0, 200, s.c_str ());
  s = "a = " + std::to_string (a);
  painter.drawText (0, 260, s.c_str ());
  s = "b = " + std::to_string (b);
  painter.drawText (0, 320, s.c_str ());
  s = "p = " + std::to_string (p);
  painter.drawText (0, 380, s.c_str ());
  if (line_up == 3)
    {
      printf ("residual_spline : %.6e\n", residual_spline);
      char buffer[32];
      memset(buffer, 0, sizeof(buffer));
      snprintf(buffer, sizeof(buffer), "%g", residual_spline);
      std::string str(buffer);
      s = "residual cube spline = " + str;
      painter.drawText (0, 440, s.c_str ());
      if (n <= 50)
        {
          printf ("residual_newton : %.6e\n", residual_newton);
          memset(buffer, 0, sizeof(buffer));
          snprintf(buffer, sizeof(buffer), "%g", residual_newton);
          std::string str2(buffer);
          s = "residual newton = " + str2;
          painter.drawText (0, 500, s.c_str ());
        }
    }
  a /= scale;
  b /= scale;
}
