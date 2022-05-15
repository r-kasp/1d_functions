
#ifndef WINDOW_H
#define WINDOW_H

#include <QtWidgets/QtWidgets>


void set_points (int n, double a, double b, double *x, double *y, double (*f) (double), double &f_max);
void calculate_sum (int n, int k, double *x, double *y, double *fx);
int bin_search (int n, double *x, double point);
double newton (int n, int k, double *x, double *y, double *fx, double point);
double cube_spline (int n, double *x, double *fx, double point, double (*d) (double));
void set_and_solve_system (int n, int k, double *a, double *b, double *x, double *fx, double (*d) (double));
double cube_spline (int n, int k, double *x, double *y, double *fx, double point, double (*d) (double));

class Window : public QWidget
{
  Q_OBJECT

private:
  int func_id;
  int line_up;
  const char *f_name;
  double a;
  double b;
  int n;
  int k;
  int p;
  double scale;
  double (*f) (double);
  double (*d) (double);
  double *x = nullptr;
  double *y = nullptr;
  double *fx = nullptr;
  double *A = nullptr;
  double *B = nullptr;
  double f_max;
  
public:
  Window (QWidget *parent, double a, double b, int n, int k);
  ~Window ();

  QSize minimumSizeHint () const;
  QSize sizeHint () const;
  void redeclare_memory ();

public slots:
  void set_func ();
  void change_line_up ();
  void increase_n ();
  void decrease_n ();
  void increase_scale ();
  void decrease_scale ();
  void increase_p ();
  void decrease_p ();

protected:
  void paintEvent (QPaintEvent *event);
};

#endif
