
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include <iostream>

#include "window.h"

int main (int argc, char *argv[])
{
  double a, b;
  int n, func_id;
  
  if (   argc != 5 
      || sscanf (argv[1], "%lf", &a) != 1
      || sscanf (argv[2], "%lf", &b) != 1
      || b - a < 1.e-6
      || sscanf (argv[3], "%d", &n) != 1
      || n <= 0
      || sscanf (argv[4], "%d", &func_id) != 1)
    {
      printf ("Wrong input, right : %s a b n k\n", argv[0]);
      return -1;
    }

  QApplication app (argc, argv);

  QMainWindow *window = new QMainWindow;
  QMenuBar *tool_bar = new QMenuBar (window);
  Window *graph_area = new Window (window, a, b, n, func_id);
  QAction *action;
    
  action = tool_bar->addAction ("&Change function", graph_area, SLOT (set_func ()));
  action->setShortcut (QString ("0"));
  
  action = tool_bar->addAction ("&Change line-up", graph_area, SLOT (change_line_up ()));
  action->setShortcut (QString ("1"));
  
  action = tool_bar->addAction ("&Increase scale", graph_area, SLOT (increase_scale ()));
  action->setShortcut (QString ("2"));
  
  action = tool_bar->addAction ("&Decrease scale", graph_area, SLOT (decrease_scale ()));
  action->setShortcut (QString ("3"));
  
  action = tool_bar->addAction ("&Increase n", graph_area, SLOT (increase_n ()));
  action->setShortcut (QString ("4"));
  
  action = tool_bar->addAction ("&Decrease n", graph_area, SLOT (decrease_n ()));
  action->setShortcut (QString ("5"));

  action = tool_bar->addAction ("&Increase p", graph_area, SLOT (increase_p ()));
  action->setShortcut (QString ("6"));
  
  action = tool_bar->addAction ("&Decrease p", graph_area, SLOT (decrease_p ()));
  action->setShortcut (QString ("7"));  

  action = tool_bar->addAction ("E&xit", window, SLOT (close ()));
  action->setShortcut (QString ("Ctrl+X"));

  tool_bar->setMaximumHeight (30);

  window->setMenuBar (tool_bar);
  window->setCentralWidget (graph_area);
  window->setWindowTitle ("Graph");

  window->show ();
  app.exec ();
  delete window;
  return 0;
}
