/******************************************************************************
 * Copyright (C) 2024 ZHJ. All Rights Reserved.
 *****************************************************************************/
#include "dijkstra.h"
#include "matplotlibcpp.h"
#include <iostream>
namespace plt = matplotlibcpp;
int main() {
  std::vector<double> start{-5, 10}, goal{50, 50};
  double grid_size = 2.0;
  double robot_radius = 1.0;

  std::vector<int> obstacle_x_list;
  std::vector<int> obstacle_y_list;

  // 障碍物边缘
  for (double i = -10; i < 60; i++) {
    obstacle_x_list.push_back(i);
    obstacle_y_list.push_back(-10.0);
  }
  for (double i = -10; i < 60; i++) {
    obstacle_x_list.push_back(60.0);
    obstacle_y_list.push_back(i);
  }
  for (double i = -10; i < 61; i++) {
    obstacle_x_list.push_back(i);
    obstacle_y_list.push_back(60.0);
  }
  for (double i = -10; i < 61; i++) {
    obstacle_x_list.push_back(-10.0);
    obstacle_y_list.push_back(i);
  }
  for (double i = -10; i < 40; i++) {
    obstacle_x_list.push_back(20.0);
    obstacle_y_list.push_back(i);
  }
  for (double i = -10; i < 50; i++) {
    obstacle_x_list.push_back(10.0);
    obstacle_y_list.push_back(i);
  }
  for (double i = -10; i < 50; i++) {
    obstacle_x_list.push_back(4.0);
    obstacle_y_list.push_back(i);
  }
  for (double i = -0; i < 50; i++) {
    obstacle_x_list.push_back(30.0);
    obstacle_y_list.push_back(60 - i);
  }
  for (double i = 0; i < 40; i++) {
    obstacle_x_list.push_back(40.0);
    obstacle_y_list.push_back(60.0 - i);
  }
  
  Dijkstra dijkstra_planner(robot_radius, grid_size);
  dijkstra_planner.Init(obstacle_x_list, obstacle_y_list);
  dijkstra_planner.Planning(start[0], start[1], goal[0], goal[1]);
  return 0;
}