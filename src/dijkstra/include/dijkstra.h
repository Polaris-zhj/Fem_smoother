/******************************************************************************
 * Copyright (C) 2024 ZHJ. All Rights Reserved.
 *****************************************************************************/
#pragma once

#include "config/config.h"
#include "fem_smoother.h"
#include <array>
#include <istream>
#include <memory>
#include <stdlib.h>
#include <time.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

using namespace Eigen;
class FemSmoother;
class Dijkstra {
public:
  struct Node {
    int x_index;
    int y_index;
    double cost;
    int parent_index;
    Node(const int &temp_x_index, const int &temp_y_index,
         const double &temp_cost, const int &temp_parent_index);
    std::string Debug() {
      std::ostringstream output;
      output << "x_index: " << x_index << ",y_index: " << y_index
             << ",cost: " << cost << ",parent_index: " << parent_index;
      return output.str();
    }
  };

  Dijkstra(const double &robot_radius, const double &resolution);
  void Init(const std::vector<int> &obstacle_x_list,
            const std::vector<int> &obstacle_y_list);
  void Planning(const int &start_x, const int &start_y, const int &goal_x,
                const int &goal_y);

private:
  void InitParam();

  /**
   * @brief 根据世界坐标系下的x或y坐标，计算的到在栅格地图中的坐标
   * @param real_xy
   * @param min_xy
   * @param xy_index
   * @return
   */
  int CalcXyIndex(const int &real_xy, const int &min_xy);
  int CalcIndex(const Node &node);
  int CalcPosition(const int &xy_index, const int &min_xy);
  void NextStepMotions();
  void CreateSearchMap(std::vector<int> obstacle_x_list,
                       std::vector<int> obstacle_y_list);
  bool VerifyNode(const Node &node);
  void RetrivePath(const Node *node, std::map<int, Node *> &closed_set,
                   std::vector<int> *const path_x_list,
                   std::vector<int> *const path_y_list);

  void DiscretePointBSplineSmoother(std::vector<int> &path_x_list,
                             std::vector<int> &path_y_list);
  void DiscretePointFemSmoother(std::vector<double> path_x_list,
                                std::vector<double> path_y_list);

private:
  int min_x_, max_x_, min_y_, max_y_;
  int min_x_index_, min_y_index_, max_x_index_, max_y_index_;
  std::vector<int> obstacle_x_list_, obstacle_y_list_;
  int x_width_;
  int y_width_;
  double resolution_;
  double robot_radius_;
  std::vector<std::array<double, 3>> next_step_motions_;
  std::vector<std::vector<bool>> obstacle_map_;
  std::vector<std::pair<int, int>> path_;
  std::shared_ptr<FemSmoother> fem_smoother_;
  ProjectConfig config_;

  std::vector<double> result_x_list_, result_y_list_;
};
