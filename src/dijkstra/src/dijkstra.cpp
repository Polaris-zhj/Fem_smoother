/******************************************************************************
 * Copyright (C) 2024 ZHJ. All Rights Reserved.
 *****************************************************************************/
#include "dijkstra.h"
#include "Eigen/src/Core/GlobalFunctions.h"
#include "include/fem_smoother.h"
#include "matplotlibcpp.h"
#include "tinysplinecpp.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <math.h>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <yaml-cpp/node/node.h>
#include <yaml-cpp/yaml.h>
namespace plt = matplotlibcpp;

#define kMathMinVal 1e-4

Dijkstra::Node::Node(const int &temp_x_index, const int &temp_y_index,
                     const double &temp_cost, const int &temp_parent_index) {
  x_index = temp_x_index;
  y_index = temp_y_index;
  cost = temp_cost;
  parent_index = temp_parent_index;
}

Dijkstra::Dijkstra(const double &robot_radius, const double &resolution) {
  robot_radius_ = robot_radius;
  resolution_ = resolution;
}

void Dijkstra::Init(const std::vector<int> &obstacle_x_list,
                    const std::vector<int> &obstacle_y_list) {
  CreateSearchMap(obstacle_x_list, obstacle_y_list);
  NextStepMotions();
  InitParam();
}

void Dijkstra::InitParam() {
  std::string project_path = (std::string)PROJECT_PATH;
  std::cout << "project_path: " << project_path << std::endl;
  std::string fix_path = project_path + "/config/static_config/config.yaml";
  YAML::Node parameters = YAML::LoadFile(fix_path);
  config_.is_use_fem_smooth = parameters["is_use_fem_smooth"].as<bool>();
  config_.is_use_tiny_spline_smooth =
      parameters["is_use_tiny_spline_smooth"].as<bool>();
  config_.weight_config.weight_smooth = parameters["weight_smooth"].as<double>();
  config_.weight_config.weight_length = parameters["weight_length"].as<double>();
  config_.weight_config.weight_ref = parameters["weight_ref"].as<double>();
}

void Dijkstra::Planning(const int &start_x, const int &start_y,
                        const int &goal_x, const int &goal_y) {
  Node *start_node = new Node(CalcXyIndex(start_x, min_x_),
                              CalcXyIndex(start_y, min_y_), 0.0, -1);
  Node *goal_node = new Node(CalcXyIndex(goal_x, min_x_),
                             CalcXyIndex(goal_y, min_y_), 0.0, -1);
  std::map<int, Node *> open_set, close_set;
  open_set[CalcIndex(*start_node)] = start_node;
  std::cout << "start planning " << std::endl;
  std::cout << "start node info: " << start_node->Debug() << std::endl;
  std::cout << "CalcIndex(*start_node): " << CalcIndex(*start_node)
            << std::endl;

  while (true) {
    int index;
    double cost = std::numeric_limits<double>::max();
    for (const auto &iter : open_set) {
      if (iter.second->cost < cost) {
        cost = iter.second->cost;
        index = iter.first;
      }
    }

    Node *current_node = open_set[index];
    open_set.erase(index);
    close_set[index] = current_node;
    if (abs(current_node->x_index - goal_node->x_index) < kMathMinVal &&
        abs(current_node->y_index - goal_node->y_index) < kMathMinVal) {
      // goal_node->parent_index = index;
      // goal_node->cost = current_node->cost;
      goal_node = current_node;
      std::cout << "goal cost: " << goal_node->cost << std::endl;
      break;
    }

    for (const auto &next_step_motion : next_step_motions_) {
      Node *tem_node =
          new Node(current_node->x_index + next_step_motion[0],
                   current_node->y_index + next_step_motion[1],
                   current_node->cost + next_step_motion[2], index);
      int tem_index = CalcIndex(*tem_node);
      if (VerifyNode(*tem_node))
        continue;
      if (close_set.find(tem_index) != close_set.end())
        continue;
      if (open_set.find(tem_index) != open_set.end()) {
        auto &iter = open_set[tem_index];
        if (tem_node->cost < iter->cost) {
          iter->cost = tem_node->cost;
          iter->parent_index = tem_node->parent_index;
        }
      } else {
        open_set[tem_index] = tem_node;
      }
    }
  }

  std::cout << "planning end" << std::endl;
  std::vector<int> path_x_list, path_y_list;
  std::vector<double> tem_path_x_list, tem_path_y_list;

  RetrivePath(goal_node, close_set, &path_x_list, &path_y_list);
  for (int i = 0; i < path_x_list.size(); ++i) {
    tem_path_x_list.emplace_back(path_x_list[i] + 0.01);
    tem_path_y_list.emplace_back(path_y_list[i] + 0.01);
    // std::cout << "xx: " << tem_path_x_list[i] << ", yy: " <<
    // tem_path_y_list[i]
    //           << std::endl;
  }
  if (config_.is_use_tiny_spline_smooth) {
    DiscretePointBSplineSmoother(path_x_list, path_y_list);

    plt::plot(result_x_list_, result_y_list_, ".");
  }
  if (config_.is_use_fem_smooth) {
    DiscretePointFemSmoother(tem_path_x_list, tem_path_y_list);
    plt::plot(fem_smoother_->GetSmootherXlist(),
              fem_smoother_->GetSmootherYlist(), "r");
  }
  plt::plot(tem_path_x_list, tem_path_y_list);
  plt::plot(obstacle_x_list_, obstacle_y_list_, "*");
  plt::show();
}

void Dijkstra::DiscretePointFemSmoother(std::vector<double> path_x_list,
                                        std::vector<double> path_y_list) {
  fem_smoother_ =
      std::make_shared<FemSmoother>(path_x_list, path_y_list, config_);
  fem_smoother_->Solve();
}

void Dijkstra::DiscretePointBSplineSmoother(std::vector<int> &path_x_list,
                                            std::vector<int> &path_y_list) {
  double accumulate_length = 0.0;
  for (int i = 0; i < path_x_list.size() - 1; ++i) {
    double dx = path_x_list[i + 1] - path_x_list[i];
    double dy = path_y_list[i + 1] - path_y_list[i];
    accumulate_length += std::hypot(dx, dy);
  }
  int degree = 3;
  double average_length = accumulate_length / (path_x_list.size() - 1);
  if (average_length > 10)
    degree = 3;
  else if (accumulate_length > 5)
    degree = 4;
  else
    degree = 5;
  tinyspline::BSpline b_spline_raw(path_x_list.size(), 2, 3);
  std::cout << "degree: " << b_spline_raw.degree() << std::endl;
  std::vector<tinyspline::real> ctrlp_raw = b_spline_raw.controlPoints();
  for (size_t i = 0; i != path_x_list.size(); ++i) {
    ctrlp_raw[2 * i] = path_x_list[i];
    ctrlp_raw[2 * i + 1] = path_y_list[i];
  }
  b_spline_raw.setControlPoints(ctrlp_raw);
  double delta_t = 1.0 / (2 * accumulate_length);
  double tem_t = 0.0;
  // std::vector<double> result_x_list, result_y_list;
  while (tem_t <= 1) {
    auto result = b_spline_raw.eval(tem_t).result();
    result_x_list_.emplace_back(result[0]);
    result_y_list_.emplace_back(result[1]);
    tem_t += delta_t;
  }
}

void Dijkstra::RetrivePath(const Node *node, std::map<int, Node *> &closed_set,
                           std::vector<int> *const path_x_list,
                           std::vector<int> *const path_y_list) {
  const Node *temp_node = node;
  while (temp_node->parent_index != -1) {
    int temp_x_index = temp_node->x_index;
    int x = CalcPosition(temp_x_index, min_x_);
    path_x_list->push_back(x);
    int temp_y_index = temp_node->y_index;
    int y = CalcPosition(temp_node->y_index, min_y_);
    path_y_list->push_back(y);
    temp_node = closed_set[temp_node->parent_index];
    path_.push_back(std::make_pair(x, y));
  }
}

int Dijkstra::CalcXyIndex(const int &real_xy, const int &min_xy) {
  return std::round((real_xy - min_xy) / resolution_);
}

int Dijkstra::CalcIndex(const Node &node) {
  return ((node.x_index) * y_width_ + (node.y_index));
}

int Dijkstra::CalcPosition(const int &xy_index, const int &min_xy) {
  return (xy_index * resolution_ + min_xy);
}

void Dijkstra::NextStepMotions() {
  next_step_motions_ = {{1, 0, 1},
                        {-1, 0, 1},
                        {0, 1, 1},
                        {0, -1, 1},
                        {1, 1, std::sqrt(2)},
                        {1, -1, std::sqrt(2)},
                        {-1, 1, std::sqrt(2)},
                        {-1, -1, std::sqrt(2)}};
}

void Dijkstra::CreateSearchMap(std::vector<int> obstacle_x_list,
                               std::vector<int> obstacle_y_list) {
  obstacle_x_list_ = obstacle_x_list;
  obstacle_y_list_ = obstacle_y_list;
  min_x_ =
      round(*std::min_element(obstacle_x_list.begin(), obstacle_x_list.end()));
  max_x_ =
      round(*std::max_element(obstacle_x_list.begin(), obstacle_x_list.end()));
  min_y_ =
      round(*std::min_element(obstacle_y_list.begin(), obstacle_y_list.end()));
  max_y_ =
      round(*std::max_element(obstacle_y_list.begin(), obstacle_y_list.end()));
  std::cout << "min_x: " << min_x_ << ", max_x:" << max_x_ << std::endl;
  std::cout << "min_y: " << min_y_ << ", max_y:" << max_y_ << std::endl;
  min_x_index_ = min_x_ / resolution_;
  max_x_index_ = max_x_ / resolution_;
  min_y_index_ = min_y_ / resolution_;
  max_y_index_ = max_y_ / resolution_;
  x_width_ = round((max_x_ - min_x_) / resolution_ + 1);
  y_width_ = round((max_y_ - min_y_) / resolution_ + 1);
  std::cout << "x_width: " << x_width_ << ", y_width: " << y_width_
            << std::endl;

  obstacle_map_ = std::move(std::vector<std::vector<bool>>(
      x_width_, std::vector<bool>(y_width_, false)));
  for (int i = 0; i < x_width_; ++i) {
    double temp_x = CalcPosition(i, min_x_);
    for (int j = 0; j < y_width_; ++j) {
      double temp_y = CalcPosition(j, min_y_);
      for (int k = 0; k < obstacle_x_list.size(); ++k) {
        double temp_d = std::hypot(temp_x - obstacle_x_list[k],
                                   temp_y - obstacle_y_list[k]);
        if (temp_d < robot_radius_) {
          obstacle_map_[i][j] = true;
          break;
        }
      }
    }
  }
}

bool Dijkstra::VerifyNode(const Node &node) {
  int x = CalcPosition(node.x_index, min_x_);
  int y = CalcPosition(node.y_index, min_y_);
  if (x < min_x_)
    return true;
  if (x > max_x_)
    return true;
  if (y < min_y_)
    return true;
  if (y > max_y_)
    return true;
  if (obstacle_map_[node.x_index][node.y_index])
    return true;
  return false;
}
