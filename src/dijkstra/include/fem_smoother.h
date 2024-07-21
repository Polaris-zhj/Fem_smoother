/******************************************************************************
 * Copyright (C) 2024 ZHJ. All Rights Reserved.
 *****************************************************************************/
#pragma once
#include "OsqpEigen/OsqpEigen.h"
#include "config/config.h"
#include <cstddef>
#include <vector>

class FemSmoother {
public:
  FemSmoother(const std::vector<double> &referenceline_x_init,
              const std::vector<double> &referenceline_y_init,
              const ProjectConfig &config);
  void Solve();
  void GetHessianMatrix(Eigen::SparseMatrix<double> *const hessian);
  void GetGradient(Eigen::VectorXd *const gradient,
                   Eigen::SparseMatrix<double> *const linearMatrix,
                   Eigen::VectorXd *const lb, Eigen::VectorXd *const ub);

  void GetGradient(Eigen::VectorXd *const gradient);

  void GetConstraintMatrix(Eigen::SparseMatrix<double> *const linearMatrix,
                           Eigen::VectorXd *const lb,
                           Eigen::VectorXd *const ub);
  void SmootherResult();

  std::vector<double> GetSmootherXlist() { return smooth_x_list_; }
  std::vector<double> GetSmootherYlist() { return smooth_y_list_; }

private:
  OsqpEigen::Solver solver_;
  std::size_t dispersed_point_size_;
  std::vector<double> referenceline_x_init_;
  std::vector<double> referenceline_y_init_;

  std::vector<double> smooth_x_list_;
  std::vector<double> smooth_y_list_;

  ProjectConfig config_;
};