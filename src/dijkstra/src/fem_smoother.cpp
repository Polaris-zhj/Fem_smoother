/******************************************************************************
 * Copyright (C) 2024 ZHJ. All Rights Reserved.
 *****************************************************************************/
#include "include/fem_smoother.h"

FemSmoother::FemSmoother(const std::vector<double> &referenceline_x_init,
                         const std::vector<double> &referenceline_y_init,
                         const ProjectConfig &config) {
  referenceline_x_init_ = referenceline_x_init;
  referenceline_y_init_ = referenceline_y_init;
  dispersed_point_size_ = referenceline_x_init.size();
  config_ = config;

  // for (int i = 0; i < dispersed_point_size_; ++i) {
  //   std::cout << "x: " << referenceline_x_init_[i]
  //             << ",y: " << referenceline_y_init_[i] << std::endl;
  // }
}

void FemSmoother::Solve() {

  size_t variable_num = dispersed_point_size_;
  solver_.settings()->setVerbosity(true);
  solver_.settings()->setWarmStart(true);
  solver_.data()->setNumberOfVariables(2 * variable_num);
  solver_.data()->setNumberOfConstraints(2 * variable_num);
  // solver_.settings()->setMaxIteration(1);
  Eigen::SparseMatrix<double> hessian;
  Eigen::VectorXd gradient;
  Eigen::SparseMatrix<double> linearMatrix;
  Eigen::VectorXd lowerBound;
  Eigen::VectorXd upperBound;
  GetHessianMatrix(&hessian);
  GetGradient(&gradient);
  GetConstraintMatrix(&linearMatrix, &lowerBound, &upperBound);

  // Input to solver.
  if (!solver_.data()->setHessianMatrix(hessian)) {
    std::cout << "初始化H矩阵失败" << std::endl;
  }

  if (!solver_.data()->setGradient(gradient)) {
    std::cout << "初始化f矩阵失败" << std::endl;
  }

  if (!solver_.data()->setLinearConstraintsMatrix(linearMatrix)) {
    std::cout << "初始化约束矩阵失败" << std::endl;
  }

  if (!solver_.data()->setLowerBound(lowerBound)) {
    std::cout << "初始化约下边界失败" << std::endl;
  }

  if (!solver_.data()->setUpperBound(upperBound)) {
    std::cout << "初始化约上边界失败" << std::endl;
  }

  // Solve.
  if (!solver_.initSolver()) {
    std::cout << "初始化求解器失败" << std::endl;
    return;
  }
  solver_.solveProblem();
  SmootherResult();
}

void FemSmoother::SmootherResult() {
  const auto &QPSolution{solver_.getSolution()};
  for (int i = 0; i < dispersed_point_size_; ++i) {
    double tem_val_x = QPSolution(2 * i);
    double tem_val_y = QPSolution(2 * i + 1);
    smooth_x_list_.emplace_back(tem_val_x);
    smooth_y_list_.emplace_back(tem_val_y);
    std::cout << "x: " << tem_val_x << ", y: " << tem_val_y << std::endl;
  }
}

void FemSmoother::GetHessianMatrix(Eigen::SparseMatrix<double> *const hessian) {

  Eigen::MatrixXd A = Eigen::MatrixXd::Constant(2 * dispersed_point_size_,
                                                2 * dispersed_point_size_, 0);
  Eigen::MatrixXd A1 = Eigen::MatrixXd::Constant(2 * dispersed_point_size_ - 4,
                                                 2 * dispersed_point_size_, 0);
  Eigen::MatrixXd A2 = Eigen::MatrixXd::Constant(2 * dispersed_point_size_ - 2,
                                                 2 * dispersed_point_size_, 0);
  Eigen::MatrixXd A3 = Eigen::MatrixXd::Constant(2 * dispersed_point_size_,
                                                 2 * dispersed_point_size_, 0);
  for (int i = 0; i < (2 * dispersed_point_size_ - 5); i = i + 2) {
    A1(i, i) = 1;
    A1(i, i + 2) = -2;
    A1(i, i + 4) = 1;
    A1(i + 1, i + 1) = 1;
    A1(i + 1, i + 3) = -2;
    A1(i + 1, i + 5) = 1;
  }

  for (int i = 0; i < (2 * dispersed_point_size_ - 3); i = i + 2) {
    A2(i, i) = 1;
    A2(i, i + 2) = -1;
    A2(i + 1, i + 1) = 1;
    A2(i + 1, i + 3) = -1;
  }

  for (int i = 0; i < dispersed_point_size_; ++i) {
    A3(2 * i, 2 * i) = 1;
    A3(2 * i + 1, 2 * i + 1) = 1;
  }

  A = config_.weight_config.weight_smooth * (A1.transpose() * A1) +
      config_.weight_config.weight_length * (A2.transpose() * A2) +
      config_.weight_config.weight_ref * A3;

  *hessian = A.sparseView();
}

void FemSmoother::GetGradient(Eigen::VectorXd *const gradient) {
  *gradient = Eigen::VectorXd::Constant(2 * dispersed_point_size_, 0);
  for (int i = 0; i < dispersed_point_size_; ++i) {
    (*gradient)(2 * i, 0) = -2 * referenceline_x_init_[i];
    (*gradient)(2 * i + 1, 0) = -2 * referenceline_y_init_[i];
  }
}

void DebugMrtrix(Eigen::MatrixXd temp_Matrix) {
  for (int i = 0; i < temp_Matrix.rows(); ++i) {
    for (int j = 0; j < temp_Matrix.cols(); ++j) {
      std::cout << temp_Matrix(i, j) << ", ";
    }
    std::cout << std::endl;
  }
}

void FemSmoother::GetConstraintMatrix(
    Eigen::SparseMatrix<double> *const linearMatrix, Eigen::VectorXd *const lb,
    Eigen::VectorXd *const ub) {
  *lb = Eigen::VectorXd::Constant(2 * dispersed_point_size_, 1, 0);
  *ub = Eigen::VectorXd::Constant(2 * dispersed_point_size_, 1, 0);
  Eigen::MatrixXd temp_linearMatrix = Eigen::MatrixXd::Zero(
      2 * dispersed_point_size_, 2 * dispersed_point_size_);
  double y_down_value = -1, y_rise_value = 1;
  double x_down_value = 0, x_rise_value = 0;
  std::cout << "gradient: " << std::endl;
  // 起始点约束
  (*lb)(0, 0) = referenceline_x_init_[0];
  (*ub)(0, 0) = referenceline_x_init_[0];

  (*lb)(1, 0) = referenceline_y_init_[0];
  (*ub)(1, 0) = referenceline_y_init_[0];
  // 终点约束
  (*lb)(2 * (dispersed_point_size_ - 1), 0) =
      referenceline_x_init_[dispersed_point_size_ - 1];
  (*ub)(2 * (dispersed_point_size_ - 1), 0) =
      referenceline_x_init_[dispersed_point_size_ - 1];

  (*lb)(2 * (dispersed_point_size_ - 1) + 1, 0) =
      referenceline_y_init_[dispersed_point_size_ - 1];
  (*ub)(2 * (dispersed_point_size_ - 1) + 1, 0) =
      referenceline_y_init_[dispersed_point_size_ - 1];
  // 中间点约束
  for (int i = 0; i < dispersed_point_size_; ++i) {
    if (i != 0 && i != dispersed_point_size_ - 1) {
      (*lb)(2 * i, 0) = referenceline_x_init_[i] + x_down_value;
      (*ub)(2 * i, 0) = referenceline_x_init_[i] + x_rise_value;

      (*lb)(2 * i + 1, 0) = referenceline_y_init_[i] + y_down_value;
      (*ub)(2 * i + 1, 0) = referenceline_y_init_[i] + y_rise_value;
    }
    temp_linearMatrix(2 * i, 2 * i) = 1;
    temp_linearMatrix(2 * i + 1, 2 * i + 1) = 1;
  }
  *linearMatrix = temp_linearMatrix.sparseView();
}
