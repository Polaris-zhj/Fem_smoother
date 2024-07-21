/******************************************************************************
 * Copyright (C) 2024 ZHJ. All Rights Reserved.
 *****************************************************************************/
#pragma once

struct CostConfig {
  double weight_smooth = 30.4;
  double weight_length = 2.3;
  double weight_ref = 1.2;
};

struct ProjectConfig {
  CostConfig weight_config;

  bool is_use_tiny_spline_smooth;
  bool is_use_fem_smooth;
};