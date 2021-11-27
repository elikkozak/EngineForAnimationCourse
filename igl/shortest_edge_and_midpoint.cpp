// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "shortest_edge_and_midpoint.h"

#include <Eigen/src/LU/Inverse.h>
#include <Eigen/Dense>

IGL_INLINE void igl::shortest_edge_and_midpoint(
  const int e,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & /*F*/,
  const Eigen::MatrixXi & E,
  const Eigen::VectorXi & /*EMAP*/,
  const Eigen::MatrixXi & /*EF*/,
  const Eigen::MatrixXi & /*EI*/,
  double & cost,
  Eigen::RowVectorXd & p)
{
  cost = (V.row(E(e,0))-V.row(E(e,1))).norm();
  p = 0.5*(V.row(E(e,0))+V.row(E(e,1)));
}

IGL_INLINE void igl::shortest_edge_and_midpoint_quadric_error_metrics(
	const int e,
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& /*F*/,
	const Eigen::MatrixXi& E,
	const Eigen::VectorXi& /*EMAP*/,
	const Eigen::MatrixXi& /*EF*/,
	const Eigen::MatrixXi& /*EI*/,
	const std::vector<Eigen::Matrix4d>& q_vertices,
	double& cost,
	Eigen::RowVectorXd& p)
{
	Eigen::Vector4d vec1(0, 0, 0, 1);
	// vec1 << 0, 0, 0, 1;
	Eigen::Matrix4d q_edge;
	
	q_edge = q_vertices[E(e, 0)] + q_vertices[E(e, 1)];
	
	Eigen::Matrix4d opt_location = Eigen::Matrix4d::Zero();
	opt_location.row(0) << q_edge.row(0).col(0), q_edge.row(0).col(1), q_edge.row(0).col(2), q_edge.row(0).col(3);
	opt_location.row(1) << q_edge.row(0).col(1), q_edge.row(1).col(1), q_edge.row(1).col(2), q_edge.row(1).col(3);
	opt_location.row(2) << q_edge.row(0).col(2), q_edge.row(1).col(2), q_edge.row(2).col(2), q_edge.row(2).col(3);
	opt_location.row(3) << 0, 0, 0, 1;
	// Eigen::Matrix4d opt_location_inverse = ;
	Eigen::Vector4d new_loc=  opt_location.colPivHouseholderQr().solve(vec1);
	//auto v_bla = opt_location.inverse().m_matrix* (vec1);
	cost = new_loc.transpose() * q_edge * new_loc;
	Eigen::RowVector3d new_loc_t;
	new_loc_t << new_loc[0], new_loc[1], new_loc[2];
	p = new_loc_t;
}
