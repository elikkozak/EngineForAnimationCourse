// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "Viewer.h"

//#include <chrono>
#include <thread>

#include <Eigen/LU>


#include <cmath>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <cassert>

#include <igl/project.h>
//#include <igl/get_seconds.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/adjacency_list.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/massmatrix.h>
#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>
#include <igl/quat_mult.h>
#include <igl/axis_angle_to_quat.h>
#include <igl/collapse_edge.h>
#include <igl/trackball.h>
#include <igl/two_axis_valuator_fixed_up.h>
#include <igl/snap_to_canonical_view_quat.h>
#include <igl/unproject.h>
#include <igl/serialize.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl\edge_flaps.cpp>

#include "Build/my_collapse_edge.h"
//#include <igl\opengl\ViewerData.cpp>
//#include <igl\WindingNumberAABB.h>
//#include <tutorial\sandBox\sandBox.h>

// Internal global variables used for glfw event handling
//static igl::opengl::glfw::Viewer * __viewer;
static double highdpi = 1;
static double scroll_x = 0;
static double scroll_y = 0;


namespace igl
{
namespace opengl
{
namespace glfw
{

  void Viewer::Init(const std::string config)
  {
	  

  }

  IGL_INLINE Viewer::Viewer():
    data_list(1),
    selected_data_index(0),
    next_data_id(1),
	isPicked(false),
	isActive(false),
	  OV(20),
	   OF(20),

	  // Prepare array-based edge data structures and priority queue
	   EMAP(20),
	   E(20), EF(20), EI(20),
	   Q(20),
	   Qit(20),
	  // If an edge were collapsed, we'd collapse it to these points:
	   C(20),
		q_vertices(20),
	   num_collapsed(20)
  {
    data_list.front().id = 0;

	

    // Temporary variables initialization
   // down = false;
  //  hack_never_moved = true;
    scroll_position = 0.0f;

    // Per face
    data().set_face_based(false);

    
#ifndef IGL_VIEWER_VIEWER_QUIET
    const std::string usage(R"(igl::opengl::glfw::Viewer usage:
  [drag]  Rotate scene
  A,a     Toggle animation (tight draw loop)
  F,f     Toggle face based
  I,i     Toggle invert normals
  L,l     Toggle wireframe
  O,o     Toggle orthographic/perspective projection
  T,t     Toggle filled faces
  [,]     Toggle between cameras
  1,2     Toggle between models
  ;       Toggle vertex labels
  :       Toggle face labels)"
);
    std::cout<<usage<<std::endl;
#endif
  }

  IGL_INLINE Viewer::~Viewer()
  {
  }

  IGL_INLINE bool Viewer::load_mesh_from_file(
      const std::string & mesh_file_name_string)
  {

    // Create new data slot and set to selected
    if(!(data().F.rows() == 0  && data().V.rows() == 0))
    {
      append_mesh();
    }
    data().clear();

    size_t last_dot = mesh_file_name_string.rfind('.');
    if (last_dot == std::string::npos)
    {
      std::cerr<<"Error: No file extension found in "<<
        mesh_file_name_string<<std::endl;
      return false;
    }

    std::string extension = mesh_file_name_string.substr(last_dot+1);

    if (extension == "off" || extension =="OFF")
    {
      Eigen::MatrixXd V;
      Eigen::MatrixXi F;
      if (!igl::readOFF(mesh_file_name_string, V, F))
        return false;
      data().set_mesh(V,F);
    }
    else if (extension == "obj" || extension =="OBJ")
    {
      Eigen::MatrixXd corner_normals;
      Eigen::MatrixXi fNormIndices;

      Eigen::MatrixXd UV_V;
      Eigen::MatrixXi UV_F;
      Eigen::MatrixXd V;
      Eigen::MatrixXi F;

      if (!(
            igl::readOBJ(
              mesh_file_name_string,
              V, UV_V, corner_normals, F, UV_F, fNormIndices)))
      {
        return false;
      }

      data().set_mesh(V,F);
      if (UV_V.rows() > 0)
      {
          data().set_uv(UV_V, UV_F);
      }

    }
    else
    {
      // unrecognized file type
      printf("Error: %s is not a recognized file type.\n",extension.c_str());
      return false;
    }

    data().compute_normals();
    data().uniform_colors(Eigen::Vector3d(51.0/255.0,43.0/255.0,33.3/255.0),
                   Eigen::Vector3d(255.0/255.0,228.0/255.0,58.0/255.0),
                   Eigen::Vector3d(255.0/255.0,235.0/255.0,80.0/255.0));

    // Alec: why?
    if (data().V_uv.rows() == 0)
    {
      data().grid_texture();
    }
    
    OF[selected_data_index] = data().F;
    OV[selected_data_index] = data().V;
    reset();

    
    //for (unsigned int i = 0; i<plugins.size(); ++i)
    //  if (plugins[i]->post_load())
    //    return true;

    return true;
  }

  IGL_INLINE bool Viewer::save_mesh_to_file(
      const std::string & mesh_file_name_string)
  {
    // first try to load it with a plugin
    //for (unsigned int i = 0; i<plugins.size(); ++i)
    //  if (plugins[i]->save(mesh_file_name_string))
    //    return true;

    size_t last_dot = mesh_file_name_string.rfind('.');
    if (last_dot == std::string::npos)
    {
      // No file type determined
      std::cerr<<"Error: No file extension found in "<<
        mesh_file_name_string<<std::endl;
      return false;
    }
    std::string extension = mesh_file_name_string.substr(last_dot+1);
    if (extension == "off" || extension =="OFF")
    {
      return igl::writeOFF(
        mesh_file_name_string,data().V,data().F);
    }
    else if (extension == "obj" || extension =="OBJ")
    {
      Eigen::MatrixXd corner_normals;
      Eigen::MatrixXi fNormIndices;

      Eigen::MatrixXd UV_V;
      Eigen::MatrixXi UV_F;

      return igl::writeOBJ(mesh_file_name_string,
          data().V,
          data().F,
          corner_normals, fNormIndices, UV_V, UV_F);
    }
    else
    {
      // unrecognized file type
      printf("Error: %s is not a recognized file type.\n",extension.c_str());
      return false;
    }
    return true;
  }
 
  IGL_INLINE bool Viewer::load_scene()
  {
    std::string fname = igl::file_dialog_open();
    if(fname.length() == 0)
      return false;
    return load_scene(fname);
  }

  IGL_INLINE bool Viewer::load_scene(std::string fname)
  {
   // igl::deserialize(core(),"Core",fname.c_str());
    igl::deserialize(data(),"Data",fname.c_str());
    return true;
  }

  IGL_INLINE bool Viewer::save_scene()
  {
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
      return false;
    return save_scene(fname);
  }

  IGL_INLINE bool Viewer::save_scene(std::string fname)
  {
    //igl::serialize(core(),"Core",fname.c_str(),true);
    igl::serialize(data(),"Data",fname.c_str());

    return true;
  }

  IGL_INLINE void Viewer::open_dialog_load_mesh()
  {
    std::string fname = igl::file_dialog_open();

    if (fname.length() == 0)
      return;
    
    this->load_mesh_from_file(fname.c_str());
  }

  IGL_INLINE void Viewer::open_dialog_save_mesh()
  {
    std::string fname = igl::file_dialog_save();

    if(fname.length() == 0)
      return;

    this->save_mesh_to_file(fname.c_str());
  }

  IGL_INLINE ViewerData& Viewer::data(int mesh_id /*= -1*/)
  {
    assert(!data_list.empty() && "data_list should never be empty");
    int index;
    if (mesh_id == -1)
      index = selected_data_index;
    else
      index = mesh_index(mesh_id);

    assert((index >= 0 && index < data_list.size()) &&
      "selected_data_index or mesh_id should be in bounds");
    return data_list[index];
  }

  IGL_INLINE const ViewerData& Viewer::data(int mesh_id /*= -1*/) const
  {
    assert(!data_list.empty() && "data_list should never be empty");
    int index;
    if (mesh_id == -1)
      index = selected_data_index;
    else
      index = mesh_index(mesh_id);

    assert((index >= 0 && index < data_list.size()) &&
      "selected_data_index or mesh_id should be in bounds");
    return data_list[index];
  }

  IGL_INLINE int Viewer::append_mesh(bool visible /*= true*/)
  {
    assert(data_list.size() >= 1);

    data_list.emplace_back();
    selected_data_index = data_list.size()-1;
    data_list.back().id = next_data_id++;
    //if (visible)
    //    for (int i = 0; i < core_list.size(); i++)
    //        data_list.back().set_visible(true, core_list[i].id);
    //else
    //    data_list.back().is_visible = 0;
    return data_list.back().id;
  }

  IGL_INLINE bool Viewer::erase_mesh(const size_t index)
  {
    assert((index >= 0 && index < data_list.size()) && "index should be in bounds");
    assert(data_list.size() >= 1);
    if(data_list.size() == 1)
    {
      // Cannot remove last mesh
      return false;
    }
    data_list[index].meshgl.free();
    data_list.erase(data_list.begin() + index);
    if(selected_data_index >= index && selected_data_index > 0)
    {
      selected_data_index--;
    }

    return true;
  }

  IGL_INLINE size_t Viewer::mesh_index(const int id) const {
    for (size_t i = 0; i < data_list.size(); ++i)
    {
      if (data_list[i].id == id)
        return i;
    }
    return 0;
  }

  Eigen::Matrix4d Viewer::CalcParentsTrans(int indx) 
  {
	  Eigen::Matrix4d prevTrans = Eigen::Matrix4d::Identity();

	  for (int i = indx; parents[i] >= 0; i = parents[i])
	  {
		  //std::cout << "parent matrix:\n" << scn->data_list[scn->parents[i]].MakeTrans() << std::endl;
		  prevTrans = data_list[parents[i]].MakeTransd() * prevTrans;
	  }

	  return prevTrans;
  }

  void Viewer::reset()
  {
      data().F = OF[selected_data_index];
      data().V = OV[selected_data_index];
      edge_flaps(data().F, E[selected_data_index], EMAP[selected_data_index], EF[selected_data_index], EI[selected_data_index]);

      Qit[selected_data_index].resize(E[selected_data_index].rows());
      q_vertices[selected_data_index].resize(OV[selected_data_index].rows());
      for (int i = 0; i < OV[selected_data_index].rows(); ++i)
      {
          q_vertices[selected_data_index][i] = Eigen::Matrix4d::Zero(); 
      }
      
      C[selected_data_index].resize(E[selected_data_index].rows(), data().V.cols());
      Eigen::VectorXd costs(E[selected_data_index].rows());
      Q[selected_data_index].clear();

  	 calcVertexCost();

      for (int e = 0; e < E[selected_data_index].rows(); e++)
      {
          double cost = e;
          Eigen::RowVectorXd p(1, 3);
      	//shortest_edge_and_midpoint(e, data().V, data().F, E[selected_data_index], EMAP[selected_data_index], EF[selected_data_index], EI[selected_data_index], cost, p);

          shortest_edge_and_midpoint_quadric_error_metrics(e, data().V, data().F, E[selected_data_index], EMAP[selected_data_index], EF[selected_data_index], EI[selected_data_index],q_vertices[selected_data_index], cost, p);

          C[selected_data_index].row(e) = p;

          Qit[selected_data_index][e] = Q[selected_data_index].insert(std::pair<double, int>(cost, e)).first;

      }
      num_collapsed[selected_data_index] = 0;
      data().set_mesh(data().V, data().F);

  }

  bool Viewer::preDraw()
  {
      if ( !Q[selected_data_index].empty())
      {
          bool something_collapsed = false;
          // collapse edge
          const int max_iter = std::ceil(0.05 * Q[selected_data_index].size());
          for (int j = 0; j < max_iter; j++)
          {
              if (!collapse_edge(
                  shortest_edge_and_midpoint, data().V, data().F, E[selected_data_index], EMAP[selected_data_index], EF[selected_data_index], EI[selected_data_index], Q[selected_data_index], Qit[selected_data_index], C[selected_data_index]))
              {
                  break;
              }
              something_collapsed = true;
              num_collapsed[selected_data_index]++;
          }

          if (something_collapsed)
          {
             
              data().set_mesh(data().V, data().F);
          }
      }
  	  return false;
  }

  bool Viewer::myPreDraw()
  {
      if (!Q[selected_data_index].empty())
      {
          bool something_collapsed = false;
          // collapse edge
          const int max_iter = std::ceil(0.05 * Q[selected_data_index].size());
          for (int j = 0; j < max_iter; j++)
          {
              if (!my_collapse_edge(
                  shortest_edge_and_midpoint_quadric_error_metrics, data().V, data().F, E[selected_data_index], EMAP[selected_data_index], EF[selected_data_index], EI[selected_data_index], q_vertices[selected_data_index], Q[selected_data_index], Qit[selected_data_index], C[selected_data_index]))
              {
                  break;
              }
              something_collapsed = true;
              num_collapsed[selected_data_index]++;
          }
          calcVertexCost();


          if (something_collapsed)
          {

              data().set_mesh(data().V, data().F);
          }
      }
      return false;
  }
  
  void Viewer::calcVertexCost()
  {

      // q_vertices[selected_data_index].clear();
      for (auto j = 0; j < data().F.rows(); j++)
      {
          Eigen::Matrix4d kp_matrix = Eigen::Matrix4d::Zero();
          Eigen::RowVector3d p_vector = data().F_normals.row(j).normalized();
          Eigen::Vector3d point = data().V.row(data().F.row(j)[0]);
          double d = -(p_vector[0] * point[0] + p_vector[1] * point[1] + p_vector[2] * point[2]);
          Eigen::Vector4d plane_vector;
          plane_vector << p_vector[0], p_vector[1], p_vector[2], d;
          kp_matrix += plane_vector * (plane_vector.transpose());
          q_vertices[selected_data_index][data().F.row(j)[0]] += kp_matrix;
          q_vertices[selected_data_index][data().F.row(j)[1]] += kp_matrix;
          q_vertices[selected_data_index][data().F.row(j)[2]] += kp_matrix;
      }

  }
  
 
} // end namespace
} // end namespace
}
