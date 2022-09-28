// This file is part of fTetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2019 Yixin Hu <yixin.hu@nyu.edu>
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//

#include <CLI/CLI.hpp>

#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

#include <Eigen/Dense>

#include <igl/write_triangle_mesh.h>

#ifdef LIBIGL_WITH_TETGEN
#include <igl/copyleft/tetgen/tetrahedralize.h>
#endif

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/mesh/mesh.h>
#include <bitset>

using namespace Eigen;

#include <geogram/basic/common.h>

//#include "CLI11.hpp" //TODO: just retrieve the header file
#include "ftetwild.hpp"
#include "ftetwild_io.hpp"

struct CLOptions
{
    bool         skip_simplify = false;
    bool         nobinary      = false;
    bool         nocolor       = false;
    bool         export_raw    = false;
    int          boolean_op    = -1;
    std::string  background_mesh;
    unsigned int max_threads = std::numeric_limits<unsigned int>::max();
};

struct InputMesh {
    std::vector<floatTetWild::Vector3> vertices;
    std::vector<Vector3i>              faces;
    std::vector<int>                   tags;
};

int main(int argc, char** argv) {
#ifndef WIN32
    setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);
#endif

    GEO::initialize();

    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("pre");
    GEO::CmdLine::import_arg_group("algo");

    CLOptions          cl_options;
    floatTetWild::Mesh mesh;

    CLI::App command_line {"float-tetwild"};
    command_line
      .add_option("-i,--input",
                  mesh.params.input_path,
                  "Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)")
      ->check(CLI::ExistingFile);
    command_line.add_option("-o,--output",
                            mesh.params.output_path,
                            "Output tetmesh OUTPUT in .msh format. (string, optional, default: "
                            "input_file+postfix+'.msh')");
    command_line
      .add_option("--tag", mesh.params.tag_path, "Tag input faces for Boolean operation.")
      ->check(CLI::ExistingFile);
    std::string csg_file;
    command_line.add_option("--op",
                            cl_options.boolean_op,
                            "Boolean operation: 0: union, 1: intersection, 2: difference.");
    command_line.add_option(
      "-l,--lr",
      mesh.params.ideal_edge_length_rel,
      "ideal_edge_length = diag_of_bbox * L. (double, optional, default: 0.05)");
    command_line.add_option("-e,--epsr",
                            mesh.params.eps_rel,
                            "epsilon = diag_of_bbox * EPS. (double, optional, default: 1e-3)");

    command_line.add_option("--max-its", mesh.params.max_its, "(for debugging usage only)");
    command_line.add_option("--stop-energy",
                            mesh.params.stop_energy,
                            "Stop optimization when max energy is lower than this.");
    command_line.add_option("--stage", mesh.params.stage, "(for debugging usage only)");
    command_line.add_option("--stop-p", mesh.params.stop_p, "(for debugging usage only)");
    command_line.add_option("--postfix", mesh.params.postfix, "(for debugging usage only)");
    command_line.add_flag("-q,--is-quiet", mesh.params.is_quiet, "Mute console output. (optional)");
    command_line.add_flag("--skip-simplify", cl_options.skip_simplify, "skip preprocessing.");
    command_line.add_flag("--no-binary", cl_options.nobinary, "export meshes as ascii");
    command_line.add_flag("--no-color", cl_options.nocolor, "don't export color");
    command_line.add_flag(
      "--not-sort-input", mesh.params.not_sort_input, "(for debugging usage only)");
    command_line.add_flag("--correct-surface-orientation",
                          mesh.params.correct_surface_orientation,
                          "(for debugging usage only)");
    command_line.add_flag(
      "--smooth-open-boundary", mesh.params.smooth_open_boundary, "Smooth the open boundary.");
    command_line.add_flag("--export-raw", cl_options.export_raw, "Export raw output.");
    command_line.add_flag(
      "--manifold-surface", mesh.params.manifold_surface, "Force the output to be manifold.");
    command_line.add_flag(
      "--coarsen", mesh.params.coarsen, "Coarsen the output as much as possible.");
    command_line.add_option("--csg", csg_file, "json file containg a csg tree")
      ->check(CLI::ExistingFile);
    command_line.add_flag("--disable-filtering",
                          mesh.params.disable_filtering,
                          "Disable filtering out outside elements.");
    command_line.add_flag(
      "--use-floodfill", mesh.params.use_floodfill, "Use flood-fill to extract interior volume.");
    command_line.add_flag(
      "--use-general-wn", mesh.params.use_general_wn, "Use general winding number.");
    command_line.add_flag(
      "--use-input-for-wn", mesh.params.use_input_for_wn, "Use input surface for winding number.");
    command_line
      .add_option(
        "--bg-mesh", cl_options.background_mesh, "Background mesh for sizing field (.msh file).")
      ->check(CLI::ExistingFile);

#ifdef NEW_ENVELOPE
    std::string epsr_tags;
    command_line
      .add_option("--epsr-tags", epsr_tags, "List of envelope size for each input faces.")
      ->check(CLI::ExistingFile);
#endif

#ifdef LIBIGL_WITH_TETGEN
    command_line.add_flag("--tetgen", run_tet_gen, "run tetgen too. (optional)");
#endif

#ifdef FLOAT_TETWILD_USE_TBB
    command_line.add_option("--max-threads", cl_options.max_threads, "Maximum number of threads used");
#endif

    try {
        command_line.parse(argc, argv);
    }
    catch (const CLI::ParseError& e) {
        return command_line.exit(e);
    }

#ifdef FLOAT_TETWILD_USE_TBB
    const size_t MB          = 1024 * 1024;
    const size_t stack_size  = 64 * MB;
    unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
    num_threads              = std::min(cl_options.max_threads, num_threads);
    mesh.params.num_threads       = num_threads;
    std::cout << "TBB threads " << num_threads << std::endl;
    tbb::task_scheduler_init scheduler(num_threads, stack_size);
#endif

    // Sanitize output_path (generate if not specified)
    if (mesh.params.output_path.empty())
        mesh.params.output_path = mesh.params.input_path;
    std::string output_mesh_name;
    if (mesh.params.output_path.size() > 3 &&
        mesh.params.output_path.substr(mesh.params.output_path.size() - 3, mesh.params.output_path.size()) ==
          "msh")
        output_mesh_name = mesh.params.output_path;
    else if (mesh.params.output_path.size() > 4 &&
             mesh.params.output_path.substr(mesh.params.output_path.size() - 4, mesh.params.output_path.size()) ==
               "mesh")
        output_mesh_name = mesh.params.output_path;
    else
        output_mesh_name = mesh.params.output_path + "_" + mesh.params.postfix + ".msh";

    // Setup the background mesh
    Eigen::VectorXd V_in;
    Eigen::VectorXi T_in;
    Eigen::VectorXd values;
    if (!cl_options.background_mesh.empty()) {
        PyMesh::MshLoader mshLoader(cl_options.background_mesh);
        V_in   = mshLoader.get_nodes();
        T_in   = mshLoader.get_elements();
        values = mshLoader.get_node_field("values");
    }

    if (V_in.rows() != 0 && T_in.rows() != 0 && values.rows() != 0) {
        mesh.params.apply_sizing_field = true;

        mesh.params.V_sizing_field      = V_in;
        mesh.params.T_sizing_field      = T_in;
        mesh.params.values_sizing_field = values;
    }

    // Input mesh
    InputMesh inputMesh;
    if (!mesh.params.tag_path.empty()) {
        inputMesh.tags.reserve(inputMesh.faces.size());
        std::string   line;
        std::ifstream fin(mesh.params.tag_path);
        if (fin.is_open()) {
            while (getline(fin, line)) {
                inputMesh.tags.push_back(std::stoi(line));
            }
            fin.close();
        }
    }
#ifdef NEW_ENVELOPE
    if (!epsr_tags.empty()) {
        std::ifstream fin(epsr_tags);
        std::string   line;
        while (std::getline(fin, line)) {
            mesh.params.input_epsr_tags.push_back(std::stod(line));
        }
        fin.close();
    }
#endif

    GEO::Mesh                sf_mesh;
    nlohmann::json           tree_with_ids;
    std::vector<std::string> meshes;
    if (!csg_file.empty()) {
        nlohmann::json csg_tree = nlohmann::json({});
        std::ifstream  file(csg_file);

        if (file.is_open())
            file >> csg_tree;
        else {
            return -1;
        }
        file.close();

        floatTetWild::CSGTreeParser::get_meshes(csg_tree, meshes, tree_with_ids);

        if (!floatTetWild::CSGTreeParser::load_and_merge(
              meshes, inputMesh.vertices, inputMesh.faces, sf_mesh, inputMesh.tags))
            return -1;
    }
    else {
        
#ifdef NEW_ENVELOPE
        if (!MeshIO::load_mesh(params.input_path,
                               input_vertices,
                               input_faces,
                               sf_mesh,
                               input_tags,
                               mesh.params.input_epsr_tags)) {
#else
        if (!floatTetWild::MeshIO::load_mesh(
              mesh.params.input_path, inputMesh.vertices, inputMesh.faces, sf_mesh, inputMesh.tags)) {
#endif
            
            floatTetWild::MeshIO::write_mesh(output_mesh_name, mesh, false);
            return -1;
        }
        else if (inputMesh.vertices.empty() || inputMesh.faces.empty()) {
            floatTetWild::MeshIO::write_mesh(output_mesh_name, mesh, false);
            return -1;
        }

        if (inputMesh.tags.size() != inputMesh.faces.size()) {
            inputMesh.tags.resize(inputMesh.faces.size());
            std::fill(inputMesh.tags.begin(), inputMesh.tags.end(), 0);
        }
    }
    floatTetWild::AABBWrapper tree(sf_mesh);
    if (!mesh.params.init(tree.get_sf_diag())) {
        return -1;
    }

#ifdef NEW_ENVELOPE
    if (!epsr_tags.empty())
        tree.init_sf_tree(
          input_vertices, input_faces, mesh.params.input_epsr_tags, mesh.params.bbox_diag_length);
    else
        tree.init_sf_tree(input_vertices, input_faces, mesh.params.eps);
#endif

#ifdef LIBIGL_WITH_TETGEN
    if (run_tet_gen) {
        Eigen::MatrixXd tetgen_pts(input_vertices.size(), 3);
        Eigen::MatrixXi tetgen_faces(input_faces.size(), 3);

        for (size_t i = 0; i < input_vertices.size(); ++i) {
            tetgen_pts.row(i) = input_vertices[i].cast<double>();
        }

        for (size_t i = 0; i < input_faces.size(); ++i) {
            tetgen_faces.row(i) = input_faces[i];
        }

        std::stringstream buf;
        buf.precision(100);
        buf.setf(std::ios::fixed, std::ios::floatfield);
        buf << "Qpq2.0a"
            << mesh.params.ideal_edge_length * mesh.params.ideal_edge_length * mesh.params.ideal_edge_length *
                 sqrt(2.) / 12.;

        Eigen::MatrixXi tetgen_generated_tets;
        Eigen::MatrixXd tetgen_generated_points;
        Eigen::MatrixXi tetgen_generated_faces;

        timer.start();
        igl::copyleft::tetgen::tetrahedralize(tetgen_pts,
                                              tetgen_faces,
                                              buf.str(),
                                              tetgen_generated_points,
                                              tetgen_generated_tets,
                                              tetgen_generated_faces);
        timer.stop();
        logger().info("Tetgen time {}s", timer.getElapsedTimeInSec());
        stats().record(StateInfo::tetgen_id,
                       timer.getElapsedTimeInSec(),
                       tetgen_generated_points.rows(),
                       tetgen_generated_tets.rows(),
                       0,
                       0);
    }
#endif

    // PREPROCESSING
    simplify(inputMesh.vertices, inputMesh.faces, inputMesh.tags, tree, mesh.params, cl_options.skip_simplify);

    tree.init_b_mesh_and_tree(inputMesh.vertices, inputMesh.faces, mesh);
    if (mesh.params.log_level <= 1)
        floatTetWild::output_component(inputMesh.vertices, inputMesh.faces, inputMesh.tags);
    std::vector<bool> is_face_inserted(inputMesh.faces.size(), false);
    floatTetWild::FloatTetDelaunay::tetrahedralize(inputMesh.vertices, inputMesh.faces, tree, mesh, is_face_inserted);
    insert_triangles(inputMesh.vertices, inputMesh.faces, inputMesh.tags, mesh, is_face_inserted, tree, false);
    optimization(
      inputMesh.vertices, inputMesh.faces, inputMesh.tags, is_face_inserted, mesh, tree, {{1, 1, 1, 1}});
    correct_tracked_surface_orientation(mesh, tree);
    if (cl_options.export_raw) {
        Eigen::Matrix<floatTetWild::Scalar, Eigen::Dynamic, 3> Vt;
        Eigen::Matrix<int, Eigen::Dynamic, 3>    Ft;

        if (!csg_file.empty()) {
            int max_id = floatTetWild::CSGTreeParser::get_max_id(tree_with_ids);

            for (int i = 0; i <= max_id; ++i) {
                get_tracked_surface(mesh, Vt, Ft, i);
                igl::write_triangle_mesh(
                  mesh.params.output_path + "_" + mesh.params.postfix + "_" + std::to_string(i) + "_all.obj",
                  Vt,
                  Ft);
            }
        }
        else {
            get_tracked_surface(mesh, Vt, Ft);
            igl::write_triangle_mesh(
              mesh.params.output_path + "_" + mesh.params.postfix + "_all.obj", Vt, Ft);
        }
        floatTetWild::MeshIO::write_mesh(mesh.params.output_path + "_" + mesh.params.postfix + "_all.msh", mesh, false);
    }
    if (!csg_file.empty())
        boolean_operation(mesh, tree_with_ids, meshes);
    else if (cl_options.boolean_op >= 0)
        boolean_operation(mesh, cl_options.boolean_op);
    else {
        if (mesh.params.smooth_open_boundary) {
            smooth_open_boundary(mesh, tree);
            for (auto& t : mesh.tets) {
                if (t.is_outside)
                    t.is_removed = true;
            }
        }
        else {
            if (!mesh.params.disable_filtering) {
                if (mesh.params.use_floodfill) {
                    filter_outside_floodfill(mesh);
                }
                else if (mesh.params.use_input_for_wn) {
                    filter_outside(mesh, inputMesh.vertices, inputMesh.faces);
                }
                else
                    filter_outside(mesh);
            }
        }
    }
    Eigen::MatrixXd V_sf;
    Eigen::MatrixXi F_sf;
    if (mesh.params.manifold_surface) {
        manifold_surface(mesh, V_sf, F_sf);
    }
    else {
        get_surface(mesh, V_sf, F_sf);
    }
    std::vector<floatTetWild::Scalar> colors;
    if (!cl_options.nocolor) {
        colors.resize(mesh.tets.size(), -1);
        for (int i = 0; i < mesh.tets.size(); i++) {
            if (mesh.tets[i].is_removed)
                continue;
            colors[i] = mesh.tets[i].quality;
        }
    }
    floatTetWild::MeshIO::write_mesh(output_mesh_name, mesh, false, colors, !cl_options.nobinary, !csg_file.empty());
    igl::write_triangle_mesh(mesh.params.output_path + "_" + mesh.params.postfix + "_sf.obj", V_sf, F_sf);
    return 0;
}