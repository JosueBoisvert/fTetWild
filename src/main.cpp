#include <CLI11.hpp>

#include <floattetwild/AABBWrapper.h>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/Statistics.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/CSGTreeParser.hpp>
#include <floattetwild/Mesh.hpp>
#include <floattetwild/MeshIO.hpp>

#include <Eigen/Dense>
#include <floattetwild/Logger.hpp>

#include <igl/write_triangle_mesh.h>

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/mesh/mesh.h>
#include <bitset>

using namespace floatTetWild;
using namespace Eigen;


#include <geogram/basic/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/numeric.h>
#include <floattetwild/Predicates.hpp>

#include <floattetwild/MshLoader.h>
#include <geogram/mesh/mesh_AABB.h>


int main(int argc, char** argv)
{
    GEO::initialize();

    std::vector<int> indices(20);
    std::iota(std::begin(indices), std::end(indices), 0);
    floatTetWild::Random::shuffle(indices);
    for (int a : indices)
        std::cout << a << " ";
    std::cout << std::endl;

    // Import standard command line arguments, and custom ones
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("pre");
    GEO::CmdLine::import_arg_group("algo");

    bool run_tet_gen   = false;
    bool skip_simplify = false;
    bool nobinary      = false;
    bool nocolor       = false;
    bool export_raw    = false;

    Mesh        mesh;
    Parameters& params = mesh.params;
    int         boolean_op = -1;
    std::string csg_file = "";
    std::string background_mesh = "";

    // Parse command line arguments
    CLI::App command_line {"float-tetwild"};
    command_line.add_option("-i,--input", params.input_path, "Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)")->check(CLI::ExistingFile);
    command_line.add_option("-o,--output", params.output_path, "Output tetmesh OUTPUT in .msh format. (string, optional, default: input_file+postfix+'.msh')");
    command_line.add_option("--tag", params.tag_path, "Tag input faces for Boolean operation.")->check(CLI::ExistingFile);
    command_line.add_option("--op", boolean_op, "Boolean operation: 0: union, 1: intersection, 2: difference.");
    command_line.add_option("-l,--lr", params.ideal_edge_length_rel, "ideal_edge_length = diag_of_bbox * L. (double, optional, default: 0.05)");
    command_line.add_option("-e,--epsr", params.eps_rel, "epsilon = diag_of_bbox * EPS. (double, optional, default: 1e-3)");
    command_line.add_option("--max-its", params.max_its, "(for debugging usage only)");
    command_line.add_option("--stop-energy", params.stop_energy, "Stop optimization when max energy is lower than this.");
    command_line.add_option("--stage", params.stage, "(for debugging usage only)");
    command_line.add_option("--stop-p", params.stop_p, "(for debugging usage only)");
    command_line.add_option("--postfix", params.postfix, "(for debugging usage only)");
    command_line.add_option("--log", params.log_path, "Log info to given file.");
    command_line.add_option("--level", params.log_level, "Log level (0 = most verbose, 6 = off).");
    command_line.add_flag("-q,--is-quiet", params.is_quiet, "Mute console output. (optional)");
    command_line.add_flag("--skip-simplify", skip_simplify, "skip preprocessing.");
    command_line.add_flag("--no-binary", nobinary, "export meshes as ascii");
    command_line.add_flag("--no-color", nocolor, "don't export color");
    command_line.add_flag("--not-sort-input", params.not_sort_input, "(for debugging usage only)");
    command_line.add_flag("--correct-surface-orientation", params.correct_surface_orientation, "(for debugging usage only)");
    command_line.add_option("--envelope-log", params.envelope_log, "(for debugging usage only)");
    command_line.add_flag("--smooth-open-boundary", params.smooth_open_boundary, "Smooth the open boundary.");
    command_line.add_flag("--export-raw", export_raw, "Export raw output.");
    command_line.add_flag("--manifold-surface", params.manifold_surface, "Force the output to be manifold.");
    command_line.add_flag("--coarsen", params.coarsen, "Coarsen the output as much as possible.");
    command_line.add_option("--csg", csg_file, "json file containg a csg tree")->check(CLI::ExistingFile);
    command_line.add_flag("--use-old-energy", floatTetWild::use_old_energy, "(for debugging usage only)");  // tmp
    command_line.add_flag("--disable-filtering", params.disable_filtering, "Disable filtering out outside elements.");
    command_line.add_flag("--use-floodfill", params.use_floodfill, "Use flood-fill to extract interior volume.");
    command_line.add_flag("--use-general-wn", params.use_general_wn, "Use general winding number.");
    command_line.add_flag("--use-input-for-wn", params.use_input_for_wn, "Use input surface for winding number.");
    command_line.add_option("--bg-mesh", background_mesh, "Background mesh for sizing field (.msh file).")->check(CLI::ExistingFile);

    try {
        command_line.parse(argc, argv);
    }
    catch (const CLI::ParseError& e) {
        return command_line.exit(e);
    }



    // Logger::init(!params.is_quiet, params.log_path);
    // params.log_level = std::max(0, std::min(6, params.log_level));
    // spdlog::set_level(static_cast<spdlog::level::level_enum>(params.log_level));
    // spdlog::flush_every(std::chrono::seconds(3));

    // GEO::Logger* geo_logger = GEO::Logger::instance();
    // geo_logger->unregister_all_clients();
    // geo_logger->register_client(new GeoLoggerForward(logger().clone("geogram")));
    // geo_logger->set_pretty(false);

    if (params.output_path.empty())
        params.output_path = params.input_path;
    if (params.log_path.empty())
        params.log_path = params.output_path;

    // simple output file name sanitization
    // std::string output_mesh_name = params.output_path;
    // if (params.output_path.size() > 3 && params.output_path.substr(params.output_path.size() - 3, params.output_path.size()) == "msh")
    //     output_mesh_name = params.output_path;
    // else if (params.output_path.size() > 4 && params.output_path.substr(params.output_path.size() - 4, params.output_path.size()) == "mesh")
    //     output_mesh_name = params.output_path;
    // else
    //     output_mesh_name = params.output_path + "_" + params.postfix + ".msh";
    std::string output_mesh_name = "test_output_mesh.msh";

    // set sizing field
    Eigen::VectorXd V_in;
    Eigen::VectorXi T_in;
    Eigen::VectorXd values;

    // Load background .msh file
    if (!background_mesh.empty()) {
        PyMesh::MshLoader mshLoader(background_mesh);
        V_in   = mshLoader.get_nodes();
        T_in   = mshLoader.get_elements();
        values = mshLoader.get_node_field("values");

        if (V_in.rows() != 0 && T_in.rows() != 0 && values.rows() != 0) {
            params.apply_sizing_field = true;

            params.V_sizing_field      = V_in;
            params.T_sizing_field      = T_in;
            params.values_sizing_field = values;
        }
    }


    // Set input tags
    std::vector<Vector3>  input_vertices;
    std::vector<Vector3i> input_faces;
    std::vector<int>      input_tags;

    if (!params.tag_path.empty()) {
        input_tags.reserve(input_faces.size());
        std::string   line;
        std::ifstream fin(params.tag_path);
        if (fin.is_open()) {
            while (getline(fin, line)) {
                input_tags.push_back(std::stoi(line));
            }
            fin.close();
        }
    }

    // set envelope
    GEO::Mesh                sf_mesh;
    json                     tree_with_ids;
    std::vector<std::string> meshes;
    
    if (!csg_file.empty()) {
        json          csg_tree = json({});
        std::ifstream file(csg_file);

        if (file.is_open())
            file >> csg_tree;
        else {
            return EXIT_FAILURE;
        }
        file.close();

        CSGTreeParser::get_meshes(csg_tree, meshes, tree_with_ids);

        if (!CSGTreeParser::load_and_merge(
              meshes, input_vertices, input_faces, sf_mesh, input_tags))
            return EXIT_FAILURE;
    }
    else {
        // Load the mesh
        if (!MeshIO::load_mesh(params.input_path, input_vertices, input_faces, sf_mesh, input_tags)) {
            MeshIO::write_mesh(output_mesh_name, mesh, false);
            return -1;
        }
        else if (input_vertices.empty() || input_faces.empty()) {
            MeshIO::write_mesh(output_mesh_name, mesh, false);
            return -1;
        }

        // Ensure that the tag vector matches the face vector length. Pad with zeros
        if (input_tags.size() != input_faces.size()) {
            input_tags.resize(input_faces.size());
            std::fill(input_tags.begin(), input_tags.end(), 0);
        }
    }

    AABBWrapper tree(sf_mesh);
    if (!params.init(tree.get_sf_diag())) {
        return -1;
    }

    simplify(input_vertices, input_faces, input_tags, tree, params, skip_simplify);
    tree.init_b_mesh_and_tree(input_vertices, input_faces, mesh);

    if (params.log_level <= 1)
        output_component(input_vertices, input_faces, input_tags);

    std::vector<bool> is_face_inserted(input_faces.size(), false);
    FloatTetDelaunay::tetrahedralize(input_vertices, input_faces, tree, mesh, is_face_inserted);

    insert_triangles(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, false);

    optimization(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, {{1, 1, 1, 1}});

    correct_tracked_surface_orientation(mesh, tree);

    if (export_raw) {
        Eigen::Matrix<Scalar, Eigen::Dynamic, 3> Vt;
        Eigen::Matrix<int, Eigen::Dynamic, 3>    Ft;

        if (!csg_file.empty()) {
            int max_id = CSGTreeParser::get_max_id(tree_with_ids);

            for (int i = 0; i <= max_id; ++i) {
                get_tracked_surface(mesh, Vt, Ft, i);
                igl::write_triangle_mesh(
                  params.output_path + "_" + params.postfix + "_" + std::to_string(i) + "_all.obj",
                  Vt,
                  Ft);
            }
        }
        else {
            get_tracked_surface(mesh, Vt, Ft);
            igl::write_triangle_mesh(
              params.output_path + "_" + params.postfix + "_all.obj", Vt, Ft);
        }
        MeshIO::write_mesh(params.output_path + "_" + params.postfix + "_all.msh", mesh, false);
    }

    if (!csg_file.empty())
        boolean_operation(mesh, tree_with_ids, meshes);
    else if (boolean_op >= 0)
        boolean_operation(mesh, boolean_op);
    else {
        if (params.smooth_open_boundary) {
            smooth_open_boundary(mesh, tree);
            for (auto& t : mesh.tets) {
                if (t.is_outside)
                    t.is_removed = true;
            }
        }
        else {
            if (!params.disable_filtering) {
                if (params.use_floodfill) {
                    filter_outside_floodfill(mesh);
                }
                else if (params.use_input_for_wn) {
                    filter_outside(mesh, input_vertices, input_faces);
                }
                else
                    filter_outside(mesh);
            }
        }
    }
    Eigen::MatrixXd V_sf;
    Eigen::MatrixXi F_sf;
    if (params.manifold_surface) {
        manifold_surface(mesh, V_sf, F_sf);
    }
    else {
        get_surface(mesh, V_sf, F_sf);
    }

    // fortest
    std::vector<Scalar> colors;
    if (!nocolor) {
        colors.resize(mesh.tets.size(), -1);
        for (int i = 0; i < mesh.tets.size(); i++) {
            if (mesh.tets[i].is_removed)
                continue;
            colors[i] = mesh.tets[i].quality;
        }
    }
    // fortest
    MeshIO::write_mesh(output_mesh_name, mesh, false, colors, !nobinary, !csg_file.empty());
    igl::write_triangle_mesh(params.output_path + "_" + params.postfix + "_sf.obj", V_sf, F_sf);

    return 0;
}