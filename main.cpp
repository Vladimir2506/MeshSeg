#include "MinCut.h"
#include "GeoGraph.h"
#include "FuzzyCluster.h"

#include "argparse.hpp"
#include <filesystem>

int main(int argc, char* argv[])
{
	argparse::ArgumentParser parser("Mesh Segmentation with K way decomposition");
	
	parser.add_argument("input")
		.help("Path to input model.");
	parser.add_argument("--K")
		.help("Number of components.")
		.action([](const std::string& value) { return std::stoi(value); })
		.default_value(2);
	parser.add_argument("--M")
		.help("Max iterations for fuzzy clustering.")
		.action([](const std::string& value) { return std::stoi(value); })
		.default_value(2);
	parser.add_argument("--delta")
		.help("Balance factor between Geo. and Ang. distances.")
		.action([](const std::string& value) { return std::stof(value); })
		.default_value(0.8);
	parser.add_argument("--eps")
		.help("Difference between the maximum and the second maximum prob.")
		.action([](const std::string& value) { return std::stof(value); })
		.default_value(0.05);

	try 
	{
		parser.parse_args(argc, argv);
	}
	catch (const std::runtime_error& err) 
	{
		std::cerr << err.what() << std::endl;
		std::cerr << parser;
		return 2;
	}

	auto inputFilename = parser.get<std::string>("input");

	auto maxIterations = parser.get<int>("--M");
	auto numCategories = parser.get<int>("--K");
	auto delta = parser.get<float>("--delta");
	auto eps = parser.get<float>("--eps");

	PaletteHSV pl(numCategories);

	MeshType mesh;

	if (!OpenMesh::IO::read_mesh(mesh, inputFilename))
	{
		std::cerr << "> Importing mesh failed." << std::endl;
		return 1;
	}

	std::cout << "> Constructing dual graph..." << std::endl;
	GeoGraph gg(mesh, delta);
	std::cout << "> Dual graph constructed." << std::endl;
	std::cout << "> Running APSP on CUDA..." << std::endl;
	auto&& distanceMatrix = gg.RunShortestPath();
	std::cout << "> APSP finished." << std::endl;	
	FuzzyCluster ff(maxIterations, numCategories, distanceMatrix);
	auto gVal = ff.GetGValue();
	std::cout << "> Iterating fuzzy cluster..." << std::endl;
	ff.RunCluster();
	std::cout << "> Cluster finished." << std::endl;
	auto&& probs = ff.GetProbabilities();
	std::cout << "> Constructing flow graph..." << std::endl;
	MinCut mc(probs, eps, mesh);
	std::cout << "> Flow graph constructed." << std::endl;
	std::cout << "> Running Mincut..." << std::endl;
	mc.Execute();
	auto&& labels = mc.GetLabels();
	std::cout << "> Mincut refinement finished." << std::endl;

	mesh.request_face_colors();
	for (auto& face : mesh.faces())
	{
		mesh.set_color(face, pl.colors[labels[face.idx()]]);
	}

	if (!std::filesystem::exists("results"))
	{
		std::filesystem::create_directories("results");
	}

	std::string stem = std::filesystem::path(inputFilename).stem().string();
	std::ostringstream formatter;
	formatter 
		<< "_K" << numCategories 
		<< "_d" << std::setprecision(2) << delta 
		<< "_e" << std::setprecision(2) << eps 
		<< "_G" << std::setprecision(4) << gVal;
	std::string outputFilename = "results/" + stem + formatter.str() + ".obj";

	OpenMesh::IO::write_mesh(mesh, outputFilename, OpenMesh::IO::Options::FaceColor);
	
	std::cout << "> Process done." << std::endl;
	std::cout << "> G(" << numCategories << ") = " << std::setprecision(4) << gVal << std::endl;

	return 0;
}