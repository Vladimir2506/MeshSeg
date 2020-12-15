#include "MinCut.h"
#include "GeoGraph.h"
#include "FuzzyCluster.h"

#include "argparse.hpp"
#include <filesystem>
#include <chrono>

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
		.default_value(0.8f);
	parser.add_argument("--eta")
		.help("Amplification for concaves in Ang. distance.")
		.action([](const std::string& value) { return std::stof(value); })
		.default_value(0.2f);
	parser.add_argument("--eps")
		.help("Difference between the maximum and the second maximum prob.")
		.action([](const std::string& value) { return std::stof(value); })
		.default_value(0.05f);
	parser.add_argument("--output")
		.help("Folder name of output models.")
		.default_value(std::string("results"));
	parser.add_argument("--debug")
		.help("Draw fuzzy results for only binary case.")
		.default_value(false)
		.implicit_value(true);
	parser.add_argument("--benchmark")
		.help("Benchmark timing for every procedure.")
		.default_value(false)
		.implicit_value(true);
	parser.add_argument("--no_cuda")
		.help("Do not use CUDA to accelerate APSP.")
		.default_value(false)
		.implicit_value(true);

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
	auto outputFoldername = parser.get<std::string>("--output");
	auto maxIterations = parser.get<int>("--M");
	auto numCategories = parser.get<int>("--K");
	auto delta = parser.get<float>("--delta");
	auto eta = parser.get<float>("--eta");
	auto eps = parser.get<float>("--eps");
	auto dbg = parser.get<bool>("--debug");
	auto tic = parser.get<bool>("--benchmark");
	auto ncuda = parser.get<bool>("--no_cuda");

	numCategories = dbg ? 2 : numCategories;

	PaletteHSV pl(numCategories);

	MeshType mesh;

	if (!OpenMesh::IO::read_mesh(mesh, inputFilename))
	{
		std::cout << "> Importing mesh failed." << std::endl;
		return 1;
	}
	auto t1 = std::chrono::system_clock::now();
	std::cout << "> Constructing dual graph..." << std::endl;
	GeoGraph gg(mesh, delta, eta);
	std::cout << "> Dual graph constructed." << std::endl;
	std::cout << "> Running APSP on CUDA..." << std::endl;
	auto&& distanceMatrix = gg.RunShortestPath(!ncuda);
	std::cout << "> APSP finished." << std::endl;
	auto t2 = std::chrono::system_clock::now();
	auto dt21 = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	if (tic)
	{
		std::cout << "Elapsed " << dt21 << "ms / " << dt21 << "ms." << std::endl;
	}

	auto t3 = std::chrono::system_clock::now();
	std::cout << "> Initiating representatives..." << std::endl;
	FuzzyCluster ff(maxIterations, numCategories, distanceMatrix);
	auto gVal = ff.GetGValue();
	std::cout << "> Representatives initiated." << std::endl;
	std::cout << "> Iterating fuzzy cluster..." << std::endl;
	ff.RunCluster();
	auto&& probs = ff.GetProbabilities();
	std::cout << "> Cluster finished." << std::endl;
	auto t4 = std::chrono::system_clock::now();
	auto dt43 = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count(), dt41 = std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t1).count();
	if (tic)
	{
		std::cout << "Elapsed " << dt43 << "ms / " << dt41 << "ms." << std::endl;
	}
	
	auto t5 = std::chrono::system_clock::now();
	std::cout << "> Constructing flow graph..." << std::endl;
	MinCut mc(probs, eps, mesh, eta);
	std::cout << "> Flow graph constructed." << std::endl;
	std::cout << "> Running Mincut..." << std::endl;
	mc.Execute();
	auto&& labels = mc.GetLabels();
	std::cout << "> Mincut refinement finished." << std::endl;
	auto t6 = std::chrono::system_clock::now();
	auto dt65 = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t5).count(), dt61 = std::chrono::duration_cast<std::chrono::milliseconds>(t6 - t1).count();
	if (tic)
	{
		std::cout << "Elapsed " << dt65 << "ms / " << dt61 << "ms." << std::endl;
	}

	mesh.request_face_colors();
	for (auto& face : mesh.faces())
	{
		mesh.set_color(face, pl.colors[labels[face.idx()]]);
	}

	if (dbg)
	{
		for (auto& face : mesh.faces())
		{
			if (mc.vcax.contains(face.idx()))
			{
				mesh.set_color(face, ColorType{ 255,63,63 });
			}
			if (mc.vcbx.contains(face.idx()))
			{
				mesh.set_color(face, ColorType{ 63,255,255 });
			}
			if (mc.vcx.contains(face.idx()))
			{
				mesh.set_color(face, ColorType{ 224,224,224 });
			}
		}
	}

	if (!std::filesystem::exists(outputFoldername))
	{
		std::filesystem::create_directories(outputFoldername);
	}

	std::string stem = std::filesystem::path(inputFilename).stem().string();
	std::ostringstream formatter;
	formatter 
		<< "_K" << numCategories 
		<< "_d" << std::setprecision(2) << delta 
		<< "_et" << std::setprecision(2) << eta
		<< "_e" << std::setprecision(2) << eps 
		<< "_G" << std::setprecision(4) << gVal;
	if (dbg)
	{
		formatter << "_dbg";
	}
	std::string outputFilename = outputFoldername + "/" + stem + formatter.str() + ".obj";

	OpenMesh::IO::write_mesh(mesh, outputFilename, OpenMesh::IO::Options::FaceColor);
	
	std::cout << "> Process done." << std::endl;
	std::cout << "> G(" << numCategories << ") = " << std::setprecision(4) << gVal << std::endl;

	return 0;
}