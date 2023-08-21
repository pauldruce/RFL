//
// Created by Paul Druce on 10/12/2022.
//
#define ARMA_USE_HDF5
#include "BarrettGlaser/Action.hpp"
#include "BarrettGlaser/Metropolis.hpp"
#include "DiracOperator.hpp"
#include "GslRng.hpp"
#include "Simulation.hpp"

#include <cstdlib>
#include <filesystem>
#include <iomanip>

using namespace arma;
namespace fs = std::filesystem;

class Recorder {
public:
  void recordEigenvalues(const DiracOperator& dirac, double g2, int diracId) {
    auto type = dirac.getType();
    auto p = type.first;
    auto q = type.second;
    auto matrixSize = dirac.getMatrixDimension();
    struct diracData data = {p, q, matrixSize};

    std::cout << "Generating and saving eigenvalues for "
              << "{p,q} = {" << p << ", " << q << "}"
              << "and N = " << matrixSize << "\n";

    auto eigenvalues = dirac.getEigenvalues();

    auto outputRootPathEnv = std::getenv("RFL_OUTPUT_DIR");
    std::string outputRootPath = outputRootPathEnv ? outputRootPathEnv : "";
    if (outputRootPath.empty()) {
      outputRootPath = "/tmp/RFL";
    }
    if (!fs::exists(outputRootPath)) {
      fs::create_directory(outputRootPath);
    }

    std::string timeString = getCurrentDateTime();
    auto outputDirectoryPath = create_output_dir(data, timeString, outputRootPath);
    auto eigenvaluesOutputPath = create_eigenvalues_dir(data, outputDirectoryPath);
    auto hdf5Filename = create_eigenvalues_hdf_filename(data, g2);
    auto datasetName = create_dirac_dataset(diracId);

    std::string hdf_filepath = eigenvaluesOutputPath + hdf5Filename;
    std::cout << "Saving eigenvalues to " << hdf_filepath << '\n';
    eigenvalues.save(arma::hdf5_name(hdf_filepath, datasetName, arma::hdf5_opts::append));
  }

private:
  struct diracData {
    int p;
    int q;
    int matrixSize;
  };

  std::string getCurrentDateTime() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::stringstream ss;
    ss << std::put_time(&tm, "%Y%m%d%H%M%S");
    return ss.str();
  }

  std::string create_dirac_dataset(int diracId) {
    std::ostringstream ss;
    ss << "/eigenvalues_" << diracId;
    return ss.str();
  }

  std::string create_eigenvalues_hdf_filename(const struct diracData& data, const double g2) {
    std::ostringstream ss;
    ss << "Eigenvalues_"
       << data.p << "_" << data.q
       << "_N_" << data.matrixSize
       << std::fixed << std::setprecision(3) << "_g2_" << -1.0 * g2
       << ".h5";
    auto eigenvalues_hdf_filename = ss.str();
    return eigenvalues_hdf_filename;
  }

  std::string create_eigenvalues_dir(const diracData& data, const std::string& outputDirectoryPath) {
    std::ostringstream ss;
    ss << "eigenvalues_" << data.p << "_" << data.q << "_N_" << data.matrixSize << "/";
    auto eigenvalues_folder = outputDirectoryPath + ss.str();
    if (!fs::exists(eigenvalues_folder)) {
      fs::create_directory(eigenvalues_folder);
    }
    return eigenvalues_folder;
  }
  std::string create_output_dir(const diracData& data, const std::string& timeString, const std::string& outputRootPath) {
    std::ostringstream ss;

    // Create output dir if needed.
    auto output_path = outputRootPath + "/output/";
    if (!fs::exists(output_path)) {
      fs::create_directory(output_path);
    }

    // Create type folder if needed.
    ss.str("");
    ss.clear();
    ss << output_path << "Type" << data.p << data.q << "/";
    std::string diracTypeDir = ss.str();
    if (!fs::exists(diracTypeDir)) {
      fs::create_directory(diracTypeDir);
    }

    // Create folder for simulation data.
    ss.str("");
    ss.clear();
    ss << diracTypeDir
       << "Simulation_" << data.q << "_" << data.q << "_N_" << data.matrixSize
       << "_" << timeString
       << "/";
    output_path = ss.str();
    if (!fs::exists(output_path)) {
      fs::create_directory(output_path);
    }
    return ss.str();
  }
};

int main() {
  double metropolisScale = 0.2;
  int iter = 10;
  auto rng = std::make_unique<GslRng>();

  auto dirac = std::make_unique<DiracOperator>(1, 3, 10);

  auto g2 = -2.7;
  auto g4 = 1.0;
  auto action = std::make_unique<Action>(g2, g4);

  auto metropolis = std::make_unique<Metropolis>(
      std::move(action),
      metropolisScale,
      iter,
      std::move(rng));

  auto simulation = Simulation(std::move(dirac), std::move(metropolis));

  for (int i = 0; i < 10; i++) {
    simulation.run();
    const auto& dirac2 = simulation.getDiracOperator();

    Recorder recorder;
    recorder.recordEigenvalues(dirac2, g2, i);
  }

  return 0;
}
