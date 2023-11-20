//
// Created by Paul Druce on 30/12/2023.
//

#include "EigenvalueRecorder.hpp"
#include <filesystem>
#include <iomanip>

namespace fs = std::filesystem;

struct diracData {
  int p;
  int q;
  int matrixSize;
};

static std::string create_eigenvalues_dataset(const int diracId) {
  std::ostringstream ss;
  ss << "/eigenvalues_" << diracId;
  return ss.str();
}

static std::string create_simulation_group(const std::string& timeString) {
  std::ostringstream ss;
  ss << "Simulation_" << timeString;
  return ss.str();
}

static std::string create_eigenvalues_hdf_filename(const diracData& data, const double g2) {
  std::ostringstream ss;

  ss << -1 * g2;
  std::string g2_string = ss.str();
  ss.str("");
  ss.clear();

  // Replace decimal with underscore in g2 string.
  std::replace(g2_string.begin(), g2_string.end(), '.', '_');

  ss << "Eigenvalues_"
     << data.p << "_" << data.q
     << "_N_" << data.matrixSize
     << std::fixed << std::setprecision(3) << "_g2_" << g2_string
     << ".h5";
  auto eigenvalues_hdf_filename = ss.str();
  return eigenvalues_hdf_filename;
}

static std::string create_eigenvalues_dir(const diracData& data, const std::string& outputDirectoryPath) {
  std::ostringstream ss;
  ss << "eigenvalues_" << data.p << "_" << data.q << "_N_" << data.matrixSize << "/";
  auto eigenvalues_folder = outputDirectoryPath + ss.str();
  if (!fs::exists(eigenvalues_folder)) {
    fs::create_directory(eigenvalues_folder);
  }
  return eigenvalues_folder;
}

static std::string create_output_dir(const diracData& data, const std::string& outputRootPath) {
  std::ostringstream ss;

  // Create output dir if needed.
  const auto output_path = outputRootPath + "/output/";
  if (!fs::exists(output_path)) {
    fs::create_directory(output_path);
  }

  // Create type folder if needed.
  ss.str("");
  ss.clear();
  ss << output_path << "Type" << data.p << data.q << "/";
  const std::string diracTypeDir = ss.str();
  if (!fs::exists(diracTypeDir)) {
    fs::create_directory(diracTypeDir);
  }
  return ss.str();
}

void EigenvalueRecorder::recordEigenvalues(int diracId) const {
  auto [p, q] = m_dirac.getType();
  auto matrixSize = m_dirac.getMatrixDimension();
  struct diracData data = {p, q, matrixSize};

  std::cout << "Generating and saving eigenvalues for "
            << "{p,q} = {" << p << ", " << q << "} "
            << "and N = " << matrixSize << "\n";

  auto eigenvalues = m_dirac.getEigenvalues();

  auto outputRootPathEnv = std::getenv("RFL_OUTPUT_DIR");
  std::string outputRootPath = outputRootPathEnv ? outputRootPathEnv : "";
  if (outputRootPath.empty()) {
    outputRootPath = "/tmp/RFL";
  }
  if (!fs::exists(outputRootPath)) {
    fs::create_directory(outputRootPath);
  }
  auto outputDirectoryPath = create_output_dir(data, outputRootPath);
  auto eigenvaluesOutputPath = create_eigenvalues_dir(data, outputDirectoryPath);
  auto hdf5Filename = create_eigenvalues_hdf_filename(data, m_g2);
  auto simulation_hdf5_group = create_simulation_group(m_simulationId);
  auto datasetName = create_eigenvalues_dataset(diracId);

  std::string hdf_filepath = eigenvaluesOutputPath + hdf5Filename;
  std::cout << "Saving eigenvalues to " << hdf_filepath << '\n';
  eigenvalues.save(arma::hdf5_name(hdf_filepath, simulation_hdf5_group + "/" + datasetName, arma::hdf5_opts::append));
}
