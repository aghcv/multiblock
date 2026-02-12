#pragma once

#include "fastvessels/env_info.hpp"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace fastvessels {

struct SolverConfig {
	std::string mode = "run"; // run|preprocess
	std::string feature = "default"; // default|cpu_accelerated|gpu_accelerated|cpu_gpu_accelerated

	std::string config_path = "output/solver_input.txt";
	std::string results_output = "output/results.txt";
	std::string obj_path = "raw/geometry.vtp";

	std::string default_feature = "";
	int default_cpu_threads = 0;
	int default_num_gpus = -1;

	std::string use_mpi = "auto"; // auto|true|false

	int cpu_threads = 0; // 0=auto, else fixed per-rank
	int min_cpu_for_gpu = 1;
	int num_gpus = -1; // -1=auto-detect, 0=force CPU-only
	double gpu_work_fraction = 0.75; // fraction of work assigned to "GPU" workers when enabled

	std::uint64_t work_items = 5'000'000;
	std::uint32_t inner_iters = 64;
	int repetitions = 3;
	bool unify_walls = false;
	int max_centerline_xlets = 25;

	int voxel_base_resolution = 32;
	int voxel_max_depth = 4;
	int voxel_inside_refine_dist = 2;
};

struct FeatureSpec {
	std::string name;
	bool uses_cpu = false;
	bool uses_gpu = false;
};

const std::vector<FeatureSpec>& GetFeatureSpecs();
const FeatureSpec& ResolveFeatureSpec(const std::string& name);
std::string JoinFeatureNames();

std::unordered_map<std::string, std::string> ReadConfigFile(const std::string& path);
void ApplyConfigOverrides(SolverConfig& config, const std::unordered_map<std::string, std::string>& values);

SolverConfig LoadConfig(const std::string& path, const SolverConfig& defaults);
void WriteConfigFile(const std::string& path,
	const SolverConfig& config,
	const CpuInfo& cpu,
	const GpuInfo& gpu,
	const MpiInfo& mpi,
	const std::string& title);

bool ShouldUseMpi(const SolverConfig& config, const MpiInfo& mpiEnv);

} // namespace fastvessels
