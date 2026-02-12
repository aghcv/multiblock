#include "fastvessels/config.hpp"

#include "fastvessels/common.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>

namespace fastvessels {

const std::vector<FeatureSpec>& GetFeatureSpecs() {
	static const std::vector<FeatureSpec> specs = {
		{"default", false, false},
		{"cpu_accelerated", true, false},
		{"gpu_accelerated", false, true},
		{"cpu_gpu_accelerated", true, true},
	};
	return specs;
}

static const FeatureSpec& DefaultFeatureSpec() {
	static const FeatureSpec fallback{"default", false, false};
	return fallback;
}

const FeatureSpec& ResolveFeatureSpec(const std::string& name) {
	std::string lowered = ToLower(name);
	for (const auto& spec : GetFeatureSpecs()) {
		if (spec.name == lowered) {
			return spec;
		}
	}
	return DefaultFeatureSpec();
}

std::string JoinFeatureNames() {
	std::string joined;
	for (const auto& spec : GetFeatureSpecs()) {
		if (!joined.empty()) {
			joined += ", ";
		}
		joined += spec.name;
	}
	return joined;
}

std::unordered_map<std::string, std::string> ReadConfigFile(const std::string& path) {
	std::unordered_map<std::string, std::string> values;
	std::ifstream in(path);
	if (!in) {
		return values;
	}
	std::string line;
	while (std::getline(in, line)) {
		line = Trim(line);
		if (line.empty()) {
			continue;
		}
		if (line.rfind("#", 0) == 0 || line.rfind("//", 0) == 0) {
			continue;
		}
		auto eq = line.find('=');
		if (eq == std::string::npos) {
			continue;
		}
		std::string key = Trim(line.substr(0, eq));
		std::string value = Trim(line.substr(eq + 1));
		values[ToLower(key)] = value;
	}
	return values;
}

void ApplyConfigOverrides(SolverConfig& config, const std::unordered_map<std::string, std::string>& values) {
	auto it = values.find("mode");
	if (it != values.end() && !it->second.empty()) {
		config.mode = ToLower(it->second);
	}
	it = values.find("feature");
	if (it != values.end() && !it->second.empty()) {
		config.feature = ToLower(it->second);
		if (config.feature == "alt_hello" || config.feature == "hello" || config.feature == "hello_alt") {
			config.feature = "cpu_accelerated";
		}
	}
	it = values.find("results_output");
	if (it != values.end() && !it->second.empty()) {
		config.results_output = it->second;
	}
	it = values.find("obj_path");
	if (it != values.end() && !it->second.empty()) {
		config.obj_path = it->second;
	}
	it = values.find("default_feature");
	if (it != values.end() && !it->second.empty()) {
		config.default_feature = ToLower(it->second);
	}
	it = values.find("default_cpu_threads");
	if (it != values.end() && !it->second.empty()) {
		config.default_cpu_threads = static_cast<int>(ParseInt(it->second, config.default_cpu_threads));
	}
	it = values.find("default_num_gpus");
	if (it != values.end() && !it->second.empty()) {
		config.default_num_gpus = static_cast<int>(ParseInt(it->second, config.default_num_gpus));
	}
	it = values.find("use_mpi");
	if (it != values.end() && !it->second.empty()) {
		config.use_mpi = ToLower(it->second);
	}

	it = values.find("cpu_threads");
	if (it != values.end() && !it->second.empty()) {
		config.cpu_threads = static_cast<int>(ParseInt(it->second, config.cpu_threads));
	}
	it = values.find("min_cpu_for_gpu");
	if (it != values.end() && !it->second.empty()) {
		config.min_cpu_for_gpu = static_cast<int>(ParseInt(it->second, config.min_cpu_for_gpu));
	}
	it = values.find("num_gpus");
	if (it != values.end() && !it->second.empty()) {
		config.num_gpus = static_cast<int>(ParseInt(it->second, config.num_gpus));
	}
	it = values.find("gpu_work_fraction");
	if (it != values.end() && !it->second.empty()) {
		config.gpu_work_fraction = ParseDouble(it->second, config.gpu_work_fraction);
	}

	it = values.find("work_items");
	if (it != values.end() && !it->second.empty()) {
		config.work_items = static_cast<std::uint64_t>(std::max(1LL, ParseInt(it->second, static_cast<long long>(config.work_items))));
	}
	it = values.find("inner_iters");
	if (it != values.end() && !it->second.empty()) {
		config.inner_iters = static_cast<std::uint32_t>(std::max(1LL, ParseInt(it->second, static_cast<long long>(config.inner_iters))));
	}
	it = values.find("repetitions");
	if (it != values.end() && !it->second.empty()) {
		config.repetitions = static_cast<int>(std::max(1LL, ParseInt(it->second, config.repetitions)));
	}
	it = values.find("unify_walls");
	if (it != values.end() && !it->second.empty()) {
		config.unify_walls = ParseBool(it->second, config.unify_walls);
	}
	it = values.find("max_centerline_xlets");
	if (it != values.end() && !it->second.empty()) {
		config.max_centerline_xlets = static_cast<int>(std::max(1LL, ParseInt(it->second, config.max_centerline_xlets)));
	}
	it = values.find("voxel_base_resolution");
	if (it != values.end() && !it->second.empty()) {
		config.voxel_base_resolution = static_cast<int>(std::max(4LL, ParseInt(it->second, config.voxel_base_resolution)));
	}
	it = values.find("voxel_max_depth");
	if (it != values.end() && !it->second.empty()) {
		config.voxel_max_depth = static_cast<int>(std::max(0LL, ParseInt(it->second, config.voxel_max_depth)));
	}
	it = values.find("voxel_inside_refine_dist");
	if (it != values.end() && !it->second.empty()) {
		config.voxel_inside_refine_dist = static_cast<int>(std::max(1LL, ParseInt(it->second, config.voxel_inside_refine_dist)));
	}
}

SolverConfig LoadConfig(const std::string& path, const SolverConfig& defaults) {
	SolverConfig config = defaults;
	auto values = ReadConfigFile(path);
	ApplyConfigOverrides(config, values);
	config.config_path = path;
	if (config.results_output.empty()) {
		config.results_output = defaults.results_output;
	}
	return config;
}

void WriteConfigFile(const std::string& path,
	const SolverConfig& config,
	const CpuInfo& cpu,
	const GpuInfo& gpu,
	const MpiInfo& mpi,
	const std::string& title) {
	std::filesystem::path outPath(path);
	if (outPath.has_parent_path()) {
		std::filesystem::create_directories(outPath.parent_path());
	}
	std::ofstream out(path);
	if (!out) {
		return;
	}

	out << "# " << title << "\n";
	out << "# Summary: minimal, user-facing configuration\n";
	out << "configuration=cpu_logical=" << cpu.logical;
	if (cpu.physical > 0) {
		out << ", cpu_physical=" << cpu.physical;
	}
	out << ", gpus=" << gpu.count;
	if (gpu.cores > 0) {
		out << ", gpu_cores=" << gpu.cores;
	}
	out << "\n\n";

	out << "mode=" << config.mode << "\n";
	out << "feature=" << config.feature << "\n";
	out << "use_mpi=" << config.use_mpi << "\n";
	out << "results_output=" << config.results_output << "\n\n";
	out << "obj_path=" << config.obj_path << "\n\n";
	out << "default_feature=" << config.default_feature << "\n";
	out << "default_cpu_threads=" << config.default_cpu_threads << "\n";
	out << "default_num_gpus=" << config.default_num_gpus << "\n\n";

	out << "work_items=" << config.work_items << "\n";
	out << "inner_iters=" << config.inner_iters << "\n";
	out << "repetitions=" << config.repetitions << "\n\n";
	out << "unify_walls=" << (config.unify_walls ? "true" : "false") << "\n\n";
	out << "max_centerline_xlets=" << config.max_centerline_xlets << "\n\n";
	out << "voxel_base_resolution=" << config.voxel_base_resolution << "\n";
	out << "voxel_max_depth=" << config.voxel_max_depth << "\n";
	out << "voxel_inside_refine_dist=" << config.voxel_inside_refine_dist << "\n\n";

	out << "cpu_threads=" << config.cpu_threads << "\n";
	out << "min_cpu_for_gpu=" << config.min_cpu_for_gpu << "\n";
	out << "num_gpus=" << config.num_gpus << "\n";
	out << "gpu_work_fraction=" << std::fixed << std::setprecision(3) << config.gpu_work_fraction << "\n\n";

	out << "# Examples\n";
	out << "# 1) Preprocess (auto-benchmark & write optimized file)\n";
	out << "#    mode=preprocess\n";
	out << "# 2) CPU-only accelerated run\n";
	out << "#    feature=cpu_accelerated\n";
	out << "#    cpu_threads=" << std::min<unsigned int>(cpu.logical, 8) << "\n";
	out << "# 3) GPU-only run (if available)\n";
	out << "#    feature=gpu_accelerated\n";
	out << "#    num_gpus=1\n";

	if (mpi.active) {
		out << "# MPI detected: size=" << mpi.size << " rank=" << mpi.rank << "\n";
	}
}

bool ShouldUseMpi(const SolverConfig& config, const MpiInfo& mpiEnv) {
	std::string value = ToLower(config.use_mpi);
	if (value == "true" || value == "1" || value == "yes") {
		return true;
	}
	if (value == "false" || value == "0" || value == "no") {
		return false;
	}
	return mpiEnv.active;
}

} // namespace fastvessels
