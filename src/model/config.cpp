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

static std::string NormalizeReportLevel(const std::string& value, const std::string& fallback) {
	if (value.empty()) {
		return fallback;
	}
	const std::string lowered = ToLower(value);
	if (lowered == "long" || lowered == "detailed" || lowered == "full") {
		return "long";
	}
	if (lowered == "short" || lowered == "brief" || lowered == "minimal") {
		return "short";
	}
	return fallback;
}

static std::string NormalizeWallDetectionMode(const std::string& value, const std::string& fallback) {
	if (value.empty()) {
		return fallback;
	}
	const std::string lowered = ToLower(value);
	if (lowered == "ranked" || lowered == "rank" || lowered == "area_flat_connect") {
		return "ranked";
	}
	if (lowered == "flatness" || lowered == "flat" || lowered == "angle") {
		return "flatness";
	}
	return fallback;
}

static std::string NormalizeVoxelRefineMode(const std::string& value, const std::string& fallback) {
	if (value.empty()) {
		return fallback;
	}
	const std::string lowered = ToLower(value);
	if (lowered == "linear" || lowered == "uniform" || lowered == "default") {
		return "linear";
	}
	if (lowered == "distance_jump" || lowered == "distance" || lowered == "jump") {
		return "distance_jump";
	}
	return fallback;
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
	it = values.find("flat_angle_rad");
	if (it != values.end() && !it->second.empty()) {
		config.flat_angle_rad = std::max(0.0, ParseDouble(it->second, config.flat_angle_rad));
	}
	it = values.find("wall_detection_mode");
	if (it != values.end() && !it->second.empty()) {
		config.wall_detection_mode = NormalizeWallDetectionMode(it->second, config.wall_detection_mode);
	}
	it = values.find("wall_rank_area_weight");
	if (it != values.end() && !it->second.empty()) {
		config.wall_rank_area_weight = std::max(0.0, ParseDouble(it->second, config.wall_rank_area_weight));
	}
	it = values.find("wall_rank_flatness_weight");
	if (it != values.end() && !it->second.empty()) {
		config.wall_rank_flatness_weight = std::max(0.0, ParseDouble(it->second, config.wall_rank_flatness_weight));
	}
	it = values.find("wall_rank_connect_weight");
	if (it != values.end() && !it->second.empty()) {
		config.wall_rank_connect_weight = std::max(0.0, ParseDouble(it->second, config.wall_rank_connect_weight));
	}
	it = values.find("global_htg");
	if (it != values.end() && !it->second.empty()) {
		config.global_htg = ParseBool(it->second, config.global_htg);
	}
	it = values.find("global_htg_output");
	if (it != values.end() && !it->second.empty()) {
		config.global_htg_output = it->second;
	}
	it = values.find("global_htg_force_single_label");
	if (it != values.end() && !it->second.empty()) {
		config.global_htg_force_single_label = ParseBool(it->second, config.global_htg_force_single_label);
	}
	it = values.find("region_htg");
	if (it != values.end() && !it->second.empty()) {
		config.region_htg = ParseBool(it->second, config.region_htg);
	}
	it = values.find("report_level");
	if (it != values.end() && !it->second.empty()) {
		config.report_level = NormalizeReportLevel(it->second, config.report_level);
	}
	it = values.find("report_table_rows");
	if (it != values.end() && !it->second.empty()) {
		config.report_table_rows = static_cast<int>(std::max(1LL, ParseInt(it->second, config.report_table_rows)));
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
	it = values.find("voxel_refine_mode");
	if (it != values.end() && !it->second.empty()) {
		config.voxel_refine_mode = NormalizeVoxelRefineMode(it->second, config.voxel_refine_mode);
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
	out << "flat_angle_rad=" << std::fixed << std::setprecision(3) << config.flat_angle_rad << "\n";
	out << "wall_detection_mode=" << config.wall_detection_mode << "\n";
	out << "wall_rank_area_weight=" << std::fixed << std::setprecision(2) << config.wall_rank_area_weight << "\n";
	out << "wall_rank_flatness_weight=" << std::fixed << std::setprecision(2) << config.wall_rank_flatness_weight << "\n";
	out << "wall_rank_connect_weight=" << std::fixed << std::setprecision(2) << config.wall_rank_connect_weight << "\n";
	out << "report_level=" << config.report_level << "\n";
	out << "report_table_rows=" << config.report_table_rows << "\n\n";
	out << "global_htg=" << (config.global_htg ? "true" : "false") << "\n";
	out << "global_htg_output=" << config.global_htg_output << "\n";
	out << "global_htg_force_single_label=" << (config.global_htg_force_single_label ? "true" : "false") << "\n\n";
	out << "region_htg=" << (config.region_htg ? "true" : "false") << "\n\n";
	out << "voxel_base_resolution=" << config.voxel_base_resolution << "\n";
	out << "voxel_max_depth=" << config.voxel_max_depth << "\n";
	out << "voxel_inside_refine_dist=" << config.voxel_inside_refine_dist << "\n";
	out << "voxel_refine_mode=" << config.voxel_refine_mode << "\n\n";

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
