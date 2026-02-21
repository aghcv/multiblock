#include "fastvessels/solver.hpp"

#include "fastvessels/common.hpp"
#include "fastvessels/io_utils.hpp"
#include "fastvessels/obj_pipeline.hpp"

#include <algorithm>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <limits>
#include <thread>
#include <vector>

namespace fastvessels {

RuntimeConfig ResolveRuntimeConfig(const SolverConfig& config, const CpuInfo& cpu, const GpuInfo& gpu) {
	RuntimeConfig runtime;
	const auto& featureSpec = ResolveFeatureSpec(config.feature);
	int detectedGpus = gpu.count;
	int requestedGpus = config.num_gpus;
	int useGpus = requestedGpus;
	if (useGpus < 0) {
		if (config.default_num_gpus >= 0) {
			useGpus = config.default_num_gpus;
		} else {
			useGpus = detectedGpus;
		}
	}
	useGpus = std::max(0, useGpus);

	int cpuThreads = config.cpu_threads;
	if (cpuThreads <= 0) {
		if (config.default_cpu_threads > 0) {
			cpuThreads = config.default_cpu_threads;
		} else {
			cpuThreads = static_cast<int>(std::max(1u, cpu.logical));
		}
	}
	cpuThreads = std::max(1, cpuThreads);

	const auto& resolvedFeature = (featureSpec.name == "default" && !config.default_feature.empty())
		? ResolveFeatureSpec(config.default_feature)
		: featureSpec;
	bool wantsGpu = resolvedFeature.uses_gpu;
	bool wantsCpuAccel = resolvedFeature.uses_cpu;
	if (resolvedFeature.name == "default") {
		cpuThreads = 1;
		useGpus = 0;
	}
	if (!wantsCpuAccel) {
		cpuThreads = 1;
	}
	if (!wantsGpu) {
		useGpus = 0;
	}
	if (wantsGpu) {
		cpuThreads = std::max(cpuThreads, config.min_cpu_for_gpu);
	}

	runtime.cpu_threads = cpuThreads;
	runtime.gpu_workers = useGpus;
	runtime.resolved_feature = resolvedFeature.name;
	return runtime;
}

static BenchmarkResult RunBenchmarkOnce(const SolverConfig& config, const CpuInfo& cpu, const GpuInfo& gpu, const MpiInfo& mpi) {
	BenchmarkResult result;
	result.feature = config.feature;
	result.work_items = config.work_items;
	result.inner_iters = config.inner_iters;

	RuntimeConfig runtime = ResolveRuntimeConfig(config, cpu, gpu);
	result.cpu_threads = runtime.cpu_threads;
	result.gpu_workers = runtime.gpu_workers;

	auto start = std::chrono::steady_clock::now();

	vtkSmartPointer<vtkMultiBlockDataSet> refined;

	std::cout << "\n=== Reading geometry ===" << std::endl;
	std::cout << "Geometry path: " << config.obj_path << std::endl;
	if (std::filesystem::is_directory(config.obj_path)) {
		std::cout << "Mode: STL directory" << std::endl;
	} else {
		std::cout << "Mode: single file" << std::endl;
	}

	if (std::filesystem::is_directory(config.obj_path)) {
		auto blocks = ReadStlDirectory_AsMultiBlock(config.obj_path);
		refined = BuildRegionsFromSurfaceBlocks(blocks);
	} else {
		auto blocks = ReadGeometry_AsMultiBlock(config.obj_path);
		refined = BuildRegionSurfaceHierarchy(blocks, "GroupId", "RegionId", true);
	}
	AnalyzeRegionGroupSurfaces(refined, runtime.cpu_threads, config.unify_walls,
		config.flat_angle_rad, config.wall_detection_mode,
		config.wall_rank_area_weight, config.wall_rank_flatness_weight,
		config.wall_rank_connect_weight, config.report_level, config.report_table_rows);
	CenterlineBase(refined, config.report_level, config.report_table_rows);
	VoxelizeRegionsBase(refined, config.voxel_base_resolution, config.voxel_max_depth,
		config.voxel_inside_refine_dist, config.region_htg, config.global_htg,
		config.global_htg_output, config.global_htg_force_single_label);

	std::filesystem::create_directories("output");
	WriteMultiBlock(refined, "output/geometry_multiblock.vtm");

	ObjPipelineStats stats = AnalyzeClosedSurfaces(refined, runtime.cpu_threads);

	auto stop = std::chrono::steady_clock::now();
	result.seconds = std::chrono::duration<double>(stop - start).count();
	result.block_count = stats.block_count;
	result.closed_surfaces = stats.closed_surfaces;
	result.blocks_with_openings = stats.blocks_with_openings;
	result.boundary_edge_cells = stats.boundary_edge_cells;
	result.checksum = static_cast<double>(stats.closed_surfaces);
	return result;
}

BenchmarkResult RunBenchmark(const SolverConfig& config, const CpuInfo& cpu, const GpuInfo& gpu, const MpiInfo& mpi) {
	if (ToLower(config.mode) == "run") {
		return RunBenchmarkOnce(config, cpu, gpu, mpi);
	}

	BenchmarkResult best;
	best.seconds = std::numeric_limits<double>::infinity();
	for (int r = 0; r < config.repetitions; ++r) {
		auto res = RunBenchmarkOnce(config, cpu, gpu, mpi);
		if (res.seconds < best.seconds) {
			best = res;
		}
	}
	return best;
}

double ThroughputItemsPerSecond(const BenchmarkResult& r) {
	if (r.seconds <= 0.0) {
		return 0.0;
	}
	if (r.block_count > 0) {
		return static_cast<double>(r.block_count) / r.seconds;
	}
	return static_cast<double>(r.work_items) / r.seconds;
}

} // namespace fastvessels
