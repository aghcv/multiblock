#pragma once

#include <string>

#include <vtkMultiBlockDataSet.h>
#include <vtkSmartPointer.h>

namespace fastvessels {

vtkSmartPointer<vtkMultiBlockDataSet> ReadGeometry_AsMultiBlock(const std::string& path);

vtkSmartPointer<vtkMultiBlockDataSet> ReadStlDirectory_AsMultiBlock(const std::string& dirPath);

vtkSmartPointer<vtkMultiBlockDataSet> BuildRegionsFromSurfaceBlocks(vtkMultiBlockDataSet* blocks);

struct ObjPipelineStats {
	int block_count = 0;
	int closed_surfaces = 0;
	int blocks_with_openings = 0;
	long long boundary_edge_cells = 0;
};

ObjPipelineStats AnalyzeClosedSurfaces(vtkMultiBlockDataSet* blocks, int cpuThreads);

void AnalyzeRegionGroupSurfaces(vtkMultiBlockDataSet* regions,
	int maxRegionThreads,
	bool unifyWalls,
	double flatAngleRad,
	const std::string& wallDetectionMode,
	double wallRankAreaWeight,
	double wallRankFlatnessWeight,
	double wallRankConnectWeight,
	const std::string& reportLevel,
	int reportTableRows);

void CenterlineBase(vtkMultiBlockDataSet* regions,
	const std::string& reportLevel,
	int reportTableRows);

void VoxelizeRegionsBase(vtkMultiBlockDataSet* regions, int baseResolution, int maxDepth, int insideRefineDist);
void VoxelizeRegionsBase(vtkMultiBlockDataSet* regions,
	int baseResolution,
	int maxDepth,
	int insideRefineDist,
	bool writeRegion,
	bool writeGlobal,
	const std::string& globalOutputPath,
	bool forceSingleLabel,
	const std::string& refineMode);

vtkSmartPointer<vtkMultiBlockDataSet> BuildRegionSurfaceHierarchy(
	vtkMultiBlockDataSet* inputMb,
	const std::string& groupArrayName = "GroupId",
	const std::string& regionArrayName = "RegionId",
	bool addOriginalBlockId = true);

} // namespace fastvessels
