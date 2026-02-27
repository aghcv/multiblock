#include "fastvessels/obj_pipeline.hpp"

#include "fastvessels/common.hpp"

#include <vtkAppendPolyData.h>
#include <vtkAbstractPointLocator.h>
#include <vtkCellData.h>
#include <vtkCompositeDataSet.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkHyperTreeGrid.h>
#include <vtkHyperTreeGridNonOrientedGeometryCursor.h>
#include <vtkHyperTreeGridToUnstructuredGrid.h>
#include <vtkIdList.h>
#include <vtkImageData.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkMassProperties.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPolyData.h>
#include <vtkPointSet.h>
#include <vtkStaticPointLocator.h>
#include <vtkSmartPointer.h>
#include <vtkTriangleFilter.h>
#include <vtkUnsignedCharArray.h>
#include <vtkXMLHyperTreeGridWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace fastvessels {

namespace {

std::string GetBlockName(vtkMultiBlockDataSet* mb, unsigned int idx) {
	if (!mb) return "";
	vtkInformation* info = mb->GetMetaData(idx);
	if (!info) return "";
	const char* name = info->Get(vtkCompositeDataSet::NAME());
	return name ? std::string(name) : "";
}

int GetSurfaceTypeField(vtkPolyData* pd) {
	if (!pd) return -1;
	vtkFieldData* field = pd->GetFieldData();
	if (!field) return -1;
	auto* arr = vtkIntArray::SafeDownCast(field->GetAbstractArray("SurfaceType"));
	if (!arr || arr->GetNumberOfTuples() < 1) return -1;
	return arr->GetValue(0);
}

int ParseGroupIdFromName(const std::string& name) {
	if (name.empty()) return -1;
	const std::string prefix = "Group_";
	if (name.rfind(prefix, 0) != 0) return -1;
	const std::string idStr = name.substr(prefix.size());
	if (idStr.empty()) return -1;
	for (char c : idStr) {
		if (c < '0' || c > '9') return -1;
	}
	return std::stoi(idStr);
}

int GetGroupId(vtkPolyData* pd, const std::string& name) {
	const int fromName = ParseGroupIdFromName(name);
	if (fromName >= 0) return fromName;
	if (!pd) return -1;
	vtkDataArray* arr = pd->GetCellData()->GetArray("GroupId");
	if (!arr || arr->GetNumberOfTuples() < 1) return -1;
	return static_cast<int>(arr->GetComponent(0, 0));
}

double ComputeSurfaceArea(vtkPolyData* pd) {
	if (!pd || pd->GetNumberOfCells() == 0) return 0.0;
	auto tri = vtkSmartPointer<vtkTriangleFilter>::New();
	tri->SetInputData(pd);
	tri->Update();

	auto mass = vtkSmartPointer<vtkMassProperties>::New();
	mass->SetInputData(tri->GetOutput());
	mass->Update();

	return mass->GetSurfaceArea();
}

vtkSmartPointer<vtkPolyData> CombineRegionSurfaces(vtkMultiBlockDataSet* regionMb) {
	if (!regionMb) return nullptr;
	auto appender = vtkSmartPointer<vtkAppendPolyData>::New();
	const int groupCount = static_cast<int>(regionMb->GetNumberOfBlocks());
	for (int g = 0; g < groupCount; ++g) {
		auto* pd = vtkPolyData::SafeDownCast(regionMb->GetBlock(static_cast<unsigned int>(g)));
		if (!pd || pd->GetNumberOfCells() == 0) continue;
		appender->AddInputData(pd);
	}
	appender->Update();

	auto tri = vtkSmartPointer<vtkTriangleFilter>::New();
	tri->SetInputData(appender->GetOutput());
	tri->Update();

	auto merged = vtkSmartPointer<vtkPolyData>::New();
	merged->ShallowCopy(tri->GetOutput());
	return merged;
}


struct VoxelizationParams {
	int baseResolution = 64;
	int maxDepth = 6;
	int marginVoxels = 2;
	int insideRefineDistanceVoxels = 3;
	std::string refineMode = "linear";
};

struct GroupSurfaceEntry {
	vtkSmartPointer<vtkPolyData> surface;
	vtkSmartPointer<vtkImplicitPolyDataDistance> implicit;
	int groupId = -1;
	int surfaceType = -1;
};

struct RegionSurfaceEntry {
	int regionId = -1;
	vtkSmartPointer<vtkPolyData> combined;
	vtkSmartPointer<vtkImplicitPolyDataDistance> implicit;
	std::vector<GroupSurfaceEntry> groups;
};

enum class VoxelState : unsigned char {
	Outside = 0,
	Inside = 1,
	On = 2
};

VoxelState ClassifyLeaf(double dist, double cellSize) {
	const double halfDiag = 0.5 * std::sqrt(3.0) * cellSize;
	if (dist > halfDiag) return VoxelState::Outside;
	if (dist < -halfDiag) return VoxelState::Inside;
	return VoxelState::On;
}

void EnsureArraySize(vtkDataArray* arr, vtkIdType idx, double defaultValue = 0.0) {
	if (!arr) return;
	const vtkIdType oldSize = arr->GetNumberOfTuples();
	if (idx < oldSize) return;
	arr->SetNumberOfTuples(idx + 1);
	for (vtkIdType i = oldSize; i <= idx; ++i) {
		arr->SetComponent(i, 0, defaultValue);
	}
}

double RangeOrNaN(int value, int minValue, int maxValue) {
	if (value < minValue || value > maxValue) {
		return std::numeric_limits<double>::quiet_NaN();
	}
	return static_cast<double>(value);
}

double SurfaceTypeOrNaN(int value) {
	if (value == 0 || value == 1) {
		return static_cast<double>(value);
	}
	return std::numeric_limits<double>::quiet_NaN();
}

vtkSmartPointer<vtkStaticPointLocator> BuildPointLocator(vtkPolyData* surface) {
	if (!surface || surface->GetNumberOfPoints() == 0) return nullptr;
	auto locator = vtkSmartPointer<vtkStaticPointLocator>::New();
	locator->SetDataSet(surface);
	locator->BuildLocator();
	return locator;
}

bool LeafHasSurfacePoint(vtkAbstractPointLocator* locator, const double bounds[6]) {
	if (!locator) return false;
	auto* pointSet = vtkPointSet::SafeDownCast(locator->GetDataSet());
	if (!pointSet) return false;
	const double center[3] = {
		0.5 * (bounds[0] + bounds[1]),
		0.5 * (bounds[2] + bounds[3]),
		0.5 * (bounds[4] + bounds[5])
	};
	const double dx = bounds[1] - bounds[0];
	const double dy = bounds[3] - bounds[2];
	const double dz = bounds[5] - bounds[4];
	const double halfDiag = 0.5 * std::sqrt(dx * dx + dy * dy + dz * dz);

	auto ids = vtkSmartPointer<vtkIdList>::New();
	locator->FindPointsWithinRadius(halfDiag, center, ids);
	for (vtkIdType i = 0; i < ids->GetNumberOfIds(); ++i) {
		double p[3] = {0.0, 0.0, 0.0};
		pointSet->GetPoint(ids->GetId(i), p);
		if (p[0] >= bounds[0] && p[0] <= bounds[1]
			&& p[1] >= bounds[2] && p[1] <= bounds[3]
			&& p[2] >= bounds[4] && p[2] <= bounds[5]) {
			return true;
		}
	}
	return false;
}

void SanitizeIntegerRangeArray(vtkDoubleArray* arr, int minValue, int maxValue) {
	if (!arr) return;
	const double nanValue = std::numeric_limits<double>::quiet_NaN();
	const vtkIdType count = arr->GetNumberOfTuples();
	for (vtkIdType i = 0; i < count; ++i) {
		const double value = arr->GetValue(i);
		if (!std::isfinite(value)) {
			arr->SetValue(i, nanValue);
			continue;
		}
		if (value < static_cast<double>(minValue)
			|| value > static_cast<double>(maxValue)) {
			arr->SetValue(i, nanValue);
			continue;
		}
		const double rounded = std::round(value);
		if (std::abs(value - rounded) > 1e-6) {
			arr->SetValue(i, nanValue);
			continue;
		}
		arr->SetValue(i, rounded);
	}
}

void SanitizeBinaryArray(vtkDoubleArray* arr) {
	SanitizeIntegerRangeArray(arr, 0, 1);
}

void SanitizeRefineArray(vtkDoubleArray* arr, int maxDepth) {
	SanitizeIntegerRangeArray(arr, 0, maxDepth);
}

int ComputeTargetDepth(double distAbs, double rootCellSize, int maxDepth) {
	if (maxDepth <= 0 || rootCellSize <= 0.0) return 0;
	const double minCellSize = rootCellSize / static_cast<double>(1 << std::max(0, maxDepth));
	const double safeDist = std::max(distAbs, minCellSize);
	const double ratio = rootCellSize / safeDist;
	if (ratio <= 1.0) return 0;
	int depth = static_cast<int>(std::floor(std::log2(ratio)));
	if (depth < 0) depth = 0;
	if (depth > maxDepth) depth = maxDepth;
	return depth;
}

std::vector<RegionSurfaceEntry> BuildRegionSurfaceEntries(vtkMultiBlockDataSet* regions) {
	std::vector<RegionSurfaceEntry> entries;
	if (!regions) return entries;
	const int regionCount = static_cast<int>(regions->GetNumberOfBlocks());
	entries.reserve(static_cast<size_t>(std::max(0, regionCount)));
	for (int r = 0; r < regionCount; ++r) {
		auto regionMb = vtkMultiBlockDataSet::SafeDownCast(regions->GetBlock(static_cast<unsigned int>(r)));
		if (!regionMb) continue;
		auto combined = CombineRegionSurfaces(regionMb);
		if (!combined || combined->GetNumberOfCells() == 0) continue;

		RegionSurfaceEntry entry;
		entry.regionId = r;
		entry.combined = combined;
		entry.implicit = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
		entry.implicit->SetInput(combined);

		const int groupCount = static_cast<int>(regionMb->GetNumberOfBlocks());
		entry.groups.reserve(static_cast<size_t>(std::max(0, groupCount)));
		for (int g = 0; g < groupCount; ++g) {
			vtkPolyData* pd = vtkPolyData::SafeDownCast(
				regionMb->GetBlock(static_cast<unsigned int>(g)));
			if (!pd || pd->GetNumberOfCells() == 0) continue;
			const std::string groupName = GetBlockName(regionMb, static_cast<unsigned int>(g));
			GroupSurfaceEntry group;
			group.surface = pd;
			group.implicit = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
			group.implicit->SetInput(pd);
			group.groupId = GetGroupId(pd, groupName);
			group.surfaceType = GetSurfaceTypeField(pd);
			entry.groups.push_back(group);
		}
		entries.push_back(entry);
	}
	return entries;
}

void RefineGlobalHyperTreeCell(
	vtkHyperTreeGridNonOrientedGeometryCursor* cursor,
	const std::vector<RegionSurfaceEntry>& regions,
	int depth,
	int maxDepth,
	double rootCellSize,
	const std::string& refineMode,
	vtkUnsignedCharArray* stateArray,
	vtkDoubleArray* levelArray,
	vtkUnsignedCharArray* activeArray,
	vtkDoubleArray* maskArray,
	vtkDoubleArray* refineArray,
	vtkIntArray* regionArray,
	vtkIntArray* groupArray,
	vtkDoubleArray* surfaceArray,
	vtkUnsignedCharArray* onMaskArray,
	vtkAbstractPointLocator* surfaceLocator,
	bool forceSingleLabel) {
	if (!cursor || regions.empty()) return;

	double bounds[6] = {0, 0, 0, 0, 0, 0};
	cursor->GetBounds(bounds);
	const double cellSize = bounds[1] - bounds[0];
	const double center[3] = {
		0.5 * (bounds[0] + bounds[1]),
		0.5 * (bounds[2] + bounds[3]),
		0.5 * (bounds[4] + bounds[5])
	};
	const double halfDiag = 0.5 * std::sqrt(3.0) * cellSize;

	int intersectRegions = 0;
	int bestRegionIndex = -1;
	double bestAbsDist = std::numeric_limits<double>::infinity();
	double bestRegionDist = 0.0;

	for (size_t i = 0; i < regions.size(); ++i) {
		const double dist = regions[i].implicit->EvaluateFunction(center[0], center[1], center[2]);
		const double distAbs = std::abs(dist);
		if (distAbs <= halfDiag) {
			++intersectRegions;
		}
		if (distAbs < bestAbsDist) {
			bestAbsDist = distAbs;
			bestRegionIndex = static_cast<int>(i);
			bestRegionDist = dist;
		}
	}

	if (intersectRegions > 1 && depth < maxDepth) {
		cursor->SubdivideLeaf();
		for (int child = 0; child < 8; ++child) {
			cursor->ToChild(child);
			RefineGlobalHyperTreeCell(cursor, regions, depth + 1, maxDepth,
				rootCellSize, refineMode, stateArray, levelArray, activeArray,
				maskArray, refineArray, regionArray, groupArray, surfaceArray,
				onMaskArray, surfaceLocator, forceSingleLabel);
			cursor->ToParent();
		}
		return;
	}

	if (bestRegionIndex < 0) {
		return;
	}

	const RegionSurfaceEntry& region = regions[static_cast<size_t>(bestRegionIndex)];
	int intersectGroups = 0;
	int bestGroupIndex = -1;
	double bestGroupAbsDist = std::numeric_limits<double>::infinity();
	for (size_t g = 0; g < region.groups.size(); ++g) {
		const double dist = region.groups[g].implicit->EvaluateFunction(center[0], center[1], center[2]);
		const double distAbs = std::abs(dist);
		if (distAbs <= halfDiag) {
			++intersectGroups;
		}
		if (distAbs < bestGroupAbsDist) {
			bestGroupAbsDist = distAbs;
			bestGroupIndex = static_cast<int>(g);
		}
	}

	if (intersectGroups > 1 && depth < maxDepth) {
		cursor->SubdivideLeaf();
		for (int child = 0; child < 8; ++child) {
			cursor->ToChild(child);
			RefineGlobalHyperTreeCell(cursor, regions, depth + 1, maxDepth,
				rootCellSize, refineMode, stateArray, levelArray, activeArray,
				maskArray, refineArray, regionArray, groupArray, surfaceArray,
				onMaskArray, surfaceLocator, forceSingleLabel);
			cursor->ToParent();
		}
		return;
	}

	if (intersectGroups > 1 && !forceSingleLabel) {
		return;
	}

	int regionId = region.regionId;
	int groupId = -1;
	int surfaceType = -1;
	if (bestGroupIndex >= 0 && bestGroupIndex < static_cast<int>(region.groups.size())) {
		const auto& group = region.groups[static_cast<size_t>(bestGroupIndex)];
		groupId = group.groupId;
		surfaceType = group.surfaceType;
	}

	const std::string mode = ToLower(refineMode);
	if (mode == "distance_jump" && depth < maxDepth) {
		const int targetDepth = ComputeTargetDepth(bestAbsDist, rootCellSize, maxDepth);
		if (targetDepth > depth) {
			cursor->SubdivideLeaf();
			for (int child = 0; child < 8; ++child) {
				cursor->ToChild(child);
				RefineGlobalHyperTreeCell(cursor, regions, depth + 1, maxDepth,
					rootCellSize, refineMode, stateArray, levelArray, activeArray,
					maskArray, refineArray, regionArray, groupArray, surfaceArray,
					onMaskArray, surfaceLocator, forceSingleLabel);
				cursor->ToParent();
			}
			return;
		}
	}

	const VoxelState state = ClassifyLeaf(bestRegionDist, cellSize);
	const vtkIdType nodeId = cursor->GetGlobalNodeIndex();
	if (nodeId == vtkHyperTreeGrid::InvalidIndex || nodeId < 0) {
		return;
	}
	const double nanValue = std::numeric_limits<double>::quiet_NaN();
	EnsureArraySize(stateArray, nodeId, static_cast<double>(VoxelState::Outside));
	EnsureArraySize(levelArray, nodeId, nanValue);
	EnsureArraySize(activeArray, nodeId, 0.0);
	EnsureArraySize(maskArray, nodeId, nanValue);
	EnsureArraySize(refineArray, nodeId, nanValue);
	EnsureArraySize(regionArray, nodeId, -1.0);
	EnsureArraySize(groupArray, nodeId, -1.0);
	EnsureArraySize(surfaceArray, nodeId, nanValue);
	EnsureArraySize(onMaskArray, nodeId, 0.0);
	stateArray->SetValue(nodeId, static_cast<unsigned char>(state));
	levelArray->SetValue(nodeId, RangeOrNaN(depth, 0, maxDepth));
	activeArray->SetValue(nodeId, cursor->IsLeaf() ? 1 : 0);
	maskArray->SetValue(nodeId, state == VoxelState::Outside ? 0.0 : 1.0);
	refineArray->SetValue(nodeId, RangeOrNaN(depth, 0, maxDepth));
	regionArray->SetValue(nodeId, regionId);
	groupArray->SetValue(nodeId, groupId);
	surfaceArray->SetValue(nodeId, SurfaceTypeOrNaN(surfaceType));
	onMaskArray->SetValue(nodeId, LeafHasSurfacePoint(surfaceLocator, bounds) ? 1 : 0);
}

void RefineHyperTreeCell(
	vtkHyperTreeGridNonOrientedGeometryCursor* cursor,
	vtkImplicitPolyDataDistance* implicit,
	int depth,
	int maxDepth,
	double rootCellSize,
	const std::string& refineMode,
	vtkUnsignedCharArray* stateArray,
	vtkDoubleArray* levelArray,
	vtkUnsignedCharArray* activeArray,
	vtkDoubleArray* maskArray,
	vtkDoubleArray* refineArray,
	vtkUnsignedCharArray* onMaskArray,
	vtkAbstractPointLocator* surfaceLocator) {
	if (!cursor || !implicit) return;

	double bounds[6] = {0, 0, 0, 0, 0, 0};
	cursor->GetBounds(bounds);
	const double cellSize = bounds[1] - bounds[0];
	const double center[3] = {
		0.5 * (bounds[0] + bounds[1]),
		0.5 * (bounds[2] + bounds[3]),
		0.5 * (bounds[4] + bounds[5])
	};
	const double dist = implicit->EvaluateFunction(center[0], center[1], center[2]);
	const double distAbs = std::abs(dist);
	const double halfDiag = 0.5 * std::sqrt(3.0) * cellSize;

	const std::string mode = ToLower(refineMode);
	if (mode == "distance_jump") {
		const int targetDepth = ComputeTargetDepth(distAbs, rootCellSize, maxDepth);
		if (targetDepth > depth && depth < maxDepth) {
			cursor->SubdivideLeaf();
			for (int child = 0; child < 8; ++child) {
				cursor->ToChild(child);
				RefineHyperTreeCell(cursor, implicit, depth + 1, maxDepth,
					rootCellSize, refineMode, stateArray, levelArray, activeArray,
					maskArray, refineArray, onMaskArray, surfaceLocator);
				cursor->ToParent();
			}
			return;
		}
	} else {
		const bool intersectsSurface = distAbs <= halfDiag;
		if (intersectsSurface && depth < maxDepth) {
			cursor->SubdivideLeaf();
			for (int child = 0; child < 8; ++child) {
				cursor->ToChild(child);
				RefineHyperTreeCell(cursor, implicit, depth + 1, maxDepth,
					rootCellSize, refineMode, stateArray, levelArray, activeArray,
					maskArray, refineArray, onMaskArray, surfaceLocator);
				cursor->ToParent();
			}
			return;
		}
	}

	const vtkIdType nodeId = cursor->GetGlobalNodeIndex();
	if (nodeId == vtkHyperTreeGrid::InvalidIndex || nodeId < 0) {
		return;
	}
	const VoxelState state = ClassifyLeaf(dist, cellSize);
	const double nanValue = std::numeric_limits<double>::quiet_NaN();
	EnsureArraySize(stateArray, nodeId, static_cast<double>(VoxelState::Outside));
	EnsureArraySize(levelArray, nodeId, nanValue);
	EnsureArraySize(activeArray, nodeId, 0.0);
	EnsureArraySize(maskArray, nodeId, nanValue);
	EnsureArraySize(refineArray, nodeId, nanValue);
	EnsureArraySize(onMaskArray, nodeId, 0.0);
	stateArray->SetValue(nodeId, static_cast<unsigned char>(state));
	levelArray->SetValue(nodeId, RangeOrNaN(depth, 0, maxDepth));
	activeArray->SetValue(nodeId, cursor->IsLeaf() ? 1 : 0);
	maskArray->SetValue(nodeId, state == VoxelState::Outside ? 0.0 : 1.0);
	refineArray->SetValue(nodeId, RangeOrNaN(depth, 0, maxDepth));
	onMaskArray->SetValue(nodeId, LeafHasSurfacePoint(surfaceLocator, bounds) ? 1 : 0);
}

vtkSmartPointer<vtkHyperTreeGrid> VoxelizeRegionAdaptiveHTG(
	vtkPolyData* surface,
	const VoxelizationParams& params) {
	if (!surface || surface->GetNumberOfCells() == 0) return nullptr;


	double bounds[6] = {0, 0, 0, 0, 0, 0};
	surface->GetBounds(bounds);

	const double dx = bounds[1] - bounds[0];
	const double dy = bounds[3] - bounds[2];
	const double dz = bounds[5] - bounds[4];
	const double diag = std::sqrt(dx * dx + dy * dy + dz * dz);
	const double voxelSize = diag / static_cast<double>(std::max(1, params.baseResolution));
	const double margin = voxelSize * static_cast<double>(std::max(0, params.marginVoxels));

	double minB[3] = {bounds[0] - margin, bounds[2] - margin, bounds[4] - margin};
	double maxB[3] = {bounds[1] + margin, bounds[3] + margin, bounds[5] + margin};

	const int nx = std::max(1, static_cast<int>(std::ceil((maxB[0] - minB[0]) / voxelSize)));
	const int ny = std::max(1, static_cast<int>(std::ceil((maxB[1] - minB[1]) / voxelSize)));
	const int nz = std::max(1, static_cast<int>(std::ceil((maxB[2] - minB[2]) / voxelSize)));

	auto grid = vtkSmartPointer<vtkHyperTreeGrid>::New();
	grid->SetDimensions(nx + 1, ny + 1, nz + 1);
	grid->SetBranchFactor(2);

	auto xCoords = vtkSmartPointer<vtkDoubleArray>::New();
	auto yCoords = vtkSmartPointer<vtkDoubleArray>::New();
	auto zCoords = vtkSmartPointer<vtkDoubleArray>::New();
	xCoords->SetNumberOfTuples(nx + 1);
	yCoords->SetNumberOfTuples(ny + 1);
	zCoords->SetNumberOfTuples(nz + 1);
	for (int i = 0; i <= nx; ++i) xCoords->SetValue(i, minB[0] + i * voxelSize);
	for (int j = 0; j <= ny; ++j) yCoords->SetValue(j, minB[1] + j * voxelSize);
	for (int k = 0; k <= nz; ++k) zCoords->SetValue(k, minB[2] + k * voxelSize);
	grid->SetXCoordinates(xCoords);
	grid->SetYCoordinates(yCoords);
	grid->SetZCoordinates(zCoords);

	auto stateArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
	stateArray->SetName("VoxelState");
	stateArray->SetNumberOfComponents(1);

	auto levelArray = vtkSmartPointer<vtkDoubleArray>::New();
	levelArray->SetName("RefineLevel");
	levelArray->SetNumberOfComponents(1);

	auto activeArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
	activeArray->SetName("HTGActive");
	activeArray->SetNumberOfComponents(1);

	auto maskArray = vtkSmartPointer<vtkDoubleArray>::New();
	maskArray->SetName("HTGMask");
	maskArray->SetNumberOfComponents(1);

	auto refineArray = vtkSmartPointer<vtkDoubleArray>::New();
	refineArray->SetName("HTGRefine");
	refineArray->SetNumberOfComponents(1);

	auto onMaskArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
	onMaskArray->SetName("HTGOnMask");
	onMaskArray->SetNumberOfComponents(1);

	auto implicit = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
	implicit->SetInput(surface);
	auto surfaceLocator = BuildPointLocator(surface);

	const vtkIdType treeCount = grid->GetMaxNumberOfTrees();
	for (vtkIdType treeId = 0; treeId < treeCount; ++treeId) {
		auto cursor = vtkSmartPointer<vtkHyperTreeGridNonOrientedGeometryCursor>::New();
		grid->InitializeNonOrientedGeometryCursor(cursor, treeId, true);
		if (!cursor->HasTree()) continue;
		double rootBounds[6] = {0, 0, 0, 0, 0, 0};
		cursor->GetBounds(rootBounds);
		const double rootCellSize = rootBounds[1] - rootBounds[0];
		RefineHyperTreeCell(cursor, implicit, 0, params.maxDepth,
			rootCellSize, params.refineMode, stateArray, levelArray, activeArray,
			maskArray, refineArray, onMaskArray, surfaceLocator);
	}

	SanitizeRefineArray(refineArray, params.maxDepth);
	SanitizeIntegerRangeArray(levelArray, 0, params.maxDepth);
	SanitizeBinaryArray(maskArray);


	grid->GetCellData()->AddArray(stateArray);
	grid->GetCellData()->AddArray(levelArray);
	grid->GetCellData()->AddArray(activeArray);
	grid->GetCellData()->AddArray(maskArray);
	grid->GetCellData()->AddArray(refineArray);
	grid->GetCellData()->AddArray(onMaskArray);
	return grid;
}

vtkSmartPointer<vtkHyperTreeGrid> VoxelizeGlobalAdaptiveHTG(
	vtkMultiBlockDataSet* regions,
	const VoxelizationParams& params,
	bool forceSingleLabel) {
	if (!regions) return nullptr;

	auto regionEntries = BuildRegionSurfaceEntries(regions);
	if (regionEntries.empty()) return nullptr;

	double bounds[6] = {0, 0, 0, 0, 0, 0};
	bool boundsInit = false;
	for (const auto& entry : regionEntries) {
		double b[6] = {0, 0, 0, 0, 0, 0};
		entry.combined->GetBounds(b);
		if (!boundsInit) {
			for (int i = 0; i < 6; ++i) bounds[i] = b[i];
			boundsInit = true;
		} else {
			bounds[0] = std::min(bounds[0], b[0]);
			bounds[1] = std::max(bounds[1], b[1]);
			bounds[2] = std::min(bounds[2], b[2]);
			bounds[3] = std::max(bounds[3], b[3]);
			bounds[4] = std::min(bounds[4], b[4]);
			bounds[5] = std::max(bounds[5], b[5]);
		}
	}
	if (!boundsInit) return nullptr;

	const double dx = bounds[1] - bounds[0];
	const double dy = bounds[3] - bounds[2];
	const double dz = bounds[5] - bounds[4];
	const double diag = std::sqrt(dx * dx + dy * dy + dz * dz);
	const double voxelSize = diag / static_cast<double>(std::max(1, params.baseResolution));
	const double margin = voxelSize * static_cast<double>(std::max(0, params.marginVoxels));

	double minB[3] = {bounds[0] - margin, bounds[2] - margin, bounds[4] - margin};
	double maxB[3] = {bounds[1] + margin, bounds[3] + margin, bounds[5] + margin};

	const int nx = std::max(1, static_cast<int>(std::ceil((maxB[0] - minB[0]) / voxelSize)));
	const int ny = std::max(1, static_cast<int>(std::ceil((maxB[1] - minB[1]) / voxelSize)));
	const int nz = std::max(1, static_cast<int>(std::ceil((maxB[2] - minB[2]) / voxelSize)));

	auto grid = vtkSmartPointer<vtkHyperTreeGrid>::New();
	grid->SetDimensions(nx + 1, ny + 1, nz + 1);
	grid->SetBranchFactor(2);

	auto xCoords = vtkSmartPointer<vtkDoubleArray>::New();
	auto yCoords = vtkSmartPointer<vtkDoubleArray>::New();
	auto zCoords = vtkSmartPointer<vtkDoubleArray>::New();
	xCoords->SetNumberOfTuples(nx + 1);
	yCoords->SetNumberOfTuples(ny + 1);
	zCoords->SetNumberOfTuples(nz + 1);
	for (int i = 0; i <= nx; ++i) xCoords->SetValue(i, minB[0] + i * voxelSize);
	for (int j = 0; j <= ny; ++j) yCoords->SetValue(j, minB[1] + j * voxelSize);
	for (int k = 0; k <= nz; ++k) zCoords->SetValue(k, minB[2] + k * voxelSize);
	grid->SetXCoordinates(xCoords);
	grid->SetYCoordinates(yCoords);
	grid->SetZCoordinates(zCoords);

	auto stateArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
	stateArray->SetName("VoxelState");
	stateArray->SetNumberOfComponents(1);

	auto levelArray = vtkSmartPointer<vtkDoubleArray>::New();
	levelArray->SetName("RefineLevel");
	levelArray->SetNumberOfComponents(1);

	auto activeArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
	activeArray->SetName("HTGActive");
	activeArray->SetNumberOfComponents(1);

	auto maskArray = vtkSmartPointer<vtkDoubleArray>::New();
	maskArray->SetName("HTGMask");
	maskArray->SetNumberOfComponents(1);

	auto refineArray = vtkSmartPointer<vtkDoubleArray>::New();
	refineArray->SetName("HTGRefine");
	refineArray->SetNumberOfComponents(1);

	auto regionArray = vtkSmartPointer<vtkIntArray>::New();
	regionArray->SetName("RegionId");
	regionArray->SetNumberOfComponents(1);

	auto groupArray = vtkSmartPointer<vtkIntArray>::New();
	groupArray->SetName("GroupId");
	groupArray->SetNumberOfComponents(1);

	auto surfaceArray = vtkSmartPointer<vtkDoubleArray>::New();
	surfaceArray->SetName("SurfaceType");
	surfaceArray->SetNumberOfComponents(1);

	auto onMaskArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
	onMaskArray->SetName("HTGOnMask");
	onMaskArray->SetNumberOfComponents(1);

	auto allSurfaces = vtkSmartPointer<vtkAppendPolyData>::New();
	for (const auto& entry : regionEntries) {
		if (entry.combined && entry.combined->GetNumberOfPoints() > 0) {
			allSurfaces->AddInputData(entry.combined);
		}
	}
	allSurfaces->Update();
	auto allTri = vtkSmartPointer<vtkTriangleFilter>::New();
	allTri->SetInputData(allSurfaces->GetOutput());
	allTri->Update();
	auto allMerged = vtkSmartPointer<vtkPolyData>::New();
	allMerged->ShallowCopy(allTri->GetOutput());
	auto surfaceLocator = BuildPointLocator(allMerged);

	const vtkIdType treeCount = grid->GetMaxNumberOfTrees();
	for (vtkIdType treeId = 0; treeId < treeCount; ++treeId) {
		auto cursor = vtkSmartPointer<vtkHyperTreeGridNonOrientedGeometryCursor>::New();
		grid->InitializeNonOrientedGeometryCursor(cursor, treeId, true);
		if (!cursor->HasTree()) continue;
		double rootBounds[6] = {0, 0, 0, 0, 0, 0};
		cursor->GetBounds(rootBounds);
		const double rootCellSize = rootBounds[1] - rootBounds[0];
		RefineGlobalHyperTreeCell(cursor, regionEntries, 0, params.maxDepth,
			rootCellSize, params.refineMode, stateArray, levelArray, activeArray,
			maskArray, refineArray, regionArray, groupArray, surfaceArray,
			onMaskArray, surfaceLocator, forceSingleLabel);
	}

	SanitizeRefineArray(refineArray, params.maxDepth);
	SanitizeIntegerRangeArray(levelArray, 0, params.maxDepth);
	SanitizeBinaryArray(maskArray);
	SanitizeBinaryArray(surfaceArray);

	grid->GetCellData()->AddArray(stateArray);
	grid->GetCellData()->AddArray(levelArray);
	grid->GetCellData()->AddArray(activeArray);
	grid->GetCellData()->AddArray(maskArray);
	grid->GetCellData()->AddArray(refineArray);
	grid->GetCellData()->AddArray(regionArray);
	grid->GetCellData()->AddArray(groupArray);
	grid->GetCellData()->AddArray(surfaceArray);
	grid->GetCellData()->AddArray(onMaskArray);
	return grid;
}

vtkSmartPointer<vtkImageData> VoxelizeRegionAdaptive(
	vtkPolyData* surface,
	const VoxelizationParams& params) {
	if (!surface || surface->GetNumberOfCells() == 0) return nullptr;

	double bounds[6] = {0, 0, 0, 0, 0, 0};
	surface->GetBounds(bounds);

	const double dx = bounds[1] - bounds[0];
	const double dy = bounds[3] - bounds[2];
	const double dz = bounds[5] - bounds[4];
	const double diag = std::sqrt(dx * dx + dy * dy + dz * dz);
	const double voxelSize = diag / static_cast<double>(std::max(1, params.baseResolution));
	const double margin = voxelSize * static_cast<double>(std::max(0, params.marginVoxels));

	double minB[3] = {bounds[0] - margin, bounds[2] - margin, bounds[4] - margin};
	double maxB[3] = {bounds[1] + margin, bounds[3] + margin, bounds[5] + margin};

	const int nx = std::max(1, static_cast<int>(std::ceil((maxB[0] - minB[0]) / voxelSize)));
	const int ny = std::max(1, static_cast<int>(std::ceil((maxB[1] - minB[1]) / voxelSize)));
	const int nz = std::max(1, static_cast<int>(std::ceil((maxB[2] - minB[2]) / voxelSize)));

	auto image = vtkSmartPointer<vtkImageData>::New();
	image->SetOrigin(minB[0], minB[1], minB[2]);
	image->SetSpacing(voxelSize, voxelSize, voxelSize);
	image->SetDimensions(nx + 1, ny + 1, nz + 1);

	auto stateArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
	stateArray->SetName("VoxelState");
	stateArray->SetNumberOfComponents(1);
	stateArray->SetNumberOfTuples(static_cast<vtkIdType>(nx) * ny * nz);

	auto levelArray = vtkSmartPointer<vtkIntArray>::New();
	levelArray->SetName("RefineLevel");
	levelArray->SetNumberOfComponents(1);
	levelArray->SetNumberOfTuples(static_cast<vtkIdType>(nx) * ny * nz);

	auto implicit = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
	implicit->SetInput(surface);

	vtkIdType idx = 0;
	const double maxInsideDist = voxelSize * static_cast<double>(std::max(1, params.insideRefineDistanceVoxels));
	for (int k = 0; k < nz; ++k) {
		const double zc = minB[2] + (k + 0.5) * voxelSize;
		for (int j = 0; j < ny; ++j) {
			const double yc = minB[1] + (j + 0.5) * voxelSize;
			for (int i = 0; i < nx; ++i) {
				const double xc = minB[0] + (i + 0.5) * voxelSize;
				double center[3] = {xc, yc, zc};
				const double dist = implicit->EvaluateFunction(center[0], center[1], center[2]);
				const VoxelState state = ClassifyLeaf(dist, voxelSize);
				stateArray->SetValue(idx, static_cast<unsigned char>(state));

				const double distAbs = std::abs(dist);
				const double halfDiag = 0.5 * std::sqrt(3.0) * voxelSize;
				int targetDepth = 0;
				if (distAbs <= halfDiag) {
					targetDepth = params.maxDepth;
				} else if (dist < 0.0 && distAbs <= maxInsideDist) {
					const double t = 1.0 - std::min(distAbs / maxInsideDist, 1.0);
					targetDepth = std::max(targetDepth, static_cast<int>(std::round(t * params.maxDepth)));
				}
				levelArray->SetValue(idx, targetDepth);
				++idx;
			}
		}
	}

	image->GetCellData()->SetScalars(stateArray);
	image->GetCellData()->AddArray(levelArray);
	return image;
}

vtkSmartPointer<vtkImageData> BuildRefinedImage(
	vtkPolyData* surface,
	const VoxelizationParams& params,
	const double minB[3],
	const double voxelSize,
	int nx,
	int ny,
	int nz,
	vtkUnsignedCharArray* baseState,
	vtkIntArray* baseLevel) {
	const int maxScale = 1 << std::max(0, params.maxDepth);
	const int rx = std::max(1, nx * maxScale);
	const int ry = std::max(1, ny * maxScale);
	const int rz = std::max(1, nz * maxScale);
	const double fineSize = voxelSize / static_cast<double>(maxScale);

	auto refined = vtkSmartPointer<vtkImageData>::New();
	refined->SetOrigin(minB[0], minB[1], minB[2]);
	refined->SetSpacing(fineSize, fineSize, fineSize);
	refined->SetDimensions(rx + 1, ry + 1, rz + 1);

	auto refinedState = vtkSmartPointer<vtkUnsignedCharArray>::New();
	refinedState->SetName("VoxelState");
	refinedState->SetNumberOfComponents(1);
	refinedState->SetNumberOfTuples(static_cast<vtkIdType>(rx) * ry * rz);

	auto refineSource = vtkSmartPointer<vtkIntArray>::New();
	refineSource->SetName("RefineSourceLevel");
	refineSource->SetNumberOfComponents(1);
	refineSource->SetNumberOfTuples(static_cast<vtkIdType>(rx) * ry * rz);

	auto implicit = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
	implicit->SetInput(surface);

	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i) {
				const vtkIdType baseIdx =
					static_cast<vtkIdType>((k * ny + j) * nx + i);
				const int level = baseLevel->GetValue(baseIdx);
				const unsigned char baseStateValue = baseState->GetValue(baseIdx);
				const int effectiveLevel = std::max(0, std::min(level, params.maxDepth));
				const int subScale = 1 << effectiveLevel;
				const int fillScale = maxScale / subScale;
				const double subSize = voxelSize / static_cast<double>(subScale);

				for (int sk = 0; sk < subScale; ++sk) {
					const double subZ0 = minB[2] + (k + sk / static_cast<double>(subScale)) * voxelSize;
					for (int sj = 0; sj < subScale; ++sj) {
						const double subY0 = minB[1] + (j + sj / static_cast<double>(subScale)) * voxelSize;
						for (int si = 0; si < subScale; ++si) {
							const double subX0 = minB[0] + (i + si / static_cast<double>(subScale)) * voxelSize;
							unsigned char refinedStateValue = baseStateValue;
							if (effectiveLevel > 0) {
								const double center[3] = {
									subX0 + 0.5 * subSize,
									subY0 + 0.5 * subSize,
									subZ0 + 0.5 * subSize
								};
								const double dist = implicit->EvaluateFunction(
									center[0], center[1], center[2]);
								refinedStateValue = static_cast<unsigned char>(
									ClassifyLeaf(dist, subSize));
							}

							const int fi0 = (i * maxScale) + (si * fillScale);
							const int fj0 = (j * maxScale) + (sj * fillScale);
							const int fk0 = (k * maxScale) + (sk * fillScale);
							for (int fk = 0; fk < fillScale; ++fk) {
								for (int fj = 0; fj < fillScale; ++fj) {
									vtkIdType baseOut = static_cast<vtkIdType>(
										((fk0 + fk) * ry + (fj0 + fj)) * rx + fi0);
									for (int fi = 0; fi < fillScale; ++fi) {
										refinedState->SetValue(baseOut + fi, refinedStateValue);
										refineSource->SetValue(baseOut + fi, effectiveLevel);
									}
								}
							}
						}
					}
				}
			}
		}
	}

	refined->GetCellData()->SetScalars(refinedState);
	refined->GetCellData()->AddArray(refineSource);
	return refined;
}

void WriteVoxelImage(vtkImageData* image, const std::string& path) {
	if (!image || path.empty()) return;
	auto writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetFileName(path.c_str());
	writer->SetInputData(image);
	writer->Write();
}

void WriteHyperTreeGrid(vtkHyperTreeGrid* grid, const std::string& path) {
	if (!grid || path.empty()) return;
	auto writer = vtkSmartPointer<vtkXMLHyperTreeGridWriter>::New();
	writer->SetFileName(path.c_str());
	writer->SetInputData(grid);
	writer->Write();
}

void WriteHyperTreeGridAsUnstructuredGrid(vtkHyperTreeGrid* grid, const std::string& path) {
	if (!grid || path.empty()) return;
	auto converter = vtkSmartPointer<vtkHyperTreeGridToUnstructuredGrid>::New();
	converter->SetInputData(grid);
	converter->Update();

	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(path.c_str());
	writer->SetInputData(converter->GetOutput());
	writer->Write();
}


} // namespace

void CenterlineBase(vtkMultiBlockDataSet* regions,
	const std::string& reportLevel,
	int reportTableRows) {
	if (!regions) {
		std::cerr << "CenterlineBase: regions is null." << std::endl;
		return;
	}

	const int regionCount = static_cast<int>(regions->GetNumberOfBlocks());
	if (regionCount <= 0) {
		std::cout << "CenterlineBase: no regions found." << std::endl;
		return;
	}

	const unsigned int hwThreads = std::max(1u, std::thread::hardware_concurrency());
	const bool detailed = (ToLower(reportLevel) == "long");
	const int maxRows = std::max(1, reportTableRows);
	int totalXlets = 0;
	int totalWalls = 0;
	int totalUnknown = 0;

	struct SurfaceRow {
		std::string region;
		int groupId = -1;
		int cells = 0;
		double area = 0.0;
	};
	std::vector<SurfaceRow> xletRows;
	std::vector<SurfaceRow> wallRows;

	for (int r = 0; r < regionCount; ++r) {
		auto regionMb = vtkMultiBlockDataSet::SafeDownCast(regions->GetBlock(static_cast<unsigned int>(r)));
		if (!regionMb) {
			continue;
		}

		const int groupCount = static_cast<int>(regionMb->GetNumberOfBlocks());
		int regionXlets = 0;
		int regionWalls = 0;
		int regionUnknown = 0;
		std::vector<std::pair<int, double>> xletAreas;
		xletAreas.reserve(static_cast<size_t>(groupCount));
		const std::string regionName = GetBlockName(regions, static_cast<unsigned int>(r));
		const std::string regionLabel = regionName.empty()
			? std::to_string(r)
			: regionName;

		for (int g = 0; g < groupCount; ++g) {
			vtkPolyData* pd = vtkPolyData::SafeDownCast(
				regionMb->GetBlock(static_cast<unsigned int>(g)));
			if (!pd) {
				continue;
			}
			const int surfaceType = GetSurfaceTypeField(pd);
			const std::string groupName = GetBlockName(regionMb, static_cast<unsigned int>(g));
			const int groupId = GetGroupId(pd, groupName);
			const int cellCount = static_cast<int>(pd->GetNumberOfCells());
			const double area = ComputeSurfaceArea(pd);
			if (surfaceType == 1) {
				++regionXlets;
				xletAreas.emplace_back(groupId, area);
				if (detailed) {
					xletRows.push_back({regionLabel, groupId, cellCount, area});
				}
			} else if (surfaceType == 0) {
				++regionWalls;
				if (detailed) {
					wallRows.push_back({regionLabel, groupId, cellCount, area});
				}
			} else {
				++regionUnknown;
			}
		}

		const unsigned int cpuIndex = static_cast<unsigned int>(r) % hwThreads;
		std::cout << "CenterlineBase: region "
				  << (regionName.empty() ? std::to_string(r) : regionName)
				  << " assigned_cpu=" << cpuIndex
				  << " groups=" << groupCount
				  << " xlets=" << regionXlets
				  << " walls=" << regionWalls
				  << std::endl;

		std::sort(xletAreas.begin(), xletAreas.end(), [](const auto& a, const auto& b) {
			return a.second > b.second;
		});

		const size_t topCount = std::min<size_t>(3, xletAreas.size());
		if (topCount > 0) {
			std::cout << "CenterlineBase: top_xlets=";
			for (size_t i = 0; i < topCount; ++i) {
				const int gid = xletAreas[i].first;
				const double area = xletAreas[i].second;
				std::cout << (i == 0 ? " " : ", ")
						  << "GroupId=" << gid
						  << " area=" << std::fixed << std::setprecision(3) << area;
			}
			std::cout << std::endl;
		}

		if (regionXlets > 100 && !xletAreas.empty()) {
			const int rootGroupId = xletAreas.front().first;
			const double rootArea = xletAreas.front().second;
			std::cout << "CenterlineBase: root_xlet=GroupId=" << rootGroupId
					  << " area=" << std::fixed << std::setprecision(3) << rootArea
					  << " (xlets>100)"
					  << std::endl;
		}

		totalXlets += regionXlets;
		totalWalls += regionWalls;
		totalUnknown += regionUnknown;
	}

	if (detailed) {
		auto byArea = [](const SurfaceRow& a, const SurfaceRow& b) {
			return a.area > b.area;
		};
		std::sort(xletRows.begin(), xletRows.end(), byArea);
		std::sort(wallRows.begin(), wallRows.end(), byArea);

		auto formatArea = [](double value) {
			std::ostringstream oss;
			oss << std::fixed << std::setprecision(3) << value;
			return oss.str();
		};
		auto buildRows = [&](const std::vector<SurfaceRow>& surfaces) {
			std::vector<std::vector<std::string>> rows;
			rows.reserve(surfaces.size());
			for (const auto& surface : surfaces) {
				rows.push_back({
					surface.region,
					std::to_string(surface.groupId),
					std::to_string(surface.cells),
					formatArea(surface.area)
				});
			}
			return rows;
		};

		const std::vector<std::string> headers = {"Region", "GroupId", "Cells", "Area"};
		const std::vector<size_t> widths = {12, 10, 8, 12};

		if (!xletRows.empty()) {
			auto rows = buildRows(xletRows);
			const size_t shown = std::min<size_t>(rows.size(), static_cast<size_t>(maxRows));
			std::cout << "CenterlineBase: xlet table (showing " << shown
					  << " of " << rows.size() << ")\n";
			std::cout << BuildTextTable(headers, rows, widths, static_cast<size_t>(maxRows));
		}
		if (!wallRows.empty()) {
			auto rows = buildRows(wallRows);
			const size_t shown = std::min<size_t>(rows.size(), static_cast<size_t>(maxRows));
			std::cout << "CenterlineBase: wall table (showing " << shown
					  << " of " << rows.size() << ")\n";
			std::cout << BuildTextTable(headers, rows, widths, static_cast<size_t>(maxRows));
		}
	}

	std::cout << "CenterlineBase: regions=" << regionCount
			  << " total_xlets=" << totalXlets
			  << " total_walls=" << totalWalls
			  << " unknown_surfaces=" << totalUnknown
			  << std::endl;
}

void VoxelizeRegionsBase(vtkMultiBlockDataSet* regions, int baseResolution, int maxDepth, int insideRefineDist) {
	VoxelizeRegionsBase(regions, baseResolution, maxDepth, insideRefineDist, true, false, "", true, "linear");
}

void VoxelizeRegionsBase(vtkMultiBlockDataSet* regions,
	int baseResolution,
	int maxDepth,
	int insideRefineDist,
	bool writeRegion,
	bool writeGlobal,
	const std::string& globalOutputPath,
	bool forceSingleLabel,
	const std::string& refineMode) {
	if (!regions) {
		std::cerr << "VoxelizeRegionsBase: regions is null." << std::endl;
		return;
	}

	const int regionCount = static_cast<int>(regions->GetNumberOfBlocks());
	if (regionCount <= 0) {
		std::cout << "VoxelizeRegionsBase: no regions found." << std::endl;
		return;
	}

	std::filesystem::create_directories("output/centerline_voxels");
	VoxelizationParams params;
	params.baseResolution = std::max(4, baseResolution);
	params.maxDepth = std::max(0, maxDepth);
	params.insideRefineDistanceVoxels = std::max(1, insideRefineDist);
	params.refineMode = ToLower(refineMode);

	if (writeRegion) {
		for (int r = 0; r < regionCount; ++r) {
			auto regionMb = vtkMultiBlockDataSet::SafeDownCast(regions->GetBlock(static_cast<unsigned int>(r)));
			if (!regionMb) continue;

			auto combined = CombineRegionSurfaces(regionMb);
			if (!combined || combined->GetNumberOfCells() == 0) continue;

			auto htg = VoxelizeRegionAdaptiveHTG(combined, params);
			if (!htg) continue;

			std::string regionName = GetBlockName(regions, static_cast<unsigned int>(r));
			if (regionName.empty()) {
				regionName = "Region_" + std::to_string(r);
			}
			std::string outPath = "output/centerline_voxels/" + regionName + ".vth";
			WriteHyperTreeGrid(htg, outPath);
			std::filesystem::path meshPath(outPath);
			meshPath.replace_extension(".vtu");
			WriteHyperTreeGridAsUnstructuredGrid(htg, meshPath.string());
			std::cout << "VoxelizeRegionsBase: wrote " << outPath << std::endl;
		}
	}

	if (writeGlobal) {
		std::string outPath = globalOutputPath.empty()
			? "output/centerline_voxels/regions_global.vth"
			: globalOutputPath;
		auto globalHtg = VoxelizeGlobalAdaptiveHTG(regions, params, forceSingleLabel);
		if (globalHtg) {
			WriteHyperTreeGrid(globalHtg, outPath);
			std::filesystem::path meshPath(outPath);
			meshPath.replace_extension(".vtu");
			WriteHyperTreeGridAsUnstructuredGrid(globalHtg, meshPath.string());
			std::cout << "VoxelizeRegionsBase: wrote " << outPath << std::endl;
		}
	}
}

} // namespace fastvessels
