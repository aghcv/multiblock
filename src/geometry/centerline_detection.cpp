#include "fastvessels/obj_pipeline.hpp"

#include <vtkAppendPolyData.h>
#include <vtkCellData.h>
#include <vtkCompositeDataSet.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkHyperTreeGrid.h>
#include <vtkHyperTreeGridNonOrientedGeometryCursor.h>
#include <vtkImageData.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkMassProperties.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTriangleFilter.h>
#include <vtkUnsignedCharArray.h>
#include <vtkXMLHyperTreeGridWriter.h>
#include <vtkXMLImageDataWriter.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
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

void RefineHyperTreeCell(
	vtkHyperTreeGridNonOrientedGeometryCursor* cursor,
	vtkImplicitPolyDataDistance* implicit,
	int depth,
	int maxDepth,
	vtkUnsignedCharArray* stateArray,
	vtkIntArray* levelArray) {
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

	const bool intersectsSurface = distAbs <= halfDiag;
	if (intersectsSurface && depth < maxDepth) {
		cursor->SubdivideLeaf();
		for (int child = 0; child < 8; ++child) {
			cursor->ToChild(child);
			RefineHyperTreeCell(cursor, implicit, depth + 1, maxDepth, stateArray, levelArray);
			cursor->ToParent();
		}
		return;
	}

	const vtkIdType nodeId = cursor->GetGlobalNodeIndex();
	if (nodeId == vtkHyperTreeGrid::InvalidIndex || nodeId < 0) {
		return;
	}
	const VoxelState state = ClassifyLeaf(dist, cellSize);
	EnsureArraySize(stateArray, nodeId, static_cast<double>(VoxelState::Outside));
	EnsureArraySize(levelArray, nodeId, 0.0);
	stateArray->SetValue(nodeId, static_cast<unsigned char>(state));
	levelArray->SetValue(nodeId, depth);
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

	auto levelArray = vtkSmartPointer<vtkIntArray>::New();
	levelArray->SetName("RefineLevel");
	levelArray->SetNumberOfComponents(1);

	auto implicit = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
	implicit->SetInput(surface);

	const vtkIdType treeCount = grid->GetMaxNumberOfTrees();
	for (vtkIdType treeId = 0; treeId < treeCount; ++treeId) {
		auto cursor = vtkSmartPointer<vtkHyperTreeGridNonOrientedGeometryCursor>::New();
		grid->InitializeNonOrientedGeometryCursor(cursor, treeId, true);
		if (!cursor->HasTree()) continue;
		RefineHyperTreeCell(cursor, implicit, 0, params.maxDepth, stateArray, levelArray);
	}


	grid->GetCellData()->AddArray(stateArray);
	grid->GetCellData()->AddArray(levelArray);
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


} // namespace

void CenterlineBase(vtkMultiBlockDataSet* regions) {
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
	int totalXlets = 0;
	int totalWalls = 0;
	int totalUnknown = 0;

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

		for (int g = 0; g < groupCount; ++g) {
			vtkPolyData* pd = vtkPolyData::SafeDownCast(
				regionMb->GetBlock(static_cast<unsigned int>(g)));
			if (!pd) {
				continue;
			}
			const int surfaceType = GetSurfaceTypeField(pd);
			if (surfaceType == 1) {
				++regionXlets;
				const std::string groupName = GetBlockName(regionMb, static_cast<unsigned int>(g));
				const int groupId = GetGroupId(pd, groupName);
				const double area = ComputeSurfaceArea(pd);
				xletAreas.emplace_back(groupId, area);
			} else if (surfaceType == 0) {
				++regionWalls;
			} else {
				++regionUnknown;
			}
		}

		const unsigned int cpuIndex = static_cast<unsigned int>(r) % hwThreads;
		const std::string regionName = GetBlockName(regions, static_cast<unsigned int>(r));
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

	std::cout << "CenterlineBase: regions=" << regionCount
			  << " total_xlets=" << totalXlets
			  << " total_walls=" << totalWalls
			  << " unknown_surfaces=" << totalUnknown
			  << std::endl;
}

void VoxelizeRegionsBase(vtkMultiBlockDataSet* regions, int baseResolution, int maxDepth, int insideRefineDist) {
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
		std::cout << "VoxelizeRegionsBase: wrote " << outPath << std::endl;
	}
}

} // namespace fastvessels
