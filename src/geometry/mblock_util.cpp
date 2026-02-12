#include "fastvessels/obj_pipeline.hpp"

#include "fastvessels/common.hpp"

#include <vtkActor.h>
#include <vtkActorCollection.h>
#include <vtkAppendPolyData.h>
#include <vtkAbstractArray.h>
#include <vtkDataArray.h>
#include <vtkFieldData.h>
#include <vtkCompositeDataIterator.h>
#include <vtkCompositeDataSet.h>
#include <vtkConnectivityFilter.h>
#include <vtkCellData.h>
#include <vtkExtractCells.h>
#include <vtkFeatureEdges.h>
#include <vtkGeometryFilter.h>
#include <vtkInformation.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkMapper.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkOBJImporter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataNormals.h>

#include <cctype>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <unordered_map>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>
#include <vtkPointData.h>
#include <vtkPointDataToCellData.h>
#include <vector>

namespace fastvessels {

static std::string GuessMTLPathFromOBJ(const std::string& objPath) {
	namespace fs = std::filesystem;
	fs::path p(objPath);
	fs::path mtl = p;
	mtl.replace_extension(".mtl");
	if (fs::exists(mtl)) return mtl.string();
	return "";
}

static std::string JoinArrayNames(vtkDataSetAttributes* attrs);
static std::string JoinFieldArrayNames(vtkFieldData* attrs);
static void LogDataArrays(const std::string& label, vtkPolyData* pd);
static bool PromoteGroupArrayFromFieldData(
	vtkPolyData* pd,
	const std::vector<std::string>& candidates,
	std::string& resolvedName);
static void SetBlockName(vtkMultiBlockDataSet* mb, unsigned int idx, const std::string& name);

vtkSmartPointer<vtkMultiBlockDataSet> ReadGeometry_AsMultiBlock(const std::string& path) {
	namespace fs = std::filesystem;

	if (!fs::exists(path)) {
		throw std::runtime_error("Geometry file not found: " + path);
	}

	std::string ext = fs::path(path).extension().string();
	for (auto& c : ext) c = static_cast<char>(std::tolower(c));

	if (ext == ".vtp") {
		auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
		reader->SetFileName(path.c_str());
		reader->Update();

		vtkPolyData* pd = reader->GetOutput();
		if (!pd || pd->GetNumberOfPoints() == 0) {
			throw std::runtime_error("Failed to read VTP: " + path);
		}

		LogDataArrays("VTP read", pd);
		std::string resolvedName;
		std::vector<std::string> candidates = {
			"GroupId",
			"GroupIds",
			"groupIds",
			"group_id",
			"groupId",
			"group_ids",
			"GroupID",
			"groupID"
		};
		if (PromoteGroupArrayFromFieldData(pd, candidates, resolvedName)) {
			std::cerr << "Info: promoted Group array '" << resolvedName
					  << "' from FieldData to DataSet attributes.\n";
			LogDataArrays("VTP after field promotion", pd);
		}

		auto blocks = vtkSmartPointer<vtkMultiBlockDataSet>::New();
		blocks->SetNumberOfBlocks(1);
		blocks->SetBlock(0, pd);
		blocks->GetMetaData(0u)->Set(vtkCompositeDataSet::NAME(), "geometry");
		return blocks;
	}

	if (ext != ".obj") {
		throw std::runtime_error("Unsupported geometry format: " + ext);
	}

	auto importer = vtkSmartPointer<vtkOBJImporter>::New();
	importer->SetFileName(path.c_str());

	std::string mtlPath = GuessMTLPathFromOBJ(path);
	if (!mtlPath.empty()) {
		importer->SetFileNameMTL(mtlPath.c_str());
		fs::path texDir = fs::path(path).parent_path();
		importer->SetTexturePath(texDir.string().c_str());
	}

	importer->Update();

	vtkRenderer* ren = importer->GetRenderer();
	if (!ren) {
		throw std::runtime_error("OBJ importer did not produce a renderer.");
	}

	auto blocks = vtkSmartPointer<vtkMultiBlockDataSet>::New();

	vtkActorCollection* actors = ren->GetActors();
	actors->InitTraversal();

	unsigned int blockIdx = 0;
	for (vtkActor* a = actors->GetNextActor(); a != nullptr; a = actors->GetNextActor()) {
		vtkMapper* mapper = a->GetMapper();
		if (!mapper) continue;
		vtkPolyData* pd = vtkPolyData::SafeDownCast(mapper->GetInput());
		if (!pd || pd->GetNumberOfPoints() == 0) continue;

		auto copy = vtkSmartPointer<vtkPolyData>::New();
		copy->ShallowCopy(pd);

		blocks->SetBlock(blockIdx, copy);

		std::string name = a->GetObjectName();
		if (name.empty()) {
			name = "part_" + std::to_string(blockIdx);
		}
		blocks->GetMetaData(blockIdx)->Set(vtkCompositeDataSet::NAME(), name.c_str());

		++blockIdx;
	}

	if (blockIdx == 0) {
		throw std::runtime_error("No polydata blocks extracted from OBJ. Check OBJ/MTL integrity.");
	}

	blocks->SetNumberOfBlocks(blockIdx);
	return blocks;
}

static long long CountBoundaryEdges(vtkPolyData* pd) {
	auto edges = vtkSmartPointer<vtkFeatureEdges>::New();
	edges->SetInputData(pd);
	edges->BoundaryEdgesOn();
	edges->FeatureEdgesOff();
	edges->ManifoldEdgesOff();
	edges->NonManifoldEdgesOff();
	edges->Update();
	return static_cast<long long>(edges->GetOutput()->GetNumberOfCells());
}

static std::string GetBlockName(vtkMultiBlockDataSet* mb, unsigned int idx) {
	if (!mb) return "";
	vtkInformation* info = mb->GetMetaData(idx);
	if (!info) return "";
	const char* name = info->Get(vtkCompositeDataSet::NAME());
	return name ? std::string(name) : "";
}

static double MeanCellNormalAngle(vtkPolyData* pd) {
	if (!pd || pd->GetNumberOfCells() == 0) {
		return 3.141592653589793;
	}

	auto normalsFilter = vtkSmartPointer<vtkPolyDataNormals>::New();
	normalsFilter->SetInputData(pd);
	normalsFilter->ComputeCellNormalsOn();
	normalsFilter->ComputePointNormalsOff();
	normalsFilter->SplittingOff();
	normalsFilter->ConsistencyOn();
	normalsFilter->AutoOrientNormalsOn();
	normalsFilter->Update();

	vtkPolyData* out = normalsFilter->GetOutput();
	if (!out) {
		return 3.141592653589793;
	}

	vtkDataArray* normals = out->GetCellData()->GetNormals();
	if (!normals || normals->GetNumberOfTuples() == 0) {
		return 3.141592653589793;
	}

	double sum[3] = {0.0, 0.0, 0.0};
	const vtkIdType n = normals->GetNumberOfTuples();
	for (vtkIdType i = 0; i < n; ++i) {
		double v[3] = {0.0, 0.0, 0.0};
		normals->GetTuple(i, v);
		sum[0] += v[0];
		sum[1] += v[1];
		sum[2] += v[2];
	}
	const double len = std::sqrt(sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]);
	if (len <= 1e-8) {
		return 3.141592653589793;
	}
	sum[0] /= len;
	sum[1] /= len;
	sum[2] /= len;

	double angleSum = 0.0;
	for (vtkIdType i = 0; i < n; ++i) {
		double v[3] = {0.0, 0.0, 0.0};
		normals->GetTuple(i, v);
		double dot = v[0] * sum[0] + v[1] * sum[1] + v[2] * sum[2];
		dot = std::max(-1.0, std::min(1.0, dot));
		angleSum += std::acos(dot);
	}
	return angleSum / static_cast<double>(n);
}

static void SetSurfaceTypeField(vtkPolyData* pd, bool isXlet) {
	if (!pd) return;
	vtkFieldData* field = pd->GetFieldData();
	if (!field) return;
	field->RemoveArray("SurfaceType");

	auto arr = vtkSmartPointer<vtkIntArray>::New();
	arr->SetName("SurfaceType");
	arr->SetNumberOfComponents(1);
	arr->SetNumberOfTuples(1);
	arr->SetValue(0, isXlet ? 1 : 0);
	field->AddArray(arr);
}


static void SetRegionCountField(vtkPolyData* pd, const char* name, int value) {
	if (!pd) return;
	vtkFieldData* field = pd->GetFieldData();
	if (!field) return;
	field->RemoveArray(name);

	auto arr = vtkSmartPointer<vtkIntArray>::New();
	arr->SetName(name);
	arr->SetNumberOfComponents(1);
	arr->SetNumberOfTuples(1);
	arr->SetValue(0, value);
	field->AddArray(arr);
}

static int GetSurfaceTypeField(vtkPolyData* pd) {
	if (!pd) return -1;
	vtkFieldData* field = pd->GetFieldData();
	if (!field) return -1;
	auto* arr = vtkIntArray::SafeDownCast(field->GetAbstractArray("SurfaceType"));
	if (!arr || arr->GetNumberOfTuples() < 1) return -1;
	return arr->GetValue(0);
}


ObjPipelineStats AnalyzeClosedSurfaces(vtkMultiBlockDataSet* blocks, int cpuThreads) {
	ObjPipelineStats stats;
	if (!blocks) {
		return stats;
	}
	int blockCount = static_cast<int>(blocks->GetNumberOfBlocks());
	stats.block_count = blockCount;
	if (blockCount <= 0) {
		return stats;
	}

	int threads = std::max(1, cpuThreads);
	threads = std::min(threads, blockCount);
	std::vector<int> localClosed(static_cast<size_t>(threads), 0);
	std::vector<int> localBlocksWithOpenings(static_cast<size_t>(threads), 0);
	std::vector<long long> localBoundaryCells(static_cast<size_t>(threads), 0);
	std::vector<std::thread> workers;
	workers.reserve(static_cast<size_t>(threads));

	auto workerFn = [&](int workerIndex, int begin, int end) {
		int closed = 0;
		int openingsBlocks = 0;
		long long boundaryCells = 0;
		for (int i = begin; i < end; ++i) {
			vtkPolyData* pd = vtkPolyData::SafeDownCast(blocks->GetBlock(static_cast<unsigned int>(i)));
			if (!pd) continue;
			long long b = CountBoundaryEdges(pd);
			if (b == 0) {
				closed++;
			} else {
				openingsBlocks++;
				boundaryCells += b;
			}
		}
		localClosed[static_cast<size_t>(workerIndex)] = closed;
		localBlocksWithOpenings[static_cast<size_t>(workerIndex)] = openingsBlocks;
		localBoundaryCells[static_cast<size_t>(workerIndex)] = boundaryCells;
	};

	int base = blockCount / threads;
	int rem = blockCount % threads;
	int start = 0;
	for (int t = 0; t < threads; ++t) {
		int count = base + (t < rem ? 1 : 0);
		int end = start + count;
		workers.emplace_back(workerFn, t, start, end);
		start = end;
	}

	for (auto& w : workers) {
		w.join();
	}

	int closedTotal = 0;
	int openingsBlocksTotal = 0;
	long long boundaryCellsTotal = 0;
	for (size_t i = 0; i < localClosed.size(); ++i) {
		closedTotal += localClosed[i];
		openingsBlocksTotal += localBlocksWithOpenings[i];
		boundaryCellsTotal += localBoundaryCells[i];
	}
	stats.closed_surfaces = closedTotal;
	stats.blocks_with_openings = openingsBlocksTotal;
	stats.boundary_edge_cells = boundaryCellsTotal;
	return stats;
}

void AnalyzeRegionGroupSurfaces(vtkMultiBlockDataSet* regions, int maxRegionThreads, bool unifyWalls) {
	if (!regions) {
		return;
	}

	const int regionCount = static_cast<int>(regions->GetNumberOfBlocks());
	if (regionCount <= 0) {
		return;
	}

	const int regionThreads = std::max(1, std::min(regionCount, maxRegionThreads));
	std::vector<int> localGroups(static_cast<size_t>(regionThreads), 0);
	std::vector<int> localXlets(static_cast<size_t>(regionThreads), 0);
	std::vector<int> localWalls(static_cast<size_t>(regionThreads), 0);
	std::vector<std::thread> workers;
	workers.reserve(static_cast<size_t>(regionThreads));

	const double flatAngleRad = 0.20;

	auto workerFn = [&](int workerIndex, int begin, int end) {
		int groups = 0;
		int xlets = 0;
		int walls = 0;
		for (int r = begin; r < end; ++r) {
			auto regionMb = vtkMultiBlockDataSet::SafeDownCast(regions->GetBlock(static_cast<unsigned int>(r)));
			if (!regionMb) continue;
			const int groupCount = static_cast<int>(regionMb->GetNumberOfBlocks());
			if (groupCount <= 0) continue;

			int innerThreads = std::max(1, std::min(groupCount, 4));
			std::vector<int> innerXlets(static_cast<size_t>(innerThreads), 0);
			std::vector<int> innerWalls(static_cast<size_t>(innerThreads), 0);
			std::vector<std::thread> innerWorkers;
			innerWorkers.reserve(static_cast<size_t>(innerThreads));

			int base = groupCount / innerThreads;
			int rem = groupCount % innerThreads;
			int start = 0;
			for (int t = 0; t < innerThreads; ++t) {
				int count = base + (t < rem ? 1 : 0);
				int gEnd = start + count;
				innerWorkers.emplace_back([&, t, start, gEnd]() {
					int xletCount = 0;
					int wallCount = 0;
					for (int g = start; g < gEnd; ++g) {
						vtkPolyData* pd = vtkPolyData::SafeDownCast(
							regionMb->GetBlock(static_cast<unsigned int>(g)));
						if (!pd) continue;
						const double meanAngle = MeanCellNormalAngle(pd);
						const bool isXlet = meanAngle <= flatAngleRad;
						SetSurfaceTypeField(pd, isXlet);
						if (isXlet) {
							++xletCount;
						} else {
							++wallCount;
						}
					}
					innerXlets[static_cast<size_t>(t)] = xletCount;
					innerWalls[static_cast<size_t>(t)] = wallCount;
				});
				start = gEnd;
			}

			for (auto& w : innerWorkers) {
				w.join();
			}

			int xletTotal = 0;
			int wallTotal = 0;
			for (int t = 0; t < innerThreads; ++t) {
				xletTotal += innerXlets[static_cast<size_t>(t)];
				wallTotal += innerWalls[static_cast<size_t>(t)];
			}

			groups += groupCount;
			xlets += xletTotal;
			walls += wallTotal;

			const std::string regionName = GetBlockName(regions, static_cast<unsigned int>(r));
			std::cout << "Region " << (regionName.empty() ? std::to_string(r) : regionName)
					  << ": groups=" << groupCount
					  << " xlets=" << xletTotal
					  << " walls=" << (unifyWalls ? (wallTotal > 0 ? 1 : 0) : wallTotal)
					  << std::endl;
		}
		localGroups[static_cast<size_t>(workerIndex)] = groups;
		localXlets[static_cast<size_t>(workerIndex)] = xlets;
		localWalls[static_cast<size_t>(workerIndex)] = walls;
	};

	int base = regionCount / regionThreads;
	int rem = regionCount % regionThreads;
	int start = 0;
	for (int t = 0; t < regionThreads; ++t) {
		int count = base + (t < rem ? 1 : 0);
		int end = start + count;
		workers.emplace_back(workerFn, t, start, end);
		start = end;
	}

	for (auto& w : workers) {
		w.join();
	}

	int totalGroups = 0;
	int totalXlets = 0;
	int totalWalls = 0;
	for (int t = 0; t < regionThreads; ++t) {
		totalGroups += localGroups[static_cast<size_t>(t)];
		totalXlets += localXlets[static_cast<size_t>(t)];
		totalWalls += localWalls[static_cast<size_t>(t)];
	}

	int unifiedWallRegions = 0;
	if (unifyWalls) {
		for (int r = 0; r < regionCount; ++r) {
			auto regionMb = vtkMultiBlockDataSet::SafeDownCast(regions->GetBlock(static_cast<unsigned int>(r)));
			if (!regionMb) continue;
			const int groupCount = static_cast<int>(regionMb->GetNumberOfBlocks());
			if (groupCount <= 0) continue;

			std::vector<vtkSmartPointer<vtkPolyData>> xletBlocks;
			std::vector<std::string> xletNames;
			auto wallAppender = vtkSmartPointer<vtkAppendPolyData>::New();
			int wallCount = 0;
			int xletCount = 0;

			for (int g = 0; g < groupCount; ++g) {
				auto* pd = vtkPolyData::SafeDownCast(regionMb->GetBlock(static_cast<unsigned int>(g)));
				if (!pd) continue;
				const int surfaceType = GetSurfaceTypeField(pd);
				if (surfaceType == 1) {
					xletBlocks.emplace_back(pd);
					xletNames.emplace_back(GetBlockName(regionMb, static_cast<unsigned int>(g)));
					++xletCount;
				} else {
					wallAppender->AddInputData(pd);
					++wallCount;
				}
			}

			unsigned int outCount = static_cast<unsigned int>(xletBlocks.size() + (wallCount > 0 ? 1 : 0));
			auto regionOut = vtkSmartPointer<vtkMultiBlockDataSet>::New();
			regionOut->SetNumberOfBlocks(outCount);

			unsigned int outIdx = 0;
			for (size_t i = 0; i < xletBlocks.size(); ++i) {
				auto xletCopy = vtkSmartPointer<vtkPolyData>::New();
				xletCopy->ShallowCopy(xletBlocks[i]);
				SetRegionCountField(xletCopy, "RegionXletCount", xletCount);
				SetRegionCountField(xletCopy, "RegionWallGroupCount", wallCount);
				regionOut->SetBlock(outIdx, xletCopy);
				std::string name = xletNames[i];
				if (name.empty()) {
					name = "Xlet_" + std::to_string(outIdx);
				}
				SetBlockName(regionOut, outIdx, name);
				++outIdx;
			}

			if (wallCount > 0) {
				wallAppender->Update();
				auto wallCombined = vtkSmartPointer<vtkPolyData>::New();
				wallCombined->ShallowCopy(wallAppender->GetOutput());
				SetSurfaceTypeField(wallCombined, false);
				SetRegionCountField(wallCombined, "RegionXletCount", xletCount);
				SetRegionCountField(wallCombined, "RegionWallGroupCount", wallCount);
				regionOut->SetBlock(outIdx, wallCombined);
				SetBlockName(regionOut, outIdx, "Walls");
				++unifiedWallRegions;
			}

			const std::string regionName = GetBlockName(regions, static_cast<unsigned int>(r));
			regions->SetBlock(static_cast<unsigned int>(r), regionOut);
			if (!regionName.empty()) {
				SetBlockName(regions, static_cast<unsigned int>(r), regionName);
			}
		}
	}

	std::cout << "Xlet/Wall summary: groups=" << totalGroups
			  << " xlets=" << totalXlets
			  << " walls=" << (unifyWalls ? unifiedWallRegions : totalWalls)
			  << std::endl;
}

static void SetBlockName(vtkMultiBlockDataSet* mb, unsigned int idx, const std::string& name) {
	if (!mb) return;
	vtkInformation* info = mb->GetMetaData(idx);
	info->Set(vtkCompositeDataSet::NAME(), name.c_str());
}

static inline int GetCellInt(vtkDataArray* arr, vtkIdType cellId) {
	return static_cast<int>(arr->GetComponent(cellId, 0));
}

static std::string JoinArrayNames(vtkDataSetAttributes* attrs) {
	if (!attrs) return "";
	std::ostringstream oss;
	for (int i = 0; i < attrs->GetNumberOfArrays(); ++i) {
		if (i > 0) oss << ", ";
		const char* name = attrs->GetArrayName(i);
		oss << (name ? name : "<unnamed>");
	}
	return oss.str();
}

static std::string JoinFieldArrayNames(vtkFieldData* attrs) {
	if (!attrs) return "";
	std::ostringstream oss;
	for (int i = 0; i < attrs->GetNumberOfArrays(); ++i) {
		if (i > 0) oss << ", ";
		const char* name = attrs->GetArrayName(i);
		oss << (name ? name : "<unnamed>");
	}
	return oss.str();
}

static void LogDataArrays(const std::string& label, vtkPolyData* pd) {
	if (!pd) return;
	std::cerr << "[" << label << "] CellData arrays: ["
			  << JoinArrayNames(pd->GetCellData()) << "] PointData arrays: ["
			  << JoinArrayNames(pd->GetPointData()) << "] FieldData arrays: ["
			  << JoinFieldArrayNames(pd->GetFieldData()) << "]\n";
}

static vtkSmartPointer<vtkIntArray> ConvertStringArrayToInt(vtkStringArray* strArr) {
	if (!strArr) return nullptr;
	auto out = vtkSmartPointer<vtkIntArray>::New();
	out->SetNumberOfComponents(1);
	out->SetNumberOfTuples(strArr->GetNumberOfTuples());

	std::unordered_map<std::string, int> mapping;
	int nextId = 0;
	for (vtkIdType i = 0; i < strArr->GetNumberOfTuples(); ++i) {
		const std::string value = strArr->GetValue(i);
		auto it = mapping.find(value);
		if (it == mapping.end()) {
			mapping.emplace(value, nextId);
			out->SetValue(i, nextId);
			++nextId;
		} else {
			out->SetValue(i, it->second);
		}
	}

	return out;
}

static vtkDataArray* ResolveGroupArrayOnCellData(
	vtkPolyData* pd,
	const std::vector<std::string>& candidates,
	std::string& resolvedName) {
	if (!pd) return nullptr;
	for (const auto& name : candidates) {
		if (name.empty()) continue;
		vtkAbstractArray* abs = pd->GetCellData()->GetAbstractArray(name.c_str());
		if (!abs) continue;
		if (auto* dataArr = vtkDataArray::SafeDownCast(abs)) {
			resolvedName = name;
			return dataArr;
		}
		if (auto* strArr = vtkStringArray::SafeDownCast(abs)) {
			auto conv = ConvertStringArrayToInt(strArr);
			if (conv) {
				conv->SetName(name.c_str());
				pd->GetCellData()->AddArray(conv);
				resolvedName = name;
				return conv;
			}
		}
	}
	return nullptr;
}

static vtkDataArray* ResolveGroupArrayOnPointData(
	vtkPolyData* pd,
	const std::vector<std::string>& candidates,
	std::string& resolvedName) {
	if (!pd) return nullptr;
	for (const auto& name : candidates) {
		if (name.empty()) continue;
		vtkAbstractArray* abs = pd->GetPointData()->GetAbstractArray(name.c_str());
		if (!abs) continue;
		if (auto* dataArr = vtkDataArray::SafeDownCast(abs)) {
			resolvedName = name;
			return dataArr;
		}
		if (auto* strArr = vtkStringArray::SafeDownCast(abs)) {
			auto conv = ConvertStringArrayToInt(strArr);
			if (conv) {
				conv->SetName(name.c_str());
				pd->GetPointData()->AddArray(conv);
				resolvedName = name;
				return conv;
			}
		}
	}
	return nullptr;
}

static bool PromoteGroupArrayFromFieldData(
	vtkPolyData* pd,
	const std::vector<std::string>& candidates,
	std::string& resolvedName) {
	if (!pd) return false;
	vtkFieldData* field = pd->GetFieldData();
	if (!field) return false;

	for (const auto& name : candidates) {
		if (name.empty()) continue;
		vtkAbstractArray* abs = field->GetAbstractArray(name.c_str());
		if (!abs) continue;
		vtkDataArray* dataArr = vtkDataArray::SafeDownCast(abs);
		if (!dataArr) continue;

		const vtkIdType tuples = dataArr->GetNumberOfTuples();
		if (tuples == pd->GetNumberOfCells()) {
			pd->GetCellData()->AddArray(dataArr);
			resolvedName = name;
			return true;
		}
		if (tuples == pd->GetNumberOfPoints()) {
			pd->GetPointData()->AddArray(dataArr);
			resolvedName = name;
			return true;
		}
	}

	return false;
}

static vtkDataArray* FindGroupArray(
	vtkPolyData* pd,
	const std::string& requestedName,
	std::string& resolvedName,
	bool& fromPointData) {
	if (!pd) return nullptr;

	std::vector<std::string> candidates = {
		requestedName,
		"GroupIds",
		"groupIds",
		"group_id",
		"groupId",
		"group_ids",
		"GroupID",
		"groupID"
	};

	if (auto arr = ResolveGroupArrayOnCellData(pd, candidates, resolvedName)) {
		fromPointData = false;
		return arr;
	}

	if (auto arr = ResolveGroupArrayOnPointData(pd, candidates, resolvedName)) {
		fromPointData = true;
		return arr;
	}

	resolvedName.clear();
	fromPointData = false;
	return nullptr;
}


vtkSmartPointer<vtkMultiBlockDataSet> BuildRegionSurfaceHierarchy(
	vtkMultiBlockDataSet* inputMb,
	const std::string& groupArrayName,
	const std::string& regionArrayName,
	bool addOriginalBlockId) {
	if (!inputMb) {
		throw std::runtime_error("BuildRegionSurfaceHierarchy: inputMb is null.");
	}

	auto appender = vtkSmartPointer<vtkAppendPolyData>::New();
	auto it = vtkSmartPointer<vtkCompositeDataIterator>::Take(inputMb->NewIterator());
	it->SkipEmptyNodesOn();

	int blockCounter = 0;
	for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextItem()) {
		vtkDataObject* dobj = it->GetCurrentDataObject();
		vtkPolyData* pd = vtkPolyData::SafeDownCast(dobj);
		if (!pd) continue;

		std::vector<std::string> candidates = {
			groupArrayName,
			"GroupIds",
			"groupIds",
			"group_id",
			"groupId",
			"group_ids",
			"GroupID",
			"groupID"
		};
		std::string promotedName;
		if (PromoteGroupArrayFromFieldData(pd, candidates, promotedName)) {
			std::cerr << "Info: promoted Group array '" << promotedName
					  << "' from FieldData on input block.\n";
		}

		auto tri = vtkSmartPointer<vtkTriangleFilter>::New();
		tri->SetInputData(pd);
		tri->Update();

		vtkPolyData* triPd = tri->GetOutput();
		if (!triPd || triPd->GetNumberOfCells() == 0) continue;

		const std::string cellNamesBefore = JoinArrayNames(triPd->GetCellData());
		const std::string pointNamesBefore = JoinArrayNames(triPd->GetPointData());
		const std::string fieldNamesBefore = JoinFieldArrayNames(triPd->GetFieldData());

		std::string resolvedGroupName;
		bool fromPointData = false;
		vtkDataArray* gArr = FindGroupArray(triPd, groupArrayName, resolvedGroupName, fromPointData);
		if (gArr && fromPointData) {
			auto p2c = vtkSmartPointer<vtkPointDataToCellData>::New();
			p2c->SetInputData(triPd);
			p2c->PassPointDataOn();
			p2c->Update();
			triPd = vtkPolyData::SafeDownCast(p2c->GetOutput());
			gArr = triPd ? triPd->GetCellData()->GetArray(resolvedGroupName.c_str()) : nullptr;
			fromPointData = false;
		}
		if (!gArr) {
			auto fallback = vtkSmartPointer<vtkIntArray>::New();
			fallback->SetName(groupArrayName.c_str());
			fallback->SetNumberOfComponents(1);
			fallback->SetNumberOfTuples(triPd->GetNumberOfCells());
			for (vtkIdType c = 0; c < triPd->GetNumberOfCells(); ++c) {
				fallback->SetValue(c, 0);
			}
			triPd->GetCellData()->AddArray(fallback);
			gArr = fallback;
			std::cerr << "Warning: missing GroupId array; defaulting all cells to GroupId=0. "
					  << "Available CellData arrays: [" << cellNamesBefore << "] "
					  << "Available PointData arrays: [" << pointNamesBefore << "] "
					  << "Available FieldData arrays: [" << fieldNamesBefore << "]\n";
		}
		if (gArr && !resolvedGroupName.empty() && resolvedGroupName != groupArrayName) {
			auto renamed = vtkSmartPointer<vtkDataArray>::Take(gArr->NewInstance());
			renamed->DeepCopy(gArr);
			renamed->SetName(groupArrayName.c_str());
			triPd->GetCellData()->AddArray(renamed);
			gArr = renamed;
		}
		if (!resolvedGroupName.empty() && resolvedGroupName != groupArrayName) {
			std::cerr << "Info: using GroupId array name '" << resolvedGroupName
					  << "' for block " << blockCounter << ".\n";
		}

		if (addOriginalBlockId) {
			auto ob = vtkSmartPointer<vtkIntArray>::New();
			ob->SetName("OriginalBlockId");
			ob->SetNumberOfComponents(1);
			ob->SetNumberOfTuples(triPd->GetNumberOfCells());
			for (vtkIdType c = 0; c < triPd->GetNumberOfCells(); ++c) {
				ob->SetValue(c, blockCounter);
			}
			triPd->GetCellData()->AddArray(ob);
		}

		appender->AddInputData(triPd);
		++blockCounter;
	}

	appender->Update();

	vtkPolyData* flatPd = appender->GetOutput();
	if (!flatPd || flatPd->GetNumberOfCells() == 0) {
		return vtkSmartPointer<vtkMultiBlockDataSet>::New();
	}

	auto conn = vtkSmartPointer<vtkConnectivityFilter>::New();
	conn->SetInputData(flatPd);
	conn->SetExtractionModeToAllRegions();
	conn->ColorRegionsOn();
	conn->Update();

	vtkPolyData* labeledPd = vtkPolyData::SafeDownCast(conn->GetOutput());
	if (!labeledPd || labeledPd->GetNumberOfCells() == 0) {
		return vtkSmartPointer<vtkMultiBlockDataSet>::New();
	}

	vtkDataArray* groupArr = labeledPd->GetCellData()->GetArray(groupArrayName.c_str());
	if (!groupArr) {
		std::ostringstream oss;
		oss << "After connectivity, missing CellData array '" << groupArrayName << "'.";
		throw std::runtime_error(oss.str());
	}

	vtkDataArray* regionArr = labeledPd->GetCellData()->GetArray(regionArrayName.c_str());
	if (!regionArr) {
		std::ostringstream oss;
		oss << "Missing CellData array '" << regionArrayName << "' produced by vtkConnectivityFilter.";
		throw std::runtime_error(oss.str());
	}

	std::map<int, std::map<int, vtkSmartPointer<vtkIdList>>> regionGroupToCells;

	const vtkIdType nCells = labeledPd->GetNumberOfCells();
	for (vtkIdType c = 0; c < nCells; ++c) {
		const int r = GetCellInt(regionArr, c);
		const int g = GetCellInt(groupArr, c);

		auto& groupMap = regionGroupToCells[r];
		auto& idList = groupMap[g];

		if (!idList) {
			idList = vtkSmartPointer<vtkIdList>::New();
		}
		idList->InsertNextId(c);
	}

	auto out = vtkSmartPointer<vtkMultiBlockDataSet>::New();
	out->SetNumberOfBlocks(static_cast<unsigned int>(regionGroupToCells.size()));

	unsigned int regionIdx = 0;
	for (const auto& regionPair : regionGroupToCells) {
		const int regionId = regionPair.first;
		const auto& groupMap = regionPair.second;

		auto regionMb = vtkSmartPointer<vtkMultiBlockDataSet>::New();
		regionMb->SetNumberOfBlocks(static_cast<unsigned int>(groupMap.size()));

		std::ostringstream regionName;
		regionName << "Region_" << regionId;
		SetBlockName(out, regionIdx, regionName.str());

		unsigned int groupIdx = 0;
		for (const auto& groupPair : groupMap) {
			const int groupId = groupPair.first;
			const auto& cellIds = groupPair.second;

			auto extract = vtkSmartPointer<vtkExtractCells>::New();
			extract->SetInputData(labeledPd);
			extract->SetCellList(cellIds);
			extract->Update();

			auto geom = vtkSmartPointer<vtkGeometryFilter>::New();
			geom->SetInputConnection(extract->GetOutputPort());
			geom->Update();

			vtkPolyData* leaf = geom->GetOutput();
			auto leafCopy = vtkSmartPointer<vtkPolyData>::New();
			leafCopy->ShallowCopy(leaf);

			regionMb->SetBlock(groupIdx, leafCopy);

			std::ostringstream groupName;
			groupName << "Group_" << groupId;
			SetBlockName(regionMb, groupIdx, groupName.str());

			++groupIdx;
		}

		out->SetBlock(regionIdx, regionMb);
		++regionIdx;
	}

	return out;
}

} // namespace fastvessels
