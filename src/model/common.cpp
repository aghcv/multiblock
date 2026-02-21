#include "fastvessels/common.hpp"

#include <algorithm>
#include <cctype>
#include <numeric>
#include <sstream>

namespace fastvessels {

std::string ToLower(std::string value) {
	std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
		return static_cast<char>(std::tolower(c));
	});
	return value;
}

std::string Trim(std::string value) {
	const char* whitespace = " \t\n\r";
	auto start = value.find_first_not_of(whitespace);
	if (start == std::string::npos) {
		return "";
	}
	auto end = value.find_last_not_of(whitespace);
	return value.substr(start, end - start + 1);
}

bool ParseBool(const std::string& value, bool defaultValue) {
	std::string lowered = ToLower(value);
	if (lowered == "true" || lowered == "1" || lowered == "yes") {
		return true;
	}
	if (lowered == "false" || lowered == "0" || lowered == "no") {
		return false;
	}
	return defaultValue;
}

long long ParseInt(const std::string& value, long long defaultValue) {
	try {
		size_t idx = 0;
		long long v = std::stoll(value, &idx);
		if (idx != value.size()) {
			return defaultValue;
		}
		return v;
	} catch (...) {
		return defaultValue;
	}
}

double ParseDouble(const std::string& value, double defaultValue) {
	try {
		size_t idx = 0;
		double v = std::stod(value, &idx);
		if (idx != value.size()) {
			return defaultValue;
		}
		return v;
	} catch (...) {
		return defaultValue;
	}
}

static std::string TruncateOrPad(const std::string& value, size_t width) {
	if (width == 0) {
		return "";
	}
	if (value.size() > width) {
		if (width <= 3) {
			return value.substr(0, width);
		}
		return value.substr(0, width - 3) + "...";
	}
	if (value.size() < width) {
		return value + std::string(width - value.size(), ' ');
	}
	return value;
}

std::string BuildTextTable(const std::vector<std::string>& headers,
	const std::vector<std::vector<std::string>>& rows,
	const std::vector<size_t>& columnWidths,
	size_t maxRows) {
	if (columnWidths.empty() || headers.empty() || headers.size() != columnWidths.size()) {
		return "";
	}
	const size_t colCount = columnWidths.size();
	const size_t rowCount = rows.size();
	const size_t showRows = std::min(maxRows, rowCount);
	const size_t rowWidth = std::accumulate(columnWidths.begin(), columnWidths.end(), size_t{0}) + (colCount - 1);
	std::string output;
	output.reserve((showRows + 2) * (rowWidth + 1));

	auto appendRow = [&](const std::vector<std::string>& cols) {
		for (size_t i = 0; i < colCount; ++i) {
			const std::string value = (i < cols.size() ? cols[i] : "");
			output += TruncateOrPad(value, columnWidths[i]);
			if (i + 1 < colCount) {
				output.push_back(' ');
			}
		}
		output.push_back('\n');
	};

	appendRow(headers);
	std::vector<std::string> sep(colCount);
	for (size_t i = 0; i < colCount; ++i) {
		sep[i] = std::string(columnWidths[i], '-');
	}
	appendRow(sep);

	for (size_t r = 0; r < showRows; ++r) {
		appendRow(rows[r]);
	}

	return output;
}

std::uint64_t SplitWorkBegin(std::uint64_t total, int partIndex, int partCount) {
	std::uint64_t base = total / static_cast<std::uint64_t>(partCount);
	std::uint64_t rem = total % static_cast<std::uint64_t>(partCount);
	std::uint64_t begin = static_cast<std::uint64_t>(partIndex) * base + std::min<std::uint64_t>(static_cast<std::uint64_t>(partIndex), rem);
	return begin;
}

std::uint64_t SplitWorkEnd(std::uint64_t total, int partIndex, int partCount) {
	std::uint64_t begin = SplitWorkBegin(total, partIndex, partCount);
	std::uint64_t base = total / static_cast<std::uint64_t>(partCount);
	std::uint64_t rem = total % static_cast<std::uint64_t>(partCount);
	std::uint64_t size = base + (static_cast<std::uint64_t>(partIndex) < rem ? 1 : 0);
	return begin + size;
}

double ComputeKernel(std::uint64_t seed, std::uint32_t innerIters) {
	double x = 0.123456789 + (static_cast<double>(seed % 1024) * 0.0009765625);
	double y = 0.987654321 - (static_cast<double>((seed / 1024) % 1024) * 0.00048828125);
	for (std::uint32_t i = 0; i < innerIters; ++i) {
		x = x * 1.0000001192092896 + y * 0.9999999403953552 + 0.0000003;
		y = y * 1.0000002384185791 - x * 0.9999998211860657 + 0.0000001;
		if (x > 1e6) x *= 1e-6;
		if (y < -1e6) y *= 1e-6;
	}
	return x + y;
}

} // namespace fastvessels
