#include "wsd.h"

#include <iostream>
#include <string>
#include <vector>

namespace {

DemInfo segment(int template_index, int start, int end,
                const std::string &strand = "+")
{
    return DemInfo(template_index, "template" + std::to_string(template_index),
                   strand, 1.0F, start, end);
}

bool ordered_without_overlap(const std::vector<DemInfo> &segments)
{
    for (std::size_t i = 1; i < segments.size(); ++i)
    {
        if (segments[i - 1].end_pos > segments[i].st_pos)
            return false;
    }
    return true;
}

bool exact_boundary()
{
    std::vector<std::vector<DemInfo>> chunks{
        {segment(0, 0, 100), segment(9, 100, 110)},
        {segment(1, 100, 200)}};
    std::vector<DemInfo> merged;
    combin_dem_res(chunks, merged);
    return merged.size() == 2 && merged[0].st_pos == 0 &&
           merged[0].end_pos == 100 && merged[1].st_pos == 100 &&
           merged[1].end_pos == 200 && ordered_without_overlap(merged);
}

bool one_base_overlap()
{
    std::vector<std::vector<DemInfo>> chunks{
        {segment(0, 0, 100), segment(9, 100, 110)},
        {segment(1, 99, 199)}};
    std::vector<DemInfo> merged;
    combin_dem_res(chunks, merged);
    return merged.size() == 2 && merged[0].end_pos == 99 &&
           merged[1].st_pos == 99 && merged[1].end_pos == 199 &&
           ordered_without_overlap(merged);
}

bool small_valid_overlap()
{
    std::vector<std::vector<DemInfo>> chunks{
        {segment(0, 0, 100), segment(9, 100, 110)},
        {segment(1, 95, 195)}};
    std::vector<DemInfo> merged;
    combin_dem_res(chunks, merged);
    return merged.size() == 2 && merged[0].end_pos == 95 &&
           merged[1].st_pos == 95 && merged[1].end_pos == 195 &&
           ordered_without_overlap(merged);
}

bool duplicate_segment()
{
    std::vector<std::vector<DemInfo>> chunks{
        {segment(0, 0, 100), segment(9, 100, 110)},
        {segment(0, 0, 100), segment(1, 100, 200)}};
    std::vector<DemInfo> merged;
    combin_dem_res(chunks, merged);
    return merged.size() == 2 && merged[0].tem_idx == 0 &&
           merged[1].tem_idx == 1 && ordered_without_overlap(merged);
}

bool conflicting_overlap()
{
    std::vector<std::vector<DemInfo>> chunks{
        {segment(0, 0, 100), segment(9, 100, 110)},
        {segment(2, 0, 150)}};
    std::vector<DemInfo> merged;
    combin_dem_res(chunks, merged);
    return merged.size() == 1 && merged[0].tem_idx == 0 &&
           merged[0].st_pos == 0 && merged[0].end_pos == 100;
}

bool reverse_strand_overlap()
{
    std::vector<std::vector<DemInfo>> chunks{
        {segment(0, 0, 100, "-"), segment(9, 100, 110, "-")},
        {segment(1, 99, 199, "-")}};
    std::vector<DemInfo> merged;
    combin_dem_res(chunks, merged);
    return merged.size() == 2 && merged[1].strand == "-" &&
           merged[0].end_pos == 99 && merged[1].st_pos == 99 &&
           ordered_without_overlap(merged);
}

bool multiple_consecutive_seams()
{
    std::vector<std::vector<DemInfo>> chunks{
        {segment(0, 0, 100), segment(9, 100, 110)},
        {segment(1, 99, 199), segment(9, 199, 210)},
        {segment(2, 198, 298)}};
    std::vector<DemInfo> merged;
    combin_dem_res(chunks, merged);
    return merged.size() == 3 && merged[0].end_pos == 99 &&
           merged[1].st_pos == 99 && merged[1].end_pos == 198 &&
           merged[2].st_pos == 198 && merged[2].end_pos == 298 &&
           ordered_without_overlap(merged);
}

} // namespace

int main()
{
    struct TestCase
    {
        const char *name;
        bool (*run)();
    };

    const std::vector<TestCase> tests{
        {"exact_boundary", exact_boundary},
        {"one_base_overlap", one_base_overlap},
        {"small_valid_overlap", small_valid_overlap},
        {"duplicate_segment", duplicate_segment},
        {"conflicting_overlap", conflicting_overlap},
        {"reverse_strand_overlap", reverse_strand_overlap},
        {"multiple_consecutive_seams", multiple_consecutive_seams},
    };

    int failures = 0;
    for (const TestCase &test : tests)
    {
        const bool passed = test.run();
        std::cout << test.name << '\t' << (passed ? "PASS" : "FAIL") << '\n';
        failures += passed ? 0 : 1;
    }
    return failures == 0 ? 0 : 1;
}
