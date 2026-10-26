#pragma once

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace mlipper::util {

inline size_t checked_add_size(size_t lhs, size_t rhs, const char* context)
{
    if (rhs > std::numeric_limits<size_t>::max() - lhs) {
        throw std::length_error(std::string(context) + ": size addition overflow");
    }
    return lhs + rhs;
}

inline size_t checked_mul_size(size_t lhs, size_t rhs, const char* context)
{
    // Check with division before multiplying: unsigned overflow is defined to
    // wrap, so inspecting the product afterward cannot recover the true size.
    if (lhs != 0 && rhs > std::numeric_limits<size_t>::max() / lhs) {
        throw std::length_error(std::string(context) + ": size multiplication overflow");
    }
    return lhs * rhs;
}

template <typename... Factors>
inline size_t checked_product(const char* context, size_t first, Factors... rest)
{
    // Validate each intermediate product; checking only the final factor would
    // miss an earlier wraparound.
    size_t result = first;
    ((result = checked_mul_size(result, static_cast<size_t>(rest), context)), ...);
    return result;
}

template <typename T>
inline size_t checked_allocation_bytes(size_t count, const char* context)
{
    return checked_mul_size(count, sizeof(T), context);
}

} // namespace mlipper::util
