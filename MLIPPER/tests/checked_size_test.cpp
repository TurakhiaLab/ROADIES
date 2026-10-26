#include <cassert>
#include <limits>
#include <stdexcept>
#include <string>

#include "util/checked_size.hpp"

template <typename Function>
void expect_length_error(Function&& function, const char* context)
{
    try {
        function();
    } catch (const std::length_error& error) {
        assert(std::string(error.what()).find(context) != std::string::npos);
        return;
    }
    assert(false && "Expected checked size arithmetic to fail");
}

int main()
{
    using mlipper::util::checked_add_size;
    using mlipper::util::checked_allocation_bytes;
    using mlipper::util::checked_mul_size;
    using mlipper::util::checked_product;

    assert(checked_add_size(20, 22, "add") == 42);
    assert(checked_mul_size(6, 7, "multiply") == 42);
    assert(checked_product("product", 2, 3, 7) == 42);
    assert(checked_product("zero", 0, std::numeric_limits<size_t>::max()) == 0);

    expect_length_error([] {
        (void)checked_add_size(std::numeric_limits<size_t>::max(), 1, "add guard");
    }, "add guard");
    expect_length_error([] {
        (void)checked_mul_size(std::numeric_limits<size_t>::max(), 2, "multiply guard");
    }, "multiply guard");
    expect_length_error([] {
        (void)checked_product("product guard",
            std::numeric_limits<size_t>::max() / 2 + 1, 2, 3);
    }, "product guard");
    expect_length_error([] {
        (void)checked_allocation_bytes<double>(
            std::numeric_limits<size_t>::max(), "byte guard");
    }, "byte guard");
}
