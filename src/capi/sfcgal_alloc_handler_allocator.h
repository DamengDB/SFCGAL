#pragma once
#include <new>
#include <SFCGAL/capi/sfcgal_c.h>

namespace SFCGAL
{

template<class T>
struct sfcgal_alloc_handler_allocator
{
    typedef T value_type;
    sfcgal_alloc_handler_allocator() noexcept {}
    T* allocate(const size_t n) const;
    void deallocate(T* const p, size_t) const noexcept;
};

template <class T>
T* sfcgal_alloc_handler_allocator<T>::allocate(const size_t n) const
{
    void* const pv = sfcgal_alloc_handler(n * sizeof(T));
    if (!pv) { throw std::bad_alloc(); }
    return static_cast<T*>(pv);
}

template <class T>
void sfcgal_alloc_handler_allocator<T>::deallocate(T* const p,size_t) const noexcept
{
    sfcgal_free_handler(p);
}

}
