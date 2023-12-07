#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
#  define EXTERNC       extern "C"
#  define EXTERNC_BEGIN extern "C" {
#  define EXTERNC_END   }
#else
#  define EXTERNC
#  define EXTERNC_BEGIN
#  define EXTERNC_END
#endif

#define SUCCESS 0
#define ERROR                     (-1)
#define ERROR_ILLEGAL_PARAMETER   (-2)
#define ERROR_NTT_TABLE_EMPTY     (-3)
#define ERROR_MEMORY_ALLOCATION   (-4)
#define ERROR_OUT_OF_NTT_RANGE    (-5)

#define GUARD(func)         \
  {                         \
    size_t ret = (func);    \
    if(SUCCESS != ret) {     \
      return ret;           \
    }                       \
  }

#define GUARD_MSG(func, msg) \
  {                          \
    if(SUCCESS != (func)) {  \
      printf(msg);           \
      return ERROR;          \
    }                        \
  }

#if defined(__GNUC__) || defined(__clang__)
#  define UNUSED __attribute__((unused))
#else
#  define UNUSED
#endif

#define MASK 0xffffffffffffffff

#define ALIGN(n) __attribute__((aligned(n)))

#define MAX_POLY_NUMS         1024
#define MAX_NTT_NUM           2048
#define PINGPONG_BUFFER_SIZE  4096

