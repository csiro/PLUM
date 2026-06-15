/**
 * @file builddate.h
 * @brief Build timestamp definition for PLUM metabolic model gap-filling framework
 *
 * This file defines a compile-time constant containing the build date and time.
 * It is intentionally touched during each compilation to ensure the timestamp
 * reflects the actual build time of the executable.
 */
#pragma once

/*
  This file is touched each time the main is compiled,
  so that compilation date and time are updated
*/

/**
 * @brief Compile-time build timestamp string
 *
 * This constant contains the date and time when the main executable was compiled,
 * formatted as "MMM DD YYYY HH:MM:SS". The timestamp is generated using the
 * preprocessor macros __DATE__ and __TIME__ at compilation time.
 *
 * @note This file is touched before each main compilation to ensure the timestamp
 *       accurately reflects the current build time.
 */
constexpr const char* _build_date = __DATE__ " " __TIME__;
