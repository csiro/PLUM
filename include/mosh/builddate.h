#pragma once

/*
  This file is touched each time the main is compiled,
  so that compilation date and time are updated
*/

constexpr const char* _build_date = __DATE__ " " __TIME__;
