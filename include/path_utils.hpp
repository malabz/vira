#pragma once
#include "common.hpp"

#include <string>

std::string join_path(const std::string& a, const std::string& b);
std::string path_basename(const std::string& path);
std::string path_dirname(const std::string& path);
std::string path_without_extension(const std::string& path);
bool file_exists(const std::string& path);
std::string absolute_path(const std::string& path);
std::string prepare_temp_dir(const Args& args, const std::string& final_out);
std::string output_prefix(const std::string& final_out);
std::string shell_quote(const std::string& s);
