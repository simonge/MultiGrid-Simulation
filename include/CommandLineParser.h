#ifndef MG_COMMAND_LINE_PARSER_H
#define MG_COMMAND_LINE_PARSER_H

#include <charconv>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

namespace MG {

struct CommandLineOptions {
  std::string macro;
  std::string session;
  std::string geometryFile;
  std::string gdmlOutputFile;
  std::string outputFile;
  std::string preset;
  long seed = 0;
  long events = -1;
  bool hasOutputFile = false;
  bool help = false;
};

inline std::string RequireValue(int& index, int argc, char** argv,
                                std::string_view option,
                                std::string_view inlineValue,
                                bool hasInlineValue)
{
  if (hasInlineValue) {
    if (inlineValue.empty())
      throw std::invalid_argument("empty value for " + std::string(option));
    return std::string(inlineValue);
  }
  if (index + 1 >= argc || argv[index + 1][0] == '-')
    throw std::invalid_argument("missing value for " + std::string(option));
  return argv[++index];
}

inline long ParseLong(std::string_view value, std::string_view option)
{
  long result = 0;
  const auto* begin = value.data();
  const auto* end = begin + value.size();
  const auto parsed = std::from_chars(begin, end, result);
  if (parsed.ec != std::errc{} || parsed.ptr != end)
    throw std::invalid_argument("invalid integer for " + std::string(option)
                                + ": " + std::string(value));
  return result;
}

inline CommandLineOptions ParseCommandLine(int argc, char** argv)
{
  CommandLineOptions options;
  bool hasMacro = false;
  bool hasGeometry = false;

  for (int i = 1; i < argc; ++i) {
    std::string_view argument(argv[i]);
    const auto equals = argument.find('=');
    const bool hasInlineValue = equals != std::string_view::npos;
    const std::string_view option = argument.substr(0, equals);
    const std::string_view inlineValue = equals == std::string_view::npos
        ? std::string_view{} : argument.substr(equals + 1);

    if (option == "-h" || option == "--help") {
      options.help = true;
    } else if (option == "-m" || option == "--macro") {
      options.macro = RequireValue(i, argc, argv, option, inlineValue, hasInlineValue);
      hasMacro = true;
    } else if (option == "-f" || option == "--geometry") {
      options.geometryFile = RequireValue(i, argc, argv, option, inlineValue, hasInlineValue);
      hasGeometry = true;
    } else if (option == "--export-gdml") {
      options.gdmlOutputFile = RequireValue(i, argc, argv, option, inlineValue, hasInlineValue);
    } else if (option == "-o" || option == "--output") {
      options.outputFile = RequireValue(i, argc, argv, option, inlineValue, hasInlineValue);
      options.hasOutputFile = true;
    } else if (option == "-r" || option == "--seed") {
      options.seed = ParseLong(RequireValue(i, argc, argv, option, inlineValue, hasInlineValue), option);
      if (options.seed < 0)
        throw std::invalid_argument("seed must be non-negative");
    } else if (option == "-u" || option == "--session") {
      options.session = RequireValue(i, argc, argv, option, inlineValue, hasInlineValue);
    } else if (option == "-n" || option == "--events") {
      options.events = ParseLong(RequireValue(i, argc, argv, option, inlineValue, hasInlineValue), option);
      if (options.events < 0)
        throw std::invalid_argument("event count must be non-negative");
    } else if (option == "-mf" || option == "--preset") {
      options.preset = RequireValue(i, argc, argv, option, inlineValue, hasInlineValue);
    } else {
      throw std::invalid_argument("unknown option: " + std::string(argument));
    }
  }

  if (!options.preset.empty()) {
    if (!hasMacro && options.gdmlOutputFile.empty())
      options.macro = "macros/" + options.preset + ".mac";
    if (!hasGeometry)
      options.geometryFile = "macros/" + options.preset + ".dat";
  }
  if (!options.gdmlOutputFile.empty()
      && (hasMacro || options.events >= 0
          || options.hasOutputFile || !options.session.empty()))
    throw std::invalid_argument("--export-gdml cannot be combined with a macro, event count, ROOT output, or UI session");
  if (!options.macro.empty() && options.events >= 0)
    throw std::invalid_argument("--events cannot be combined with --macro; event counts belong in the macro");
  return options;
}

inline void PrintUsage(std::ostream& output, const char* program)
{
  output << "Usage: " << program << " [options]\n"
         << "  -m, --macro FILE       Execute a Geant4 macro (batch mode)\n"
         << "  -o, --output FILE      Write ROOT output to FILE (overrides macro setting)\n"
         << "  -f, --geometry FILE    Read detector geometry parameters from FILE\n"
         << "      --export-gdml FILE Export initialized geometry to a GDML file\n"
         << "  -r, --seed INTEGER     Set random seed (0 selects a time-based seed)\n"
         << "  -n, --events INTEGER   Run this many events (without a macro)\n"
         << "  -u, --session TYPE     Select the Geant4 interactive UI session\n"
         << "      --preset NAME      Use macros/NAME.mac and macros/NAME.dat by default\n"
         << "  -h, --help             Show this help\n"
         << "\nLegacy -mf NAME is an alias for --preset NAME.\n";
}

} // namespace MG

#endif
