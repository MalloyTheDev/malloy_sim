#include <test_check.hpp>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

// Issue #13: the module graph is enforced at LINK time and nowhere else.
//
// Every module declares `target_include_directories(<mod> PUBLIC include)`, so
// the whole `include/` tree is on every consumer's search path. A module can
// therefore include another module's header without linking it, and the build
// succeeds. Nothing catches it, and the rule it breaks is the one the whole
// architecture rests on: domains share only malloy::math and the sim_core
// vocabulary, and no domain knows that any other domain exists (ADR 0006).
//
// This checks the source against the build files, so the two cannot drift.
//
// Since M41 it also keeps the graph pointing one way (every library on a
// layer, linking only lower ones), keeps third-party code out of every
// library (a library links only malloy:: targets and includes only malloy and
// standard headers, and no graphics header at all), and keeps the domains
// from knowing the scenario and render layers built over them. The 3D
// viewer, which is the one target that may use a graphics library, is an
// executable, so it is not a module here.
namespace
{
namespace fs = std::filesystem;

std::string read_file(const fs::path& path)
{
    std::ifstream in(path);
    return std::string((std::istreambuf_iterator<char>(in)),
                       std::istreambuf_iterator<char>());
}

// The text between a `target_link_libraries(` and its matching close paren.
// Both spellings appear in this project: INTERFACE targets put it on one line,
// STATIC ones spread it over several.
std::vector<std::pair<std::string, std::string>> link_blocks(const std::string& text)
{
    std::vector<std::pair<std::string, std::string>> blocks;
    const std::string call = "target_link_libraries(";
    for (std::size_t at = text.find(call); at != std::string::npos;
         at = text.find(call, at + 1))
    {
        std::size_t open = at + call.size();
        std::size_t depth = 1;
        std::size_t scan = open;
        while (scan < text.size() && depth > 0)
        {
            if (text[scan] == '(')
            {
                ++depth;
            }
            else if (text[scan] == ')')
            {
                --depth;
            }
            ++scan;
        }
        const std::string body = text.substr(open, scan - open - 1);

        // The target name is the first token.
        std::size_t end = body.find_first_of(" \t\r\n");
        if (end == std::string::npos)
        {
            continue;
        }
        blocks.emplace_back(body.substr(0, end), body);
    }
    return blocks;
}

// Every `malloy::<name>` in a block, minus the two that carry only flags.
std::set<std::string> aliases_in(const std::string& body)
{
    std::set<std::string> names;
    const std::string prefix = "malloy::";
    for (std::size_t at = body.find(prefix); at != std::string::npos;
         at = body.find(prefix, at + 1))
    {
        std::size_t start = at + prefix.size();
        std::size_t end = start;
        while (end < body.size() &&
               (std::isalnum(static_cast<unsigned char>(body[end])) != 0 ||
                body[end] == '_'))
        {
            ++end;
        }
        const std::string name = body.substr(start, end - start);
        if (name != "project_options" && name != "project_warnings" && !name.empty())
        {
            names.insert(name);
        }
    }
    return names;
}

// Every `#include <malloy/<module>/...>` in one file.
std::set<std::string> malloy_includes(const std::string& text)
{
    std::set<std::string> modules;
    const std::string prefix = "#include <malloy/";
    for (std::size_t at = text.find(prefix); at != std::string::npos;
         at = text.find(prefix, at + 1))
    {
        const std::size_t start = at + prefix.size();
        const std::size_t slash = text.find('/', start);
        if (slash != std::string::npos)
        {
            modules.insert(text.substr(start, slash - start));
        }
    }
    return modules;
}

std::vector<fs::path> sources_of(const fs::path& root, const std::string& module)
{
    std::vector<fs::path> files;
    for (const fs::path& dir :
         {root / "include" / "malloy" / module, root / "src" / module})
    {
        if (!fs::exists(dir))
        {
            continue;
        }
        for (const auto& entry : fs::directory_iterator(dir))
        {
            if (entry.is_regular_file())
            {
                files.push_back(entry.path());
            }
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}
// The libraries: the targets made by add_library(), read from the build file
// rather than guessed from target names (M41). An executable such as the 3D
// viewer, a test or a fuzzer is not a module, whatever it is called, and the
// two INTERFACE targets that carry only flags are not modules either.
std::set<std::string> library_names(const std::string& cmake)
{
    std::set<std::string> names;
    const std::string call = "add_library(";
    const std::string prefix = "malloy_";
    for (std::size_t at = cmake.find(call); at != std::string::npos;
         at = cmake.find(call, at + 1))
    {
        const std::size_t close = cmake.find(')', at);
        const std::string body = cmake.substr(at + call.size(), close - at - call.size());
        const std::size_t start = body.find_first_not_of(" \t\r\n");
        if (start == std::string::npos)
        {
            continue;
        }
        const std::size_t end = body.find_first_of(" \t\r\n", start);
        const std::string name = body.substr(start, end - start);
        if (name.rfind(prefix, 0) == 0 && body.find("ALIAS") == std::string::npos &&
            name != "malloy_project_options" && name != "malloy_project_warnings")
        {
            names.insert(name.substr(prefix.size()));
        }
    }
    return names;
}

// The tokens of a link block after the target name, with CMake comments
// dropped.
std::vector<std::string> link_tokens(const std::string& body)
{
    std::vector<std::string> tokens;
    std::size_t scan = body.find_first_of(" \t\r\n"); // past the target name
    while (scan != std::string::npos && scan < body.size())
    {
        const std::size_t start = body.find_first_not_of(" \t\r\n", scan);
        if (start == std::string::npos)
        {
            break;
        }
        if (body[start] == '#')
        {
            scan = body.find('\n', start);
            continue;
        }
        const std::size_t stop = body.find_first_of(" \t\r\n", start);
        tokens.push_back(body.substr(start, stop - start));
        scan = stop;
    }
    return tokens;
}

// Every header a file includes, as written between the brackets or quotes,
// with whether it was quoted.
std::vector<std::pair<std::string, bool>> includes_of(const std::string& text)
{
    std::vector<std::pair<std::string, bool>> headers;
    const std::string directive = "#include";
    for (std::size_t at = text.find(directive); at != std::string::npos;
         at = text.find(directive, at + 1))
    {
        const std::size_t open = text.find_first_not_of(" \t", at + directive.size());
        if (open == std::string::npos || (text[open] != '<' && text[open] != '"'))
        {
            continue;
        }
        const char closing = text[open] == '<' ? '>' : '"';
        const std::size_t close = text.find(closing, open + 1);
        if (close == std::string::npos)
        {
            continue;
        }
        headers.emplace_back(text.substr(open + 1, close - open - 1), text[open] == '"');
    }
    return headers;
}

// Every file under a directory, at any depth, in a fixed order.
std::vector<fs::path> files_under(const fs::path& dir)
{
    std::vector<fs::path> files;
    if (fs::exists(dir))
    {
        for (const auto& entry : fs::recursive_directory_iterator(dir))
        {
            if (entry.is_regular_file())
            {
                files.push_back(entry.path());
            }
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

// A graphics or windowing header: raylib's (raylib.h, rlgl.h, raymath.h) or
// any OpenGL or GLFW header, wherever it is spelled from.
bool is_graphics_header(const std::string& header)
{
    const std::size_t slash = header.find_last_of('/');
    const std::string name = slash == std::string::npos ? header : header.substr(slash + 1);
    return name == "raylib.h" || name == "rlgl.h" || name == "raymath.h" ||
           header.rfind("GL/", 0) == 0 || header.rfind("GLFW/", 0) == 0 ||
           header.find("/GL/") != std::string::npos ||
           header.find("/GLFW/") != std::string::npos;
}
} // namespace

int main()
{
    const fs::path root{MALLOY_REPO_ROOT};
    const std::string cmake = read_file(root / "CMakeLists.txt");
    MALLOY_CHECK_TRUE(!cmake.empty());

    // --- The libraries, from add_library(). A graph this small going empty
    //     would make everything below vacuous. ---
    const std::set<std::string> libraries = library_names(cmake);
    MALLOY_CHECK_TRUE(libraries.size() >= 11);
    MALLOY_CHECK_TRUE(libraries.count("charges") == 1);
    MALLOY_CHECK_TRUE(libraries.count("scenario") == 1);
    MALLOY_CHECK_TRUE(libraries.count("nbody_terminal") == 0);

    // --- The declared graph, read out of the build files rather than kept in
    //     a table here that could drift from them.
    //
    //     M41: a library links nothing but other malloy libraries. Graphics,
    //     windowing and any other third-party code belong to an executable
    //     (the 3D viewer), never to a library another project could link
    //     (issues #80 and #82). ---
    std::map<std::string, std::set<std::string>> declared;
    for (const auto& block : link_blocks(cmake))
    {
        const std::string& target = block.first;
        const std::string prefix = "malloy_";
        if (target.rfind(prefix, 0) != 0 ||
            libraries.count(target.substr(prefix.size())) == 0)
        {
            continue; // library targets only
        }
        declared[target.substr(prefix.size())] = aliases_in(block.second);

        for (const std::string& token : link_tokens(block.second))
        {
            if (token == "PUBLIC" || token == "PRIVATE" || token == "INTERFACE" ||
                token.rfind("malloy::", 0) == 0)
            {
                continue;
            }
            std::cerr << target << " links " << token
                      << "; a library may link only malloy:: libraries" << '\n';
            return 1;
        }
    }
    // Every library declares its links, so none escapes the checks below.
    MALLOY_CHECK_EQ(declared.size(), libraries.size());

    // --- Every module's includes must be covered by its declared links.
    //
    //     This is the check the link step cannot make: the whole include tree
    //     is on every consumer's path, so an undeclared include compiles. ---
    for (const auto& entry : declared)
    {
        const std::string& module = entry.first;
        const std::set<std::string>& allowed = entry.second;

        for (const fs::path& file : sources_of(root, module))
        {
            for (const std::string& included : malloy_includes(read_file(file)))
            {
                if (included == module || allowed.count(included) == 1)
                {
                    continue;
                }
                std::cerr << "malloy_" << module << " includes <malloy/" << included
                          << "/...> in " << file.filename().string()
                          << " but does not link malloy::" << included << '\n';
                return 1;
            }
        }
    }

    // --- ADR 0006's central rule, stated directly: no physics domain may
    //     reference another. The check above would allow it if someone added
    //     the link, so this one is about the architecture rather than the
    //     build files agreeing with themselves.
    //
    //     M41: nor may a domain reference the layers built over the domains,
    //     the scenario loader and its run layer, or the 3D render model. They
    //     know every domain; a domain that knew them would know them all. ---
    const std::vector<std::string> domains = {"nbody", "particles", "rigid",
                                              "springs", "charges"};
    for (const std::string& module : domains)
    {
        MALLOY_CHECK_TRUE(declared.count(module) == 1);
        std::vector<std::string> forbidden = {"scenario", "view3d"};
        for (const std::string& other : domains)
        {
            if (other != module)
            {
                forbidden.push_back(other);
            }
        }
        for (const std::string& other : forbidden)
        {
            if (declared[module].count(other) == 1)
            {
                std::cerr << "malloy_" << module << " links malloy::" << other
                          << ", but no domain may know another module that "
                             "knows every domain, nor another domain (ADR 0006)"
                          << '\n';
                return 1;
            }
            for (const fs::path& file : sources_of(root, module))
            {
                if (malloy_includes(read_file(file)).count(other) == 1)
                {
                    std::cerr << "malloy_" << module << " includes <malloy/" << other
                              << "/...> in " << file.filename().string()
                              << ", but no domain may know another module that "
                                 "knows every domain, nor another domain (ADR 0006)"
                              << '\n';
                    return 1;
                }
            }
        }

        // And each domain really does rest on the shared vocabulary, so the
        // loop above is checking something rather than passing on empty sets.
        MALLOY_CHECK_TRUE(declared[module].count("math") == 1);
        MALLOY_CHECK_TRUE(declared[module].count("sim_core") == 1);
    }

    // --- malloy_sim_core must stay tiny, which CLAUDE.md rule 13 states as
    //     "adding a domain must not grow it". It may depend on math and
    //     nothing else. ---
    MALLOY_CHECK_TRUE(declared.count("sim_core") == 1);
    for (const std::string& dep : declared["sim_core"])
    {
        if (dep != "math")
        {
            std::cerr << "malloy_sim_core depends on malloy::" << dep
                      << "; it may depend on malloy::math and nothing else" << '\n';
            return 1;
        }
    }

    // --- M41: the dependency direction points one way, with nothing pointing
    //     upward (ADR 0006). Every library sits on a layer and links only
    //     libraries on lower layers, so no two libraries on one layer know each
    //     other and nothing reaches back up. A new library has to be given a
    //     layer here before the build passes, which is the moment to decide
    //     where it belongs.
    //
    //     0: the vocabulary. 1: the support libraries. 2: the five domains.
    //     3: the layers that know every domain, the scenario loader with its
    //     run layer and the 3D render model, which do not know each other. ---
    const std::map<std::string, int> layer = {
        {"math", 0},
        {"time", 1},     {"sim_core", 1},  {"collide", 1}, {"ascii", 1},
        {"nbody", 2},    {"particles", 2}, {"rigid", 2},   {"springs", 2}, {"charges", 2},
        {"scenario", 3}, {"view3d", 3},
    };
    for (const auto& entry : declared)
    {
        const auto own = layer.find(entry.first);
        if (own == layer.end())
        {
            std::cerr << "malloy_" << entry.first
                      << " has no layer; give it one in layering_tests.cpp" << '\n';
            return 1;
        }
        for (const std::string& dep : entry.second)
        {
            const auto other = layer.find(dep);
            if (other == layer.end() || other->second >= own->second)
            {
                std::cerr << "malloy_" << entry.first << " (layer " << own->second
                          << ") links malloy::" << dep
                          << ", which is not on a lower layer" << '\n';
                return 1;
            }
        }
    }

    // --- M41: no library source includes a graphics or windowing header, or
    //     any header that is neither a malloy header nor a standard one (no
    //     '.' and no '/' in its name), and none uses a quoted include. Every
    //     file under include/ and src/ is read, at any depth, so a header that
    //     belongs to no listed module is covered too. Rendering lives only in
    //     the viewer executable (issue #80). ---
    std::size_t scanned = 0;
    for (const char* tree : {"include", "src"})
    {
        for (const fs::path& file : files_under(root / tree))
        {
            ++scanned;
            const std::string name = fs::relative(file, root).generic_string();
            for (const auto& [header, quoted] : includes_of(read_file(file)))
            {
                if (is_graphics_header(header))
                {
                    std::cerr << name << " includes " << header
                              << "; a library never includes a graphics header, "
                                 "only the viewer does" << '\n';
                    return 1;
                }
                if (quoted)
                {
                    std::cerr << name << " uses a quoted include of " << header
                              << "; a library includes <malloy/...> and standard "
                                 "headers only" << '\n';
                    return 1;
                }
                if (header.rfind("malloy/", 0) != 0 &&
                    header.find_first_of("./") != std::string::npos)
                {
                    std::cerr << name << " includes <" << header
                              << ">; a library includes <malloy/...> and standard "
                                 "headers only" << '\n';
                    return 1;
                }
            }
        }
    }
    MALLOY_CHECK_TRUE(scanned >= 40);

    std::cout << "malloy_layering_tests passed\n";
    return 0;
}
